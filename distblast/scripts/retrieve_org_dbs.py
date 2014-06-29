#!/usr/bin/env python
"""Retrieve full genome databases, preparing them for BLAST analysis.

Usage:
    retrieve_org_dbs.py <YAML config file>

Requires:
    - NCBI's blast+ -- for preparing the organism databases
      ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    - Biopython libraries
"""
import os
import sys
import csv
import glob
import ftplib
import subprocess
import contextlib
import urllib2
import socket
import time

import yaml

from Bio import Entrez

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    Entrez.email = config.get('email', 'test@example.com')
    socket.setdefaulttimeout(config['url_timeout'])
    local_get = LocalRetrieval(config)
    ncbi_get = NcbiEntrezRetrieval(config)
    ensembl_get = EnsemblFtpRetrieval(config)
    organisms = read_org_list(config['org_file'])
    db_dir = config['db_dir']
    ensembl_db_dir = os.path.join(db_dir, "ensembl")
    for check_dir in [db_dir, ensembl_db_dir]:
        if not os.path.exists(check_dir):
            os.makedirs(check_dir)
    org_files = []
    for org in organisms:
        check_glob = os.path.join(config["db_dir"], "custom", "%s*" % org)
        print "Preparing organism:", org
        check_custom = [x for x in glob.glob(check_glob)
                        if not x.endswith((".phr", ".pin", ".psq"))]
        if org in config.get('problem_orgs', []):
            db_file = ''
        elif len(check_custom) == 1:
            db_file = local_get.retrieve_db(org, check_custom[0], db_dir)
        else:
            print("Did not find single pre-downloaded FASTA file in '%s'\n"
                  "Instead Found %s\n"
                  "Attempting to download from Ensembl or NCBI" % (check_glob, check_custom))
            db_file = ensembl_get.retrieve_db(org, ensembl_db_dir)
            if db_file:
                print "Ensembl"
                db_file = os.path.join(os.path.basename(ensembl_db_dir), db_file)
            else:
                print "NCBI"
                db_file = ncbi_get.retrieve_db(org, db_dir)
        org_files.append((org, db_file))
    with open(os.path.join(db_dir, "organism_dbs.txt"), "w") as out_handle:
        for org, fname in org_files:
            out_handle.write("%s\t%s\n" % (org, fname))

def read_org_list(in_file):
    with open(in_file, 'rU') as in_handle:
        reader = csv.reader(in_handle)
        orgs = [r[-1] for r in reader]
    return orgs

class _BaseRetrieval:
    def _make_blast_db(self, db_dir, final_file, db_name, organism):
        with _chdir(db_dir):
            if not os.path.exists("%s.pin" % db_name):
                cmd = self._config.get("blastdb_cmd", "makeblastdb")
                cl = [cmd, "-in", os.path.basename(final_file),
                      "-dbtype", "prot",
                      "-out", db_name,
                      "-title", organism]
                subprocess.check_call(cl)

class LocalRetrieval(_BaseRetrieval):
    """Prepare a database file from a local FASTA ref.
    """
    def __init__(self, config):
        self._config = config

    def retrieve_db(self, org, fname, db_dir):
        self._make_blast_db(os.path.dirname(fname), os.path.basename(fname),
                            os.path.splitext(os.path.basename(fname))[0], org)
        return os.path.splitext(fname.replace("%s/" % db_dir, ""))[0]

class NcbiEntrezRetrieval(_BaseRetrieval):
    """Pull down fasta protein genome sequences using NCBI Entrez.
    """
    def __init__(self, config):
        self._max_tries = 5
        self._config = config

    def retrieve_db(self, organism, db_dir):
        genome_ids = self._query_for_ids(organism)
        out_file = os.path.join(db_dir, "%s-entrez.fa" %
                organism.replace(" ", "_"))
        db_name = os.path.splitext(os.path.basename(out_file))[0]
        if not os.path.exists(out_file):
            num_tries = 1
            while 1:
                try:
                    self._download_and_error_out(out_file, genome_ids)
                    break
                except urllib2.URLError:
                    print "Timeout error"
                    time.sleep(5)
                    if num_tries > self._max_tries:
                        raise
                    else:
                        num_tries += 1
        self._make_blast_db(db_dir, os.path.basename(out_file), db_name,
                organism)
        return db_name

    def _download_and_error_out(self, out_file, genome_ids):
        """Do the full genome downloading, raising timeout errors to be handled.
        """
        with open(out_file, "w") as out_handle:
            for genome_id in genome_ids:
                print "Downloading", genome_id
                self._download_to_file(genome_id, out_handle)

    def _download_to_file(self, genome_id, out_handle):
        entrez_url = "http://www.ncbi.nlm.nih.gov/sites/entrez?Db=genome&" \
                     "Cmd=File&dopt=Protein+FASTA&list_uids=%s" % genome_id
        download_handle = urllib2.urlopen(entrez_url)
        # read off garbage at the beginning of the file related to the genome
        while 1:
            line = download_handle.readline()
            if line.startswith(">"):
                out_handle.write(line)
                break
            if not line:
                break
            print line
        for line in download_handle:
            out_handle.write(line)
        download_handle.close()
        # be sure output has trailing newlines. Who knows what could be there.
        out_handle.write("\n")

    def _query_for_ids(self, organism):
        handle = Entrez.esearch(db="genome", term="%s[Organism]" % organism)
        record = Entrez.read(handle)
        return record['IdList']

class EnsemblFtpRetrieval(_BaseRetrieval):
    """Handle obtaining a reference genome from Ensembl
    """
    def __init__(self, config):
        self._main_ftp = "ftp://ftp.ensembl.org/pub/current_fasta/"
        self._genome_ftp = "ftp://ftp.ensemblgenomes.org/pub/%s/current/fasta/"
        self._genome_dbs = ["bacteria", "protists", "metazoa", "fungi",
                "plants"]
        self._initialized = False
        self._config = config

    def _initialize(self):
        if not self._initialized:
            urls = [self._genome_ftp % d for d in self._genome_dbs] + \
                   [self._main_ftp]
            self._org_to_urls = dict()
            for url in urls:
                orgs = self._files_at_url(url)
                for org in orgs:
                    self._org_to_urls[org] = url
            self._initialized = True

    def _files_at_url(self, url):
        """Add organisms available at the provided FTP url.
        """
        parts = url.replace("ftp://", "").split("/")
        ftp = ftplib.FTP(parts[0])
        ftp.login()
        orgs = ftp.nlst("/".join(parts[1:]))
        return [o.split("/")[-1] for o in orgs]

    def retrieve_db(self, organism, db_dir):
        self._initialize()
        ftp_url = self._get_ftp_url(organism)
        if ftp_url is None:
            return ""
        file_name = ftp_url.split("/")[-1]
        final_file = os.path.join(db_dir, file_name.replace(".gz", ""))
        db_name = os.path.splitext(os.path.basename(final_file))[0]
        if not os.path.exists(final_file):
            with _chdir(db_dir):
                cl = ["wget", ftp_url]
                subprocess.check_call(cl)
                cl = ["gunzip", file_name]
                subprocess.check_call(cl)
        self._make_blast_db(db_dir, final_file, db_name, organism)
        return db_name

    def _get_ftp_url(self, organism):
        """Retrieve the protein database link for a given organism.
        """
        ftp_url = None
        org_parts = organism.split()
        for check_org in [organism.replace(" ", "_").lower(),
                "_".join([org_parts[0][0], org_parts[1]]).lower()]:
            try:
                ftp_url = self._org_to_urls[check_org]
                break
            except KeyError:
                pass
        if ftp_url:
            ftp_url = ftp_url + check_org + "/pep/"
            files = self._files_at_url(ftp_url)
            for f in files:
                if f.endswith("pep.all.fa.gz"):
                    ftp_url = ftp_url + f
                    break
        return ftp_url

@contextlib.contextmanager
def _chdir(new_dir):
    orig_dir = os.getcwd()
    try:
        os.chdir(new_dir)
        yield
    finally:
        os.chdir(orig_dir)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        main(*sys.argv[1:])
    else:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
