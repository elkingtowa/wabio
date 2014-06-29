#!/usr/bin/env python
"""Setup script for distributed phylogenetic analysis by BLAST.
"""
from setuptools import setup, find_packages

setup(name = "bcbio-phyloblast",
      version = "0.1.2",
      author = "Brad Chapman",
      author_email = "chapmanb@50mail.com",
      description = "Distributed BLAST for gene phylogeny estimation",
      license = "MIT",
      url = "http://bcbio.wordpress.com",
      namespace_packages = ["bcbio"],
      packages = find_packages(),
      scripts = ["scripts/blast_cross_orgs.py",
                 "scripts/retrieve_org_dbs.py"],
      install_requires=[
          "biopython >= 1.62",
          "PyYAML >= 3.10"]
      )
