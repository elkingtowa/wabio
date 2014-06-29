"""Next-gen alignments with BWA (http://bio-bwa.sourceforge.net/)
"""
import os
import subprocess

from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction

galaxy_location_file = "bwa_index.loc"

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
          rg_name=None):
    """Perform a BWA alignment, generating a SAM file.
    """
    sai1_file = os.path.join(align_dir, "%s_1.sai" % out_base)
    sai2_file = (os.path.join(align_dir, "%s_2.sai" % out_base)
                 if pair_file else None)
    sam_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not file_exists(sam_file):
        if not file_exists(sai1_file):
            with file_transaction(sai1_file) as tx_sai1_file:
                _run_bwa_align(fastq_file, ref_file, tx_sai1_file, config)
        if sai2_file and not file_exists(sai2_file):
            with file_transaction(sai2_file) as tx_sai2_file:
                _run_bwa_align(pair_file, ref_file, tx_sai2_file, config)
        align_type = "sampe" if sai2_file else "samse"
        sam_cl = [config_utils.get_program("bwa", config), align_type, ref_file, sai1_file]
        if sai2_file:
            sam_cl.append(sai2_file)
        sam_cl.append(fastq_file)
        if sai2_file:
            sam_cl.append(pair_file)
        with file_transaction(sam_file) as tx_sam_file:
            with open(tx_sam_file, "w") as out_handle:
                logger.info(" ".join(sam_cl))
                subprocess.check_call(sam_cl, stdout=out_handle)
    return sam_file

def _bwa_args_from_config(config):
    cores = config.get("resources", {}).get("bwa", {}).get("cores", None)
    core_flags = ["-t", str(cores)] if cores else []
    qual_format = config["algorithm"].get("quality_format", "").lower()
    qual_flags = ["-I"] if qual_format == "illumina" else []
    return core_flags + qual_flags

def _run_bwa_align(fastq_file, ref_file, out_file, config):
    aln_cl = [config_utils.get_program("bwa", config), "aln",
              "-n %s" % config["algorithm"]["max_errors"],
              "-k %s" % config["algorithm"]["max_errors"]]
    aln_cl += _bwa_args_from_config(config)
    aln_cl += [ref_file, fastq_file]
    with open(out_file, "w") as out_handle:
        logger.info(" ".join(aln_cl))
        subprocess.check_call(aln_cl, stdout=out_handle)

