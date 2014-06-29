"""Task definitions for the Celery message queue (http://celeryproject.org/).
"""
import time

from celery.task import task

from bcbio.pipeline import sample, lane, toplevel, storage, shared, variation
from bcbio.variation import realign, genotype, ensemble, recalibrate, multi

# Global configuration for tasks in the main celeryconfig module
import celeryconfig

@task(ignore_results=True, queue="toplevel")
def analyze_and_upload(*args):
    """Run full analysis and upload results to Galaxy instance.

    Workers need to run on the machine with Galaxy installed for upload,
    but the actual processing can be distributed to multiple nodes.
    """
    config_file = celeryconfig.BCBIO_CONFIG_FILE
    remote_info = args[0]
    toplevel.analyze_and_upload(remote_info, config_file)

@task(ignore_results=True, queue="storage")
def long_term_storage(*args):
    config_file = celeryconfig.BCBIO_CONFIG_FILE
    remote_info = args[0]
    storage.long_term_storage(remote_info, config_file)

@task
def process_lane(*args):
    return lane.process_lane(*args)

@task
def process_alignment(*args):
    return lane.process_alignment(*args)

@task
def merge_sample(*args):
    return sample.merge_sample(*args)

@task
def prep_recal(*args):
    return recalibrate.prep_recal(*args)

@task
def write_recal_bam(*args):
    return recalibrate.write_recal_bam(*args)

@task
def realign_sample(*args):
    return realign.realign_sample(*args)

@task
def process_sample(*args):
    return sample.process_sample(*args)

@task
def split_variants_by_sample(*args):
    return multi.split_variants_by_sample(*args)

@task
def postprocess_variants(*args):
    return sample.postprocess_variants(*args)

@task
def generate_bigwig(*args):
    return sample.generate_bigwig(*args)

@task
def combine_bam(*args):
    return shared.combine_bam(*args)

@task
def variantcall_sample(*args):
    return genotype.variantcall_sample(*args)

@task
def combine_variant_files(*args):
    return genotype.combine_variant_files(*args)

@task
def detect_sv(*args):
    return variation.detect_sv(*args)

@task
def combine_calls(*args):
    return ensemble.combine_calls(*args)

@task
def test(x):
    print x
    time.sleep(5)
    return x
