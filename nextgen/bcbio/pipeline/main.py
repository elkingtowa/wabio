"""Main entry point for distributed next-gen sequencing pipelines.

Handles running the full pipeline based on instructions
"""
import os
import sys
import math
import argparse

from bcbio.solexa.flowcell import get_fastq_dir
from bcbio import utils
from bcbio.log import setup_logging
from bcbio.distributed.messaging import parallel_runner
from bcbio.pipeline.run_info import get_run_info
from bcbio.pipeline.demultiplex import add_multiplex_across_lanes
from bcbio.pipeline.merge import organize_samples
from bcbio.pipeline.qcsummary import write_metrics, write_project_summary
from bcbio.variation.realign import parallel_realign_sample
from bcbio.variation.genotype import parallel_variantcall, combine_multiple_callers
from bcbio.variation import ensemble, recalibrate

def run_main(config, config_file, work_dir, parallel,
             fc_dir=None, run_info_yaml=None):
    """Run toplevel analysis, processing a set of input files.

    config_file -- Main YAML configuration file with system parameters
    fc_dir -- Directory of fastq files to process
    run_info_yaml -- YAML configuration file specifying inputs to process
    """
    setup_logging(config)
    fc_name, fc_date, run_info = get_run_info(fc_dir, config, run_info_yaml)
    fastq_dir, galaxy_dir, config_dir = _get_full_paths(get_fastq_dir(fc_dir) if fc_dir else None,
                                                        config, config_file)
    config_file = os.path.join(config_dir, os.path.basename(config_file))
    dirs = {"fastq": fastq_dir, "galaxy": galaxy_dir,
            "work": work_dir, "flowcell": fc_dir, "config": config_dir}
    config = _set_resources(parallel, config)
    run_parallel = parallel_runner(parallel, dirs, config, config_file)

    # process each flowcell lane
    run_items = add_multiplex_across_lanes(run_info["details"], dirs["fastq"], fc_name)
    lanes = ((info, fc_name, fc_date, dirs, config) for info in run_items)
    lane_items = run_parallel("process_lane", lanes)
    align_items = run_parallel("process_alignment", lane_items)
    # process samples, potentially multiplexed across multiple lanes
    samples = organize_samples(align_items, dirs, config_file)
    samples = run_parallel("merge_sample", samples)
    samples = run_parallel("prep_recal", samples)
    samples = recalibrate.parallel_write_recal_bam(samples, run_parallel)
    samples = parallel_realign_sample(samples, run_parallel)
    samples = parallel_variantcall(samples, run_parallel)
    samples = run_parallel("postprocess_variants", samples)
    samples = combine_multiple_callers(samples)
    samples = run_parallel("detect_sv", samples)
    samples = run_parallel("combine_calls", samples)
    run_parallel("process_sample", samples)
    run_parallel("generate_bigwig", samples, {"programs": ["ucsc_bigwig"]})
    write_project_summary(samples)
    write_metrics(run_info, fc_name, fc_date, dirs)

def _set_resources(parallel, config):
    """Set resource availability for programs, downsizing to local runs.
    """
    for program in ["gatk", "novoalign"]:
        if not config["resources"].has_key(program):
            config["resources"][program] = {}
        if parallel["type"] == "local":
            import multiprocessing
            cores = min(parallel["cores"], multiprocessing.cpu_count())
            config["resources"][program]["cores"] = cores
    return config

# ## Utility functions

def parse_cl_args(in_args):
    """Parse input commandline arguments, handling multiple cases.

    Returns the main config file and set of kwargs.
    """
    parser = argparse.ArgumentParser(
        description="Best-practice pipelines for fully automated high throughput sequencing analysis")
    parser.add_argument("inputs", nargs="*")
    parser.add_argument("-n", "--numcores", type=int, default=0)
    parser.add_argument("-t", "--paralleltype")
    parser.add_argument("-s", "--scheduler")
    parser.add_argument("-q", "--queue")
    parser.add_argument("-u", "--upgrade", help="Perform an upgrade of bcbio_nextgen in place.",
                        choices = ["stable", "development", "system"])

    args = parser.parse_args(in_args)
    config_file = args.inputs[0] if len(args.inputs) > 0 else None
    kwargs = {"numcores": args.numcores if args.numcores > 0 else None,
              "paralleltype": args.paralleltype,
              "scheduler": args.scheduler,
              "queue": args.queue,
              "upgrade": args.upgrade}
    if len(args.inputs) == 3:
        kwargs["fc_dir"] = args.inputs[1]
        kwargs["run_info_yaml"] = args.inputs[2]
    elif len(args.inputs) == 2:
        extra = args.inputs[1]
        if os.path.isfile(extra):
            kwargs["run_info_yaml"] = extra
        else:
            kwargs["fc_dir"] = extra
    elif args.upgrade is None:
        parser.print_help()
        sys.exit()
    return config_file, kwargs

def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    if fastq_dir:
        fastq_dir = utils.add_full_path(fastq_dir)
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config.get("galaxy_config", "universe_wsgi.ini"),
                                             config_dir)
    return fastq_dir, os.path.dirname(galaxy_config_file), config_dir
