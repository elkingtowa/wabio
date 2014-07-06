#!/usr/bin/env python
"""Merge bed files from MACS and NPS, making a unified set of ChIP predictions.

This first overlaps replicate nucleosome predictions for each set of marks, and
then combines these predictions into one final set of overlap regions that are
used as the reference set of nucleosomes in the experiment.

This reference set is then overlapped against each of the MACS gene calls,
producing a set of nucleosome specific MACS calls. These are output as a BED
file for further downstream analysis.

http://liulab.dfci.harvard.edu/MACS/
http://liulab.dfci.harvard.edu/NPS/

Usage:
    merge_nps_macs.py <config_file> <output_dir>

config_file is a YAML configuration file with entries for each experiment:

- name: exp_name
  macs:
  - file1
  - file2
  nps:
  - file1
  - file2
"""
import sys
import os
import csv
import collections
import re
import string

import yaml
from bx.intervals.intersection import IntervalTree
from bx.intervals.cluster import ClusterTree

PERCENT_OL = 0.75

def main(config_file, output_dir):
    with open(config_file) as config_handle:
        config = yaml.load(config_handle)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    base_nucleosomes = list(combine_regions(
            (intersect_bed(exp['nps'], PERCENT_OL) for exp in config)))
    for exp in config:
        out_file = os.path.join(output_dir, "%s.bed" % _sanitize(exp['name']))
        final_regions = combine_regions(
                (intersect_regions_and_bed(base_nucleosomes, macs, PERCENT_OL)
                    for macs in exp['macs']),
                len(exp['macs']))
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            [writer.writerow(p) for p in final_regions]

def combine_regions(all_regions, required_regions=1):
    """Generate the combination of a set of chrom, start, end regions.

    If required_regions is 1 then this is a union combination. Otherwise
    it is an intersection.
    """
    clusters = collections.defaultdict(lambda: ClusterTree(0, required_regions))
    i = 0
    for region_gen in all_regions:
        for chrom, start, end in region_gen:
            clusters[chrom].insert(start, end, i)
            i += 1
    for chrom, cluster in clusters.iteritems():
        for (s, e, _) in cluster.getregions():
            yield chrom, s, e

def intersect_regions_and_bed(regions, bed_file, percent_ol):
    """Intersect a set of regions with a BED file.
    """
    do_interesect = _intersect_by_percent(percent_ol)
    return do_interesect(regions, bed_reader(bed_file))

def intersect_bed(bed_files, percent_ol):
    """Intersect a set of replicate bed files.
    """
    return reduce(_intersect_by_percent(percent_ol),
                  (bed_reader(f) for f in bed_files))

def _intersect_by_percent(percent_ol):
    """Intersect two sets of intervals, requiring a percentage overlap.
    """
    def _intersect_two(regions1, regions2):
        """Intersect two regions, generating chromosome, start, end intervals.
        """
        int_trees = collections.defaultdict(lambda: IntervalTree())
        for chrom, start, end in regions2:
            info = dict(start=start, end=end)
            int_trees[chrom].insert(start, end, info)
        int_trees = dict(int_trees)
        for chrom, start, end in regions1:
            try:
                overlaps = int_trees[chrom].find(start, end)
            except KeyError:
                overlaps = []
            for ol in overlaps:
                cur_ol = (float(len(set(range(start, end)) &
                                    set(range(ol['start'], ol['end'])))) /
                          float(end - start))
                if cur_ol >= percent_ol:
                    yield chrom, start, end
                    break
    return _intersect_two

def bed_reader(bed_file):
    """Generator producing chrom, start, end coordinates for a BED file.
    """
    with open(bed_file) as in_handle:
        reader = csv.reader(in_handle, dialect="excel-tab")
        for parts in reader:
            try:
                yield parts[0], int(parts[1]), int(parts[2])
            except ValueError:
                # ignore header lines without integer starts and ends
                pass

to_remove = re.compile("[%s%s]" % (string.whitespace, string.punctuation))
def _sanitize(name):
    return to_remove.sub("_", name)

if __name__ == "__main__":
    main(*sys.argv[1:])
