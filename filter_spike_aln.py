import click
import pysam
import os
import pandas as pd
from collections import defaultdict


def filter_and_count(bamfile):
    """
    Filters and counts reads in one bam file.
    Returns dictionary with the counts
    key -- spike
    value -- count of reads aligned to that spike
    """
    samfile = pysam.Samfile(bamfile, "rb")
    # get names of size spikein that reads have been align to
    size_spike_refs = [r for r in samfile.references if r.startswith("SSPK_")]

    # make a dictionary to lookup min lenght of aligned reads
    min_read_size = {}
    for spike in size_spike_refs:
        size = int(spike.split("_")[2])
        if size <= 40:
            min_size = round(0.8 * size)
        else:
            min_size = 32
        min_read_size[spike] = min_size

    # iterate over references and count reads that have min_read_size
    spike_counts = defaultdict(int)
    for spike in size_spike_refs:
        alignments = samfile.fetch(spike)
        for aln in alignments:
            if aln.rlen >= min_read_size[spike]:
                spike_counts[spike] += 1
    return(spike_counts)


# create a click group to use subcommands
@click.group()
def cli():
    pass


# subcommand used for filtering and counting spikes
@cli.command()
@click.option('--basename', help='basename of the bam files with alignments')
@click.option('--bam-dir',
              type=click.Path(exists=True, resolve_path=True),
              help='directory with bam files')
@click.option('--count-dir',
              type=click.Path(),
              help='directory to save count table in')
def count_spikes(basename, bam_dir, count_dir):
    """
    Applies filter rules and counts reads witch adhere to ther rules.
    """
    # get the list of files that match the basname
    bams = [f for f in os.listdir(bam_dir) if basename in f]
    bams = [f for f in bams if f.endswith(".bam")]
    # extract sample names from bam file names
    samples = ["".join(b.split(basename)).strip() for b in bams]

    count_dicts = []
    for b in bams:
        bamfile = os.path.join(bam_dir, b)
        counts = filter_and_count(bamfile)
        count_dicts.append(counts)

    # put count dicts into a pandas datatable
    count_table = pd.DataFrame(count_dicts, index=samples)
    # that gives table with spikes as columns and we need samples coulmns
    # and spikes rows, therefore we transpoze the table
    count_table = count_table.T

    # check if output dir is available, if not create one
    if not os.path.exists(count_dir):
        os.makedirs(count_dir)

    # save count table as tab delimited file
    count_file = os.path.join(count_dir, "CountTable_spike.txt")
    with open(count_file, 'w') as fh:
        count_table.to_csv(fh, sep="\t")

if __name__ == '__main__':
    cli()