import click
import pysam
import os
import csv
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from BCBio import GFF
from itertools import groupby
from functools import reduce
from collections import Counter


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


# subcommand used for counting number of reads
@cli.command()
@click.option('--fq-dir',
              type=click.Path(exists=True, resolve_path=True),
              help='directory with fastq files')
@click.option('--count-dir',
              type=click.Path(),
              help='directory to save count table in')
def count_fq_reads(fq_dir, count_dir):
    """
    Count reads in fastq files and save csv file with counts per sample.
    """
    fastqs = [f for f in os.listdir(fq_dir) if f.endswith(".fastq")]
    with open(os.path.join(count_dir, "total_reads.csv"), 'w') as fout:
        csv_writer = csv.writer(fout, delimiter="\t")
        for f in fastqs:
            fq = SeqIO.parse(os.path.join(fq_dir, f), "fastq")
            sample_name = f.split(".fastq")[0]
            count = sum([1 for i in fq])
            csv_writer.writerow([sample_name, count])


@cli.command()
@click.option('--input',
              help="fastq file to tream",
              type=click.File('rt'))
@click.option('--output',
              help='filtered reads file name',
              type=click.File('wt'))
@click.option('--min-len',
              type=int,
              help='min lenght of the read')
def filter_short_reads(input, output, min_len):
    """
    Removes reads shorter than min-len from fastq file
    """
    reads = [r for r in SeqIO.parse(input, 'fastq') if len(r) >= min_len]
    SeqIO.write(reads, output, format='fastq')


@cli.command()
@click.option('--dat',
              help="mirna annotation in embl format",
              type=click.File('rt'))
@click.option('--org',
              help="organism (ie. dre, dme, cel)",
              type=str)
@click.option('--gff',
              help="output gff file",
              type=click.File('wt'))
def embl2gff(dat, org, gff):
    """
    Parse embl file and estract mature miRNA location information.
    """
    # extract records
    dat_parser = SeqIO.parse(dat, "embl")
    # extract organism specific miRNAs
    org_mirnas = [mirna for mirna in dat_parser if mirna.name.startswith(org)]
    for mirna in org_mirnas:
        mirna.id = mirna.name
    GFF.write(org_mirnas, gff)


@cli.command()
@click.option('--dat',
              help="mirna annotation in embl format",
              type=click.File('rt'))
@click.option('--bamfile',
              help="sorted and indexed alignments in bam format",
              type=click.STRING)
@click.option('--out',
              help="output counts file",
              type=click.File('wt'))
def count_miRNAs(dat, bamfile, out):
    """
    Count mirRNAs from alignments to pre-miRNAs.

    Read alignment manipulation is done with pysam.
    Documentation of pysam can be found:
    http://pysam.readthedocs.org/en/latest/api.html
    """

    def make_group_aln_func(hairpin):
        """
        Create a function that will group alignments.

        Alignments are asigned into two groups:
        'mature': those overlaping with mature miRNA
        'pre': those mapped two hairpin outside mature annotations
        """
        mature_mirnas = hairpin.features
        mature_ranges = [(int(mirna.location.start), int(mirna.location.end))
                         for mirna in mature_mirnas]

        def group_aln(aln):
            for mr in mature_ranges:
                is_in = aln.pos >= mr[0] and aln.pos < mr[1]
                extend_3 = (aln.pos + aln.inferred_length) < (mr[1] + 3)
                if (is_in and extend_3):
                    # mature miRNA
                    if is_quality_mature(aln):
                        return('mature')
                    else:
                        return('low_qual')
                else:
                    # pre-miRNA
                    if is_quality_pre(aln):
                        return('pre')
                    else:
                        return('low_qual')
        return(group_aln)

    def is_quality_mature(aln):
        """
        Discard reads aligning to mature miRNAs
        according to specific rules.
        """
        # allow 3' softcliping but not 5'
        # sofcliping is a tupple (4,...) in a cigar tupple list
        # check if it is as a first typple of the alignment cigar
        cig = aln.cigar
        if cig[0][0] == 4:
            return(False)

        # drop alignments with more than 3 nucleotides softclipped on 3'
        if (cig[-1][0] == 4 and cig[-1][1] > 3):
            return(False)
        # allow for 10%  of mismatches + indels
        # bowtie does not log mismatches in cigar string
        # but uses tag 'NM'
        n_mismatch = [tag[1] for tag in aln.tags if tag[0] == 'NM'][0]
        # insertion has code 1 and deletion 2 in cigar string
        n_indel = sum([c[1] for c in cig if c[0] in (1, 2)])
        if (n_mismatch + n_indel > 0.1*aln.rlen):
            return(False)
        return(True)

    def is_quality_pre(aln):
        """
        Discard reads aligned to hairpin
        according to specific rules.
        """
        cig = aln.cigar
        # allow for 10%  of mismatches + indels + sofclip
        # bowtie does not log mismatches in cigar string
        # but uses tag 'NM'
        n_mismatch = [tag[1] for tag in aln.tags if tag[0] == 'NM'][0]
        # insertion has code 1 and deletion 2 in cigar string
        n_indel = sum([c[1] for c in cig if c[0] in (1, 2)])
        n_softclip = sum([c[1] for c in cig if c[0] == 4])
        if (n_mismatch + n_indel + n_softclip > 0.1*aln.rlen):
            return(False)
        return(True)

    def make_count_alns_func(hairpin):
        def count_alns(countlist, group_aln):
            mature_mirnas = [feature for feature in hairpin.features]
            group, g_aln = group_aln
            for aln in g_aln:
                start = aln.pos
                end = aln.pos + aln.inferred_length
                length = aln.inferred_length
                if group == 'pre':
                    seqid = "%s_%s_%d_%d_%d" % (hairpin.name,
                                                "pre", start, end,
                                                length)
                    countlist += [seqid]
                elif group == 'mature':
                    # check to wich mature miRNA the read aligns
                    mirna = [mirna for mirna in mature_mirnas if (start >= int(mirna.location.start) and start < int(mirna.location.end))][0]
                    name = mirna.qualifiers['product'][0]
                    seqid = "%s_%s_%d_%d_%d" % (hairpin.name,
                                                name, start, end,
                                                length)
                    countlist += [seqid]
                else:
                    continue
            return(countlist)
        return(count_alns)

    # read in miRNA annotations from embl formated dat file
    annots = SeqIO.parse(dat, 'embl')
    alignments = pysam.Samfile(bamfile, 'rb')
    references = alignments.references

    # fetch annotation for reference sequences
    # prsent in the bam file
    # for some reason name is valid seqid fo the sequence
    hairpin_annots = [a for a in annots if a.name in references]

    # loop over all SeqRecords in annotations
    counts_tupples = []
    for hairpin in hairpin_annots:
        # mature miRNAs are the features of hairpin SeqRecords
        # create grouping function which knows miRNA annotations
        group_aln = make_group_aln_func(hairpin)
        # get reads aligning to this hairpin
        hairpin_alns = alignments.fetch(hairpin.name)
        # group alignments into those mapped to mature miRNA
        # and those mapping autside (pre-miRNA)
        grouped_alns = groupby(hairpin_alns, group_aln)
        # count
        count_alns = make_count_alns_func(hairpin)
        counts_list = reduce(count_alns, grouped_alns, [])
        # cont same occurences of seqid
        counts_dict = Counter(counts_list)
        # store them as list of tupples, because it's difficult
        # to join dictionaries
        counts_tupples += [i for i in counts_dict.items()]

    # create dir for counts if needed
    counts_dir = os.path.dirname(out.name)
    if not os.path.exists(counts_dir):
        os.makedirs(counts_dir)

    # export counts as csv file
    csv_writer = csv.writer(out)
    for row in counts_tupples:
        csv_writer.writerow(row)


@cli.command()
@click.option('--dir',
              help="count tables directory",
              type=click.STRING)
@click.option('--suffix',
              help="common suffix of count files, used to extract smapleids",
              type=click.STRING)
@click.option('--out',
              help="output counts file",
              type=click.File('wt'))
def merge_count_tables(dir, suffix, out):
    """
    Merge count tables on seqid's into single table.
    """
    # list files in dir
    count_files = [f for f in os.listdir(dir) if f.endswith(".csv")]
    samples = [f.split(suffix)[0] for f in count_files]
    dfs = [pd.read_csv(os.path.join(dir, f), index_col=0) for f in count_files]
    merged_dfs = reduce(lambda df1, df2: pd.merge(df1, df2,
                                                  left_index=True,
                                                  right_index=True,
                                                  how="outer"), dfs)
    merged_dfs = merged_dfs.fillna(0)
    merged_dfs.columns = samples
    merged_dfs.to_csv(out, sep="\t")


if __name__ == '__main__':
    cli()
