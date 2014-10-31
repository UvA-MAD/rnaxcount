import os
from snakemake.utils import R
# species for which pipeline is run
SPEC = "dre"

# reference spikein sequences
SPIKES_REF = "/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/RBAB/spikes/spikes40"

MIRNA_REF = "/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/dre/RNA/miRNA/hairpin21"

MIRNA_ANNOTATIONS = "/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/RNA/mirBase21/miRNA.dat"

PIRNA_REF = "/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/dre/RNA/piRNA/piRNA"

# experiment directory in which while analysis is conducted
ANALYSIS_HOME = "./"

# directory with raw reads in fastq format
FQ_DIR = "./fq"
SAMPLES =  [s[:-6] for s in os.listdir(FQ_DIR) if s.endswith(".fastq")]

# trimmed reads dir
TRIM_DIR = "./trim"

# alignments di
BAM_DIR = "./bam"

# spike count table directory
SPIKE_COUNTS_DIR = "./spikes/counts"

# length of the reads used for alignment to spikes
SPIKE_TRIMMED_LEN = 40 

# length of the reads used for alignment to miRNA
MIRNA_TRIMMED_LEN = 40
MIRNA_MIN_LEN = '12'

# bowtie settings
# alignment score differs between 'local' and 'end-to-end' alignment
# this is why -score-min differs for 10% mm
SPIKE_BOWTIE_PARAMS_LIST = [
    "-L 6",                      # seed length
    "-i S,0,0.5",                # interval between extracted seeds
    "--ignore-quals",            # treat as all qualities would be max possible
    "--norc",                    # do not align to reverse strand
    "--score-min L,-1,-0.6",     # -1-0.6*read_length -- 10% mismatches allowed
    "-D 20",                     # consecutive seed extension attempts
    "-t",                        # print clock time
    "-p 16"                      # number of threads
]
SPIKE_BOWTIE_PARAMS = " ".join(SPIKE_BOWTIE_PARAMS_LIST)

MIRNA_BOWTIE_PARAMS_LIST = [
    "-L 6"                       # seed length
    "--ignore-quals",            # treat as all qualities would be max possible
    "--norc",                    # do not align to reverse strand
    "--local",                   # allow for softcliping
    "--score-min L,0,1.2",       # 1.2*read_length -- 10% mismatches allowed
    "-p 6"                       # number of threads
]
MIRNA_BOWTIE_PARAMS = " ". join(MIRNA_BOWTIE_PARAMS_LIST)

PIRNA_BOWTIE_PARAMS_LIST = [
    "-L 6",                      # seed length
    "--ignore-quals",            # treat as all qualities would be max possible
    "--norc",                    # do not align to reverse strand
    "--score-min L,-1,-0.6",     # -1-0.6*read_length -- 10% mismatches allowed
    "-D 20",                     # consecutive seed extension attempts
    "-t",                        # print clock time
    "-p 16"                      # number of threads
]
PIRNA_BOWTIE_PARAMS = " ".join(PIRNA_BOWTIE_PARAMS_LIST)

rule all:
    input: "./spikes/counts/CountTable_spike.txt",
           "./miRNA/counts/CountTable_mirna.txt",
           "./piRNA/counts/CountTable_pirna.txt",
           "./spikes/counts/norm_count.png"

# count piRNAs
rule count_piRNA:
    input:("./piRNA/bam/{sample}_pirna_aln_sorted.bam".format(sample=s) for s in SAMPLES)
    output: "./piRNA/counts/CountTable_pirna.txt"
    params: bam_dir="./piRNA/bam/",
            count_dir="./piRNA/counts/"
    shell: "python sRNA_tools.py count_pirna --bam-dir {params.bam_dir} --count-dir {params.count_dir}"     

#
# sam to bam, sort and index
rule sam2bam_sort_index:
    input: "./piRNA/bam/{sample}_pirna_aln.sam"
    output: "./piRNA/bam/{sample}_pirna_aln.bam",
            "./piRNA/bam/{sample}_pirna_aln_sorted.bam",
            "./piRNA/bam/{sample}_pirna_aln_sorted.bam.bai"
    params: bam="./piRNA/bam/{sample}_pirna_aln.bam",
            sorted_base="./piRNA/bam/{sample}_pirna_aln_sorted",
            sorted_bam="./piRNA/bam/{sample}_pirna_aln_sorted.bam"
    message: "Converting piRNA alignments to bam, sorting and indexing."
    shell:
        """
        samtools view -bS {input} > {params.bam}
        samtools sort {params.bam} {params.sorted_base}
        samtools index {params.sorted_bam}
        """

# piRNA
rule aln_piRNA:
    input: fq="fq/{sample}.fastq",
           dir="./piRNA/bam/"
    output: "./piRNA/bam/{sample}_pirna_aln.sam"
    message: "Aligning reads to piRNA sequences."
    shell:
        """
        bowtie2 {PIRNA_BOWTIE_PARAMS} -x {PIRNA_REF} -U {input.fq} -S {output}
        """

rule create_piRNA_dir:
    output: "./piRNA/bam/",
            "./piRNA/counts/"
    shell: "mkdir -p ./piRNA/bam ./piRNA/counts"

# mirRNA
rule merge_miRNA_counts:
    input: ("./miRNA/counts/{sample}_mirna.csv".format(sample=s) for s in SAMPLES)
    output: "./miRNA/counts/CountTable_mirna.txt"
    shell: "python sRNA_tools.py merge_count_tables --dir ./miRNA/counts --suffix _mirna.csv --out {output}"     

rule miRNA_count:
    input: "./miRNA/bam/{sample}_mirna_aln_sorted.bam"
    output: "./miRNA/counts/{sample}_mirna.csv"
    message: "Counting miRNAs"
    shell: "python sRNA_tools.py count_mirnas --dat {MIRNA_ANNOTATIONS} --bamfile {input} --out {output}"

rule miRNA_aln_index:
    input: "./miRNA/bam/{sample}_mirna_aln_sorted.bam"
    output: "./miRNA/bam/{sample}_mirna_aln_sorted.bam.bai"
    message: "Indexing miRNA alignments"
    shell: "samtools index {input}"

rule miRNA_aln_sort:
    input: "./miRNA/bam/{sample}_mirna_aln.bam"
    output: "./miRNA/bam/{sample}_mirna_aln_sorted.bam"
    params: base="miRNA/bam/{sample}_mirna_aln_sorted"
    message: "Sorting miRNA alignments"
    shell: "samtools sort {input} {params.base}"

rule miRNA_sam2bam:
    input: "./miRNA/bam/{sample}_mirna_aln.sam"
    output: "./miRNA/bam/{sample}_mirna_aln.bam"
    message: "Converting miRNA alignments from sam to bam"
    shell: "samtools view -Sb {input} > {output}"

rule aln_miRNA:
    """
    Align trimmed sequences to miRNA sequences.
    """
    input: "./miRNA/trim/{sample}_trimmed.fastq"
    output: "./miRNA/bam/{sample}_mirna_aln.sam"
    message: "Aligning reads to miRNA sequences."
    shell: 
        """
        bowtie2 {MIRNA_BOWTIE_PARAMS} -x {MIRNA_REF} -U {input} -S {output}
        """

rule filter_reads_miRNA:
    input: "./miRNA/trim/{sample}_trimmed.fastq"
    output: "./miRNA/filtered/{sample}_filtered.fastq"
    params: min_len=MIRNA_MIN_LEN
    message: "Filtering out reads shorter than %s" % MIRNA_MIN_LEN
    shell:
        """
        mkdir -p ./miRNA/filtered
        python sRNA_tools.py filter_short_reads --input {input} --output {output} --min-len {params.min_len}
        """

rule trim_reads_miRNA:
    input: fqs="fq/{sample}.fastq", dir="./miRNA"
    output: "./miRNA/trim/{sample}_trimmed.fastq"
    message: "Trimming sequencing reads fro miRNA alignment."
    shell:
        """
        fastx_trimmer -l {MIRNA_TRIMMED_LEN} -i {input.fqs} -o {output}
        """

# miRNA alignment and counting
rule create_miRNA_dir:
    output: "./miRNA"
    shell: "mkdir -p ./miRNA/trim ./miRNA/bam ./miRNA/counts"

# R command in snakemake don't take input and output
# putting them into constants
SPIKE_COUNT = "./spikes/counts/CountTable_spike.txt"
TOTAL_COUNT = os.path.join(SPIKE_COUNTS_DIR, "total_reads.csv")
COUNT_PNG = "./spikes/counts/count.png",
NORM_COUNT_PNG = "./spikes/counts/norm_count.png"
# make spike count plots
rule spike_count_plots:
    input: totalcount=os.path.join(SPIKE_COUNTS_DIR, "total_reads.csv"),
           spikecount="./spikes/counts/CountTable_spike.txt"
    output: countpng="./spikes/counts/count.png",
            normcountpng="./spikes/counts/norm_count.png"
    run: R("""
           library(faradr);
           png(filename="{COUNT_PNG}");
           plot(PlotSpikeCounts("{SPIKE_COUNT}"));
           dev.off();
           png(filename="{NORM_COUNT_PNG}");
           plot(PlotNormalSpikeCounts("{SPIKE_COUNT}", "{TOTAL_COUNT}"));
           dev.off();
           """)
             
# calculate amount of reads per sample for spike normalise plots
rule total_reads_count:
   input: "./spikes/counts/CountTable_spike.txt"
   params: fqdir=FQ_DIR,
           spikedir = SPIKE_COUNTS_DIR
   output: os.path.join(SPIKE_COUNTS_DIR, "total_reads.csv") 
   message: "Count nubmer of reads in each sample"
   shell: "python sRNA_tools.py count_fq_reads --fq-dir {params.fqdir} --count-dir {params.spikedir}"

# spike alignment and counting
rule spike_count:
    input: ("./spikes/bam/{sample}_spike_aln_sorted.bam.bai".format(sample=s) for s in SAMPLES)
    params: basename="_spike_aln_sorted.bam" ,
            bam_dir="./spikes/bam",
            count_dir="./spikes/counts"
    output: "./spikes/counts/CountTable_spike.txt"
    message: "Filtering and counting spike reads"
    shell: "python sRNA_tools.py count_spikes --basename {params.basename} --bam-dir {params.bam_dir} --count-dir {params.count_dir}" 

rule spike_aln_index:
    input: "./spikes/bam/{sample}_spike_aln_sorted.bam"
    output: "./spikes/bam/{sample}_spike_aln_sorted.bam.bai"
    message: "Indexing spike alignments"
    shell: "samtools index {input}"

rule spike_aln_sort:
    input: "./spikes/bam/{sample}_spike_aln.bam"
    output: "./spikes/bam/{sample}_spike_aln_sorted.bam"
    params: base="spikes/bam/{sample}_spike_aln_sorted"
    message: "Sorting spike alignments"
    shell: "samtools sort {input} {params.base}"

rule spike_sam2bam:
    input: "./spikes/bam/{sample}_spike_aln.sam"
    output: "./spikes/bam/{sample}_spike_aln.bam"
    message: "Converting sam to bam"
    shell: "samtools view -Sb {input} > {output}"
    
rule aln_spikes:
    """
    Align trimmed sequences to trimmed spikes.
    """
    input: "./spikes/trim/{sample}_trimmed.fastq"
    output: "./spikes/bam/{sample}_spike_aln.sam"
    message: "Aligning reads to spike sequences."
    shell: 
        """
        bowtie2 {SPIKE_BOWTIE_PARAMS} -x {SPIKES_REF} -U {input} -S {output}
        """

rule trim_reads_spike:
    input: fqs="fq/{sample}.fastq", dir="./spikes"
    output: "./spikes/trim/{sample}_trimmed.fastq"
    message: "Trimming sequencing reads for spike alignment."
    shell:
        """
        fastx_trimmer -l {SPIKE_TRIMMED_LEN} -i {input.fqs} -o {output}
        """

rule create_spike_dir:
    output: "./spikes"
    shell: "mkdir -p ./spikes/trim ./spikes/bam ./spikes/counts"
