# species for which pipeline is run
SPEC = "dre"

# reference spikein sequences
SPIKES_REF = "/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/RBAB/spikes/spikes40"

MIRNA_REF = "/zfs/datastore0/group_root/MAD-RBAB/05_Reference/external/dre/RNA/miRNA/hairpin21.fa"

# experiment directory in which while analysis is conducted
ANALYSIS_HOME = "./"

# directory with raw reads in fastq format
FQ_DIR = "./fq"
SAMPLES =  [s[:-6] for s in os.listdir(FQ_DIR) if s.endswith(".fastq")]

# trimmed reads dir
TRIM_DIR = "./trim"

# alignments dir
BAM_DIR = "./bam"

# count tables directory
COUNTS_DIR = "./counts"

# length of the reads used for alignment to spikes
SPIKE_TRIMMED_LEN = 40 

# length of the reads used for alignment to miRNA
MIRNA_TRIMMED_LEN = 40
MIRNA_MIN_LEN = '12'

# bowtie settings
BOWTIE_PARAMS_LIST = [
    "-L 6",                      # seed length
    "-i S,0,0.5",                 # interval between extracted seeds
    "--ignore-quals",            # treat as all qualities would be max possible
    "--norc",                    # do not align to reverse strand
    "--score-min L,-1,-0.6",   # linear function y=-1-0.6*x
    "-D 20",                     # consecutive seed extension attempts
    "-t",                        # print clock time
    "-p 16"                      # number of threads
]
BOWTIE_PARAMS = " ".join(BOWTIE_PARAMS_LIST)

rule all:
    input: "./spikes/counts/CountTable_spike.txt",
           ("./miRNA/filtered/{sample}_filtered.fastq".format(sample=s) for s in SAMPLES)

rule filter_reads_miRNA:
    input: "./miRNA/trim/{sample}_trimmed.fastq"
    output: "./miRNA/filtered/{sample}_filtered.fastq"
    params: min_len=MIRNA_MIN_LEN
    message: "Filtering out reads shorter than %s" % MIRNA_MIN_LEN
    shell:
        """
        mkdir -p ./miRNA/filtered
        python filter_spike_aln.py filter_short_reads --input {input} --output {output} --min-len {params.min_len}
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

# spike alignment and counting
rule spike_count:
    input: ("./spikes/bam/{sample}_spike_aln_sorted.bam.bai".format(sample=s) for s in SAMPLES)
    params: basename="_spike_aln_sorted.bam" ,
            bam_dir="./spikes/bam",
            count_dir="./spikes/counts"
    output: "./spikes/counts/CountTable_spike.txt"
    message: "Filtering and counting spike reads"
    shell: "python filter_spike_aln.py count_spikes --basename {params.basename} --bam-dir {params.bam_dir} --count-dir {params.count_dir}" 

rule spike_aln_index:
    input: "./spikes/bam/{sample}_spike_aln_sorted.bam"
    output: "./spikes/bam/{sample}_spike_aln_sorted.bam.bai"
    message: "Aligning spike alignments"
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
        bowtie2 {BOWTIE_PARAMS} -x {SPIKES_REF} -U {input} -S {output}
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
