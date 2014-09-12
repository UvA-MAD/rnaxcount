# species for which pipeline is run
SPEC = "dre"

# reference spikein sequences
SPIKES_REF = "/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/RBAB/spikes/spikes40"

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
    input: "bam/{sample}_spike_aln_sorted.bam.bai".format(sample=s) for s in SAMPLES

rule spike_aln_index:
    input: "bam/{sample}_spike_aln_sorted.bam"
    output: "bam/{sample}_spike_aln_sorted.bam.bai"
    message: "Aligning spike alignments"
    shell: "samtools index {input}"

rule spike_aln_sort:
    input: "bam/{sample}_spike_aln.bam"
    output: "bam/{sample}_spike_aln_sorted.bam"
    params: base="bam/{sample}_spike_aln_sorted"
    message: "Sorting spike alignments"
    shell: "samtools sort {input} {params.base}"

rule spike_sam2bam:
    input: "bam/{sample}_spike_aln.sam"
    output: "bam/{sample}_spike_aln.bam"
    message: "Converting sam to bam"
    shell: "samtools view -Sb {input} > {output}"
    
rule aln_spikes:
    """
    Align trimmed sequences to trimmed spikes.
    """
    input: "trim/{sample}_trimmed.fastq"
    output: "bam/{sample}_spike_aln.sam"
    message: "Aligning reads to spike sequences."
    shell: 
        """
        bowtie2 {BOWTIE_PARAMS} -x {SPIKES_REF} -U {input} -S {output}
        """

rule trim_reads:
    input: "fq/{sample}.fastq"
    output: "trim/{sample}_trimmed.fastq"
    message: "Trimming sequencing reads."
    shell:
        """
        fastx_trimmer -l {SPIKE_TRIMMED_LEN} -i {input} -o {output}
        """
