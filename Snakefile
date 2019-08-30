shell.prefix("source ~/.bash_profile; set -euo pipefail;")
from util.varsub import varsub
configfile: "config.yaml"
varsub(config)

#Wildcard rule which takes all the fastq files ending with .fq.gz

SAMPLES, = glob_wildcards(config['datadirs']['fastq'] + "/" + "{file}_1.fq.gz")
READS = ["1", "2"]


rule all:
   input:
        config['reference']['rsemgenomedir']['hg38'] + ".n2g.idx.fa",
        config['reference']['stargenomedir']['hg38'] + "/" + "SAindex",       
        expand(config['datadirs']['qc'] + "/" + "{file}_{read}_fastqc.html", file = SAMPLES, read = READS),
        expand(config['datadirs']['trim'] + "/" + "{file}_{read}_val_{read}.fq.gz", file = SAMPLES, read= READS),
        expand(config['datadirs']['fatsqctrim'] + "/" + "{file}_{read}_val_{read}_fastqc.html", file = SAMPLES, read= READS),
        expand(config['datadirs']['bam'] + "/" + "{file}_SJ.out.tab", file = SAMPLES),
        config['datadirs']['sj_files'] + "/" + "SJ.out.pass1_merged.tab",
        expand(config['datadirs']['pass2'] + "/" + "{file}_Aligned.sortedByCoord.out.bam", file = SAMPLES ),
        expand(config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam", file =SAMPLES),
        expand(config['datadirs']['rnaseq_qc'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam.metrics.tsv", file =SAMPLES),
        expand(config['datadirs']['quant'] + "/" + "{file}.genes.results", file = SAMPLES)

rule fastqc:
   input:
      f1 = config['datadirs']['fastq'] + "/" + "{file}_{read}.fq.gz"
   output: config['datadirs']['qc'] + "/" + "{file}_{read}_fastqc.html", config['datadirs']['qc'] + "/" + "{file}_{read}_fastqc.zip"
   params:
      prefix =  config['datadirs']['qc'],
   resources:
      mem_mb= 10000
   shell:"""
        fastqc  --thread 8 --outdir {params.prefix} --nogroup {input.f1}
        """
          
          
rule trim_galore_pe:
   input:
      f1 = config['datadirs']['fastq'] + "/" + "{file}_1.fq.gz",
      f2 = config['datadirs']['fastq'] + "/" + "{file}_2.fq.gz"
   output:
      fwd_pai = config['datadirs']['trim'] + "/" + "{file}_1_val_1.fq.gz",
      rev_pai = config['datadirs']['trim'] + "/" + "{file}_2_val_2.fq.gz",
   params:
      extra = " -j 8 --illumina -q 20 --phred33 --length 20",  
      prefix =  config['datadirs']['trim'],
   resources:
      mem_mb= 20000
   shell:"""
        trim_galore \
        {params.extra} \
        --paired {input.f1} {input.f2} \
        -o {params.prefix} 
         """   
         
rule fastqctrim:
   input:
      f1 = config['datadirs']['trim'] + "/" + "{file}_{read}_val_{read}.fq.gz",
      trimmedFiles = rules.trim_galore_pe.output.rev_pai 
   output: config['datadirs']['fatsqctrim'] + "/" +"{file}_{read}_val_{read}_fastqc.html"
   params:
      preefix = config['datadirs']['fatsqctrim'],
   resources:
      mem_mb= 10000
   shell:"""
       fastqc  --thread 8 --outdir {params.prefix} --nogroup {input.f1}
       """         


rule pass1:
   input:
      f1 = config['datadirs']['trim'] + "/" + "{file}_1_val_1.fq.gz",
      f2 = config['datadirs']['trim'] + "/" + "{file}_2_val_2.fq.gz",
      queue = rules.trim_galore_pe.output.rev_pai
   output: config['datadirs']['bam'] + "/" + "{file}_SJ.out.tab", config['datadirs']['bam'] + "/" + "{file}_Aligned.toTranscriptome.out.bam"
   params:
      genomedir = config['reference']['star_ref'],
      prefix =  config['datadirs']['bam'] + "/" + "{file}_"
   threads: 16
   resources:
      mem_mb= 40000
   shell: """
        STAR  \
        --runThreadN {threads} \
        --genomeDir {params.genomedir} \
        --readFilesIn {input.f1} {input.f2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype None \
        --outSAMunmapped Within \
        --quantMode TranscriptomeSAM \
        --outSAMattributes NH HI AS NM MD \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --limitBAMsortRAM 50000000000
        """


rule rsem_genome:
    input:
        fasta = config['reference']['fasta']['hg38'],
        gtf = config['reference']['gtf']['hg38'],
        genomedir = config['reference']['rsemgenomedir']['hg38']
    output:
        rsemindex = config['reference']['rsemgenomedir']['hg38'] + ".n2g.idx.fa"
    params:
        prepref = config['tools']['rsem']['prepref']
    threads: 6
    shell:"""
        {params.prepref} \
        -p {threads} \
        --gtf {input.gtf} {input.fasta} {input.genomedir}
        """

rule merge:
    input:
      sjs =  expand(config['datadirs']['bam'] + "/" + "{file}_SJ.out.tab" , file = SAMPLES )
    output:
      sjs=  config['datadirs']['sj_files'] + "/" + "SJ.out.pass1_merged.tab"
    threads: 1
    shell: """
         cat {input.sjs} | awk '$7 >= 3' | cut -f1-4 | sort -u > {output.sjs}
          """


rule star_genome:
    input:
        fasta = config['reference']['fasta']['hg38'],
        gtf = config['reference']['gtf']['hg38'],
        sjs =  config['datadirs']['sj_files'] + "/" + "SJ.out.pass1_merged.tab",
        genomedir = config['reference']['stargenomedir']['hg38'],
        queue = rules.merge.output.sjs
    output:
        starindex = config['reference']['stargenomedir']['hg38'] + "/" + "SAindex"
    params:
        overhang = 99
    threads: 12
    resources:
        mem_mb = 40000
    shell: """
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {input.genomedir} \
        --outFileNamePrefix {input.genomedir} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbFileChrStartEnd  {input.sjs} \
        --sjdbOverhang {params.overhang}
          """

rule pass2:
   input:
      f1 = config['datadirs']['trim'] + "/" + "{file}_1_val_1.fq.gz",
      f2 = config['datadirs']['trim'] + "/" + "{file}_2_val_2.fq.gz",
      line = rules.star_genome.output.starindex
   output: config['datadirs']['pass2'] + "/" + "{file}_Aligned.toTranscriptome.out.bam", config['datadirs']['pass2'] + "/" + "{file}_Aligned.sortedByCoord.out.bam"
   params:
      genomedir = config['reference']['stargenomedir']['hg38'],
      prefix =  config['datadirs']['pass2'] + "/" + "{file}_"
   threads: 16
   resources:
      mem_mb= 50000
   shell: """
        STAR  \
        --runThreadN {threads} \
        --genomeDir {params.genomedir} \
        --readFilesIn {input.f1} {input.f2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.prefix} \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterType BySJout \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --alignSoftClipAtReferenceEnds Yes \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype  BAM SortedByCoordinate  \
        --outSAMunmapped Within \
        --genomeLoad NoSharedMemory \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType WithinBAM SoftClip \
        --chimMainSegmentMultNmax 1 \
        --outSAMattributes NH HI AS nM NM ch \
        --outSAMattrRGline ID:rg1 SM:sm1
        """

  rule mark_dups:
     input:
        bam = config['datadirs']['pass2'] + "/" + "{file}_Aligned.sortedByCoord.out.bam"
     output:
        dbam = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam",
        metric = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.metrics.txt"
     params:
        picard = "java -jar $EBROOTPICARD/picard.jar"
     resources:
        mem_mb = 10000
     shell: """
         module load picard
        {params.picard} MarkDuplicates  INPUT={input.bam} OUTPUT={output.dbam} METRICS_FILE={output.metric} ASSUME_SORT_ORDER=coordinate OPTICAL_DUPLICATE_PIXEL_DISTANCE=100
          """



rule RNA_SeqC:
    input:
        gtf = config['reference']['collapsed_gtf'], # dont forget to collaspe the transcript using python script else it wont be proccessed
        fasta = config['reference']['fasta']['hg38'],
        bam = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam"
    output: config['datadirs']['rnaseq_qc'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam.metrics.tsv",
    params:
        prefix = config['datadirs']['rnaseq_qc'],
    resources:
        mem_mb= 10000
    shell:"""
         source ~/.profile
         rnaseqc {input.gtf} {input.bam} {params.prefix} --verbose
        """


rule rsem_norm:
   input:
      bam = config['datadirs']['pass2'] + "/" + "{file}_Aligned.toTranscriptome.out.bam"
   output:
      genes = config['datadirs']['quant'] + "/" + "{file}.genes.results"
   params:
      calcexp = config['tools']['rsem']['calcexp'],
      genomedir = config['reference']['rsemgenomedir']['hg38'],
      prefix = config['datadirs']['quant'] + "/" + "{file}"
   threads: 16
   resources:
      mem_mb = 15000
   shell:"""
        {params.calcexp} \
        --paired-end \
        --no-bam-output \
        --quiet \
        --no-qualities \
        -p {threads} \
        --forward-prob 0.5 \
        --seed-length 25 \
        --fragment-length-mean -1.0 \
        --bam {input.bam} {params.genomedir} {params.prefix}
        """      
