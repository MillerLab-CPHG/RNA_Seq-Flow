# RNA_Seq-Analysis-workflow

RNA-Seq analyis using STAR two pass mode for alliging the raw reads and RSEM for transcript quantification using snakemake workflow. 

Required Tools

1 - FastQC (Quality control)

2 - Trim-galore ( Adapter trimming)

3 - STAR (Allign the reads to reference genome) 

4 - Picard (Remove Duplicates)

5 - RNA-SeQC(Qualiy Control of RNA-Seq data)

6 - RSEM(Transcript Quantification)


# run the pipeline on cluster using this cammand 'modify parameters according to your cluster configuration

snakemake -j 999 --configfile config.yaml --use-conda --nolock --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition}  -N {cluster.N} -n {cluster.n}  -t {cluster.time} --mem {cluster.mem}"
