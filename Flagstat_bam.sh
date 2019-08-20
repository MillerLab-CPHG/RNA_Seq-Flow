#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=15:00:00
#SBATCH --output=output_file
#SBATCH --partition=standard
#SBATCH -A mygroup       


module load samtools
# Input folder
DIRIN= /bam                                                                        


# Output file
stat=Pass1_flagstat.txt

for file in `find $DIRIN -type f -name "*.bam" | sort`; do
                echo $file >> $FOUT;
        samtools flagstat $file >> $stat;
done
