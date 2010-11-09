#$ -N %(jobname)s
#$ -S /bin/bash
#$ -l mem_free=500M,mem_total=500M
#$ -cwd
#$ -j y
#$ -r y

echo "APBS Born profile job running on $HOSTNAME"
ulimit -c 64
module load abps/32 
apbs %(infile)s > %(outfile)s
