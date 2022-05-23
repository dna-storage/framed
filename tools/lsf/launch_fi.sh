#!/bin/sh
#wrapper script for launching fault injection campaigns across many configurations

python $LSF_FI_DIR_PATH/command_generation_fi.py --config $1 #generates the command lines from a config file


jobCounter=0
while read line; do
    bsub -W 10000 -N -C 0 -n 12 -o $HOME/dna_fi_dumps/stdout_benchmark_$jobCounter.txt -e $HOME/dna_fi_dumps/sterr_benchmark_$jobCounter.txt -q long \
	 -P DNA_FI -R "select[hc] rusage[mem=8000]" python fault_injection.py  $line
    jobCounter=$((jobCounter+1))
done < intermediate.commands

rm intermediate.commands
