#wrapper script for launching fault injection campaigns across many configurations



while read line; do

     bsub -W 100 -N -C 0 -n 12 -o $HOME/dna_fi_dumps/stdout_benchmark.txt -e $HOME/dna_fi_dumps/sterr_benchmark.txt -q standard -P DNA_FI -R "select[hc] rusage[mem=8000]"  $line

done < $1 
