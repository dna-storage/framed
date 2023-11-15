#!/bin/tcsh

cd $DNASTORAGE_HOME

source dnastorage.env

cd $DNASTORAGE_LSF

echo $PWD

python generate_seq_jobs.py --dump_dir `pwd` --params $ONT_RNA_CONFIGS/target_file_1.json --memory 16 --cores 128 --time 5 --queue standard --sequencing_path $ONT_DATA --sequencing_regex "^([1-6])___.*\.fastq" --modules "PrgEnv-intel" --conda_env_path $FRAMED_CONDA --experiment_prefix 230316_rna 

python generate_seq_jobs.py --dump_dir `pwd` --params $ONT_RNA_CONFIGS/target_file_2.json --memory 16 --cores 128 --time 5 --queue standard --sequencing_path $ONT_DATA --sequencing_regex "^([7-9]|1[0-2])___.*\.fastq" --modules "PrgEnv-intel" --conda_env_path $FRAMED_CONDA  --experiment_prefix 230316_rna 

python generate_seq_jobs.py --dump_dir `pwd` --params $ONT_RNA_CONFIGS/target_file_3.json --memory 16 --cores 128 --time 5 --queue standard --sequencing_path $ONT_DATA --sequencing_regex "^(1[3-8])___.*\.fastq" --modules "PrgEnv-intel" --conda_env_path $FRAMED_CONDA --experiment_prefix 230316_rna 
