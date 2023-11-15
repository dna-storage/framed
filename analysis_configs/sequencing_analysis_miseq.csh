#!/bin/tcsh

cd $DNASTORAGE_HOME

source dnastorage.env

cd $DNASTORAGE_LSF

echo $PWD

python generate_seq_jobs.py --dump_dir `pwd` --params $MISEQ_CONFIGS/target_file_1.json --memory 16 --cores 128 --time 5 --queue standard --sequencing_path $MISEQ_DATA --sequencing_regex "^([1-9]|[1-5][0-9]|60)___.*\.fastq" --modules "PrgEnv-intel" --conda_env_path $FRAMED_CONDA

python generate_seq_jobs.py --dump_dir `pwd` --params $MISEQ_CONFIGS/target_file_2.json --memory 16 --cores 128 --time 5 --queue standard --sequencing_path $MISEQ_DATA --sequencing_regex "^(6[1-9]|7[0-9]|8[0-1])___.*\.fastq" --modules "PrgEnv-intel" --conda_env_path $FRAMED_CONDA

python generate_seq_jobs.py --dump_dir `pwd` --params $MISEQ_CONFIGS/target_file_3.json --memory 16 --cores 128 --time 5 --queue standard --sequencing_path $MISEQ_DATA --sequencing_regex "^(8[2-9]|9[0-9]|10[0-2])___.*\.fastq" --modules "PrgEnv-intel" --conda_env_path $FRAMED_CONDA
