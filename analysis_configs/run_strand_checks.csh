#!/bin/tcsh


cd $DNASTORAGE_HOME

source dnastorage.env

cd $DNASTORAGE_LSF

echo $PWD

python generate_fi_jobs.py --params $FINAL_CONFIGS/sdc-muscle-final-strand_check.json --memory 16 --cores 1 --core_depth 1 --dump_dir $PWD --queue single_chassis --time 4 --experiment_prefix kevin_lin_strand_check --job_name strand_check --modules "PrgEnv-intel" "julia" --conda_env_path $FRAMED_CONDA 

python generate_fi_jobs.py --params $FINAL_CONFIGS/sdc-zebra-final-strand_check.json --memory 16 --cores 1 --core_depth 1 --dump_dir $PWD --queue single_chassis --time 4 --experiment_prefix kevin_lin_strand_check --job_name strand_check --modules "PrgEnv-intel" "julia" --conda_env_path $FRAMED_CONDA

python generate_fi_jobs.py --params $FINAL_CONFIGS/sdc-phos-final-strand_check.json --memory 16 --cores 1 --core_depth 1 --dump_dir $PWD --queue single_chassis --time 4 --experiment_prefix kevin_lin_strand_check --job_name strand_check --modules "PrgEnv-intel" "julia" --conda_env_path $FRAMED_CONDA
