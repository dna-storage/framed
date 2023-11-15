#!/bin/tcsh

cd $DNASTORAGE_HOME

source dnastorage.env

cd $DNASTORAGE_LSF

echo $PWD

python generate_fi_jobs.py --params $MISEQ_CONFIGS/muscle_sequencing_fi.json --cores 1 --core_depth 64 --dump_dir $PWD --queue standard --time 4 --experiment_prefix kevin_lin_230112_miseq_230619 --job_name Lin_Miseq_FI --modules "PrgEnv-intel" "julia" --conda_env_path $FRAMED_CONDA 

python generate_fi_jobs.py --params $MISEQ_CONFIGS/zebra_sequencing_fi.json  --cores 1 --core_depth 64 --dump_dir $PWD --queue standard --time 4 --experiment_prefix kevin_lin_230112_miseq_230619 --job_name Lin_Miseq_FI --modules "PrgEnv-intel" "julia" --conda_env_path $FRAMED_CONDA

python generate_fi_jobs.py --params $MISEQ_CONFIGS/phos_sequencing_fi.json --cores 1 --core_depth 64 --dump_dir $PWD --queue standard --time 4 --experiment_prefix kevin_lin_230112_miseq_230619 --job_name Lin_Miseq_FI --modules "PrgEnv-intel" "julia" --conda_env_path $FRAMED_CONDA
