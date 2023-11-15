#!/bin/tcsh

cd $DNASTORAGE_HOME

source dnastorage.env

cd $DNASTORAGE_LSF

echo $PWD

set JSON_PATH = $ONT_RNA_CONFIGS

python generate_fi_jobs.py --params $JSON_PATH/muscle_sequencing_fi.json --cores 1 --core_depth 128 --dump_dir $PWD --queue standard --time 24 --experiment_prefix kevin_lin_230322_rna --job_name Lin_rna_FI --modules "PrgEnv-intel" "julia" --conda_env_path $FRAMED_CONDA --memory 16

python generate_fi_jobs.py --params $JSON_PATH/zebra_sequencing_fi.json  --cores 1 --core_depth 128 --dump_dir $PWD --queue standard --time 24 --experiment_prefix kevin_lin_230322_rna --job_name Lin_rna_FI --modules "PrgEnv-intel" "julia" --conda_env_path $FRAMED_CONDA  --memory 16

python generate_fi_jobs.py --params $JSON_PATH/phos_sequencing_fi.json --cores 1 --core_depth 128 --dump_dir $PWD --queue standard --time 24 --experiment_prefix kevin_lin_230322_rna --job_name Lin_rna_FI --modules "PrgEnv-intel" "julia" --conda_env_path $FRAMED_CONDA --memory 16
