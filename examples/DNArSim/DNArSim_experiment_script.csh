#!/bin/tcsh

python $DNASTORAGE_LSF/generate_fi_jobs.py --params  $FRAMED_CONFIGS/DNArSim/base4_RS_ideal_cluster_nano.json --cores 128 --core_depth 1 --dump_dir $PWD --queue standard --time 24 --experiment_prefix framed_DNArSIM --job_name "DNArSim_rs_cluster_fi" --conda_env_path $FRAMED_CONDA

python $DNASTORAGE_LSF/generate_fi_jobs.py --params $FRAMED_CONFIGS/DNArSim/basic_hedges_cluster_ideal_nano.json --cores 128 --core_depth 1 --dump_dir $PWD --queue standard --time 24 --experiment_prefix framed_DNArSIM --job_name "DNArSim_hedges_cluster_fi" --conda_env_path $FRAMED_CONDA

python $DNASTORAGE_LSF/generate_fi_jobs.py --params $FRAMED_CONFIGS/DNArSim/basic_hedges_nano.json --cores 128 --core_depth 1 --dump_dir $PWD --queue standard --time 24 --experiment_prefix framed_DNArSIM  --job_name "DNArSim_hedges_fi" --conda_env_path $FRAMED_CONDA 


 
