#!/bin/tcsh

python $DNASTORAGE_LSF/generate_fi_jobs.py --params  $FRAMED_CONFIGS/iid/base4_RS_ideal_cluster.json --cores 128 --core_depth 1 --dump_dir $PWD --queue standard --time 24 --experiment_prefix framed_iid --conda_env_path $FRAMED_CONDA --job_name "rs_cluster_fi" 

python $DNASTORAGE_LSF/generate_fi_jobs.py --params  $FRAMED_CONFIGS/iid/basic_hedges_ideal_cluster.json --cores 128 --core_depth 1 --dump_dir $PWD --queue standard --time 24 --experiment_prefix framed_iid --conda_env_path $FRAMED_CONDA --job_name "hedges_cluster_fi" 

python $DNASTORAGE_LSF/generate_fi_jobs.py --params  $FRAMED_CONFIGS/iid/basic_hedges.json --cores 128 --core_depth 1 --dump_dir $PWD --queue standard --time  24 --experiment_prefix framed_iid  --conda_env_path $FRAMED_CONDA --job_name "hedges_fi" 


