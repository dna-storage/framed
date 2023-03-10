#!/bin/tcsh

python $DNASTORAGE_LSF/generate_fi_jobs.py --params $FRAMED_CONFIGS/small/base4_RS_ideal_cluster_iid.json --cores 3 --core_depth 1 --dump_dir $PWD --experiment_prefix framed_small --conda_env_path $FRAMED_CONDA --submission "shell"

python $DNASTORAGE_LSF/generate_fi_jobs.py --params  $FRAMED_CONFIGS/small/basic_hedges_iid.json --cores 3 --core_depth 1 --dump_dir $PWD --experiment_prefix framed_small --conda_env_path $FRAMED_CONDA  --submission "shell"

python $DNASTORAGE_LSF/generate_fi_jobs.py --params  $FRAMED_CONFIGS/small/basic_hedges_DNArSim.json --cores 3 --core_depth 1 --dump_dir $PWD --experiment_prefix framed_small --conda_env_path $FRAMED_CONDA  --submission "shell"

