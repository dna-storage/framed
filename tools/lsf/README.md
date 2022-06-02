This directory contains some helper scripts to test out different encoding architectures in various fault injection environments. The main script to call is done with

`python generate_fi_jobs.py --params <param_json> --memory <mem> --cores <cores> --time <time> --queue <queue> --dump_dir <dump_dir>`

`<param_json>` is the path to a json file that provides the main paramters to the overall simulation environment including the encoder, read distributions and error models. An example is given in `test_config.json`. The parameter names given in this file must much parameters of `tools/fault_injection.py`, and there must be 4 dictionaries, `encoder_params`, `header_params`,`fault_params`,`distribution_params`. These dictionaries are used to bundle related parameters together, and each of the parameter values in `fault_params` and `distribution_params` must have the same name as those in `tools/fault_injection.py` parameters. The example file provides an example for fixed fault rate with a poisson distribution. The parameters within `encoder_params`, and `header_params` should match those params expected as `kwargs` parameters in the `dnastorage/arch/builder.py` function for the given architecture. You can either specify a constant value for each parameter in the dictionary sets, or you can enter a range. A range is given as a json array `[start,stop,step]`. Alternatively mutiple options can be specificied as a list such as `["value_list",value_1,value2, ... ]` where you can have any n values. In this case, each value is an exact optoin, rather than a range. By using ranges and value lists, each option will be combined with every other possible option value to instantiate a unique encoding and simulation environment. The generate_fi_jobs.py will ensure that unique locatoins are made to dump data. All parameters related to a given instantiation will be at the end of the directory tree made for the given experiment. 

`<mem>` is memory to use on each experiment
`<cores>` is the number of cores to use for each experiment
`<time>` is the max time to run an experiment
`<queue>` is the hpc queue
`<dump_dir>` is the base directory to start the directory trees usef to store information for each experiment.


