#!/bin/csh
source dnastorage.env

git clone https://github.com/BHam-1/DNArSim/

rm -rf probEdit

mv DNArSim/simulator/probEdit $PWD

set patch_prefixes = ("channel" "functions" "loadProb")

foreach i ($patch_prefixes)
    set patch_file = "${DNArSimPath}/${i}.patch"
    set original_file = "DNArSim/simulator/${i}.jl"
    set output_file = "${DNArSimPath}/${i}.jl"
    patch -o $output_file $original_file $patch_file
end

rm -rf DNArSim

git rev-parse --git-dir >& /dev/null


set is_git_status = $?
if (${is_git_status} == 0) then
    git submodule update --init --recursive
else
    git clone https://github.com/kvolkel/schwimmbad
endif

if ( "$argv[1]" == "-no-env" ) then
    exit
endif

conda env create --prefix $FRAMED_CONDA --file $PWD/dnastorage.yml

conda activate $FRAMED_CONDA

pip install -r requirements.txt

cd schwimmbad; pip install .

cd ../

python -c "import julia; julia.install()"
