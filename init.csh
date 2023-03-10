#!/bin/csh

git clone https://github.com/BHam-1/DNArSim/

mv DNArSim/simulator/probEdit $PWD

set patch_prefixes = ("channel" "functions" "loadProb")

foreach i ($patch_prefixes)
    set patch_file = "${DNArSimPath}/${i}.patch"
    set original_file = "DNArSim/simulator/${i}.jl"
    set output_file = "${DNArSimPath}/${i}.jl"
    patch -o $output_file $original_file $patch_file
end

rm -rf DNArSim

conda env create --prefix $FRAMED_CONDA --file $PWD/dnastorage.yml

conda activate $FRAMED_CONDA

git submodule update --init --recursive

cd schwimmbad; pip install .

cd ../

python -c "import julia; julia.install()"
