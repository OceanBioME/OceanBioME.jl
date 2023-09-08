#!/bin/bash
/home/js2430/rds/hpc-work/OceanBioME-runner/_work/_temp/julia-1.9.3/bin/julia -O0 --color=yes --project test/gpu_runtests.jl > "output.test"; ec=$?
if [ $ec -eq 0 ]; then
    echo "1" > "results.test"; 
else 
    echo "0" > "results.test"; 
fi