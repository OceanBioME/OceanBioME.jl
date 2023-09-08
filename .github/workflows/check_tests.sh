#!/bin/bash

if ["$(/home/js2430/rds/hpc-work/OceanBioME-runner/_work/_temp/julia-1.9.3/bin/julia -O0 --color=yes --project test/gpu_runtests.jl)" == ""]; then 
    echo "true" > "results.test"; 
else 
    echo "false" > "results.test"; 
fi