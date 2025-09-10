import Profiling
using Oceananigans
using CUDA

backend = GPU()

@info "Precompiling PISCES benchmark case"
@time Profiling.Cases.big_PISCES(;backend, fast_kill = true)

@info "Precompiling PISCES benchmark case again"
@time Profiling.Cases.big_PISCES(;backend, fast_kill = true)

@info "Precompiling PISCES benchmark case yet again!"
@time Profiling.Cases.big_PISCES(;backend, fast_kill = true)

@info "Running the PISCES case"
@time @CUDA.profile external=true Profiling.Cases.big_PISCES(;backend, fast_kill = false)

println(prof)

