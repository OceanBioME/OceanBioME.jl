import Profiling
using Oceananigans
using CUDA
using DataFrames
using CSV
using Base.Filesystem


function timed_run(::GPU, case_func)
    CUDA.@timed case_func()
end

function timed_run(::CPU, case_func)
    @timed case_func()
end

mutable struct GPUMemoryUsageStats
    max_gpu_bytes::Union{Int,Nothing}
    min_gpu_bytes::Union{Int,Nothing}
end


"""
    map_or_else(opt::Union{T,Nothing}, mapf, else_value)

If option coontains a value applies mapf and returns the result, otherwise
returns else_value.
"""
function map_or_else(opt::Union{T,Nothing}, mapf, else_value) where {T}
    isnothing(opt) ? else_value : mapf(opt)
end

function update!(stats::GPUMemoryUsageStats, ::GPU)
    # Print info about GPU memory usage
    used_mem = CUDA.used_memory()
    cached_mem = CUDA.cached_memory()
    println(
        "GPU memory usage: $(Base.format_bytes(used_mem)) Pool size: $(Base.format_bytes(cached_mem))",
    )

    stats.max_gpu_bytes = map_or_else(stats.max_gpu_bytes, v -> max(v, used_mem), used_mem)
    stats.min_gpu_bytes = map_or_else(stats.min_gpu_bytes, v -> min(v, used_mem), used_mem)

end

function update!(stats::GPUMemoryUsageStats, ::CPU)
    # No-op for CPU
end

function run_PISCES_mem_benchmark_case(
    backend;
    enable_io::Bool,
    grid_size::Tuple{Int,Int,Int},
    runlength_scale::Float64,
)

    @info "Running PISCES memory benchmark case with grid size: $grid_size and with IO: $enable_io"
    filename = "PISCES_MEMORY_SCALING_RUN_GPU"
    @info "Using output name: $filename"

    @info "Precompiling PISCES memory benchmark case"
    @time Profiling.Cases.big_PISCES(;
        backend,
        grid_size,
        enable_io,
        fast_kill = true,
        runlength_scale,
        filename,
    )

    @info "Precompiling PISCES memory benchmark case again"
    @time Profiling.Cases.big_PISCES(;
        backend,
        grid_size,
        enable_io,
        fast_kill = true,
        runlength_scale,
        filename,
    )

    @info "Precompiling PISCES memory benchmark case yet again!"
    @time Profiling.Cases.big_PISCES(;
        backend,
        grid_size,
        enable_io,
        fast_kill = true,
        runlength_scale,
        filename,
    )

    @info "Precompiling PISCES memory benchmark case"
    @time Profiling.Cases.big_PISCES(;
        backend,
        grid_size,
        enable_io,
        fast_kill = true,
        runlength_scale,
        filename,
    )

    @info "Precompiling PISCES memory benchmark case again"
    @time Profiling.Cases.big_PISCES(;
        backend,
        grid_size,
        enable_io,
        fast_kill = true,
        runlength_scale,
        filename,
    )

    @info "Precompiling PISCES memory benchmark case yet again!"
    @time Profiling.Cases.big_PISCES(;
        backend,
        grid_size,
        enable_io,
        fast_kill = true,
        runlength_scale,
        filename,
    )

    gpu_mem_bounds = GPUMemoryUsageStats(nothing, nothing)

    function hook()
        update!(gpu_mem_bounds, backend)
    end

    GC.gc(true)
    GC.gc(true)
    GC.gc(true)
    CUDA.reclaim()

    # The selection of the timing macro below is not pretty...
    # Refactor this to use proper dispatch
    @info "Running the PISCES memory benchmark case"
    run_data = timed_run(
        backend,
        () -> Profiling.Cases.big_PISCES(;
            backend,
            grid_size,
            enable_io,
            fast_kill = false,
            runlength_scale,
            filename,
            progress_hook = hook,
        ),
    )

    if enable_io
        # Remove the output if it exists
        rm("$(filename)_bgc.jld2")
    end

    return (;
        run_data...,
        gpu_max_bytes = gpu_mem_bounds.max_gpu_bytes,
        gpu_min_bytes = gpu_mem_bounds.min_gpu_bytes,
    )
end


"""
    repackage_run_parameters(run_parameters)

Takes NamedTuple with the parameters defining case and repackages them to
field we wish to preserve in the CSV file.
"""
function repackage_run_parameters(run_parameters)
    (;
        grid_cells = prod(run_parameters.grid_size),
        Nx = run_parameters.grid_size[1],
        Ny = run_parameters.grid_size[2],
        Nz = run_parameters.grid_size[3],
        io = run_parameters.enable_io,
        runlength_scale = run_parameters.runlength_scale,
    )
end

"""
    repackage_run_data(::CPU, run_data)

Takes the row run_data collected via Base.@timed macro and flattens it into a
single NamedTuple with more explicit names for each field.

# Note
- We ignore 'lock_conflicts' number in the run_data at the moment.
- We use entries that were added in Julia 1.11!
"""
function repackage_run_data(::CPU, run_data)
    (;
        # Runtime data
        elapsed_time = run_data.time,
        bytes_allocated = run_data.bytes,
        gc_time = run_data.gctime,
        compile_time = run_data.compile_time,
        recompile_time = run_data.recompile_time,

        # Garbage collector data
        gc_bytes_allocated = run_data.gcstats.allocd,
        gc_number_of_pauses = run_data.gcstats.pause,
        gc_number_of_malloc_calls = run_data.gcstats.malloc,
        gc_number_of_ralloc_calls = run_data.gcstats.realloc,
        gc_number_of_pool_allocations = run_data.gcstats.poolalloc,
        gc_number_of_big_nonpool_allocations = run_data.gcstats.bigalloc,
        gc_number_of_free_calls = run_data.gcstats.freecall,
        gc_number_of_full_sweeps = run_data.gcstats.full_sweep,
    )
end

"""
    repackage_run_data(::GPU, run_data)

Takes the row run_data collected via CUDA.@timed macro and flattens it into a
single NamedTuple with more explicit names for each field.
"""
function repackage_run_data(::GPU, run_data)
    (;
        # Runtime data CPU
        elapsed_time = run_data.time,
        bytes_allocated = run_data.cpu_bytes,
        gc_time = run_data.cpu_gctime,
        compile_time = NaN, 
        recompile_time = NaN, 

        # Runtime data GPU
        gpu_bytes_allocated = run_data.gpu_bytes,
        gpu_memtime = run_data.gpu_memtime,

        # Data from the hook
        gpu_max_bytes = run_data.gpu_max_bytes,
        gpu_min_bytes = run_data.gpu_min_bytes,

        # CPU Garbage collector data
        gc_bytes_allocated = run_data.cpu_gcstats.allocd,
        gc_number_of_pauses = run_data.cpu_gcstats.pause,
        gc_number_of_malloc_calls = run_data.cpu_gcstats.malloc,
        gc_number_of_ralloc_calls = run_data.cpu_gcstats.realloc,
        gc_number_of_pool_allocations = run_data.cpu_gcstats.poolalloc,
        gc_number_of_big_nonpool_allocations = run_data.cpu_gcstats.bigalloc,
        gc_number_of_free_calls = run_data.cpu_gcstats.freecall,
        gc_number_of_full_sweeps = run_data.cpu_gcstats.full_sweep,

        # GPU memory stats
        gpu_bytes_allocated_stats = run_data.gpu_memstats.alloc_bytes,
        gpu_number_of_allocations = run_data.gpu_memstats.alloc_count,
        gpu_bytes_freed = run_data.gpu_memstats.free_bytes,
        gpu_number_of_frees = run_data.gpu_memstats.free_count,
    )
end

backend = GPU()

cases_PISCES_grid = [
    (grid_size = (16, 16, 8), enable_io = false, runlength_scale = 0.5),
    (grid_size = (16, 32, 8), enable_io = false, runlength_scale = 0.5),
    (grid_size = (32, 32, 8), enable_io = false, runlength_scale = 0.5), 
    (grid_size = (32, 32, 16), enable_io = false, runlength_scale = 0.5),
    (grid_size = (64, 32, 16), enable_io = false, runlength_scale = 0.5),
    (grid_size = (64, 64, 16), enable_io = false, runlength_scale = 0.5),
    (grid_size = (64, 64, 32), enable_io = false, runlength_scale = 0.5),
    (grid_size = (128, 64, 32), enable_io = false, runlength_scale = 0.5),
    (grid_size = (128, 128, 32), enable_io = false, runlength_scale = 0.5),
    (grid_size = (128, 128, 64), enable_io = false, runlength_scale = 0.5),
    (grid_size = (256, 128, 64), enable_io = false, runlength_scale = 0.5),
    (grid_size = (256, 256, 64), enable_io = false, runlength_scale = 0.5),
    (grid_size = (256, 256, 128), enable_io = false, runlength_scale = 0.5),
]

df = DataFrame()

for case in cases_PISCES_grid
    run_data = run_PISCES_mem_benchmark_case(backend; case...)

    println("Run data: $run_data")
    data_to_save =
        (; repackage_run_parameters(case)..., repackage_run_data(backend, run_data)...)


    push!(df, data_to_save)
    # Override each timestep to preserve the data in case of faliure !
    CSV.write(
        "mem_scaling_PISCES_run_gpu_without_io.csv",
        df;
        transform = (col, val) -> something(val, missing),
    )
end












