import Profiling
using Oceananigans
using CUDA
using DataFrames
using CSV
using Base.Filesystem


function timed_run(::GPU, case_func)
    CUDA.@timed case_func
end

function timed_run(::CPU, case_func)
    @timed case_func
end


function run_benchmark_case(
    backend;
    enable_io::Bool,
    grid_size::Tuple{Int,Int,Int},
    n_particles::Int,
    runlength_scale::Float64,
)
    @info "Running benchmark case with grid size: $grid_size and $n_particles particles end with IO: $enable_io"
    filename = "SCALING_RUN_GPU"
    @info "Using output name: $filename"

    @info "Precompiling benchmark case"
    @time Profiling.Cases.big_LOBSTER(;
        backend,
        grid_size,
        n_particles,
        enable_io,
        fast_kill = true,
        runlength_scale,
        filename,
    )

    @info "Precompiling benchmark case again"
    @time Profiling.Cases.big_LOBSTER(;
        backend,
        grid_size,
        n_particles,
        enable_io,
        fast_kill = true,
        runlength_scale,
        filename,
    )

    @info "Precompiling bemchmark case yet again!"
    @time Profiling.Cases.big_LOBSTER(;
        backend,
        grid_size,
        n_particles,
        enable_io,
        fast_kill = true,
        runlength_scale,
        filename,
    )

    # The selection of the timing macro below is not pretty...
    # Refactor this to use proper dispatch
    @info "Running the case"
    run_data = timed_run(
        backend,
        () -> Profiling.Cases.big_LOBSTER(;
            backend,
            grid_size,
            n_particles,
            enable_io,
            fast_kill = false,
            runlength_scale,
            filename,
        ),
    )


    # run_data = if isa(backend, CPU)
    #     @timed Profiling.Cases.big_LOBSTER(;backend, grid_size, n_particles, enable_io, fast_kill = false, runlength_scale, filename)
    # elseif isa(backend, GPU)
    #     CUDA.@timed Profiling.Cases.big_LOBSTER(;backend, grid_size, n_particles, enable_io, fast_kill = false, runlength_scale, filename)
    # else
    #     error("Unknown backend of type: $(typeof(backend))")
    # end

    if enable_io
        # Remove the output if it exists
        rm("$(filename)_bgc.jld2")
        rm("$(filename)_particles.jld2")
    end

    return run_data
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
        n_particles = run_parameters.n_particles,
        runlength_scale = run_parameters.runlength_scale,
    )
end


"""
    repackge_run_data(::CPU, run_data)

Takes the row run_data collected via Base.@timed macro and flattens it into a
single NamedTuple with more explicit names for each field.

# Note
- We ignore 'lock_conflicts' number in the run_data at the moment.
- We use entries that were added in Julia 1.11!
"""
function repackge_run_data(::CPU, run_data)
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
    repackge_run_data(::GPU, run_data)

Takes the row run_data collected via CUDA.@timed macro and flattens it into a
single NamedTuple with more explicit names for each field.
"""
function repackge_run_data(::GPU, run_data)
    (;
        # Runtime data CPU
        elapsed_time = run_data.time,
        bytes_allocated = run_data.cpu_bytes,
        gc_time = run_data.cpu_gctime,
        compile_time = NaN, # Missing in the CUDA.jl stats
        recompile_time = NaN, # Missing in the CUDA.jl stats

        # Runtime data GPU
        gpu_bytes_allocated = run_data.gpu_bytes,
        gpu_memtime = run_data.gpu_memtime,

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


function repackge_run_data(run_data, run_parameters)
    (;
        # Run parameters
        grid_cells = prod(run_parameters.grid_size),
        Nx = run_parameters.grid_size[1],
        Ny = run_parameters.grid_size[2],
        Nz = run_parameters.grid_size[3],
        io = run_parameters.enable_io,
        n_particles = run_parameters.n_particles,
        runlength_scale = run_parameters.runlength_scale,

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


backend = CPU()

# Changes the grid size by the factors of 2
cases_grid = [
# (grid_size = (16, 16, 8), n_particles = 5, enable_io = true, runlength_scale = 1.0),
# (grid_size = (16, 32, 8), n_particles = 5, enable_io = true, runlength_scale = 1.0),
    (grid_size = (32, 32, 8), n_particles = 5, enable_io = true, runlength_scale = 1.0), # Nominal Case
# (grid_size = (32, 32, 16), n_particles = 5, enable_io = true, runlength_scale = 1.0),
# (grid_size = (64, 32, 16), n_particles = 5, enable_io = true, runlength_scale = 1.0),
# (grid_size = (64, 64, 16), n_particles = 5, enable_io = true, runlength_scale = 1.0),
# (grid_size = (64, 64, 32), n_particles = 5, enable_io = true, runlength_scale = 1.0),
# (grid_size = (128, 64, 32), n_particles = 5, enable_io = true, runlength_scale = 1.0),
# (grid_size = (128, 128, 32), n_particles = 5, enable_io = true, runlength_scale = 1.0),
# (grid_size = (128, 128, 64), n_particles = 5, enable_io = true, runlength_scale = 1.0),
# (grid_size = (256, 128, 64), n_particles = 5, enable_io = true, runlength_scale = 1.0),
# (grid_size = (256, 256, 64), n_particles = 5, enable_io = true, runlength_scale = 1.0),
# (grid_size = (256, 256, 128), n_particles = 5, enable_io = true, runlength_scale = 1.0),
]

df = DataFrame()
for case in cases_grid
    run_data = run_benchmark_case(backend; case...)

    data_to_save =
        (; repackage_run_parameters(case)..., repackge_run_data(backend, run_data)...)

    push!(df, data_to_save)
    # Override each timestep to preserve the data in case of faliure !
    CSV.write("scaling_run_gpu_withio.csv", df)
end
