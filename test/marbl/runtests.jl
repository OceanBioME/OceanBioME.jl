#####
##### Run both MARBL Fortran comparisons and report once, at the end.
#####
##### Each harness is a long run (tens of minutes) which prints as it goes, so reading its output while
##### it is still going will show you the previous section's numbers, or a previous run's. This driver
##### exists so there is one place that tells you the answer, and only when there is one.
#####
##### Run: julia --check-bounds=yes --project=test/marbl test/marbl/runtests.jl
#####

include(joinpath(@__DIR__, "marbl_names.jl"))
include(joinpath(@__DIR__, "marbl_baselines.jl"))

const HARNESSES = ("compare_marbl.jl", "test_stepped_marbl.jl")

if !baselines_available()
    @warn "MARBL comparisons skipped: baselines not present under $MARBL_INPUTS"
    exit(0)
end

results = Pair{String, Bool}[]

for harness in HARNESSES
    println("\n", "="^100)
    println("running $harness")
    println("="^100)

    # a separate process each: they define overlapping names, and one failing must not stop the other
    cmd = `$(Base.julia_cmd()) --check-bounds=yes --project=$(@__DIR__) $(joinpath(@__DIR__, harness))`

    ok = success(pipeline(cmd; stdout, stderr))

    push!(results, harness => ok)
end

println("\n", "="^100)
println("MARBL Fortran comparison summary")
println("="^100)

for (harness, ok) in results
    println("  ", rpad(harness, 26), ok ? "passed" : "FAILED")
end

all(last, results) || exit(1)
