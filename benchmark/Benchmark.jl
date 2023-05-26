# These functions are copied from the brilliant Oceananigans team: https://github.com/CliMA/Oceananigans.jl/blob/main/benchmark/src/Benchmarks.jl
using BenchmarkTools, DataFrames, PrettyTables
function run_benchmarks(benchmark_fun; kwargs...)
    keys = [p.first for p in kwargs]
    vals = [p.second for p in kwargs]

    cases = Iterators.product(vals...)
    n_cases = length(cases)

    tags = string.(keys)
    suite = BenchmarkGroup(tags)
    for (n, case) in enumerate(cases)
        @info "Benchmarking $n/$n_cases: $case..."
        suite[case] = benchmark_fun(case...)
        GC.gc()
        GC.gc(true)
        GC.gc()
    end
    return suite
end

function benchmarks_dataframe(suite)
    names = Tuple(Symbol(tag) for tag in suite.tags)
    df_names = (names..., :min, :median, :mean, :max, :memory, :allocs, :samples)
    empty_cols = Tuple([] for k in df_names)
    df = DataFrame(; NamedTuple{df_names}(empty_cols)...)

    for case in keys(suite)
        trial = suite[case]
        entry = NamedTuple{names}(case) |> pairs |> Dict{Any,Any}

        entry[:min] = minimum(trial.times/1e9) |> prettytime
        entry[:median] = median(trial.times/1e9) |> prettytime
        entry[:mean] = mean(trial.times/1e9) |> prettytime
        entry[:max] = maximum(trial.times/1e9) |> prettytime
        entry[:memory] = trial.memory
        entry[:allocs] = trial.allocs
        entry[:samples] = length(trial)

        push!(df, entry)
    end

    return df
end

function benchmarks_pretty_table(df; title="")
    header = propertynames(df) .|> String

    html_filename = replace(title, ' ' => '_') * ".html"
    @info "Writing $html_filename..."
    open(html_filename, "w") do io
        html_table = pretty_table(String, df; header=header, nosubheader=true,
                                  title=title, title_alignment=:c,
                                  backend=Val(:html), tf=tf_html_simple)
        write(io, html_table)
    end

    return nothing
end