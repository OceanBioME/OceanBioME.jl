using OceanBioME, CairoMakie
# check that we cover the parameter space

cc = CarbonChemistry()

DIC = [1000.0:5.0:3000.0;]
Alk = [1000.0:5.0:3000.0;]
T = [-5:1.0:40;]
S = [0:1.0:40;]

function pCO₂_kernel!(idx, cc, DIC, Alk, T, S, pCO₂, pH)
    i, j, k, l = @inbounds CartesianIndices(pCO₂)[idx].I

    DIC = @inbounds DIC[i]
    Alk = @inbounds Alk[j]
    T = @inbounds T[k]
    S = @inbounds S[l]

    try
        pCO₂[i, j, k, l] = cc(; DIC, Alk, T, S)
        pH[i, j, k, l] = cc(; DIC, Alk, T, S, return_pH = true)
    catch
        pCO₂[i, j, k, l] = NaN
        pH[i, j, k, l] = NaN
    end
end

pCO₂ = zeros(length(DIC), length(Alk), length(T), length(S))
pH = similar(pCO₂)

backend = get_backend(pCO₂)

@info "Going for $(length(pCO₂)) samples on $(Threads.nthreads()) threads"

ts = time_ns()

Threads.@threads for idx in 1:length(pCO₂)
    pCO₂_kernel!(idx, cc, DIC, Alk, T, S, pCO₂, pH)
end

te = time_ns()

@info "Completed in $((te - ts)/1e9)s with $(sum(isnan.(pCO₂))/length(pCO₂)*100)% failure rate"