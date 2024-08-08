# # [Box model](@id box_example)
# In this example we setup a [LOBSTER](@ref LOBSTER) biogeochemical model in a single box configuration.
# This example demonstrates:
# - How to setup OceanBioME's biogeochemical models as a stand-alone box model

# ## Install dependencies
# First we check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME"
# ```

# ## Model setup
# Load the packages and setup the initial and forcing conditions
using OceanBioME, Oceananigans, Oceananigans.Units
using Oceananigans.Fields: FunctionField

const year = years = 365day
nothing #hide

grid = BoxModelGrid()
clock = Clock(time = zero(grid))

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

const z = - 10 # specify the nominal depth of the box for the PAR profile
# Modify the PAR based on the nominal depth and exponential decay
PAR_func(t) = 300.0 # Modify the PAR based on the nominal depth and exponential decay

PAR_func1(t) = 100.0
PAR_func2(t) = 100.0
PAR_func3(t)= 100.0

PAR = FunctionField{Center, Center, Center}(PAR_func, grid; clock)
PAR¹ = FunctionField{Center, Center, Center}(PAR_func1, grid; clock)
PAR² = FunctionField{Center, Center, Center}(PAR_func2, grid; clock)
PAR³ = FunctionField{Center, Center, Center}(PAR_func3, grid; clock)

nothing #hide

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = BoxModel(; biogeochemistry = PISCES(; grid, light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR¹, PAR², PAR³)), flux_feeding_rate = 2.0e-3),
                   clock)
#set!(model, P = 0.0, D = 1.0, Z = 2.0, M = 3.0, Pᶜʰˡ= 4.0, Dᶜʰˡ = 5.0, Pᶠᵉ = 6.0, Dᶠᵉ = 7.0, Dˢⁱ = 8.0, DOC = 9.0, POC = 10.0, GOC = 11.0, SFe = 12.0, BFe = 13.0, PSi = 14.0, NO₃ = 15.0, NH₄ = 16.0, PO₄ = 17.0, Fe = 18.0, Si = 19.0, CaCO₃ = 20.0, DIC = 21.0, Alk = 22.0, O₂ = 23.0, T = 14.0)
set!(model, NO₃ = 4.0, NH₄ = 0.1, P = 4.26, D = 4.26, Z = 0.426, M = 0.426, Pᶠᵉ = 3.50, Dᶠᵉ = 3.50, Pᶜʰˡ = .1, Dᶜʰˡ = .01, Dˢⁱ = 0.525, Fe = 0.00082410, O₂ = 264.0, Si = 4.557, Alk = 2360.0, PO₄ = 0.8114, DIC = 2000.0, CaCO₃ = 0.0001, T = 14.0)
#set!(model,P = 3.963367728460601, D = 3.763831823528108, Z = 0.620887631503286,  M = 0.4911996116700677, Pᶜʰˡ = 0.1263393104069646,  Dᶜʰˡ = 0.0966272698878372, Pᶠᵉ = 2.916749891527781, Dᶠᵉ = 2.6966762460922764, Dˢⁱ = 0.5250058442518801, DOC = 5.492834645446811e-5, POC = 0.00010816947467085888, GOC = 1.541376629008023, SFe = 6.94778354330689e-5, BFe = 1.3780182342394662, PSi = 0.138718322180627, NO₃ = 3.862629483089866, NH₄ = 0.10480738012675432, PO₄ = 0.8031309301476024, Fe = 0.00024547654218086575, Si = 4.413896794698411, CaCO₃ = 0.011644257272404535, DIC = 1998.9796292207268, Alk = 2360.118267032333, O₂ = 265.37453137881016, T = 14.0)

simulation = Simulation(model; Δt = 5minutes, stop_time = 40minutes)

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box.jld2", schedule = TimeInterval(5minutes), overwrite_existing = true)

prog(sim) = @info "$(prettytime(time(sim))) in $(prettytime(simulation.run_wall_time))"

function non_zero_fields!(model) 
    #for tracer in model.fields
      #  if tracer != model.fields[:T]
       #     tracer = tracer[1,1,1]
        #    getproperty(model, tracer) = max(0.0, getproperty(model, tracer))
        #end
    @inbounds for (idx, field) in enumerate(model.fields)
        if isnan(field[1,1,1])
            throw("$(keys(model.fields)[idx]) has gone NaN")
        else
            field[1, 1, 1] = max(0, field[1, 1, 1])
        end
        
    end
    #end
    return nothing

end

simulation.callbacks[:progress] = Callback(prog, IterationInterval(100))
simulation.callbacks[:non_zero_fields] = Callback(non_zero_fields!, callsite = UpdateStateCallsite())


@info "Running the model..."
run!(simulation)

# ## Load the output

times = FieldTimeSeries("box.jld2", "P").times

timeseries = NamedTuple{keys(model.fields)}(FieldTimeSeries("box.jld2", "$field")[1, 1, 1, :] for field in keys(model.fields))

# ## And plot
using CairoMakie

fig = Figure(size = (1200, 7200), fontsize = 24)

axs = []
for (name, tracer) in pairs(timeseries)
    idx = (length(axs))
    push!(axs, Axis(fig[floor(Int, idx/2), Int(idx%2)], ylabel = "$name", xlabel = "Year", xticks=(0:10)))
    lines!(axs[end], times / year, tracer, linewidth = 3)
end

fi = length(timeseries.P)

println("P = $(timeseries.P[fi]), D = $(timeseries.D[fi]), Z = $(timeseries.Z[fi]),  M = $(timeseries.M[fi]), Pᶜʰˡ = $(timeseries.Pᶜʰˡ[fi]),  Dᶜʰˡ = $(timeseries.Dᶜʰˡ[fi]), Pᶠᵉ = $(timeseries.Pᶠᵉ[fi]), Dᶠᵉ = $(timeseries.Dᶠᵉ[fi]), Dˢⁱ = $(timeseries.Dˢⁱ[fi]), DOC = $(timeseries.DOC[fi]), POC = $(timeseries.POC[fi]), GOC = $(timeseries.GOC[fi]), SFe = $(timeseries.SFe[fi]), BFe = $(timeseries.BFe[fi]), PSi = $(timeseries.PSi[fi]), NO₃ = $(timeseries.NO₃[fi]), NH₄ = $(timeseries.NH₄[fi]), PO₄ = $(timeseries.PO₄[fi]), Fe = $(timeseries.Fe[fi]), Si = $(timeseries.Si[fi]), CaCO₃ = $(timeseries.CaCO₃[fi]), DIC = $(timeseries.DIC[fi]), Alk = $(timeseries.Alk[fi]), O₂ = $(timeseries.O₂[fi]), T = 14.0")

Carbon_at_start = timeseries.P[1] + timeseries.D[1] + timeseries.Z[1] + timeseries.M[1] + timeseries.DOC[1] + timeseries.POC[1] + timeseries.GOC[1] + timeseries.DIC[1] + timeseries.CaCO₃[1]
Carbon_at_end = timeseries.P[fi] + timeseries.D[fi] + timeseries.Z[fi] + timeseries.M[fi] + timeseries.DOC[fi] + timeseries.POC[fi] + timeseries.GOC[fi] + timeseries.DIC[fi] + timeseries.CaCO₃[fi]

println("Carbon at start = ", Carbon_at_start, "Carbon at end = ", Carbon_at_end)

fig