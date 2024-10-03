grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 1))

particle_biogeochemistry = OceanBioME.Models.SugarKelpModel.SugarKelp()

particles = BiogeochemicalParticles(2; grid, 
                                       biogeochemistry = particle_biogeochemistry, 
                                       advection = nothing)

biogeochemistry = LOBSTER(; grid, 
                            particles, 
                            carbonates = true, 
                            variable_redfield = true, 
                            oxygen = true, 
                            sinking_speeds = NamedTuple())

model = NonhydrostaticModel(; grid, biogeochemistry, advection = nothing, tracers = (:T, :S))

set!(model, NO₃ = 10.0, NH₄ = 1.0, DIC = 2000, Alk = 2000, T = 10, S = 35)

particles.x .= 0.5
particles.y .= 0.5
particles.z .= -0.5

particles.fields.A .= 2
particles.fields.N .= 1
particles.fields.C .= 1