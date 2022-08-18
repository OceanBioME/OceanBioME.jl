# This script illustrates how to run OceanBioME as a box model

using OceanBioME, HDF5, Statistics, Interpolations, Plots  # load required modules

day=60*60*24  # define the length of a day in seconds
year=day*365  # define the length of a year in days

# This demonstrates how to read in a timeseries of photosynthetic available radiation (PAR) data
path="./OceanBioME_example_data/subpolar/"  # the folder where the data are stored
par_mean_timeseries=zeros(365) # create an empty array
for i in 1:365  # loop over day number  
    string_i = lpad(string(i), 3, '0')  # create string with day number, padding with zeros to get a 3 digit number, e.g. 001
    filename3=path*"V2020"*string_i*".L3b_DAY_SNPP_PAR.x.nc" 
    fid = h5open(filename3, "r") 
    par=read(fid["level-3_binned_data/par"])  # read PAR data from file
    BinList=read(fid["level-3_binned_data/BinList"])  # read information on PAR bins.  The format of BinList is (:bin_num, :nobs, :nscenes, :weights, :time_rec) 
    par_mean_timeseries[i] = mean([par[i][1]/BinList[i][4] for i in 1:length(par)])*3.99e-10*545e12/(1day)  # average PAR values in bins and convert from einstin/m^2/day to W/m^2
end

surface_PAR_itp = LinearInterpolation((0:364)*day, par_mean_timeseries) # create a function to interpolate


# Set the initial conditions
Pᵢ = 0.1
Zᵢ = 0.01            
Dᵢ = 0
DDᵢ = 0
NO₃ᵢ = 1
NH₄ᵢ = 0.1
DOMᵢ = 0 

z=-10 # specify the depth for the light level
PAR(x, y, z, t) = surface_PAR_itp(mod(t,364*day))*exp(z*0.2) # Set the PAR

# Create a list of parameters.  Here, use the default parameters from the LOBSTER model with the PAR function
params = merge(LOBSTER.defaults, (PAR=PAR, ))

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = Setup.BoxModel(:LOBSTER, params, (NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, DOM=DOMᵢ), 0.0, 1.0*year)

solution = BoxModel.run(model) # call BoxModel to timestep the biogeochemical model

# Create an array containing all model varaibles
# The dimensions of 'values' will be (time, variable number)
values = vcat(transpose.(solution.u)...)

# build a series of plots for all tracers contained in model.tracers
plts=[]
for (i, tracer) in enumerate(model.tracers)
    push!(plts, plot(solution.t[1:400:end]/day, values[1:400:end, i], ylabel=tracer, xlabel="Day", legend=false))
end

# Display the resulting figure
plot(plts...)
