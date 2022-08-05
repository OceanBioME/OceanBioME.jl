using BGC, HDF5, Statistics, Interpolations, Plots

day=60*60*24
year=day*365

path="./subpolar/"    #subtropical   #./subpolar/
par_mean_timeseries=zeros(365)
for i in 1:365    #https://discourse.julialang.org/t/leading-zeros/30450
    string_i = lpad(string(i), 3, '0')
    filename3=path*"V2020"*string_i*".L3b_DAY_SNPP_PAR.x.nc"
    fid = h5open(filename3, "r")
    par=read(fid["level-3_binned_data/par"])
    BinList=read(fid["level-3_binned_data/BinList"])  #(:bin_num, :nobs, :nscenes, :weights, :time_rec) 
    par_mean_timeseries[i] = mean([par[i][1]/BinList[i][4] for i in 1:length(par)])*3.99e-10*545e12/(1day)  #from einstin/m^2/day to W/m^2
end

surface_PAR_itp = LinearInterpolation((0:364)*day, par_mean_timeseries)

z=-10

Pᵢ = (tanh((z+250)/100)+1)/2*(0.038)+0.002          # ((tanh((z+100)/50)-1)/2*0.23+0.23)*16/106  
Zᵢ = (tanh((z+250)/100)+1)/2*(0.038)+0.008          # ((tanh((z+100)/50)-1)/2*0.3+0.3)*16/106         
Dᵢ =0
DDᵢ =0
NO₃ᵢ = (1-tanh((z+300)/150))/2*6+11.4   #  # 17.5*(1-tanh((z+100)/10))/2
NH₄ᵢ = (1-tanh((z+300)/150))/2*0.05+0.05       #1e-1*(1-tanh((z+100)/10))/2
DOMᵢ = 0 

PAR(x, y, z, t) = surface_PAR_itp(mod(t,364*day))*exp(z*0.2)#this won't adjust with the P conc.

params = merge(LOBSTER.default, (PAR=PAR, ))

model = Setup.BoxModel(:LOBSTER, params, (P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ), 0.0, 1.0*year)
solution = BoxModel.run(model)
values = vcat(transpose.(solution.u)...)
plts=[]
for (i, tracer) in enumerate(model.tracers)
    push!(plts, plot(solution.t[1:400:end]/day, values[1:400:end, i], ylabel=tracer, xlabel="Day", legend=false))
end
plot(plts...)
