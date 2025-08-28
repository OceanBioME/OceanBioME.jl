using OceanBioME
using Plots
using Measurements

#Plots.scalefontsizes(1.5) # only do this once

carbon_chemistry = CarbonChemistry()

DIC = 2100.0 ± 1.0
Alk = 2300.0 ± 1.0
T = 12 ± 0.01
S = 35 ± 0.01
# fCO2 = carbon_chemistry(DIC = DIC, Alk = Alk, T=T, S=S)
pH = carbon_chemistry(DIC = DIC, Alk = Alk, T=T, S=S, return_pH=true)
# upper_pH_bound=measurement(14.0),
# lower_pH_bound=measurement(0.0),
# atol = measurement(1e-10),
# guess = (8.0 ± 0.0))

scatter([DIC], [fCO2], markersize=8, legend = nothing)

xlims!(2090, 2110)
ylims!(380, 390)
xlabel!("DIC (μmol/kg)")
ylabel!("fCO₂ (μatm)")

pH = 8.05 ± 0.001
fCO2 = carbon_chemistry(DIC = DIC, pH = pH, T=T, S=S)
# , upper_pH_bound=measurement(14.0),
# lower_pH_bound=measurement(0.0),
# atol = measurement(1e-10),
# guess = (8.0 ± 0.0))

scatter([DIC],[fCO2],markersize=8, label=nothing, markercolor=:red)

xlims!(2090, 2110)
ylims!(380, 390)
xlabel!("DIC (μmol/kg)")
ylabel!("fCO₂ (μatm)")


DIC_measurement = 2100.0 ± 1.0
Alk_measurement = 2300.0 ± 4.0
T_measurement = 12.0 ± 0.1
S_measurement = 35.0 ± 0.1

fCO2 = carbon_chemistry(DIC=DIC_measurement, 
                        Alk=Alk_measurement, 
                        T=T_measurement, 
                        S=S_measurement)
                        # ,
                        # upper_pH_bound=measurement(14.0),
                        # lower_pH_bound=measurement(0.0),
                        # atol = measurement(1e-6),
                        # guess = (8.0 ± 0.1))

scatter([DIC_measurement], [fCO2], markersize=6, 
        label="Alk measurement")

pH_measurement = 8.05 ± 0.001

fCO2 = carbon_chemistry(DIC=DIC_measurement, 
                        pH=pH_measurement, 
                        T=T_measurement, 
                        S=S_measurement)
                        # ,
                        # upper_pH_bound=measurement(14.0),
                        # lower_pH_bound=measurement(0.0),
                        # atol = measurement(1e-10),
                        # guess = (8.0 ± 0.1))

#scatter([DIC_measurement], [fCO2], markersize=6, 
#        label="pH measurement")
xlims!(2090, 2110)
ylims!(384-20, 384+20)
xlabel!("DIC (μmol/kg)")
ylabel!("fCO₂ (μatm)")

