using CSV, DataFrames, Random, CairoMakie
using OceanBioME: CarbonChemistry

using OceanBioME.Models: teos10_density, teos10_polynomial_approximation
using OceanBioME.Models.CarbonChemistryModel: K0, K1, K2, KF



# get the glodap data - downloaded "merged" file from https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0283442/
data = CSV.read("validation/GLODAP.csv", DataFrame)

DIC_name = "G2tco2"
Alk_name = "G2talk"
T_name = "G2temperature"
S_name = "G2salinity"
P_name = "G2pressure"
Si_name = "G2silicate"
PO_name = "G2phosphate"
pH_name = "G2phtsinsitutp"
depth_name = "G2depth"
cruise_name = "G2cruise"
fCO2_name = "G2fco2"
lon_name = "G2longitude"
lat_name = "G2latitude"

DIC_name, Alk_name, T_name, S_name, P_name, Si_name, PO_name, pH_name

# remove low quality data and non QCed
# and drop a weird point from cruise with a massive errors (only 2 points)
#=data = filter(row -> any([row[name*"f"] == 2 for name in [DIC_name, Alk_name, S_name]]) && 
                     !any([row[name] == -9999 for name in [DIC_name, Alk_name, T_name, S_name]]) &&
                     #row[cruise_name] != 5005 &&
                     #row[cruise_name] != 345 &&
                     row[Alk_name] > 1500 &&
                     row["G2phtsqc"] == 1,  data)
=#
# if we don't have Si and PO₄ the error gets a lot bigger but not in the surface values (as much)
data = filter(row -> any([row[name*"f"] == 2 for name in [DIC_name, Alk_name, S_name, pH_name]]) && 
                     !any([row[name] == -9999 for name in [DIC_name, Alk_name, T_name, S_name, P_name, Si_name, PO_name, pH_name]]) &&
                     row["G2expocode"] != "325020210316" &&
                     row["G2expocode"] != "33RO20071215" &&
                     row[cruise_name] != 270 &&
                     #row[Alk_name] > 1500 &&
                     row["G2phtsqc"] == 1,
                     data)

Nrows = size(data, 1)

# setup the model

density_function = teos10_density

# pH_error = 0.00022 ± 0.01223 -> mean is zero
# pCO₂_error = -3 ± 28.8 -> mean is also zero
# with teos polynomial appriximation error increases to -12 ± 31
carbon_chemistry =
    CarbonChemistry(; # Weiss, R.F. (1974, Mar. Chem., 2, 203–215)
                    solubility = K0(-60.2409, 93.4517 * 100, 23.3585, 0.0, 0.023517, -0.023656 / 100, 0.0047036 / 100^2),
                    # Lueke, et. al (2000, Mar. Chem., 70, 105–119;)
                    carbonic_acid = (K1 = K1(constant=61.2172, inverse_T=-3633.86, log_T=-9.67770, S=0.011555, S²=-0.0001152), 
                                     K2 = K2(constant=-25.9290, inverse_T=-471.78, log_T=3.16967, S=0.01781, S²=-0.0001122)),
                    # Perez and Fraga (1987, Mar. Chem., 21, 161–168).
                    fluoride = KF(constant=-9.68, inverse_T=874.0, sqrt_S=0.111, log_S=0.0, log_S_KS=0.0),
                    density_function)
# pH_error = 0.006 ± 0.012 -> mean is zero
# pCO₂_error = 11 ± 28 -> mean is also zero
#carbon_chemistry = CarbonChemistry()

# compare
glodap_pH = zeros(Nrows)
glodap_pCO₂ = zeros(Nrows)

model_pH  = zeros(Nrows)
model_pCO₂ = zeros(Nrows)

pH_error  = zeros(Nrows)
pCO₂_error  = zeros(Nrows)

Threads.@threads for n in 1:Nrows
    T = data[n, T_name]
    S = data[n, S_name]
    P = data[n, P_name] * 0.1 # dbar to bar

    lon = data[n, lon_name]
    lat = data[n, lat_name]

    density = density_function(T, S, P, lon, lat)

    DIC = data[n, DIC_name] * density * 1e-3
    Alk = data[n, Alk_name] * density * 1e-3
    
    silicate = data[n, Si_name] * density * 1e-3
    phosphate = data[n, PO_name] * density * 1e-3

    silicate = ifelse(silicate == -9999, 0, silicate)
    phosphate = ifelse(phosphate == -9999, 0, phosphate)
    P = ifelse(P == -9999, 0, P)

    if data[n, pH_name*"f"] == 2 && data[n, pH_name] != -9999
        model_pH[n] = carbon_chemistry(; DIC, Alk, T, S, P, silicate, phosphate, lon, lat, return_pH = true)
        glodap_pH[n] = data[n, pH_name]
        pH_error[n] = model_pH[n] - glodap_pH[n]
    else
        model_pH[n] = NaN
        glodap_pH[n] = NaN
        pH_error[n] = NaN
    end

    if data[n, fCO2_name*"f"] == 2 && data[n, fCO2_name*"temp"] != -9999
        # fCO2 is reported at no pressure and at a different temp
        model_pCO₂[n] = carbon_chemistry(; DIC, Alk, data[n, fCO2_name*"temp"], S, silicate, phosphate, lon, lat)
        glodap_pCO₂[n] = data[n, fCO2_name]
        pCO₂_error[n] = model_pCO₂[n] - glodap_pCO₂[n]
    else
        model_pCO₂[n] = NaN
        glodap_pCO₂[n] = NaN
        pCO₂_error[n] = NaN
    end

end

function plot_errors(; subset = :,
                       to_plot = pH_error,
                       title = "ΔpH",
                       color_name = S_name,
                       color_marker = data[:, color_name][subset],
                       limit_y = false,
                       ylim = 100)
    fig = Figure();

    ax = Axis(fig[1, 1], title = "DIC")#, aspect = DataAspect(), title = "Local grid")
    ax2 = Axis(fig[1, 2], title = "Alk")
    ax3 = Axis(fig[1, 3], title = "pH")
    ax4 = Axis(fig[2, 1], title = "T")
    ax5 = Axis(fig[2, 2], title = "S")
    ax6 = Axis(fig[2, 3], title = "Depth")

    sc=scatter!(ax, data[:, DIC_name][subset], to_plot[subset], markersize = 3, alpha = 0.5, color = color_marker)
    scatter!(ax2, data[:, Alk_name][subset], to_plot[subset], markersize = 3, alpha = 0.5, color = color_marker)
    scatter!(ax3, data[:, pH_name][subset], to_plot[subset], markersize = 3, alpha = 0.5, color = color_marker)
    scatter!(ax4, data[:, T_name][subset], to_plot[subset], markersize = 3, alpha = 0.5, color = color_marker)
    scatter!(ax5, data[:, S_name][subset], to_plot[subset], markersize = 3, alpha = 0.5, color = color_marker)
    scatter!(ax6, data[:, depth_name][subset], to_plot[subset], markersize = 3, alpha = 0.5, color = color_marker)
    
    if limit_y
        [ylims!(a, -ylim, ylim) for a in (ax, ax2, ax3, ax4, ax5, ax6)]
    end

    Colorbar(fig[1:2, 4], sc)

    Label(fig[0, :], title)

    return fig
end