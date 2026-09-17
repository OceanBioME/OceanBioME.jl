#####
##### Heterotrophs
#####

# one predator's relationship with ONE PREY CLASS. A prey class is a group of prey — `PREY` lists its
# member names, which may be autotrophs and/or heterotrophs. The class is grazed as a single pool of
# biomass Σ prime(m) and the flux is split back over members in proportion to their biomass.
struct Grazing{PREY, FT}
       maximum_grazing_rate :: FT    # 1/s; MARBL z_umax_0
    grazing_half_saturation :: FT    # mmol C/m³; MARBL z_grz (squared internally when sigmoidal)

    fraction_to_zooplankton :: FT    # ; MARBL graze_zoo
    fraction_to_particulate :: FT    # ; MARBL graze_poc
      fraction_to_dissolved :: FT    # ; MARBL graze_doc (remainder → DIC by residual)
          detritus_fraction :: FT    # ; MARBL f_zoo_detr (food-weighted zoo-loss → POC)

                  sigmoidal :: Bool  # ; MARBL grazing_function (false = Michaelis–Menten)
end

@inline prey_snames(::Grazing{PREY}) where PREY = PREY

# one zooplankton with carbon biomass (tracer <sname>C) and its row of the grazing matrix. The P:C and
# Fe:C quotas are NOT here: they are global scalars shared by every zooplankton and live on
# `ManyPhytoZoo`. That sharing is not cosmetic — once a zooplankter eats another, a per-predator quota
# would leave an unrouted P/Fe surplus and break conservation.
struct Heterotroph{GR, FT}
                         grazing :: GR  # NamedTuple of Grazing, one per prey class this predator grazes

           linear_mortality_rate :: FT  # 1/s; MARBL z_mort_0
        quadratic_mortality_rate :: FT  # 1/s / (mmol C/m³); MARBL z_mort2_0
    loss_threshold_concentration :: FT  # mmol C/m³; MARBL loss_thres

         temperature_sensitivity :: FT  # ; MARBL Q_10
           temperature_reference :: FT  # °C; MARBL Tref
end
