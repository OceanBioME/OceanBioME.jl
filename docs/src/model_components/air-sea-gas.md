# [Air-sea gas exchange](@id air-sea-gas)

We currently have one air-sea gas exchange model implemented. The model, proposed by [Wanninkhof1992](@citet), calculates the solubility of the gas in the water dependent on the temperature and salinity, and calculates the flux depending on the solubility and mixing from the wind.

Currently, the parameters for CO₂ and oxygen are included, but it would be very straightforward to add the parameters given in the original publication for other gases (e.g. inert tracers of other nutrients such as N₂). We also currently have a very simple formulation of the gas transfer velocity which depends on an average wind speed, but it would straightforwardly be expanded to permit variable wind speed (e.g. to simulate enhanced exchange from storms).

It is straightforward to set up a boundary as an air-sea gas exchange:

```@setup gasexchange
using OceanBioME
CO₂_flux = GasExchange(; gas = :CO₂)
using Oceananigans

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

model = NonhydrostaticModel(; grid,
                              biogeochemistry = LOBSTER(; grid, carbonates = true),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ))
```


```@example gasexchange
using OceanBioME

CO₂_flux = GasExchange(; gas = :CO₂)
```

Where the symbol specifies the exchanged gas (currently `:CO₂` or `:O₂`). This can then be passed in the setup of a BGC model, for example:

```@example gasexchange
using Oceananigans

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

model = NonhydrostaticModel(; grid,
                              biogeochemistry = LOBSTER(; grid, carbonates = true),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ))
```

If the temperature and salinity are not included in the model they can be passed as functions
(or even anonymous functions):

```@example gasexchange
T_function(x, y, z, t) = 12.0

CO₂_flux = GasExchange(; gas = :CO₂, temperature = T_function, salinity = (args...) -> 35)
```

## Model equations

### Gas transfer

The gas flux is given by:

$F = k(C_w - \alpha C_a),$

where $C_w$ is the concentration in the water, $C_a$ the concentration in the air, $\alpha$ the Oswald solubility coefficient, and $k$ the gas transfer velocity. For carbon dioxide the flux is modified to:

$F = k\beta\rho_o(pCO_{2w} - pCO_{2a}),$

where $pCO_{2w}$ and $pCO_{2a}$ are the partial pressure of carbon dioxide in the water and air, $\beta$ is the Bunsen Solubility Coefficient, and $\rho_o$ is the density of the water.

The gas transfer velocity is parameterised by the wind speed and Schmidt number, which in turn is parameterised by the temperature and salinity. The gas transfer velocity is given by:

$k = 1.08\times10^{-6}u^2\left(\frac{Sc}{660}\right)^{-1/2},$

where $u$ is the winds speed 10m above the surface, and Sc is the Schmidt number parameterised as:

$Sc = A - BT + CT^2 - DT^3,$

where $T$ is temperature in Kelvin and the other parameters are dependent on the gas type and given in [Parameters](@ref parameters).

The solubilities are given by:

$\alpha = 0.00367 T \exp{A_1 + 100\frac{A_2}{T} + A_3 \ln{\frac{T}{100}} + S\left(B_1 + \frac{B_2}{T} + \frac{B_3}{T^2}\right)},$

and

$\beta = \exp{A_1 + 100\frac{A_2}{T} + A_3 \ln{\frac{T}{100}} + S\left(B_1 + \frac{B_2}{T} + \frac{B_3}{T^2}\right)},$

where $S$ is salinity in practical units and the other default parameters are given in [Parameters](@ref parameters).

### Partial pressure of carbon dioxide

We currently do not have the full OCMIP partial pressure formulation, instead we follow the simplified formulation where the partial pressure of CO$_2$ ($\mu$atm) in gas is found from:

$pCO_{2a} = f_{CO_2}P_a,$

where $f_{CO_2}$ is the air fraction (ppmv), and $P_a$ is the atmospheric pressure (atm). The ocean partial pressure ($\mu$atm) is derived from the dissolved inorganic carbon content, and alkalinity, as described by [tan1998](@citet), from the equilibrium of the following chemical system:

$\ce{CO_2(g)} \ce{<=> CO_2(aq)},\ K_0=\frac{[\ce{O_2(aq)}]}{\ce{pCO_2}},$

$\ce{CO_2(aq) + H_2O}\ce{<=> H^+ + HCO^-_3},\ K_1=\frac{\ce{[HCO^-_3][H^+]}}{\ce{[CO_2(aq)]}},$

$\ce{HCO^-_3} \ce{<=> H^+ + CO^{2-}_3},\ K_2=\frac{\ce{[CO^{2-}_3][H^+]}}{\ce{HCO^-_3}}.$

In this system DIC and Alk are defined to be:

$\ce{DIC} = \ce{[CO_2(aq)] + [HCO^-_3] + [CO_3^{2-}]},$

$\ce{Alk} = \ce{[HCO^-_3] + 2[CO^{2-}_3] + [B(OH)^-_4] - [H^+] \pm [minor species]}.$

To solve this we must find the Boric acid dissociation from:

$\ce{B(OH)_3 <=>H^+ + B(OH)^-_4},\ K_B = \frac{\ce{[B(OH)^-_4][H^+]}}{\ce{[B(OH)_3]}}.$

Finally, taking the DIC and Alkalinity in micro equivalents (i.e. scaled by $10^{-3}/\rho_o$ from mmol C/m$^3$) denoted by $\bar{DIC}$ and $\bar{Alk}$, the carbonate alkalinity is given by $Alk_C = \bar{Alk} - \frac{K_BB}{K_B + H}$, and we define the boron content, $B$, to be $1.179\times10^{-5}S$ mol/kg.

The system can be rearranged to give:

$\ce{DIC} = \left(1 + \frac{K_1}{[H^+]} + \frac{K_1K_2}{[H^+]^2}\right)[CO_2(aq)],$
    
$\ce{Alk_C} = \left(\frac{K_1}{[H^+]} + \frac{2K_1K_2}{[H^+]^2}\right)[CO_2(aq)].$

This results in a quadratic equation for $H$:

$0 = Alk_CH^2 + K_1(Alk_C - \bar{DIC})H + K_1K_2(Alk_C - \bar{DIC}).$

The equilibrium constants ($K_0$, $K_1$, and $K_2$) are given by:

$K_0  = \exp{A_0 + \frac{B_0}{T} + C_0\log(\frac{T}{100}) + S  (D_0 - E_0T + F_0\left(\frac{T}{100}\right)^2)},$

$K_1 = \exp{A_1 - \frac{B_1}{T} - C_1  \log(T) - (D_1 + \frac{E_1}{T})\sqrt{S} + F_1S - G_1S^{1.5}},$

$K_2 = \exp{A_2 - \frac{B_2}{T} - C_2 \log(T) - (D_2 + \frac{E_2}{T})\sqrt{S} + F_2S - G_1S^{1.5}},$

$K_b = \exp{\frac{A_b - B_b\sqrt{S} - C_bS + D_b  S^{1.5} - E_bS^2}{T} + F_b + G_B  \sqrt{S} + H_b  S - (I_b + J_b \sqrt{S} + S)  \log(T) + L_b  \sqrt{S}  T},$

where the default parameters are given in [Parameters](@ref parameters).

We solve these equations iteratively from an initial guess of $pH=8$ to find $H$, from which the partial pressure of $CO_2$ is calculated as:

$pCO_2 = 10^6\frac{Alk_CH^2}{K_0(K_1H + 2 K_1K_2)}$