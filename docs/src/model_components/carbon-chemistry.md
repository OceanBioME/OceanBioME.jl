# [Carbon Chemistry](@id carbon-chemistry)

Our carbon chemistry model follows the best practice described by [dickson2007](@citet), as implemented by e.g. `cbsyst` [branson2023](@citep) and `CO2SYS` [humphreys2022](@citep) in Python.

The carbon chemistry model is primarily used for diagnosing the partial pressure of carbon dioxide in the surface water for [gas exchange with the air](@ref air-sea-gas), but is capable of diagnosing other species such as carbonate concentration, useful in the calculation of calcite dissolution.
The model works by computing the pH from the total dissolved inorganic carbon and total alkalinity (and optionally silicate, phosphate, boron, sulfate, and fluoride concentration), along with temperature and salinity, which can then be used to find the concentration of different species.
We will first describe how to use the model, followed by the underlying chemistry and parameterisations.

## Using the model
To use the carbon chemistry model it is constructed by simply writing:
```@example carbon-chem
using OceanBioME

carbon_chemistry = CarbonChemistry()
```

### Computing ``fCO_2`` and ``pH``
To compute the fugacity of carbon dioxide (``fCO_2``) you call the model with the `DIC`, `Alk`alinity, `T`emperature, and `S`alinity:
```@example carbon-chem
DIC = 2145
Alk = 2448
T = 25.4
S = 36.45
carbon_chemistry(; DIC, Alk, T, S)
```
This is sufficient when computing ``fCO_2`` at the surface, but if we wanted to know ``fCO_2`` at depth where there is higher pressure we can specify the pressure in bar like:
```@example carbon-chem
carbon_chemistry(; DIC, Alk, T, S, P = 100)
```

We may also be interested in the pH so we can request that be outputted:
```@example carbon-chem
carbon_chemistry(; DIC, Alk, T, S, return_pH = true)
```

These function calls assume a constant boron, sulfate, and fluoride ratio relative to the salinity (as described below), but can be specified instead:
```@example carbon-chem
carbon_chemistry(; DIC, Alk, T, S, boron = 0)
```

And the silicate and phosphate concentrations are assumed to be zero but can similarly be specified:
```@example carbon-chem
carbon_chemistry(; DIC, Alk, T, S, silicate = 2, phosphate = 1)
```

The same code can also be used to compute ``fCO_2`` when the pH is already known by passing it in the same way:
```@example carbon-chem
carbon_chemistry(; DIC, pH = 8.1, T, S)
```

Finally, for slightly improved accuracy you may wish to specify that the seawater density calculation, is an intermediary step in the calculations above, is computed using the full TEOS-10 [feistel2008](@citet) standard by setting this as the `density_function` when you construct the carbon chemistry model:
```@example carbon-chem
using OceanBioME.Models: teos10_density
carbon_chemistry = CarbonChemistry(; density_function = teos10_density)
```
But this comes at the cost of significantly increased computational expense and GPU incompatability.

During the conversion from practical to absolute salinity units the location can then also be entered for (very slightly, in this example ~``10^{-4}\mu``atm) improved accuracy:
```@example carbon-chem
carbon_chemistry(; DIC, Alk, T, S, lon = -31.52, lat = 33.75)
```

The default uses the polynomial approximation described in [roquet2015](@citet) as provided by [`SeawaterPolynomials.jl`](https://github.com/CliMA/SeawaterPolynomials.jl/).

### Computing the carbonate concentration
So that this model can be used in calcite dissolution models it can also return the carbonate saturation by calling the function `carbonate_saturation`
```@example carbon-chem
using OceanBioME.Models.CarbonChemistryModel: carbonate_saturation

carbonate_saturation(carbon_chemistry; DIC, Alk, T, S)
```
This function takes all of the same arguments (e.g. `boron`) as `carbon_chemistry` above.

## Chemistry
### pH computation
When carbon dioxide is dissolved in seawater it dissociates into carbonate and bicarbonate species according to the equilibria:

```math
\ce{CO_2(g)} \ce{<=> CO_2(aq)},
```

```math
\ce{CO_2(aq)} + \ce{H_2O} \ce{<=> H_2CO_3(aq)},
```

```math
\ce{CO_2(aq) + H_2O}\ce{<=> H^+ HCO^-_3},
```

```math
\ce{HCO^-_3} \ce{<=> H^+ + CO^{2-}_3},
```

from which we define the total dissolved inorganic carbon (`DIC`) to be ``DIC = [\ce{CO_2(aq)}] + [\ce{HCO^-_3}] + [\ce{CO^{2-}_3}]``. 
The equilibrium constants for these equations are defined as:

```math
K_0=\frac{[\ce{CO_2(aq)}]}{\ce{pCO_2}},
```

```math
K_1=\frac{\ce{[HCO^-_3][H^+]}}{\ce{[CO_2(aq)]}},
```
and
```math
K_2=\frac{\ce{[CO^{2-}_3][H^+]}}{\ce{HCO^-_3}}.
```

These equilibria depend on the total hydrogen ion concentration ``[H^+]``, which depends on the concentration of acids and buffering of bases.

#### Alkalinity: acid-base buffering

Alkalinity is a measure of the buffering capacity of seawater and we can approximate it as:
```math
Alk = [\ce{HCO_3^-}] + 2[\ce{CO_3^{2-}}] + [\ce{B(OH)_4^-}] + [\ce{OH^-}] + [\ce{HPO_4^{2-}}] + 2[\ce{PO_4^{3-}}] + [\ce{SiO(OH)_3^-}] + [\ce{NH_3}] + [\ce{HS^-}] - [\ce{H^+}] - [\ce{HSO_4^-}] - [\ce{HF}] - [\ce{H_3PO_4}] + \text{ minor acids and bases},
```
where the "minor" species are either in sufficiently low concentrations as to be neglected (sometimes the photphate, silicate and sulfate species can also be neglected), and we shall neglect ``[\ce{NH_3}]``.

Each of these acid and base species are also in an equilibrated dissociation state given below.

#### Hydrogen sulfate
Sulfate in the form of hydrogen sulfate dissociates to form sulfate ions in the equilibrium
```math
\ce{HSO_4^-}\ce{<=> H^+ + SO_4^{2-}},
```
with an equilibrium constant given by
```math
K_S=\frac{\ce{[SO_4^{2-}][H^+]}}{\ce{HSO_4^-}}.
```

#### Boric acid
Boron in the form of boric acid equilibrates with water in the reaction
```math
\ce{B(OH)_3+H_2O}\ce{<=> H^+ + B(OH)_4^-},
```
with equilibrium constant
```math
K_B=\frac{\ce{[B(OH)_4^-][H^+]}}{\ce{B(OH)_3}}.
```

#### Hydrogen fluoride
Hydrogen fluoride dissociated in the equilibrium
```math
\ce{HF}\ce{<=> H^+ + F^-},
```
with equilibrium constant
```math
K_F=\frac{\ce{[F^-][H^+]}}{\ce{HF}}.
```

#### Phosphoric acid
Phosphoric acid undergoes a three stage dissociation given by the equilibrium
```math
\ce{H_3PO_4}\ce{<=> H^+ + H_2PO_4^-},
```

```math
\ce{H_2PO_4^-}\ce{<=> H^+ + HPO_4^{2-}},
```

```math
\ce{HPO_4^{2-}}\ce{<=> H^+ + PO_4^{3-}},
```

with equilibrium constants

```math
K_{P1} = \frac{\ce{[H^+][H_2PO_4^-]}}{\ce{H_3PO_4}},
```

```math
K_{P2} = \frac{\ce{[H^+][HPO_4^{2-}]}}{\ce{H_2PO_4^-}},
```

```math
K_{P3} = \frac{\ce{[H^+][PO_4^{3-}]}}{\ce{HPO_4^{2-}}}.
```

#### Silicic acid
Silicic acid dissociates in the equilibrium
```math
\ce{Si(OH)_4}\ce{<=> H^+ + SiO(OH)_3^-},
```

with equilibrium constant
```math
K_{Si} = \frac{\ce{[H^+][SiO(OH)_3^-]}}{\ce{Si(OH)_4}}.
```

#### Water
Finally, water dissociates in the equilibrium
```math
\ce{H_2O}\ce{<=> H^+ + OH^-},
```

with equilibrium constant
```math
K_w = \\ce{[H^+][OH^-]}.
```

#### Alkalinity equilibration
From these rate constants we can rewrite the total alkalinity (given above) in terms of only the rate constants, total hydrogen ion concentration (``\ce{H^+}``), the total dissolved inorganic carbon (``[DIC]``), boron (``B``), phosphate (``P``), silicate (``Si``), Sulfate (``Sulf``), and fluorine (``F``) content of the water, by rearranging the equations above. 
This results in a form of the total alkalinity:
```math
\begin{align}
Alk &\approx \frac{[DIC] K_1 [\ce{H^+}]}{[\ce{H^+}]^2 + K_1[\ce{H^+}] + K_1K_2}\\
    &+ \frac{2[DIC]K_1K_2}{[\ce{H^+}]^2 + K_1[\ce{H^+}] + K_1K_2}\\
    &+ \frac{B}{1+[\ce{H^+}]/K_B}\\
    &+ \frac{K_w}{[\ce{H^+}]}\\
    &+ \frac{PK_{P1}K_{P2}[\ce{H^+}]}{[\ce{H^+}]^3+K_{P1}[\ce{H^+}]^2+K_{P1}K_{P2}[\ce{H^+}] + K_{P1}K_{P2}K_{P3}}\\
    &+ \frac{2PK_{P1}K_{P2}K_{P3}}{[\ce{H^+}]^3+K_{P1}[\ce{H^+}]^2+K_{P1}K_{P2}[\ce{H^+}] + K_{P1}K_{P2}K_{P3}}\\
    &+ \frac{Si}{1+[\ce{H^+}]/K_{Si}}\\
    &- \frac{[\ce{H^+}]}{1+S/K_S}\\
    &- \frac{Sulf}{1+K_S/[\ce{H^+}]/(1+S/K_S)}\\
    &- \frac{F}{1+K_F/[\ce{H^+}]}\\
    &- \frac{[\ce{H^+}]^2}{[\ce{H^+}]^3+K_{P1}[\ce{H^+}]^2+K_{P1}K_{P2}[\ce{H^+}] + K_{P1}K_{P2}K_{P3}}.
\end{align}
```

This gives us a second equation from total alkalinity (which we must already know) with one unknown, ``[\ce{H^+}]``. 
Our model solves for ``[\ce{H^+}]`` using a bisection method to an accuracy of ``10^{-10}``, having approximated the equilibrium constants from parametrisations described below.
We can then determine the pH as, ``pH = -\log_{10}([\ce{H^+}])``.

### Carbon dioxide
Now that we know ``[DIC]``, ``Alk``, and ``[\ce{H^+}]`` we can return to the equation for total dissolved inorganic carbon to find the concentration of aqueous carbon dioxide
```math
[CO_2(aq)] = \frac{[DIC][\ce{H^+}]^2}{[\ce{H^+}]^2 + K_1[\ce{H^+}] + K_1K_2},
```

from which we determine gas phase carbon dioxide concentration as
```math
fCO_2 = \frac{[CO_2(aq)]}{K_0},
```
in atmospheres.

### Carbonate concentration and calcite saturation
Similarly we can also diagnose the calcite concentration
```math
[CO_3^{2-}] = \frac{[DIC]K_1K_2}{[\ce{H^+}]([\ce{H^+}] + K_1)+K_1K_2}.
```

This concentration is important in the dissolution of calcium carbonate which reacts according to the equilibrium
```math
\ce{CaCO_3(aq)}\ce{<=> Ca^+ + CO_3^{2-}},
```
which has an equilibrium constant
```math
K_{SP} = [\ce{Ca^{2+}}][\ce{CO_3^{2-}}].
```

The calcite saturation can then be defined as ``\Omega=\frac{[\ce{CO_3^{2-}}]}{[\ce{CO_{3, saturation}^{2-}}]}`` which can be diagnosed as:
```math
\Omega = \frac{[\ce{Ca^{2+}}][\ce{CO_3^{2-}}]}{K_{SP}}.
```

### Missing pieces
In most cases the chemistry described above requires information about more elements that is usually available.
This means that we must parameterise their concentrations, usually this results in the boron, sulfate, and fluoride concentrations being set as constant ratios of the salinity. 
Usually these ratios are:
```math
\begin{align}
\text{Boron} &= \frac{0.000232}{10.811}\frac{S}{1.80655},\\
\text{Sulfate} &= \frac{0.14}{96.06} \frac{S}{1.80655},\\
\text{Fluoride} &= \frac{0.000067}{18.9984} \frac{S}{1.80655}.\\
\end{align}
```
We use these ratios by default in the model (but they can be changed when calling the models functions described below).
Additionally, silicate and phosphate concentrations are often unavailable (both in observations and models), so are by default set to zero and their alkalinity contribution assumed to be small.

The "ionic strength" must also be parameterised for some equilibrium constants and is usually assumed to be:
```math
Is = \frac{19.924}{1000 - 1.005S},
```
but these parameters are also variable in the model.

## Model parameterisation
The chemical system described above has a large number of equilibrium constants, the constants are typically assumed to depend on the temperature, salinity, and pressure (hydrostatic pressure from depth underwater).
[//]: # (This sentence is not currently actually true, but will be by the time of PR merge)
By default, this model parameterises them based on the "best practice" guidelines of [dickson2007](@citet).
These parameterisations are:
```@example carbon-chem
using OceanBioME.Models.CarbonChemistryModel: K0, K1, K2, KB, KW, KS, KF, KP1, KP2, KP3, KSi, KSP_calcite, KSP_aragonite
K0() # Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
```
```@example carbon-chem
K1() # Millero (1995, Geochim. Cosmochim. Acta, 59, 664)
```
```@example carbon-chem
K2() # Millero (1995, Geochim. Cosmochim. Acta, 59, 664)
```
```@example carbon-chem
KB() # Dickson (1990, Deep-Sea Res., 37, 755–766)
```
```@example carbon-chem
KW() # Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677)
```
```@example carbon-chem
KS() # Dickson (1990, Chem. Thermodyn., 22, 113–127)
```
```@example carbon-chem
KF() # Dickson and Riley (1979, Mar. Chem., 7, 89–99)
```
```@example carbon-chem
KP1() # Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677)
```
```@example carbon-chem
KP2() # Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677)
```
```@example carbon-chem
KP3() # Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677)
```
```@example carbon-chem
KSi() # Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677)
```
```@example carbon-chem
KSP_calcite() # Millero, F. J. (2007, Chemical Reviews, 107(2), 308–341)
```
```@example carbon-chem
KSP_aragonite() # Millero, F. J. (2007, Chemical Reviews, 107(2), 308–341)
```
The pressure corrections from Millero, F. J. (2007, Chemical Reviews, 107(2), 308–341) are:
```@example carbon-chem
K1().pressure_correction
```
```@example carbon-chem
K2().pressure_correction
```
```@example carbon-chem
KB().pressure_correction
```
```@example carbon-chem
KW().pressure_correction
```
```@example carbon-chem
KS().pressure_correction
```
```@example carbon-chem
KF().pressure_correction
```
```@example carbon-chem
KP1().pressure_correction
```
```@example carbon-chem
KP2().pressure_correction
```
```@example carbon-chem
KP3().pressure_correction
```
```@example carbon-chem
KSP_calcite().pressure_correction
```
```@example carbon-chem
KSP_aragonite().pressure_correction
```