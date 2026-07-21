# [Library](@id library_api)

Documenting the user interface.

## OceanBioME.jl
```@autodocs
Modules = [OceanBioME]
```

## Biogeochemical Models

### Nutrients-Plankton-Detritus (NPD) framework

The [NPD framework](@ref npd_framework) constructor and its preset models
([`LOBSTER`](@ref LOBSTER), [`NPZD`](@ref NPZD), [`ImplicitBiology`](@ref ImplicitBiology)):

```@autodocs
Modules = [OceanBioME.Models.NutrientsPlanktonDetritusModels]
```

#### Nutrient components

```@autodocs
Modules = [OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels]
```

#### Plankton components

```@autodocs
Modules = [OceanBioME.Models.NutrientsPlanktonDetritusModels.PlanktonModels]
```

#### Detritus components

```@autodocs
Modules = [OceanBioME.Models.NutrientsPlanktonDetritusModels.DetritusModels]
```

#### Inorganic carbon components

```@autodocs
Modules = [OceanBioME.Models.NutrientsPlanktonDetritusModels.InorganicCarbonModels]
```

#### Oxygen component

```@autodocs
Modules = [OceanBioME.Models.NutrientsPlanktonDetritusModels.OxygenModels]
```

### Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES)

```@autodocs
Modules = [OceanBioME.Models.PISCESModel]
```

```@autodocs
Modules = [OceanBioME.Models.PISCESModel.DissolvedOrganicMatter]
```

```@autodocs
Modules = [OceanBioME.Models.PISCESModel.ParticulateOrganicMatter]
```

```@autodocs
Modules = [OceanBioME.Models.PISCESModel.Iron]
```

```@autodocs
Modules = [OceanBioME.Models.PISCESModel.InorganicCarbons]
```

```@autodocs
Modules = [OceanBioME.Models.PISCESModel.Zooplankton]
```

```@autodocs
Modules = [OceanBioME.Models.PISCESModel.Phytoplankton]
```

```@autodocs
Modules = [OceanBioME.Models.PISCESModel.Phosphates]
```

```@autodocs
Modules = [OceanBioME.Models.PISCESModel.Silicates]
```

```@autodocs
Modules = [OceanBioME.Models.PISCESModel.Nitrogen]
```

### Sugar kelp (Saccharina latissima)

```@autodocs
Modules = [OceanBioME.Models.SugarKelpModel]
```

### Carbon Chemistry 

```@autodocs
Modules = [OceanBioME.Models.CarbonChemistryModel]
```

## Light Attenuation Models

```@autodocs
Modules = [OceanBioME.Light]
```

## Sediments

```@autodocs
Modules = [OceanBioME.Models.SedimentModels]
```

## Gas exchange boundary conditions

```@autodocs
Modules = [OceanBioME.Models.GasExchangeModel, OceanBioME.Models.GasExchangeModel.ScaledGasTransferVelocity]
```

## Box Model

```@autodocs
Modules = [OceanBioME.BoxModels]
```

## Particles

```@autodocs
Modules = [OceanBioME.Particles]
```