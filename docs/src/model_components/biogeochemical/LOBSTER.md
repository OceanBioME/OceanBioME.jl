# [The Lodyc-DAMTP Ocean Biogeochemical Simulation Tools for Ecosystem and Resources (LOBSTER) model](@id LOBSTER)

LOBSTER is a medium complexity BGC model with seven core prognostic variables: phytoplankton, zooplankton, small and large detritus, nitrates, ammonia, and dissolved organic matter. LOBSTER was originally proposed by in [Levy2005](@cite) and subsequently added to by [Levy2001](@cite), [Resplandy2009](@cite), [Karleskind2011](@cite), and [Resplandy2012](@cite). 

Additionally, this implementation of LOBSTER optionally models simple carbonate chemistry (`DIC` and `Alk`alinity), `Oxy`gen, and variable redfield ratios for the now dissolved and particulate organic groups (which then allows carbon to be conserved). For details see [StrongWrightInPrep](@cite). These are activated in the model setup, for example:
```
LOBSTER(; grid, carbonates = true)
```