# [Sugar kelp (Saccharina latissima) individuals](@id SLatissima)

We have implemented a model of sugar kelp growth within this spatially infinitesimal Lagrangian particles framework originally based on the model of [Broch2012](@cite) and updated by [Broch2013](@cite), [Fossberg2018](@cite), and [Broch2019](@cite). This is the same model we passively forced in [StrongWright2022](@cite). 

The model tracks three variables, the frond area, A (dm²), carbon reserve, C (gC / gSW), and nitrate reserve, N (gN / gSW). The growth depends on the nitrate (and optionally ammonia) availability in the water, the temperature, and light availability. The minimum required coupling is with nitrates so the model can be coupled with an NPZD model, but can optionally uptake ammonia, DIC (CO₂), oxygen, and release dissolved organic matter (from exudation) and large detritus.  

Results could look something like this (from [StrongWright2022](@cite)):
![Example A, N, and C profiles from [StrongWright2022](@cite)](https://www.frontiersin.org/files/Articles/793977/fmars-08-793977-HTML/image_m/fmars-08-793977-g002.jpg)