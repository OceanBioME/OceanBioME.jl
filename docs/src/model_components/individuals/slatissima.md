# [Sugar kelp (Saccharina latissima) individuals](@id SLatissima)

We have implemented a model of sugar kelp growth within this spatially infinitesimal Lagrangian particles framework originally based on the model of [Broch2012](@citet) and updated by [Broch2013](@citet), [Fossberg2018](@citet), and [Broch2019](@citet). This is the same model passively forced by [StrongWright2022](@citet).

The model tracks three variables, the frond area, A (dm²), carbon reserve, C (gC / gSW), and nitrate reserve, N (gN / gSW). The growth depends on the nitrate (and optionally ammonia) availability in the water, the temperature, and light availability. The minimum required coupling is with nitrates so the model can be coupled with an NPZD model, but can optionally uptake ammonia, DIC (CO₂), oxygen, and release dissolved organic matter (from exudation) and large detritus.  

Results could look something like this (from [StrongWright2022](@citet)):
![Example A, N, and C profiles from [StrongWright2022](@citet)](https://www.frontiersin.org/files/Articles/793977/fmars-08-793977-HTML/image_m/fmars-08-793977-g002.jpg)