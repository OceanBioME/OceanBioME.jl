# [Implementing sediment models](@id sed-imp) 

## Notes:
- You do not need a boundary condition on the sinking particles, Oceananigans will let them fall out of the bottom if you do not bring the sinking velocity smoothly to zero at the bottom of the domain so if we let it be open, and then add to the sediment model what is falling out then we get the correct values