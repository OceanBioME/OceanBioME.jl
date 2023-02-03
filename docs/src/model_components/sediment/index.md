# [Sediment models](@id sediments)

> WARNING: we are currently overhauling sediment models but will leave what exists here for others to develop on

Sediments are an important bottom boundary. Generally sediment models receive particulate organic matter sinking from the bottom of models, and then simulate some kind of bacterial diagenesis which may store some fraction or release the matter back into the water in different forms. 

Sediment models can be 3D porous media diffusion models which possibly could be implemented alongside Oceananigans, but so far we have only implemented single layer models. These generally provide auxiliary fields with forcing functions, and boundary conditions for other tracers (which generally depend on the sediment fields).

Currently, we do not have a unified setup procedure, so each model is described separately.