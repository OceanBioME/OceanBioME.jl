# [Positivity Preservation](@id pos-preservation)

It is common in BGC models to behave badly if any tracers go bellow zero, analytically this is fine because they (usually) can not get below zero, and it is unphysical for the concentration of something to be negative. Issues arise when the inaccuracy in numerical integration making some tracer become negative, usually leading to explosions (e.g. ``\exp(-C) \to \inf``), or bounds errors (e.g. ``\log(C)``). Essentially this occurs when the local error in the numerical scheme is sufficiently large that more than the available amount of tracer is consumed.

There exists a set of numerical schemes which overcome this and guarantee positivity (provided a positivity preserving advection scheme and well-behaved diffusion scheme are used) but are complex to implement. Although we may do this in the future we have not yet done so in the meantime have provided some utilities which maintain positivity. The simplest option is to reset any negative tracers to zero, but this causes the model to gain mass. A slightly more complex version is to increase negative tracers to zero and remove the difference from other tracers with in the same conserved system.

We have found this to be a satisfactory solution (when balanced against using much smaller time steps), as it tends to cause only a small and local transient change to the solution.

See model components page for usage.