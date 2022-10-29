# [Positivity Preservation](@id pos-preservation)

It is common in BGC models to behave badly if any tracers go bellow zero, analytically this is fine because they (usually) can not get below zero, and it is unphysical for the concentration of something to be negative. Issues arise when the inaccuracy in numerical integration making some tracer go negative, usually leading to explosions (e.g. $exp(-C)\to\inf$), or bounds errors (e.g. $log(C)$). Essentially this occurs when the local error in the numerical scheme is sufficiently large that more than the available amount of tracer is consumed.

There exists a set of numerical schemes which overcome this and garantee positivity (provided a positivity preserving advection scheme and well-behaved diffusion scheme are used). 