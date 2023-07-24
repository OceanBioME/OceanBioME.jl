# [Pelagic Interactions Scheme for Carbon and Ecosystem Studies volume 2) (PISCES)](@id PISCES)

!!! warning "PISCES is not validated"
    The PISCES implementation is very early in its development and has not yet been validated (we can't
    even promise that it will run properly). If you are going to run experiments with it, we strongly recommend
    validating it against [other implementations](https://sites.nemo-ocean.io/user-guide/install.html) and
    would be very excited if you could share the results and any suggestions you have over an
    [issues](https://github.com/OceanBioME/OceanBioME.jl//issues) or a [discussion](https://github.com/OceanBioME/OceanBioME.jl/discussions).

PISCES [Aumont2015](@cite) is a more complex mixed Mondo-quota BGC model with multiple classes of phyto and zoo plankton with variable composition (C, Fe, and Chl), multiple nutrients (NO₃, NH₄, PO₄, Fe, and Si), dissolved organic matter, and the option for either a two size class or continuum size particulate model. With a total of 24 prognostic variables the complexity of PISCES is significantly higher than other models so far implemented in this package. This diagram shows the high level architecture of the PISCES model:

![Diagram of PISCES architecture](assets/PISCES_architecture.png)

PISCES is considered to be a state of the art BGC model and is widely used. It is part of [CMIP6](https://gmd.copernicus.org/articles/14/5863/2021/) and so it often used in climate modelling but with this implementation we intended to make it easier to use it for smaller scale ecosystem modelling.