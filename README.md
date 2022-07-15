# Lobster.jl (placeholder name)

If this is working the example should be runnable by going into the example folder and running `julia -i --project=.. subpolar.jl`. The first time you run it here you might have todo something for the dependencies listed in Project.toml.

(JSW) I have also changed the root finding for H to use `find_zero` from the Roots package instead of the manual way it was before.

Lobster model based on
- 2005 A four-dimensional mesoscale map of the spring bloom in the northeast Atlantic (POMME experiment): Results of a prognostic model
- 2001 Impact of sub-mesoscale physics on production and subduction of phytoplankton in an oligotrophic regime
- 2012How does dynamical spatial variability impact 234Th-derived estimates of organic export

Notes
- using flux boundary condition
- ignore aggregation term 
- annual cycle 
- add callback to diagnose pCO2
