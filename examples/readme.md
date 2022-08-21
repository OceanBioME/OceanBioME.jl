# Ocean Biogeochemial Modelling Environment

## Example Lists
This document lists several examples included in this repository for users to learn how to use OceanBioME.

### box.jl
This script illustrates how to run OceanBioME as a box model.

### kelp.jl
This script illustrates how to run OceanBioME as a 1D column model using the LOBSTER biogeochemical model using a PAR timeseries from the North Atlantic subpolar gyre. 
This script also shows how to read in data files (in netCDF format) for forcing the model. 
In this case, the script reads in the mixed layer depth from the Mercator Ocean state esimate and uses this to construct an idealized diffusivity profile. 
The script also reads in a timeseries of the photosynthetically available radiation for forcing the LOBSTER model. 

### lobster_example.jl
This script illustrates how to run OceanBioME as a 1D column model using the LOBSTER biogeochemical model. 

In this case, the script 
- uses an idealized mixed layer depth (mld) to construct an idealized diffusivity profile;
- uses an idealized surface photosynthetically available radiation (PAR) timeseries to derive the PAR profile for forcing the LOBSTER model, according to a light absorption model. 
    The PAR profile will be updated every 1 timestep as an auxiliary field; and
- shows how to interpolate a timeseries of data to access to arbitrary time. 

### lobster_example2.jl
This script illustrates how to run OceanBioME as a 1D column model using the LOBSTER biogeochemical model.

In this case, the script 
- introduces two additonal tracer types: carbonates (DIC and ALK) and oxygen (OXY);
- reads in the temperature, salinity, and mixed layer depth from the Mercator Ocean state esimate of the North Atlantic subpolar gyre;

    The temperature and salinity are needed to calculate the air-sea CO2 flux.  
    The mixed layer depth is used to construct an idealized diffusivity profile.
- reads in the surface PAR timeseries data of the same geographic region as above from the NASA's OceanColor Web; 

    Both the data downloaded are saved seperately and a small sample of them are stored in subpolar.nc for ease of use. 
    The surface PAR timeseries is used to derive the PAR profile for forcing the LOBSTER model based on a light absorption model. 
    The PAR profile will be updated every 1 timestep as an auxiliary field. 
- shows how to read *.nc files.

### subpolar_sediment.jl
This script illustrates how to run OceanBioME as a 1D column model using the LOBSTER biogeochemical model using a PAR timeseries from the North Atlantic subpolar gyre. 
This script also shows how to read in data files (in netCDF format) for forcing the model. 
In this case, the script reads in the mixed layer depth from the Mercator Ocean state esimate and uses this to construct an idealized diffusivity profile. 
The script also reads in a timeseries of the photosynthetically available radiation for forcing the LOBSTER model

### howtodownload.mp4
This video shows how to download relevant data files from https://resources.marine.copernicus.eu/products for forcing the model. 
 
https://user-images.githubusercontent.com/14152233/185803720-7f435b93-6664-40be-8fcf-17514ea36a64.mp4





