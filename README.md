# Drought-mediated extinction of an arid-land amphibian: insights from a spatially explicit dynamic occupancy model

### Erin R. Zylstra, Don E. Swann, Blake R. Hossack, Erin Muths, and Robert J. Steidl

### Ecological Applications 29(3):e01859 [10.1002/eap.1859](https://doi.org/10.1002/eap.1859)

### Code and Data DOI: [10.17605/OSF.IO/43PVJ](https://doi.org/10.17605/OSF.IO/43PVJ)

### Please contact the first author for questions about the code or data: Erin R. Zylstra (erinzylstra@gmail.com)
_______________________________________________________________________________________________________________________________________

## Abstract:
Understanding how natural and anthropogenic processes affect population dynamics of species with patchy distributions is critical to predicting their responses to environmental changes. Despite considerable evidence that demographic rates and dispersal patterns vary temporally in response to an array of biotic and abiotic processes, few applications of metapopulation theory have sought to explore factors that explain spatio-temporal variation in extinction or colonization rates. To facilitate exploring these factors, we extended a spatially explicit model of metapopulation dynamics to create a framework that requires only binary presence-absence data, makes few assumptions about the dispersal process, and accounts for imperfect detection. We apply this framework to 22 years of biannual survey data for lowland leopard frogs, *Lithobates yavapaiensis*, an amphibian that inhabits arid stream systems in the southwestern U.S. and northern Mexico. Our results highlight the importance of accounting for factors that govern temporal variation in transition probabilities, as both extinction and colonization rates varied with hydrologic conditions. Specifically, local extinctions were more frequent during drought periods, particularly at sites without reliable surface water. Colonization rates increased when larval and dispersal periods were wetter than normal, which increased the probability that potential emigrants metamorphosed and reached neighboring sites. Extirpation of frogs from all sites in one watershed during a period of severe drought demonstrated the influence of site-level features, as frogs persisted only in areas where most sites held water consistently and where the amount of sediment deposited from high-elevation wildfires was low. Application of our model provided novel insights into how climate-related processes affected the distribution and population dynamics of an arid-land amphibian. The approach we describe has application to a wide array of species that inhabit patchy environments, can improve our understanding of factors that govern metapopulation dynamics, and can inform strategies for conservation of imperiled species.
## Code 
1. [LeopardFrog_Occupancy.R](LeopardFrog_Occupancy.R): R code used to model metapopulation dynamics of lowland leopard frogs (*Lithobates yavapaiensis*) in southern Arizona based on 22 years of detection-nondetection data.  

## Data
1. [LeopardFrog_Occupancy_SurveyData.csv](Data/LeopardFrog_Occupancy_SurveyData.csv): Survey data for lowland leopard frogs and survey-level covariates.  Each row (n = 2354) represents one survey of one pool complex (i.e., site).  The columns are:
    - y: detection/nondetection data
    - site: an index for each site, 1-55
    - period: an index for each sampling period, 1-43
    - fall: an indicator for fall sampling periods
    - water: the proportion of pools with water in each complex during each survey
    - inexp.obs: an indicator for those surveys led by a less-experienced observer
    - rep: an index to differentiate surveys at a given site within a given sampling period, 1-13
2. [LeopardFrog_Occupancy_SeasonCovariates.csv](Data/LeopardFrog_Occupancy_SeasonCovariates.csv): Seasonal covariates that do not vary among sites.  Each row (n = 42) represents a summer or winter season.  The columns are:
    - season: an index for season, 1-42
    - winter: an indicator for winter seasons
3. [LeopardFrog_Occupancy_SiteXSeasonCovariates.csv](Data/LeopardFrog_Occupancy_SiteXSeasonCovariates.csv): Site-level covariates (time-invariant) and seasonal covariates that vary among sites.  Each row (n = 55) represents a site.  The columns are:
    - site: a site index, 1-55
    - watershed: a label indicating whether the site is in the North or South watershed
    - elevation: mean elevation at each site, in meters
    - n.pools: number of pools at each site
    - basin.area: the area of the catchment basin associated with the lowest-elevation pool at each site, in square kilometers
    - reliability: an index of surface-water reliability, ranging from 0 to 1
    - n.Neighbors: the number of “neighbors” for each site (here, the number of sites in the same watershed)
    - 42 columns that begin with “pdsi6mo”: Palmer Drought Severity Index (PDSI) averaged over the six months immediately prior to each winter or summer season
    - 42 columns that begin with “pdsiJA”: PDSI averaged over the preceding larval period
    - 42 columns that begin with “precip”: percent of 30-year precipitation normals during each winter or summer season
    - 42 columns that begin with “sed”: an index of sediment levels (0 = low, 1 = moderate, 2 = high) associated with each winter and summer season
    - *Note: we obtained downscaled estimates of PDSI and precipitation for each canyon and not each site.  Before running models in JAGS, we standardized values using the mean and standard deviation (SD) among canyon-level observations, which differ slightly from the mean and SD of the site-specific values provided in the .csv file. PDSI.6mo: Mean = -1.797, SD = 2.298. PDSI.JA: Mean = -1.865, SD = 2.327. Precipitation: Mean = 88.629, SD = 36.444.*
4. [LeopardFrog_Occupancy_Neighbors.csv](Data/LeopardFrog_Occupancy_Neighbors.csv): A matrix that contains indexes for the neighbors of each site where each row (n = 55) represents a site. The number of columns is equal to the maximum number of neighbors (here, 31).  Because we assumed frogs could not move between the two watersheds, each row contains the indexes for sites in the same watershed.  Because each site in the southern watershed had only 22 neighbors, NAs were entered in columns 23 through 31 for rows associated with southern sites (rows 1-23).
5. [LeopardFrog_Occupancy_Distances.csv](Data/LeopardFrog_Occupancy_Distances.csv): A matrix of distances, in kilometers, between each site and its neighbors.  Each row (n = 55) represents a site.  The number of columns is equal to the maximum number of neighbors (here, 31).  Because each site in the southern watershed had only 22 neighbors, NAs were entered in columns 23 through 31 for rows associated with southern sites (rows 1-23).
