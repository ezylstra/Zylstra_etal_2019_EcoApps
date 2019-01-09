######################################################################################
# R code used to model metapopulation dynamics of lowland leopard frogs 
# in southern Arizona based on 22 years of detection-nondetection data.
# Code adapted from Chandler et al. 2015 (Journal of Applied Ecology)
#
# Reference:
# Zylstra, E.R., D.E. Swann, B.R. Hossack, E. Muths, and R.J. Steidl. Drought-mediated
# extinction of an arid-land amphibian: insights from a spatially-explicit dynamic 
# occupancy model. Ecological Applications.
######################################################################################

#--- Load libraries -----------------------------------------------------------------#

#install.packages(c("reshape2","jagsUI"))
library(reshape2)
library(jagsUI)

#--- Load data ----------------------------------------------------------------------#

#set working directory
setwd("C:/...")

rm(list=ls(all=TRUE))

#Data stored in Open Science Framework: https://doi.org/10.17605/OSF.IO/43PVJ
surveys <- read.table("LeopardFrog_Occupancy_SurveyData.csv", header=TRUE, sep=",", na.strings=c(NA,""))
site.cov <- read.table("LeopardFrog_Occupancy_SiteXSeasonCovariates.csv", header=TRUE, sep=",", na.strings=c(NA,""))
seas.cov <- read.table("LeopardFrog_Occupancy_SeasonCovariates.csv", header=TRUE, sep=",", na.strings=c(NA,""))
neighbors <- read.table("LeopardFrog_Occupancy_Neighbors.csv", header=FALSE, sep=",", na.strings=c(NA,""))
distances <- read.table("LeopardFrog_Occupancy_Distances.csv", header=FALSE, sep=",", na.strings=c(NA,""))

#--- Format and standardize survey-level covariates ---------------------------------#

#Index of surface-water availability
  water.orig <- surveys$water
  water <- (water.orig - mean(water.orig))/sd(water.orig)
  
#Indicator for fall sampling periods
  fall <- surveys$fall
  
#Indicator for surveys led by less-experienced observers
  obs <- surveys$inexp.obs

#--- Format and standardize site-level covariates -----------------------------------#

#Index of surface-water reliability
  rel.orig <- site.cov$reliability
  rel <- (rel.orig - mean(rel.orig))/sd(rel.orig)

#Area of catchment basin
  area.orig <- site.cov$basin.area
  area <- (area.orig - mean(area.orig))/sd(area.orig)
  
#Number of pools
  npools.orig <- site.cov$n.pools
  npools <- (npools.orig - mean(npools.orig))/sd(npools.orig)

#Elevation
  elev.orig <- site.cov$elevation
  elev <- (elev.orig - mean(elev.orig))/sd(elev.orig)

#--- Format and standardize covariates that vary among sites and seasons ------------#

#Palmer Drought Severity Index (PDSI), averaged over six months prior to each winter/summer
  pdsi.orig <- site.cov[,grep('pdsi6mo',names(site.cov))]
  #standardize using the mean and SD among canyon-level observations:
  pdsi6mo.mn <- -1.797
  pdsi6mo.sd <- 2.298
  pdsi <- (pdsi.orig - pdsi6mo.mn)/pdsi6mo.sd

#Palmer Drought Severity Index (PDSI), averaged over prior larval period
  pdsiJA.orig <- site.cov[,grep('pdsiJA',names(site.cov))]
  #standardize using the mean and SD among canyon-level observations:
  pdsiJA.mn <- -1.865
  pdsiJA.sd <- 2.327
  pdsiJA <- (pdsiJA.orig - pdsiJA.mn)/pdsiJA.sd

#Percent of 30-year precipitation normals during current winter/summer season 
  precip.orig <- site.cov[,grep('precip',names(site.cov))]
  #standardize using the mean and SD among canyon-level observations:
  precip.mn <- 88.629
  precip.sd <- 36.444
  precip <- (precip.orig - precip.mn)/precip.sd

#Sediment levels during each winter/summer season 
  sed <- site.cov[,grep('sed',names(site.cov))]

#--- Neighbor index and distance matrices ------------------------------------------------#

#Indices for neighbors of each site
  nsi <- as.matrix(neighbors)

#Squared distances between a site and each of its neighbors
  distSq <- as.matrix(distances^2) 
  
#--- User-defined functions ---------------------------------------------------------#
  
#Function to create a matrix with information about known latent states, z[i,t]
  known.state.dynocc <- function(obsarray){
    state <- apply(obsarray,c(1,3),function(x) ifelse(sum(is.na(x))==length(x),NA,max(x,na.rm=T)))
    state[state==0] <- NA
    return(state)
  }
  
#Function to create a matrix with initial values for unknown latent states, z[i,t]
  inits.state.dynocc <- function(obsarray){
    state <- apply(obsarray,c(1,3),function(x) ifelse(sum(is.na(x))==length(x),NA,max(x,na.rm=T)))
    #Initial value of 1 whenever occupancy state is unknown
    state[state==1] <- 2
    state[is.na(state) | state==0] <- 1
    state[state==2] <- NA
    return(state)
  }  

#Creating arrays of observed survey data:
  y.array <- acast(surveys[,c('y','site','period','rep')], site ~ rep ~ period, value.var="y")

#--- Bundle data for JAGS -----------------------------------------------------------#

#Create matrix of detection covariates: season, water, observer
  water2 <- water*water
  cov.p <- cbind(fall, water, water2, obs)
  ncov.p <- ncol(cov.p)

#create matrix of site covariates for initial occupancy: elev, npools, area 
  elev2 <- elev*elev
  cov.psi <- cbind(elev, elev2, npools, area)
  ncov.psi <- ncol(cov.psi) 

#See LeopardFrog_Occupancy_JAGSdata.pdf for detailed decriptions of each object in jags.data
jags.data <- list(y=as.vector(surveys$y),
                  nSites=max(surveys$site),
                  nSites1=sum(site.cov$watershed=='South'),
                  nSampPeriod=max(surveys$period),
                  nObs=nrow(surveys),
                  ncov.p=ncov.p,
                  cov.p=cov.p,
                  ncov.psi=ncov.psi,
                  cov.psi=cov.psi,
                  site=as.vector(surveys$site),
                  period=as.vector(surveys$period),
                  distSq=distSq,
                  nNeighbors=as.vector(site.cov$n.Neighbors),
                  nsi=nsi,
                  winter=as.vector(seas.cov$winter),
                  pdsi=as.matrix(pdsi),
                  pdsiJA=as.matrix(pdsiJA),
                  precip=as.matrix(precip),
                  sed=as.matrix(sed),
                  rel=as.vector(rel),
                  area=as.vector(area),
                  z=known.state.dynocc(y.array))

#--- Parameters to monitor ----------------------------------------------------------#

params <- c('p0','psi0','eps0','rho0','theta',
            'beta.p0','beta.psi0','beta.eps0','beta.rho0',
            'beta.p','beta.psi',
            'beta.eps1.win','beta.eps2.pdsi','beta.eps3.winpdsi',
            'beta.eps4.sed','beta.eps5.area','beta.eps6.rel',
            'beta.rho1.win','beta.rho2.pdsiJA','beta.rho3.winpdsiJA','beta.rho4.precip',
            'beta.rho5.winprecip','beta.rho6.sed','beta.rho7.area','beta.rho8.rel',
            'PAO')

#--- Initial values -----------------------------------------------------------------#

inits <- function(){list(beta.p0=runif(1,-2,2),beta.psi0=runif(1,-2,2),
                         beta.eps0=runif(1,-2,2),beta.rho0=runif(1,-2,2),
                         theta=runif(1,0,5),
                         beta.p=runif(ncov.p,-2,2),beta.psi=runif(ncov.psi,-2,2),
                         beta.eps1.win=runif(1,-2,2),
                         beta.eps2.pdsi=runif(1,-2,2),
                         beta.eps3.winpdsi=runif(1,-2,2),
                         beta.eps4.sed=runif(1,-2,2),
                         beta.eps5.area=runif(1,-2,2),
                         beta.eps6.rel=runif(1,-2,2),
                         beta.rho1.win=runif(1,-2,2),
                         beta.rho2.pdsiJA=runif(1,-2,2),
                         beta.rho3.winpdsiJA=runif(1,-2,2),
                         beta.rho4.precip=runif(1,-2,2),
                         beta.rho5.winprecip=runif(1,-2,2),
                         beta.rho6.sed=runif(1,-2,2),
                         beta.rho7.area=runif(1,-2,2),
                         beta.rho8.rel=runif(1,-2,2),
                         z=inits.state.dynocc(y.array))}

#--- Create JAGS model file ---------------------------------------------------------#

sink("LeopardFrog_Occupancy_JAGSmodel.txt")
cat("
model{
  
  #Priors
  
   # Detection
   beta.p0 ~ dlogis(0,1)  
   for(i in 1:ncov.p){  
     beta.p[i] ~ dnorm(0,0.1) 
   } #i
   
   # Initial occupancy
   beta.psi0 ~ dlogis(0,1)
   for(i in 1:ncov.psi){
     beta.psi[i] ~ dnorm(0,0.1)
   } #i
    
   # Extinction
   beta.eps0 ~ dlogis(0,1)
   beta.eps1.win ~ dnorm(0,0.1)
   beta.eps2.pdsi ~ dnorm(0,0.1)
   beta.eps3.winpdsi ~ dnorm(0,0.1)
   beta.eps4.sed ~ dnorm(0,0.1)
   beta.eps5.area ~ dnorm(0,0.1)
   beta.eps6.rel ~ dnorm(0,0.1)
  
   # Baseline colonization
   beta.rho0 ~ dlogis(0,1)
   beta.rho1.win ~ dnorm(0,0.1)
   beta.rho2.pdsiJA ~ dnorm(0,0.1)
   beta.rho3.winpdsiJA ~ dnorm(0,0.1)
   beta.rho4.precip ~ dnorm(0,0.1)
   beta.rho5.winprecip ~ dnorm(0,0.1)
   beta.rho6.sed ~ dnorm(0,0.1)
   beta.rho7.area ~ dnorm(0,0.1)
   beta.rho8.rel ~ dnorm(0,0.1)
   
   theta ~ dunif(0,15)  
  
  #State model
    
   for(i in 1:nSites){
      
     logit(psi[i]) <- beta.psi0 + inprod(beta.psi[], cov.psi[i,])
     z[i,1] ~ dbern(psi[i])
      
     for(t in 2:nSampPeriod){
      
       # Extinction probability at site i between sampling period t-1 and t
       logit(eps[i,t-1]) <- beta.eps0 + 
                            beta.eps1.win*winter[t-1] +        
                            beta.eps2.pdsi*pdsi[i,t-1] +
                            beta.eps3.winpdsi*winter[t-1]*pdsi[i,t-1] +                       
                            beta.eps4.sed*sed[i,t-1] + 
                            beta.eps5.area*area[i] + 
                            beta.eps6.rel*rel[i]
        
       # Baseline colonization probability (probability a site is colonized by a frog 
       # from a coincident site) between sampling period t-1 and t:
       logit(rho[i,t-1]) <- beta.rho0 + 
                            beta.rho1.win*winter[t-1] + 
                            beta.rho2.pdsiJA*pdsiJA[i,t-1] +
                            beta.rho3.winpdsiJA*winter[t-1]*pdsiJA[i,t-1] + 
                            beta.rho4.precip*precip[i,t-1] +
                            beta.rho5.winprecip*winter[t-1]*precip[i,t-1] +   
                            beta.rho6.sed*sed[i,t-1] +
                            beta.rho7.area*area[i] + 
                            beta.rho8.rel*rel[i]
        
       for(n in 1:nNeighbors[i]){
        
         # Pairwise colonization probability 
         gammapair[i,n,t-1] <- rho[i,t-1]*exp(-distSq[i,n]/(2*theta^2))*z[nsi[i,n],t-1]
           # We used a Gaussian dispersal kernel, but a different kernel could be  
           # specified by replacing exp(-distSq[i,n]/(2*theta^2)) with another kernel 
           # (e.g., exponential: lambda*exp(-dist[i,n]*lambda)
         gammapairsinv[i,n,t-1] <- 1 - gammapair[i,n,t-1]
  
       } #n
        
       # Colonization probability
       gamma[i,t-1] <- 1 - prod(gammapairsinv[i,1:nNeighbors[i],t-1])
        
       # Probability of being occupied at time t
        
       # Use the following if you don't want to include a rescue effect:
       # Ez[i,t-1] <- gamma[i,t-1]*(1-z[i,t-1])+(1-eps[i,t-1])*z[i,t-1]
        
       # Use the following if you want to include a pseudo-rescue effect:
       Ez[i,t-1] <- gamma[i,t-1]*(1-z[i,t-1])+(1-eps[i,t-1]*(1-gamma[i,t-1]))*z[i,t-1]
       
       z[i,t] ~ dbern(Ez[i,t-1])
      
     } #t
   } #i
  
  #Observation model (data in long form)
  
   for(w in 1:nObs){
      
     logit(p[w]) <- beta.p0 + inprod(beta.p[], cov.p[w,])
     y[w] ~ dbern(z[site[w],period[w]]*p[w])
    
   } #w
  
  #Derived parameters

    logit(p0) <- beta.p0
    logit(psi0) <- beta.psi0
    logit(eps0) <- beta.eps0
    logit(rho0) <- beta.rho0
  
   for(t in 1:nSampPeriod){
  
     # Estimated occupancy in the southern watershed    
     PAO[1,t] <- sum(z[1:nSites1,t])/nSites1 
     # Estimated occupancy in the northern watershed
     PAO[2,t] <- sum(z[(nSites1+1):nSites,t])/(nSites-nSites1) 
  
   } #t

} #model
  ",fill=TRUE)
sink()

#--- Run model in JAGS --------------------------------------------------------------#

ni <- 83000 
nt <- 40 
nc <- 3  
nb <- 3000 
na <- 2000

out <- jags(data=jags.data, inits=inits, parameters.to.save=params,
            model.file="LeopardFrog_Occupancy_JAGSmodel.txt",
            n.chains=nc, n.adapt=na, n.iter=ni, n.burnin=nb, n.thin=nt,
            parallel=T)

print(out)
plot(out,ask=T)


