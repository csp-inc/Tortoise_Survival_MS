## ---------------------------
##
## Tortoise Survival Functions
##
## This script provides functions used in a seperate script ("Tortoise_Survival_Analysis.R")
## to replicate the analyses of Hromada et. al.,
## "An integrated model improves inferences about survival in the Mojave desert tortoise"
## 
## Author(s): Hromada, S.
##
## Date Created: 2025-08-15
## Date Last Modified: 2025-08-15
## 
## Email: stevehromada@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

rm(list=ls())
options(scipen = 6, digits = 4)

## ---------------------------

## load up packages 

#require()

## ---------------------------

## load up our functions into memory

# source(".R") 

## ---------------------------
# Function that writes out JAGS code for the known-fate survival model
WriteTelFileForJAGS <- function(){
  BUGSfilename <- "./DT_Telem.txt"
  cat("
    
    model{
    
      ###########
      # PROBABILITY OF SURVIVING
      ###########
   
      phi0 ~ dunif(0, 1)  #dbeta(15.3, 1.7)   # mean/intercept survival
      logit.phi0 <- log(phi0 / (1-phi0))
      
      phi.male.eff ~ dnorm(0, 1)          # effect of sex on mean survival (slightly regularized prior)
      phi.juv.eff ~ dnorm(0, 1)           # effect of age on survival
      phi.sub.eff ~ dnorm(0,1)            # effect of age on survival
      phi.precip.eff ~ dnorm(0, 0.01)    # effect of precipitation
      phi.precip.x.sex.eff ~ dnorm(0, 0.01) # interaction effect of precip x sex
      
      # Random site (not pixel) x year effect on survival
      phi.sitexyear.eff.sd ~ dgamma(0.1, 0.1)     
      phi.sitexyear.eff.prec <- pow(phi.sitexyear.eff.sd, -2)
      for(year in 1:nyears){
        for(site in 1:nsites){   #this is a psuedo spatial autocorrelation 
          phi.siteyear.eff[year,site] ~ dnorm(0, phi.sitexyear.eff.prec)
        }
      }
        
      #############
      # DEAL WITH MISSING DATA
      #############
      
      prob.male0 <- 0.5  
      prob.male0.logit <- log(prob.male0 / (1-prob.male0))
      
      prob.male.site.eff.prec ~ dgamma(0.1, 0.1)
      prob.male.site.eff.sd <- pow(prob.male.site.eff.prec, -2)
      for(site in 1:nsites){
        prob.male.site.eff[site] ~ dnorm(0, prob.male.site.eff.prec)
      }
      
      for(site in 1:nsites){
        for(ind in 1:ninds_telem_surv[site]){
          logit(prob.male[site,ind]) <- prob.male0.logit + prob.male.site.eff[site]    # look for variation in sex ratio among sites
          ismale_telem[site,ind] ~ dbern(prob.male[site,ind])
        }
      }
      
      for(site in 1:nsites){
        prob.age[site,1:3] ~ ddirch(c(1,1,10))   #each site has its own probabilities of starting off in each age class
      }
      # Transition from juvenile to adult
    	trans.to.sa ~ dunif(0, 1)
    	trans.to.ad ~ dunif(0, 1)
      for(site in 1:nsites){
      	for(ind in 1:ninds_telem_surv[site]){ #loop through individuals
      	  age.class_telem[site,ind,2] ~ dcat(prob.age[site,1:3]) #since year 1 is not used in modeling, start at 2
      	  is.juv_telem[site,ind,2] <- ifelse(age.class_telem[site,ind,2]==1,1,0)
      		is.sub_telem[site,ind,2] <- ifelse(age.class_telem[site,ind,2]==2,1,0)
      		for(t in 3:(ntimes_telem_surv[site,ind] + 1)){
      		  prob.trans[site,ind,t,1] <- ifelse(age.class_telem[site,ind,t-1]==1,1-trans.to.sa,0)
      		  prob.trans[site,ind,t,2] <- ifelse(age.class_telem[site,ind,t-1]==1,trans.to.sa,ifelse(age.class_telem[site,ind,t-1]==2,1-trans.to.ad,0))
      		  prob.trans[site,ind,t,3] <- ifelse(age.class_telem[site,ind,t-1]==2,trans.to.ad,ifelse(age.class_telem[site,ind,t-1]==3,1,0))
            age.class_telem[site,ind,t] ~ dcat(prob.trans[site,ind,t,1:3])
      			is.juv_telem[site,ind,t] <- ifelse(age.class_telem[site,ind,t]==1,1,0)
      			is.sub_telem[site,ind,t] <- ifelse(age.class_telem[site,ind,t]==2,1,0)
      		}
      	}
      }

    

      ##########
      # TELEMETRY BASED SURVIVAL
      ##########

      ### SURVIVAL
      for(site in 1:nsites){
        for(ind in 1:ninds_telem_surv[site]){
          for(t in 2:(ntimes_telem_surv[site,ind] + 1)){
            logit(psurviv_telem[site,ind,t]) <- logit.phi0 +
              phi.male.eff * ismale_telem[site,ind] +
              phi.juv.eff * is.juv_telem[site,ind,t] +
              phi.sub.eff * is.sub_telem[site,ind,t] +
              phi.precip.eff * precip[site,ind,t] +
              phi.precip.x.sex.eff * precip[site,ind,t] * ismale_telem[site,ind] +
              phi.siteyear.eff[year_ind_tel[site,ind,t],site_telem_surv[site,ind]]
          }
        }
      }
      
      ##############
      ### LIKELIHOOD
      ##############
      for(site in 1:nsites){
        for(ind in 1:ninds_telem_surv[site]){
          alive_telem[site,ind,1] ~ dbern(1) #first one is always alive
          for(t in 2:(ntimes_telem_surv[site,ind] + 1)){
            palive_telem[site,ind,t] <- alive_telem[site,ind,t-1] * pow(psurviv_telem[site,ind,t], 1)
            alive_telem[site,ind,t] ~ dbern(palive_telem[site,ind,t])
          }
        }
      }
          
      ##############
      ### DERIVED PARAMETERS
      ##############
      logit(phi.female) <- logit.phi0   
      logit(phi.male) <- logit.phi0 + phi.male.eff
      logit(phi.juv.female) <- logit.phi0 + phi.juv.eff  
      logit(phi.juv.male) <- logit.phi0 + phi.juv.eff + phi.male.eff  
      logit(phi.sub.male) <- logit.phi0 + phi.sub.eff + phi.male.eff
      logit(phi.sub.female) <- logit.phi0 + phi.sub.eff 
    
      for(cov in 1:length(precip_predict)){
        logit(phi.precip.male[cov]) <- logit.phi0 + phi.precip.eff * precip_predict[cov] + phi.male.eff + phi.precip.x.sex.eff * precip_predict[cov] 
        logit(phi.precip.female[cov]) <- logit.phi0 + phi.precip.eff * precip_predict[cov]
      }
    
      for(year in 1:nyears){
        logit(phi.yr.female[year]) <- logit.phi0 + mean(phi.siteyear.eff[year, 1:nsites])
        for (site in 1:nsites){
          # Female site x year survival from telemetry data (?)
          logit(phi.female.sitexyear[year,site]) <- logit.phi0 + phi.siteyear.eff[year,site]
          # Male site x year survival
          logit(phi.male.sitexyear[year,site]) <- logit.phi0 + phi.siteyear.eff[year,site] + phi.male.eff
        }
      }

    }",file=BUGSfilename)
  return(BUGSfilename)
}

# Function that writes out JAGS code for the mark-recapture survival model
WriteMRCFileForJAGS = function(){
  BUGSfilename <- "./DT_CJS.txt"
  cat("
    
    model{
    
      #############
      # PRIORS
      #############
      
      ###########
      # PROBABILITY OF CAPTURE WITHIN SURVEYED AREA
      ###########
      
      p0 ~ dunif(0, 1)                  # mean detection prob within surveyed areas
      logit.p0 <- log(p0 / (1-p0))      # convert to logit (intercept for logit linear model with standardized covariates)
      p.male.eff ~ dnorm(0, 1)          # effect of sex on mean capture probability  (slightly regularized prior)
      p.juv.eff ~ dnorm(0, 1)           # effect of age on mean capture probability  (slightly regularized prior)
      p.sub.eff ~ dnorm(0, 1)           # effect of age on mean capture probability  (slightly regularized prior)
      p.effort.eff ~ dnorm(1, 1)        # effect of increased effort on mean capture probability
      p.fall.eff ~ dnorm(1, 1)        # effect of fall survey timing on mean capture probability
    
      # Variation in capture probability among sites  # KTS made consistent with year effect
      p.site.eff.prec <- pow(p.site.eff.sd, -2)
      p.site.eff.sd ~ dunif(0, 2)
      for(site in 1:nsites_mrc){
        p.site.eff[site] ~ dnorm(0, p.site.eff.prec)
      }

      # Variation in capture probability among years   # KTS made consistent with site effect
      p.year.eff.prec <- pow(p.year.eff.sd, -2)
      p.year.eff.sd ~ dunif(0, 2)
      for(year in 1:nyears){
        p.year.eff[year] ~ dnorm(0, p.year.eff.prec)
      }

      #### Capture probability model (logit-linear)
      for(site in 1:nsites_mrc){
        for(ind in 1:ninds[site]){ # loop through individuals
          for(year in firsts[site,ind]:nyears){
            for(session in 1:nsessions[site,year]){
              logit(thisp[site,ind,year,session]) <- logit.p0 +
                p.male.eff * is.male[site,ind] +
                p.juv.eff * is.juv[site,ind,year] + 
                p.sub.eff * is.sub[site,ind,year] + 
                p.fall.eff * is.fall[site,year,session] +
                p.effort.eff * effort[site,year,session] +
                p.site.eff[site] + p.year.eff[year]
            }
          }
        }
      }
      
      #############
      # DEAL WITH MISSING DATA
      #############
      
      # Male probability
      prob.male0 <- 0.5  
      prob.male0.logit <- log(prob.male0 / (1-prob.male0))
      prob.male.site.eff.prec ~ dgamma(0.1, 0.1)
      prob.male.site.eff.sd <- pow(prob.male.site.eff.prec, -2)
      for(site in 1:nsites_mrc){
        prob.male.site.eff[site] ~ dnorm(0, prob.male.site.eff.prec)
      }
      for(site in 1:nsites_mrc){
        # Variation in sex ratio among sites
        logit(prob.male[site]) <- prob.male0.logit + prob.male.site.eff[site]  # KTS moved outside ind loop
        for(ind in 1:ninds[site]){
          is.male[site,ind] ~ dbern(prob.male[site])
        }
      }
      
      ######################
      # Age class (KTS: changed to 3 classes- juvenile=1 , subadult=2, and adult=3)
      ######################
      
      for(site in 1:nsites_mrc){
        prob.age[site,1:3] ~ ddirch(c(1,1,10))   #each site has its own probabilities of starting off in each age class
      }
      # Transition from juvenile to adult
    	trans.to.sa ~ dunif(0, 1)
    	trans.to.ad ~ dunif(0, 1)
      for(site in 1:nsites_mrc){
      	for(ind in 1:ninds[site]){ #loop through individuals
      	  age.class[site,ind,firsts[site,ind]] ~ dcat(prob.age[site,1:3])
      	  is.juv[site,ind,firsts[site,ind]] <- ifelse(age.class[site,ind,firsts[site,ind]]==1,1,0)
      		is.sub[site,ind,firsts[site,ind]] <- ifelse(age.class[site,ind,firsts[site,ind]]==2,1,0)
      		for(year in (firsts[site,ind]+1):nyears){
      		  prob.trans[site,ind,year,1] <- ifelse(age.class[site,ind,year-1]==1,1-trans.to.sa,0)
      		  prob.trans[site,ind,year,2] <- ifelse(age.class[site,ind,year-1]==1,trans.to.sa,ifelse(age.class[site,ind,year-1]==2,1-trans.to.ad,0))
      		  prob.trans[site,ind,year,3] <- ifelse(age.class[site,ind,year-1]==2,trans.to.ad,ifelse(age.class[site,ind,year-1]==3,1,0))
            age.class[site,ind,year] ~ dcat(prob.trans[site,ind,year,1:3])
      			is.juv[site,ind,year] <- ifelse(age.class[site,ind,year]==1,1,0)
      			is.sub[site,ind,year] <- ifelse(age.class[site,ind,year]==2,1,0)
      		}
      	}
      }

	    ###########
      # PROBABILITY OF SURVIVING
      ###########
   
      phi0 ~ dunif(0, 1)     # Prior for mean survival
      logit.phi0 <- log(phi0 / (1-phi0))    # Logit transformation of mean survival
      phi.male.eff ~ dnorm(0, 1)    # Fixed effect of male on mean survival (slightly regularized prior)
      phi.juv.eff ~ dnorm(0, 1)   # Fixed effect of juvenile on mean survival
      phi.sub.eff ~ dnorm(0, 1)   # Fixed effect of subadult on mean survival  # KTS added effect of subadult status on survival
      phi.precip.eff ~ dnorm(0, 1)    # Fixed effect of timelag precip
      phi.precip.x.sex.eff ~ dnorm(0, 1)   # Precip x sex interaction effect
    
      # Random site x year effect on survival
      phi.siteyear.eff.mrc.sd ~ dunif(0,2)        # Variation in survival among sites  # KTS changed to match other random effects
      phi.site.eff.prec <- pow(phi.siteyear.eff.mrc.sd, -2)     # KTS changed to match other raneffs
      for(site in 1:nsites_mrc){
        for(year in 1:nyears){
          phi.siteyear.eff[year,site] ~ dnorm(0, phi.site.eff.prec)
        }
      }

      #### Survival model (logit-linear)
      for(site in 1:nsites_mrc){
        for(ind in 1:ninds[site]){  # loop through individuals
          for(year in firsts[site,ind]:(nyears-1)){
            logit(thisphi[site,ind,year]) <- logit.phi0 +
              phi.male.eff * is.male[site,ind] +
              phi.juv.eff * is.juv[site,ind,year] +
              phi.sub.eff * is.sub[site,ind,year] +
              phi.precip.eff * precip_mrc[site,year] +
              phi.precip.x.sex.eff * precip_mrc[site,year] * is.male[site,ind] +
              phi.siteyear.eff[year,site]
          }
        }
      }
      
      ###########
      # SPATIAL COMPONENT
      ###########
    
      ###########
      # SIZE OF ACTIVITY AREA (log scale)
      ###########
      
      sigma0 ~ dunif(10, 700)    # proxy for size of activity area (m)
      logsigma0 <- log(sigma0)
      male.hr.eff ~ dnorm(0, 1)   # was unif -1 to 1
      sub.hr.eff ~ dnorm(0, 1)
      juv.hr.eff ~ dnorm(0, 1)
      fall.hr.eff ~ dnorm(0, 1)
      
      for(site in 1:nsites_mrc){
        for(ind in 1:ninds[site]){
          # Specify activity center for each individual (random from within larger area around plot)
          hrx[site,ind] ~ dunif(xlim[site,1], xlim[site,2])
          hry[site,ind] ~ dunif(ylim[site,1], ylim[site,2])
          
          for(year in firsts[site,ind]:nyears){
            log(sigma[site,ind,year]) <- logsigma0 +
              male.hr.eff * is.male[site,ind] +
              juv.hr.eff * is.juv[site,ind,year] +
              sub.hr.eff * is.sub[site,ind,year] +
              fall.hr.eff * is.fall[site,year,1] #leaving this as 1 for now, since it doesnt change across a siteyear
            tau[site,ind,year] <- pow(sigma[site,ind,year], -2)
          }
        }
      }
      
      ###########
      # SAMPLE INDIVIDUAL LOCATIONS
      ###########
      
      for(site in 1:nsites_mrc){  
        for(ind in 1:ninds[site]){
          for(year in firsts[site,ind]:nyears){
            for(session in 1:nsessions[site,year]){
              locationx[site,ind,year,session] ~ dnorm(hrx[site,ind], tau[site,ind,year])    # data node (can be missing)
              locationy[site,ind,year,session] ~ dnorm(hry[site,ind], tau[site,ind,year])		 # data node (can be missing)
              inplotx[site,ind,year,session] <- step(locationx[site,ind,year,session]-plot.xlim[site,1]) * step(plot.xlim[site,2]-locationx[site,ind,year,session])   
              inploty[site,ind,year,session] <- step(locationy[site,ind,year,session]-plot.ylim[site,1]) * step(plot.ylim[site,2]-locationy[site,ind,year,session])
              inplot[site,ind,year,session] <- inplotx[site,ind,year,session] * inploty[site,ind,year,session]
            }
          }
        }
      }
      
      #############
      # LIKELIHOOD
      #############

      for(site in 1:nsites_mrc){ 
        for(ind in 1:ninds[site]){
          alive[site,ind,firsts[site,ind]] ~ dbern(1)   # always alive on first year of capture
          for(session in (firsts2[site,ind]+1):nsessions[site,firsts[site,ind]]){
            p[site,ind,firsts[site,ind],session] <- thisp[site,ind,firsts[site,ind],session] * alive[site,ind,firsts[site,ind]] * inplot[site,ind,firsts[site,ind],session]    
            y[site,ind,firsts[site,ind],session] ~ dbern(p[site,ind,firsts[site,ind],session])
          }
          for(year in (firsts[site,ind]+1):nyears){
            expected.alive[site,ind,year] <- alive[site,ind,year-1] * thisphi[site,ind,year-1]
            alive[site,ind,year] ~ dbern(expected.alive[site,ind,year])
            for(session in 1:nsessions[site,year]){
  		        # Probability of catching an individual at a given point given it is in the population, in the study area, and real (after first PP)
              p[site,ind,year,session] <- thisp[site,ind,year,session] * alive[site,ind,year] * inplot[site,ind,year,session] #+ # if alive 
                # p0.dead * (1-alive[site,ind,year]) # KTS removed dead recovery component (will have to change data so that dead individuals can not be observed...)
              y[site,ind,year,session] ~ dbern(p[site,ind,year,session])                          
            }
          }    
        }
      }
      
      #############
      # DERIVED PARAMETERS
      #############
      logit(phi.female) <- logit.phi0                        # Female survival
      logit(phi.male) <- logit.phi0 + phi.male.eff           # Male survival
      logit(phi.juv.female) <- logit.phi0 + phi.juv.eff             # Juvenile female survival
      logit(phi.juv.male) <- logit.phi0 + phi.juv.eff + phi.male.eff  # Juvenile male survival
      logit(phi.sub.female) <- logit.phi0 + phi.sub.eff # subadult  female survival
      logit(phi.sub.male) <- logit.phi0 + phi.sub.eff + phi.male.eff # subadult male survival
    
      for(cov in 1:50){ #50 = as long as the input prediction dataset    
          logit(phi.precip.female[cov]) <- logit.phi0 + phi.precip.eff * precip_predict[cov]
          logit(phi.precip.male[cov]) <- logit.phi0 + phi.precip.eff * precip_predict[cov] + phi.male.eff + phi.precip.x.sex.eff * precip_predict[cov]
      }
    
      for(year in 2:nyears){   
        logit(phi.yr.female[year]) <- logit.phi0 + mean(phi.siteyear.eff[year, 1:nsites_mrc])
       
        for(site in 1:nsites_mrc){
          logit(phi.female.sitexyear[year,site]) <- logit.phi0 + phi.siteyear.eff[year,site] + phi.precip.eff * precip_mrc[site,year] #+ phi.pdsi.eff * pdsi_mrc[site,year-1]
          logit(phi.male.sitexyear[year,site]) <- logit.phi0 + phi.siteyear.eff[year,site] + phi.male.eff + phi.precip.eff * precip_mrc[site,year] + phi.precip.x.sex.eff * precip_mrc[site,year]  #+ phi.pdsi.eff * pdsi_mrc[site,year - 1]
        }
      }

    }   
    
    ",file=BUGSfilename)
  
  return(BUGSfilename)
}


# Function that writes out JAGS code for the integrated survival model
WriteIMFileForJags = function(){
  ##START THE FILE FOR JAGS! build function into bugs model
  
  BUGSfilename <- "./DT_IM.txt"
  
  cat("
    
    model{
    
      #############
      # PRIORS
      #############
      
      ###########
      # PROBABILITY OF CAPTURE WITHIN SURVEYED AREA
      ###########
      
      p0 ~ dunif(0, 1)                  # mean detection prob within surveyed areas
      logit.p0 <- log(p0 / (1-p0))       # convert to logit (intercept for logit linear model with standardized covariates)
      p.male.eff ~ dnorm(0, 1)          # effect of sex on mean capture probability (slightly regularized prior)
      p.juv.eff ~ dnorm(0, 1)           # effect of age on mean capture probability (slightly regularized prior)
      p.sub.eff ~ dnorm(0, 1)           # effect of age on mean capture probability  (slightly regularized prior)
      p.effort.eff ~ dnorm(0, 1)        # effect of increased effort on mean capture probability
      p.fall.eff ~ dnorm(0, 1)        # effect of fall survey timing on mean capture probability
    
      # Variation in capture probability among sites
      p.site.eff.prec <- pow(p.site.eff.sd, -2)
      p.site.eff.sd ~ dunif(0, 2)
      for(site in 1:nsites_mrc){
        p.site.eff[site] ~ dnorm(0, p.site.eff.prec)
      }

      # Variation in capture probability among years
      p.year.eff.prec <- pow(p.year.eff.sd, -2)
      p.year.eff.sd ~ dunif(0, 2) #dgamma(0.1,0.1)
      for(year in 1:nyears){
        p.year.eff[year] ~ dnorm(0, p.year.eff.prec)
      }

      #### Capture probability model (logit-linear)
      for(site in 1:nsites_mrc){
        for(ind in 1:ninds[site]){
          for(year in firsts[site,ind]:nyears){
            for(session in 1:nsessions[site,year]){
              logit(thisp[site,ind,year,session]) <- logit.p0 +
                p.male.eff * is.male[site,ind] +
                p.juv.eff * is.juv[site,ind,year] +
                p.sub.eff * is.sub[site,ind,year] + 
                p.effort.eff * effort[site,year,session] +
                p.fall.eff * is.fall[site,year,session] +
                p.site.eff[site] +
                p.year.eff[year]
            }
          }
        }
      }
      
      #############
      # DEAL WITH MISSING DATA
      #############
      
      # Male probability for mark-recapture data
      prob.male0 <- 0.5  
      prob.male0.logit <- log(prob.male0 / (1-prob.male0))
      prob.male.site.eff.prec ~ dgamma(0.1, 0.1)
      prob.male.site.eff.sd <- pow(prob.male.site.eff.prec, -2)
      for(site in 1:nsites_mrc){
        prob.male.site.eff[site] ~ dnorm(0, prob.male.site.eff.prec)
      }
      for(site in 1:nsites_mrc){
        # Variation in sex ratio among sites
        logit(prob.male[site]) <- prob.male0.logit + prob.male.site.eff[site]  # KTS moved outside ind loop
        for(ind in 1:ninds[site]){
          is.male[site,ind] ~ dbern(prob.male[site])
        }
      }
      
      # Male probability for telemetry data
      prob.male.site.eff.tel.prec ~ dgamma(0.1, 0.1)
      prob.male.site.eff.tel.sd <- pow(prob.male.site.eff.tel.prec, -2)
      for(site in 1:nsites_tel){
        prob.male.site.eff.tel[site] ~ dnorm(0, prob.male.site.eff.tel.prec)
      }
      for(site in 1:nsites_tel){
        # Variation in sex ratio among sites
        logit(prob.male.tel[site]) <- prob.male0.logit + prob.male.site.eff.tel[site]
        for(ind in 1:ninds_telem_surv[site]){
         ismale_telem[site,ind] ~ dbern(prob.male.tel[site])
        }
      }
      

     ######################
      # Age class (KTS: changed to 3 classes- juvenile=1 , subadult=2, and adult=3)
      ######################
      
      for(site in 1:nsites_mrc){
        prob.age[site,1:3] ~ ddirch(c(1,1,10))   #each site has its own probabilities of starting off in each age class
      }
      # Transition from juvenile to adult
    	trans.to.sa ~ dunif(0, 1)
    	trans.to.ad ~ dunif(0, 1)
      for(site in 1:nsites_mrc){
      	for(ind in 1:ninds[site]){ #loop through individuals
      	  age.class[site,ind,firsts[site,ind]] ~ dcat(prob.age[site,1:3])
      	  is.juv[site,ind,firsts[site,ind]] <- ifelse(age.class[site,ind,firsts[site,ind]]==1,1,0)
      		is.sub[site,ind,firsts[site,ind]] <- ifelse(age.class[site,ind,firsts[site,ind]]==2,1,0)
      		for(year in (firsts[site,ind]+1):nyears){
      		  prob.trans[site,ind,year,1] <- ifelse(age.class[site,ind,year-1]==1,1-trans.to.sa,0)
      		  prob.trans[site,ind,year,2] <- ifelse(age.class[site,ind,year-1]==1,trans.to.sa,ifelse(age.class[site,ind,year-1]==2,1-trans.to.ad,0))
      		  prob.trans[site,ind,year,3] <- ifelse(age.class[site,ind,year-1]==2,trans.to.ad,ifelse(age.class[site,ind,year-1]==3,1,0))
            age.class[site,ind,year] ~ dcat(prob.trans[site,ind,year,1:3])
      			is.juv[site,ind,year] <- ifelse(age.class[site,ind,year]==1,1,0)
      			is.sub[site,ind,year] <- ifelse(age.class[site,ind,year]==2,1,0)
      		}
      	}
      }

    
    
    #Age class imputation for telemetry data
      for(site in 1:nsites_tel){
        prob.age_telem[site,1:3] ~ ddirch(c(1,1,10))   #each site has its own probabilities of starting off in each age class
      }

      for(site in 1:nsites_tel){
      	for(ind in 1:ninds_telem_surv[site]){ #loop through individuals
      	  age.class_telem[site,ind,2] ~ dcat(prob.age_telem[site,1:3]) #since year 1 is not used in modeling, start at 2
      	  is.juv_telem[site,ind,2] <- ifelse(age.class_telem[site,ind,2]==1,1,0)
      		is.sub_telem[site,ind,2] <- ifelse(age.class_telem[site,ind,2]==2,1,0)
      		for(t in 3:(ntimes_telem_surv[site,ind] + 1)){
      		  prob.trans_tel[site,ind,t,1] <- ifelse(age.class_telem[site,ind,t-1]==1,1-trans.to.sa,0)
      		  prob.trans_tel[site,ind,t,2] <- ifelse(age.class_telem[site,ind,t-1]==1,trans.to.sa,ifelse(age.class_telem[site,ind,t-1]==2,1-trans.to.ad,0))
      		  prob.trans_tel[site,ind,t,3] <- ifelse(age.class_telem[site,ind,t-1]==2,trans.to.ad,ifelse(age.class_telem[site,ind,t-1]==3,1,0))
            age.class_telem[site,ind,t] ~ dcat(prob.trans_tel[site,ind,t,1:3])
      			is.juv_telem[site,ind,t] <- ifelse(age.class_telem[site,ind,t]==1,1,0)
      			is.sub_telem[site,ind,t] <- ifelse(age.class_telem[site,ind,t]==2,1,0)
      		}
      	}
      }
	    ###########
      # PROBABILITY OF SURVIVING
      ###########
   
      # Survival intercept
      phi0 ~ dunif(0, 1)   # prior for mean/intercept of logit survival   
      logit.phi0 <- log(phi0 / (1 - phi0))
      
      # Effects on survival
      phi.male.eff ~ dnorm(0, 1)    # effect of male on mean survival
      phi.juv.eff ~ dnorm(0, 1)   # Fixed effect of juvenile on mean survival
      phi.sub.eff ~ dnorm(0, 1)   # Fixed effect of subadult on mean survival 
      phi.precip.eff ~ dnorm(0, 1)    # effect of timelag precip
      phi.precip.x.sex.eff ~ dnorm(0,1) #effect of timelag precip x sex
    
      # Random site x year effect on survival - mark recap
      phi.siteyear.eff.mrc.sd ~ dunif(0, 2) # variation in survival among sites
      phi.sitexyear.eff.mrc.prec <- pow(phi.siteyear.eff.mrc.sd, -2)
      for(site in 1:nsites_mrc){
        for(year in 1:nyears){
          phi.siteyear.mrc.eff[year,site] ~ dnorm(0, phi.sitexyear.eff.mrc.prec)
        }
      }
      
      # Random site (not pixel) x year effect on survival - telemetry
      phi.siteyear.eff.tel.sd ~ dunif(0, 2)     # variation in survival among sites
      phi.sitexyear.eff.tel.prec <- pow(phi.siteyear.eff.tel.sd, -2)
      for(site in 1:nsites_tel){ 
        for(year in 1:nyears){
          phi.siteyear.tel.eff[year,site] ~ dnorm(0, phi.sitexyear.eff.tel.prec)
        }
      }
      
      ###########
      # PROBABILITY OF EMIGRATING
      ###########
      
      # Emigration - difference between true and apparent survival
      epsilon ~ dnorm(0, 1)

      # Effect of plot area on emigration rate
      eps.area.eff ~ dnorm(0, 1)
      
      ##########
      # MARK-RECAP-BASED SURVIVAL (logit-linear)
      ##########
      for(site in 1:nsites_mrc){
        for(ind in 1:ninds[site]){  # loop through individuals
          for(year in firsts[site,ind]:(nyears-1)){
            logit(thisphi[site,ind,year]) <- logit.phi0 +
              epsilon +
              eps.area.eff * site.area[site] +
              phi.male.eff * is.male[site,ind] +
              phi.juv.eff * is.juv[site,ind,year] +
              phi.sub.eff * is.sub[site,ind,year] +
              phi.precip.eff * precip_mrc[site,year] + 
              phi.precip.x.sex.eff * precip_mrc[site,year] * is.male[site,ind] +
              phi.siteyear.mrc.eff[year,site]
          }
        }
      }
      
      ##########
      # TELEMETRY-BASED SURVIVAL
      ##########

      for(site in 1:nsites_tel){
        for(ind in 1:ninds_telem_surv[site]){
          for(t in 2:(ntimes_telem_surv[site,ind]+1)){
            logit(psurviv_telem[site,ind,t]) <- logit.phi0 + 
              phi.male.eff * ismale_telem[site,ind] + 
              phi.juv.eff * is.juv_telem[site,ind,t] +
              phi.sub.eff * is.sub_telem[site,ind,t] +
              phi.precip.eff * precip_tel[site,ind,t] + 
              phi.precip.x.sex.eff * precip_tel[site,ind,t] * ismale_telem[site,ind] +
              phi.siteyear.tel.eff[year_ind_tel[site,ind,t],site_telem_surv[site,ind]]
          }
        }
      }

      ###########
      # SPATIAL COMPONENT
      ###########
    
      ###########
      # SIZE OF ACTIVITY AREA (log scale)
      ###########
      sigma0 ~ dunif(10, 700) # proxy for size of activity area (m)
      logsigma0 <- log(sigma0)

      male.hr.eff ~ dnorm(0, 1) # was unif -1 to 1
      juv.hr.eff ~ dnorm(0, 1)
      sub.hr.eff ~ dnorm(0, 1)
      fall.hr.eff ~ dnorm(0,1) #effect of fall survey method (3-day window)
    
      for(site in 1:nsites_mrc){
        for(ind in 1:ninds[site]){
          # Specify activity center for each individual (random from within larger area around plot)
          hrx[site,ind] ~ dunif(xlim[site,1], xlim[site,2])
          hry[site,ind] ~ dunif(ylim[site,1], ylim[site,2])
          
          for(year in firsts[site,ind]:nyears){
            log(sigma[site,ind,year]) <- logsigma0 +
              male.hr.eff * is.male[site,ind] +
              juv.hr.eff * is.juv[site,ind,year] +
              sub.hr.eff * is.sub[site,ind,year] + 
              fall.hr.eff * is.fall[site,year,1] #leaving this as 1 for now, since it doesnt change across a siteyear
            tau[site,ind,year] <- pow(sigma[site,ind,year], -2)
          }
        }
      }
      
      ###########
      # SAMPLE INDIVIDUAL LOCATIONS
      ###########
      
      for(site in 1:nsites_mrc){  
        for(ind in 1:ninds[site]){
          for(year in firsts[site,ind]:nyears){
            for(session in 1:nsessions[site,year]){
              locationx[site,ind,year,session] ~ dnorm(hrx[site,ind], tau[site,ind,year]) # data node (can be missing)
              locationy[site,ind,year,session] ~ dnorm(hry[site,ind], tau[site,ind,year]) # data node (can be missing)
              inplotx[site,ind,year,session] <- step(locationx[site,ind,year,session]-plot.xlim[site,1]) * step(plot.xlim[site,2]-locationx[site,ind,year,session])   
              inploty[site,ind,year,session] <- step(locationy[site,ind,year,session]-plot.ylim[site,1]) * step(plot.ylim[site,2]-locationy[site,ind,year,session])
              inplot[site,ind,year,session] <- inplotx[site,ind,year,session] * inploty[site,ind,year,session]
            } #session
          } #year
        } #ind
      } #site   
      
      #############
      # LIKELIHOOD
      #############
      
      ## Mark-recapture model
      for(site in 1:nsites_mrc){
        for(ind in 1:ninds[site]){
          alive[site,ind,firsts[site,ind]] ~ dbern(1)   # always alive on first year of capture
          #resident[site,ind,firsts[site,ind]] ~ dbern(1)    #always in the population at first cap
          for(session in (firsts2[site,ind]+1):nsessions[site,firsts[site,ind]]){
            p[site,ind,firsts[site,ind],session] <- thisp[site,ind,firsts[site,ind],session] * alive[site,ind,firsts[site,ind]] * inplot[site,ind,firsts[site,ind],session]   
            y[site,ind,firsts[site,ind],session] ~ dbern(p[site,ind,firsts[site,ind],session])
          }
          for(year in (firsts[site,ind]+1):nyears){
            expected.alive[site,ind,year] <- alive[site,ind,year-1] * thisphi[site,ind,year-1]
            alive[site,ind,year] ~ dbern(expected.alive[site,ind,year])
            #resident[site,ind,year] ~ dbern(disperse.prob[site]) # Chance to disperse out of pop
            for(session in 1:nsessions[site,year]){
  		        # Probability of catching an individual at a given point given it is in the population, in the study area, and real (after first PP)
              p[site,ind,year,session] <- thisp[site,ind,year,session] * alive[site,ind,year] * inplot[site,ind,year,session] #+ # if alive 
                 
              y[site,ind,year,session] ~ dbern(p[site,ind,year,session])                          
            } #session
          } #year
        } #ind
      } #site
      
      ## Telemetry model
      for(site in 1:nsites_tel){
        for(ind in 1:ninds_telem_surv[site]){
          alive_telem[site,ind,1] ~ dbern(1)
          for(t in 2:(ntimes_telem_surv[site,ind] + 1)){
            palive_telem[site,ind,t] <- alive_telem[site,ind,t-1] * pow(psurviv_telem[site,ind,t], 1)
            alive_telem[site,ind,t] ~ dbern(palive_telem[site,ind,t])
          }
        }
      }

      #############
      # DERIVED PARAMETERS
      #############
      logit(phi.female) <- logit.phi0 # True survival of females
      logit(phi.male) <- logit.phi0 + phi.male.eff # True survival of males
      logit(phi.juv.female) <- logit.phi0 + phi.juv.eff # True survival of telemetered juveniles
      logit(phi.juv.male) <- logit.phi0 + phi.juv.eff + phi.male.eff # True survival of juvenile males
      logit(phi.sub.female) <- logit.phi0 + phi.sub.eff # True survival of telemetered subadult female
      logit(phi.sub.male) <- logit.phi0 + phi.sub.eff + phi.male.eff # True survival of telemetered subadult male
    
      # Emigration rate as a consequence of plot area
      for(i in 1:20){
        logit(epsilon.area[i]) <- logit.phi0 + epsilon + eps.area.eff * ex.area[i]
      }
  
    for(cov in 1:50){ # 50 = as long as the input prediction dataset
        logit(phi.precip.female[cov]) <- logit.phi0 + phi.precip.eff * precip_predict[cov]
        logit(phi.precip.male[cov]) <- logit.phi0 + phi.precip.eff * precip_predict[cov] + phi.male.eff + phi.precip.x.sex.eff * precip_predict[cov] 
}
      for(year in 2:nyears){ #nyears - 1 as the connectivity mrc sites have an extra year to account for fall sampling
       
    for(site in 1:nsites_mrc){
          # Female site x year survival from CMR data (?)
          logit(phi.female.sitexyear.cmr[year,site]) <- logit.phi0 + phi.siteyear.mrc.eff[year,site]  + phi.precip.eff * precip_mrc[site,(year)] 
          }
      }


    }   ## end BUGS model
    
    ",file=BUGSfilename)
  
  return(BUGSfilename)
}



# Function to run a JAGS analysis

RunJags <- function(dfj, nits = 30000, burn = 10000){ #, filename ="DT_IPM_%sxxx.RData"){
  
  # system.time(
  mod <- jagsUI(
    data = dfj$data.for.bugs,
    inits = dfj$init.bugs,
    parameters.to.save = dfj$tomonitor,
    model.file = dfj$filename,
    n.chains = 3,
    n.adapt = 1000,
    n.iter = nits,   #nits,
    n.burnin = burn, 
    n.thin = 4,   # 4
    parallel = T
  )
  # )
  
  # 
  return(mod)
  
}
