Make_brms_script_PAN_FEM <- function(i,Option_df , cores = 4, iter = 2000){
  
  TP = Option_df$TP[i]
  Comp = Option_df$Comp[i]
  Err = Option_df$Err[i]
  
  Name = paste(paste0('PAN_TP.', TP)   ,
               paste0('Comp.', Comp)   ,
               paste0('Err.', Err)   ,
               sep = '_')
  
  
  ### Thermal Performance Curve
  if( TP == 'TaylorSexton'){
    TP_Priors <- c('set_prior("normal(80, 15)", nlpar = "B0", lb =10, ub = 200)',
                   'set_prior("normal(25, 50)", nlpar = "Tpk", lb = 20, ub=30)',
                   'set_prior("normal(2, 50)", nlpar = "Tmin", lb = 1, ub=25)')
    TP_function <- '(B0 * ( -(PrevTemp - Tmin)^4 + 2 * ((PrevTemp - Tmin)^2) * (Tpk - Tmin)^2 ) / (Tpk - Tmin)^4)'
    TP_NL =  paste0( c("B0" , "Tpk","Tmin"), collapse = '+')
  }
  
  if( TP == 'Briere'){
    TP_Priors <- c('set_prior("normal(80, 15)", nlpar = "B0", lb = 0.0001, ub = 200)',
                   'set_prior("normal(20, 50)", nlpar = "Tmin", lb = 10, ub=22.5)',
                   'set_prior("normal(30, 50)", nlpar = "Tmax", lb = 28.7, ub=40)')
    TP_function <- '(B0 * PrevTemp * ( PrevTemp - Tmin ) * sqrt( Tmax - PrevTemp ))'
    TP_NL =  paste0( c("B0" , "Tmin","Tmax"), collapse = '+')
  }
  
  if( TP ==   'Beta'){
    TP_Priors <- c('set_prior("normal(80, 15)", nlpar = "B0", lb = 0.001, ub = 200)',
                   'set_prior("normal(25, 50)", nlpar = "a", lb = 0.0001, ub=50)',
                   'set_prior("normal(2, 5)", nlpar = "b", lb = 0.00001, ub=5)')
    TP_function <- '(B0 * (a - (PrevTemp/10))*(PrevTemp/10)^b)'
    TP_NL =  paste0( c("B0" , "a","b"), collapse = '+')
  }

  
  ## Competition Kernel
  
  if( Comp  == 'BH'){
    Comp_Priors <-  c('set_prior("normal(0.00001, 0.5)", nlpar = "alphaii", lb = 0.00001, ub=20)',
                      'set_prior("normal(0.00001, 0.5)", nlpar = "alphaij", lb = 0.00001, ub=20)')
    Comp_function <- '(1/(1+alphaii*Prev_PAN_FEMALE_comp + alphaij*Prev_PAL_FEMALE))'
    Comp_NL =  paste0( c("alphaii", "alphaij" ), collapse = '+')
  }
  
  if( Comp  == 'BHBeta'){
    Comp_Priors <-  c('set_prior("normal(1, 5)", nlpar = "beta", lb = 0.1, ub=5)',
                      'set_prior("normal(0.00001, 0.5)", nlpar = "alphaii", lb = 0.00001, ub=20)',
                      'set_prior("normal(0.00001, 0.5)", nlpar = "alphaij", lb = 0.00001, ub=20)')
    Comp_function <- '(1/((1+alphaii*Prev_PAN_FEMALE_comp + alphaij*Prev_PAL_FEMALE)^beta))'
    Comp_NL =  paste0( c("beta", "alphaii" ,"alphaij" ), collapse = '+')
  }
  if( Comp  == 'BHTEMP'){
    Comp_Priors <-  c('set_prior("normal(0.00001, 0.5)", nlpar = "alphaii0", lb = 0.00001, ub=20)',
                      'set_prior("normal(0, 0.05)", nlpar = "alphaiiT", lb = -2, ub=2)',
                      'set_prior("normal(0.00001, 0.5)", nlpar = "alphaij0", lb = 0.00001, ub=20)',
                      'set_prior("normal(0, 0.05)", nlpar = "alphaijT", lb = -2, ub=2)')
    
    Comp_function <- '(1/((1+(alphaii0+alphaiiT*Temp20)*Prev_PAN_FEMALE_comp + (alphaij0+alphaijT*Temp20)*Prev_PAL_FEMALE)))'
    Comp_NL =  paste0( c("alphaii0" ,"alphaij0",  "alphaiiT" ,"alphaijT" ), collapse = '+')
  }
  
  if( Comp  == 'BHTEMP1'){
    Comp_Priors <-  c('set_prior("normal(0.00001, 1)", nlpar = "alphaii", lb = 0.00001, ub=20)',
                      'set_prior("normal(0.00001, 1)", nlpar = "alphaij0", lb = 0.00001, ub=20)',
                      'set_prior("normal(0, 0.05)", nlpar = "alphaijT", lb = -2, ub=2)')
    
    Comp_function <- '(1/((1+alphaii*Prev_PAN_FEMALE_comp + (alphaij0+alphaijT*Temp20)*Prev_PAL_FEMALE)))'
    Comp_NL =  paste0( c("alphaii" ,"alphaij0" ,"alphaijT" ), collapse = '+')
  }
  
  if( Comp  == 'BHTEMP1Beta'){
    Comp_Priors <-  c('set_prior("normal(1, 5)", nlpar = "beta", lb = 0.1, ub=5)',
                      'set_prior("normal(0.00001, 1)", nlpar = "alphaii", lb = 0.00001, ub=20)',
                      'set_prior("normal(0.00001, 1)", nlpar = "alphaij0", lb = 0.00001, ub=20)',
                      'set_prior("normal(0, 0.05)", nlpar = "alphaijT", lb = -2, ub=2)')
    
    Comp_function <- '(1/((1+alphaii*Prev_PAN_FEMALE_comp + (alphaij0+alphaijT*Temp20)*Prev_PAL_FEMALE)^beta))'
    Comp_NL =  paste0( c("alphaii" ,"alphaij0" ,"alphaijT" , 'beta'), collapse = '+')
  }
  
  if( Comp  == 'BHTEMPBeta'){
    Comp_Priors <-  c('set_prior("normal(1, 5)", nlpar = "beta", lb = 0.1, ub=5)',
                      'set_prior("normal(0.00001, 0.5)", nlpar = "alphaii0", lb = 0.00001, ub=20)',
                      'set_prior("normal(0, 0.05)", nlpar = "alphaiiT", lb = -2, ub=2)',
                      'set_prior("normal(0.00001, 0.5)", nlpar = "alphaij0", lb = 0.00001, ub=20)',
                      'set_prior("normal(0, 0.05)", nlpar = "alphaijT", lb = -2, ub=2)')
    
    Comp_function <- '(1/((1+(alphaii0+alphaiiT*Temp20)*Prev_PAN_FEMALE_comp + (alphaij0+alphaijT*Temp20)*Prev_PAL_FEMALE))^beta)'
    Comp_NL =  paste0( c("alphaii0" ,"alphaij0",  "alphaiiT" ,"alphaijT" , 'beta'), collapse = '+')
  }
  
  
  if( Comp  == 'BHTEMP1Log'){
    Comp_Priors <-  c('set_prior("normal(0.00001, 1)", nlpar = "alphaii", lb = 0.00001, ub=20)',
                      'set_prior("normal(0.00001, 1)", nlpar = "alphaij0", lb = 0.00001, ub=20)',
                      'set_prior("normal(0, 0.05)", nlpar = "alphaijT", lb = -2, ub=2)')
    
    Comp_function <- '(1/((1+alphaii*log(Prev_PAN_FEMALE_comp+1) + (alphaij0+alphaijT*Temp20)*log(Prev_PAL_FEMALE+1))))'
    Comp_NL =  paste0( c("alphaii" ,"alphaij0" ,"alphaijT" ), collapse = '+')
  }
  
  if( Comp  == 'BHTEMP1LogBeta'){
    Comp_Priors <-  c('set_prior("normal(1, 5)", nlpar = "beta", lb = 0.1, ub=5)',
                      'set_prior("normal(0.00001, 1)", nlpar = "alphaii", lb = 0.00001, ub=20)',
                      'set_prior("normal(0.00001, 1)", nlpar = "alphaij0", lb = 0.00001, ub=20)',
                      'set_prior("normal(0, 0.05)", nlpar = "alphaijT", lb = -2, ub=2)')
    
    Comp_function <- '(1/((1+alphaii*log(Prev_PAN_FEMALE_comp+1) + (alphaij0+alphaijT*Temp20)*log(Prev_PAL_FEMALE+1))^beta))'
    Comp_NL =  paste0( c("alphaii" ,"alphaij0" ,"alphaijT" , 'beta'), collapse = '+')
  }
  
  ##  Error Distribution
  
  if(Err == 'Norm'){
    Err_Priors = NULL
    Err_family = 'gaussian()' 
    Err_function=NULL
  }
  
  if(Err == 'Pois'){
    Err_Priors = NULL
    Err_family = 'poisson(link = "identity")' 
    Err_function=NULL
  }
  if(Err == 'NB'){
    Err_Priors = NULL
    Err_family = 'negbinomial(link = "identity")'
    Err_function=NULL
  }
  if(Err == 'ZiNB'){
    Err_Priors = NULL
    Err_family = 'zero_inflated_negbinomial(link = "identity")' 
    Err_function=NULL
  }
  if(Err == 'ZiNB_temp'){
    Err_Priors = NULL
    Err_family = 'zero_inflated_negbinomial(link = "identity")' 
    Err_function = 'zi~PrevTemp,'
  }
  if(Err == 'ZiNB_comp'){
    Err_Priors = NULL
    Err_family = 'zero_inflated_negbinomial(link = "identity")' 
    Err_function = 'zi~Prev_PAN_FEMALE,'
  }
  
  if(Err == 'ZiNB_shapetemp'){
    Err_Priors = NULL
    Err_family = 'zero_inflated_negbinomial(link = "identity")' 
    Err_function = 'shape~PrevTemp,'
  }
  if(Err == 'ZiNB_shapecomp'){
    Err_Priors = NULL
    Err_family = 'zero_inflated_negbinomial(link = "identity")' 
    Err_function = 'shape~Prev_PAN_FEMALE,'
  }
  
  Intro = paste0('### This script builds and fits the D. pandora model: ', Name,
                 '\n ## Last built on:', Sys.time(),
                 '\n ## Author Chris Terry\n ## This script currently assumes data (LaggedFemales_PAN) is in the working directory\n')  
  
  Priors <-  paste0('PRIORS = c(', paste0(
    
    c(TP_Priors,Comp_Priors,Err_Priors), 
    collapse = ',\n' ), 
    ') ## end of priors')
  
  NLs = paste0( c(TP_NL,Comp_NL),collapse = '+')
  
  
  ModelFunction <- paste0('\n\nMODEL<-bf(PAN_FEMALE~Prev_PAN_FEMALE*',TP_function, '*', Comp_function, ',\n',
                          Err_function, '\n',
                          NLs, '~1, nl = TRUE, family = ', Err_family, ')\n\n')
  
  
  ModelFit <- paste0( Name, '_fit <- brm(MODEL, data =  LaggedFemales_PAN,\n',
                      'prior = PRIORS, chains=4,cores =', cores ,
                      ', iter=',iter,')\n ',  Name,
                      '_fit <- add_criterion(',Name,
                      '_fit,criterion = c("loo", "waic"))\n' )
  
  Saving <- paste0('FIT <- ',Name,'_fit ',
                   '\nsave(', Name,'_fit , file ="../AutoModels/Fits/FEMS_', Name,'")\n print(Sys.time())')
  
  FileName = paste0('../AutoModels/Scripts/', Name, 'SCRIPT.R')
  
  write_lines(c(Intro,Priors, ModelFunction, ModelFit,Saving),
              file = FileName,
              num_threads = 1) 
  return( FileName)
  
}