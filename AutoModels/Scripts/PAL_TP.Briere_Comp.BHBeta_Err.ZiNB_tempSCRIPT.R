### This script builds and fits the model: PAL_TP.Briere_Comp.BHBeta_Err.ZiNB_temp
 ## Last built on:2024-01-12 15:29:31.306012
 ## Author Chris Terry
 ## This script currently assumes data (LaggedFemales_PAL) is in the working directory

PRIORS = c(set_prior("normal(80, 15)", nlpar = "B0", lb = 0.0001, ub = 200),
set_prior("normal(20, 50)", nlpar = "Tmin", lb = 10, ub=22.5),
set_prior("normal(30, 50)", nlpar = "Tmax", lb = 27.9, ub=35),
set_prior("normal(1, 5)", nlpar = "beta", lb =0.1, ub=5),
set_prior("normal(0.00001, 0.5)", nlpar = "alphaii", lb = 0.00001, ub=20),
set_prior("normal(0.00001, 0.5)", nlpar = "alphaij", lb = 0.00001, ub=20)) ## end of priors


MODEL<-bf(PAL_FEMALE~Prev_PAL_FEMALE*(B0 * PrevTemp * ( PrevTemp - Tmin ) * sqrt( Tmax - PrevTemp ))*(1/((1+alphaii*Prev_PAL_FEMALE_comp + alphaij*Prev_PAN_FEMALE)^beta)),
zi~PrevTemp,
B0+Tmin+Tmax+beta+alphaii+alphaij~1, nl = TRUE, family = zero_inflated_negbinomial(link = "identity"))


PAL_TP.Briere_Comp.BHBeta_Err.ZiNB_temp_fit <- brm(MODEL, data =  LaggedFemales_PAL,
prior = PRIORS, chains=4,cores =4, iter=4000)
 PAL_TP.Briere_Comp.BHBeta_Err.ZiNB_temp_fit <- add_criterion(PAL_TP.Briere_Comp.BHBeta_Err.ZiNB_temp_fit,criterion = c("loo", "waic"))

FIT <- PAL_TP.Briere_Comp.BHBeta_Err.ZiNB_temp_fit 
save(PAL_TP.Briere_Comp.BHBeta_Err.ZiNB_temp_fit , file ="../AutoModels/Fits/FEMS_PAL_TP.Briere_Comp.BHBeta_Err.ZiNB_temp")
 print(Sys.time())
