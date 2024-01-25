### This script builds and fits the D. pandora model: PAN_TP.Briere_Comp.BH_Err.ZiNB_temp
 ## Last built on:2023-12-11 17:51:25.07825
 ## Author Chris Terry
 ## This script currently assumes data (LaggedFemales_PAN) is in the working directory

PRIORS = c(set_prior("normal(80, 15)", nlpar = "B0", lb = 0.0001, ub = 200),
set_prior("normal(20, 50)", nlpar = "Tmin", lb = 10, ub=22.5),
set_prior("normal(30, 50)", nlpar = "Tmax", lb = 28.7, ub=40),
set_prior("normal(0.00001, 0.5)", nlpar = "alphaii", lb = 0.00001, ub=20),
set_prior("normal(0.00001, 0.5)", nlpar = "alphaij", lb = 0.00001, ub=20)) ## end of priors


MODEL<-bf(PAN_FEMALE~Prev_PAN_FEMALE*(B0 * PrevTemp * ( PrevTemp - Tmin ) * sqrt( Tmax - PrevTemp ))*(1/(1+alphaii*Prev_PAN_FEMALE_comp + alphaij*Prev_PAL_FEMALE)),
zi~PrevTemp,
B0+Tmin+Tmax+alphaii+alphaij~1, nl = TRUE, family = zero_inflated_negbinomial(link = "identity"))


PAN_TP.Briere_Comp.BH_Err.ZiNB_temp_fit <- brm(MODEL, data =  LaggedFemales_PAN,
prior = PRIORS, chains=4,cores =4, iter=4000)
 PAN_TP.Briere_Comp.BH_Err.ZiNB_temp_fit <- add_criterion(PAN_TP.Briere_Comp.BH_Err.ZiNB_temp_fit,criterion = c("loo", "waic"))

FIT <- PAN_TP.Briere_Comp.BH_Err.ZiNB_temp_fit 
save(PAN_TP.Briere_Comp.BH_Err.ZiNB_temp_fit , file ="../AutoModels/Fits/FEMS_PAN_TP.Briere_Comp.BH_Err.ZiNB_temp")
 print(Sys.time())
