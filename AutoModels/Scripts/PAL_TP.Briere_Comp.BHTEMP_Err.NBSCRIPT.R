### This script builds and fits the model: PAL_TP.Briere_Comp.BHTEMP_Err.NB
 ## Last built on:2024-01-12 15:29:31.340488
 ## Author Chris Terry
 ## This script currently assumes data (LaggedFemales_PAL) is in the working directory

PRIORS = c(set_prior("normal(80, 15)", nlpar = "B0", lb = 0.0001, ub = 200),
set_prior("normal(20, 50)", nlpar = "Tmin", lb = 10, ub=22.5),
set_prior("normal(30, 50)", nlpar = "Tmax", lb = 27.9, ub=35),
set_prior("normal(0.00001, 0.5)", nlpar = "alphaii0", lb = 0.00001, ub=20),
set_prior("normal(0, 0.05)", nlpar = "alphaiiT", lb = -2, ub=2),
set_prior("normal(0.00001, 0.5)", nlpar = "alphaij0", lb = 0.00001, ub=20),
set_prior("normal(0, 0.05)", nlpar = "alphaijT", lb = -2, ub=2)) ## end of priors


MODEL<-bf(PAL_FEMALE~Prev_PAL_FEMALE*(B0 * PrevTemp * ( PrevTemp - Tmin ) * sqrt( Tmax - PrevTemp ))*(1/((1+(alphaii0+alphaiiT*Temp20)*Prev_PAL_FEMALE_comp + (alphaij0+alphaijT*Temp20)*Prev_PAN_FEMALE))),

B0+Tmin+Tmax+alphaii0+alphaij0+alphaiiT+alphaijT~1, nl = TRUE, family = negbinomial(link = "identity"))


PAL_TP.Briere_Comp.BHTEMP_Err.NB_fit <- brm(MODEL, data =  LaggedFemales_PAL,
prior = PRIORS, chains=4,cores =4, iter=4000)
 PAL_TP.Briere_Comp.BHTEMP_Err.NB_fit <- add_criterion(PAL_TP.Briere_Comp.BHTEMP_Err.NB_fit,criterion = c("loo", "waic"))

FIT <- PAL_TP.Briere_Comp.BHTEMP_Err.NB_fit 
save(PAL_TP.Briere_Comp.BHTEMP_Err.NB_fit , file ="../AutoModels/Fits/FEMS_PAL_TP.Briere_Comp.BHTEMP_Err.NB")
 print(Sys.time())
