# JordanFallsWshWRR
The files here include data sets and codes for the following manuscript: Contrasting annual and summer phosphorus loading and retention rates using a hybrid Bayesian watershed model

Items in this repository:

Annual_TP.rds- This is a R data file that includes all the input data sets for the annual model. It needs to be imported into R with the appropriate command.

Summer_TP.rds- This is a R data file that includes all the input data sets for the summer model. It needs to be imported into R with the appropriate command.

RStan_codes_WRR.R- RStan code needed to run the models. The "rstan" (https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) and "rstudioapi" packages need to be installed prior to running the models. 
# Example
Step 1: Install the packages. 
Step 2: Run the "stanmodelcode_annual" for the annual model.
Step 3: Upload the annual model input data set (Annual_TP).
Step 4: Run the stan function with the following parameters:
model_annual = stan(model_code=stanmodelcode_annual, data=Annual_TP, iter=4000, 
             warmup=2000, thin=5, chains=3,cores=3,
             control = list(adapt_delta =0.99 ,max_treedepth =26 ))
