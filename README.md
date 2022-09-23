# JordanFallsWshShare
The files here include data sets and codes for the following manuscript: Contrasting annual and summer phosphorus export using a hybrid Bayesian watershed model.

Items in this repository:

Annual_TP.rds- This is a R data file that includes all the input data sets for the annual model. It needs to be imported into R with the appropriate command.

Summer_TP.rds- This is a R data file that includes all the input data sets for the summer model. It needs to be imported into R with the appropriate command.

RStan_codes_WRR.R- RStan code needed to run the models. The "rstan" (https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) and "rstudioapi" packages need to be installed prior to running the models. 
# Example
Step 1: Install "rstan" and "rstudioapi" packages. 

Step 2: Compile the "stanmodelcode_annual" annual model (line #5 in the "RStan_codes_WRR.R" code) .

Step 3: Read the annual model input data set (Annual_TP.rds), if necessary (line #168)

Step 4: Run the stan function (line #170). Suggest using the following parameters: iter=4000, warmup=2000, thin=5, chains=3,cores=3,adapt_delta =0.99 ,max_treedepth =25. The estimated runtime on desktop is about one day. Note that more iterations and warmup steps may improve model convergence.
