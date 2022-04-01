stanmodelcode_annual = '

data { 
int nl;         //subwatersheds separated by year                            
int nr;         //incremental watershed-years
int nd;         // dischargers  
vector [nr] load;       //WRTDS load at LMS
vector [nr] SD;         	//WRTDS load SD
vector [nl] chick;        //number of chickens in subwatershed
vector [nl] cow;          //number of cows in subwatershed
vector [nl] hog;          //number of hogs (swine) in subwatershed
int wsd [nr];             //count variable for LMSs
int wshed_size;           //number of LMS watersheds
vector [nr] increm_area;  //Incremental area for each loading station
vector [nr] av_prec;      //normalized precipitation
vector [nr] av_prec2;     //scaled precipitation
vector [nr] up_t_load1;   //upstream loading for nested wsds
vector [nr] up_t_load2;   //upstream loading for nested wsds
vector [nr] up_t_load3;   //upstream loading for nested wsds
vector [nr] str_loss_load1;       //stream losses of upstream loading
vector [nr] str_loss_load2;       //stream losses of upstream loading
vector [nr] str_loss_load3;       //stream losses of upstream loading
vector [nr] res_loss_load1;       //reservoir losses of upstream loading
vector [nr] res_loss_load2;       //reservoir losses of upstream loading
vector [nr] res_loss_load3;       //reservoir losses of upstream loading
vector [nd] d_loss_str;           //stream losses of dischargers to LMSs
vector [nd] d_loss_res;           //reservoir losses of dischargers to LMSs
vector [nd] d_vals;               //discharger loadings
vector [nl] l_loss_str;           //stream losses for subwatersheds
vector [nl] l_loss_res;           //reservoir losses for subwatersheds
vector [nl] ag;                   //agriculture area in subwatershed
vector [nl] devpre;               //pre-1980 urban area in subwatershed
vector [nl] devpost;		          //post-1980 urban area in subwatershed
vector [nl] wild;			            //undeveloped area in subwatershed
vector [nl] tot_l;			          //total size of subwatershed
int l_start [nr];			            //index to link subwatersheds to LMSs
int l_end [nr];			              //index to link subwatersheds to LMSs
int d_start [nr];			            //index to link dischargers to LMSs
int d_end [nr];			              //index to link dischargers to LMSs

}

transformed data{

}

parameters {  

real<lower =0> Be_a;        //Agriculture EC
real<lower =0> Be_d_pre;        //pre-1980 developed EC 
real<lower =0> Be_d_post;        //post-1980 developed EC
real<lower =0> Be_w;       //Undeveloped EC
real<lower =0, upper = 1> Be_ch;       //Chickens EC
real<lower =0, upper = 1> Be_h;       //Hogs (swine) EC
real<lower =0, upper = 1> Be_cw;       //Cows EC
real <lower =0, upper = 1> Sn;          // Stream retention
real <lower =1, upper = 60> Sn2;        //Reservoir retention
real <lower =0, upper = 0.40> PIC_q;      //PIC for retention
vector<lower =0, upper = 10> [7] pic_p;   //PIC for land classes
real <lower=0> Be_dch;      			//Point source DC 
real<lower=0, upper = 1000000> sigma_res;		//Model residual SD
real <lower=0, upper = 1000000> sigma_w;		//Random effect SD
vector [wshed_size] alpha;			//# of watersheds
real<lower = 0, upper = 2> sigma_B1;		//PIC SD
real<lower = 0, upper = 3> Bp_mean;		//PIC mean
vector [nr] ly;					// unknown true loads

}

transformed parameters {


}

model {
vector [nr] tot; 			// Sum of all loadings from all sources
vector [nr] sigma;		//SD for watershed random effects
vector [nr] y_hat;		//total loadings plus random effects minus waterbody losses
vector [nr] Dpr;			//To compile pre-1980 urban with PIC
vector [nr] Dpt;			//To compile post-1980 urban with PIC
vector [nr] A;			//To compile Agricultur with PIC 
vector [nr] W;			//To compile undeveloped urban with PIC
vector [nr] D;
vector [nr] Dch; 		 	//To compile dischargers with stream/reservoir losses
vector [nr] alpha_vals;		// Watershed indicator
vector [nr] A_lc;			//Agriculture vector
vector [nr] D_lc_pre;		//pre-1980  vector		
vector [nr] D_lc_post;		//post-1980  vector
vector [nr] W_lc;		//Undeveloped  vector
vector [nr] Disch;		//point source dischargers
vector [nr] C;   			//chickens for adding PIC
vector [nr] H;			//hogs for adding PIC
vector [nr] Cw;			//cows for adding PIC
vector [nr] C_r;   		//chickens for aggregating subwatersheds
vector [nr] H_r;			//swine for aggregating subwatersheds
vector [nr] Cw_r;		//cows  for aggregating subwatersheds
vector [nr] y;				//loading
int w;

// Loop to determine export for each watershed-year
for(i in 1:nr){
    //Looping to aggregate subwatersheds and corresponding stream and reservoir retention.
A_lc[i]= sum((ag[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
D_lc_pre[i]= sum((devpre[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
D_lc_post[i]= sum((devpost[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
W_lc[i]= sum((wild[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Disch[i]=sum((d_vals[d_start[i]:d_end[i]]) .* exp((-Sn/(1+PIC_q*av_prec[i])) * d_loss_str[d_start[i]:d_end[i]]) .* exp((-Sn2/(1+PIC_q*av_prec[i])) ./ d_loss_res[d_start[i]:d_end[i]]));
C_r[i]= sum((chick[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
H_r[i]= sum((hog[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Cw_r[i]= sum((cow[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));

    //Adding precipitation impact coefficient to land export 

A[i]=  Be_a * pow(av_prec2[i],pic_p[1]) .* A_lc[i];
Dpr[i] = Be_d_pre * pow(av_prec2[i],pic_p[2]) .* D_lc_pre[i];
Dpt[i] = Be_d_post * pow(av_prec2[i],pic_p[3]) .* D_lc_post[i];
D[i]   = Dpr[i]+Dpt[i];
W[i] = Be_w * pow(av_prec2[i],pic_p[4]) .* W_lc[i];

C[i] = Be_ch * pow(av_prec2[i],pic_p[5]) .* C_r[i];
H[i] = Be_h * pow(av_prec2[i],pic_p[6]) .* H_r[i];
Cw[i] = Be_cw * pow(av_prec2[i],pic_p[7]) .* Cw_r[i];

Dch[i] = Be_dch * Disch[i];
}


for (i in 1:nr){
  //Loop to determine random effect for each watershed
w= wsd[i];
sigma[i] = sqrt(pow(SD[i],2)+pow(sigma_res,2));   //with regular sigma values
alpha_vals[i] = alpha[w];

}
//Sum loadings from all sources
tot =  A + D + W + C + H + Cw + Dch;
//Add random effects to source loadings and subtract losses from upstream loads
y_hat = tot + alpha_vals - up_t_load1 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load1) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load1))) - up_t_load2 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load2) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load2))) - up_t_load3 .* (1-(exp((-Sn ./(1+PIC_q*av_prec)) .* str_loss_load3) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load3)));


//priors
Be_a ~ normal(100,65);  //Prior for agriculture
Be_d_pre ~ normal(100,90);  //Prior for pre-1980 development
Be_d_post ~ normal(100,90);  //Prior forpost-1980 development
Be_w ~ normal(15,5);   //Prior for undeveloped

Be_ch ~ normal(0.005,0.0025);  //Prior for chickens
Be_h ~ normal(0.02,0.01);  //Prior for hogs (swine)
Be_cw ~ normal(0,5);  //Uninformed Prior for cows

Be_dch ~ normal(1,.03);   //Prior for point source delivery
sigma_res ~ normal(0,1000000); //st error of the model
sigma_w ~ normal(0,1000000);     //st. deviation of random effect hyperdistribution
alpha ~ normal(0,sigma_w);    //watershed random effects
sigma_B1 ~ normal(0,.5);  //st. deviation of precipitation coefficient
Bp_mean ~ normal(1,.5);    //mean PIC for hyperdistribution
pic_p[1] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for ag
pic_p[2] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for pre-1980 deve
pic_p[3] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for post-1980 dev
pic_p[4] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for undev
pic_p[5] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for chicken
pic_p[6] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for swine
pic_p[7] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for cow


Sn ~ normal(.2,.08);    //Prior for stream retention rate
Sn2 ~ normal(30,8.5);   //Prior for reservoir retention rate
PIC_q ~ normal(0,1);    //prior for PIC for retention
ly ~ normal(y_hat,sigma_res);        //parameter that calibrates ly_hat with ly () 
load ~ normal(ly,SD);                // load = WRTDS estimate, SD = WRTDS sd

}

generated quantities {

}
'

stanmodelcode_summer = '

data { 
int nl;         //subwatersheds separated by year                            
int nr;         //incremental watershed-years
int nd;         // dischargers  
vector [nr] load;       //WRTDS load at LMS
vector [nr] SD;         	//WRTDS load SD
vector [nl] chick;        //number of chickens in subwatershed
vector [nl] cow;          //number of cows in subwatershed
vector [nl] hog;          //number of hogs (swine) in subwatershed
int wsd [nr];             //count variable for LMSs
int wshed_size;           //number of LMS watersheds
vector [nr] increm_area;  //Incremental area for each loading station
vector [nr] up_t_load1;   //upstream loading for nested wsds
vector [nr] up_t_load2;   //upstream loading for nested wsds
vector [nr] up_t_load3;   //upstream loading for nested wsds
vector [nr] str_loss_load1;       //stream losses of upstream loading
vector [nr] str_loss_load2;       //stream losses of upstream loading
vector [nr] str_loss_load3;       //stream losses of upstream loading
vector [nr] res_loss_load1;       //reservoir losses of upstream loading
vector [nr] res_loss_load2;       //reservoir losses of upstream loading
vector [nr] res_loss_load3;       //reservoir losses of upstream loading
vector [nd] d_loss_str;           //stream losses of dischargers to LMSs
vector [nd] d_loss_res;           //reservoir losses of dischargers to LMSs
vector [nd] d_vals;               //discharger loadings
vector [nl] l_loss_str;           //stream losses for subwatersheds
vector [nl] l_loss_res;           //reservoir losses for subwatersheds
vector [nl] ag;                   //agriculture area in subwatershed
vector [nl] devpre;               //pre-1980 urban area in subwatershed
vector [nl] devpost;		          //post-1980 urban area in subwatershed
vector [nl] wild;			            //undeveloped area in subwatershed
vector [nl] tot_l;			          //total size of subwatershed
int l_start [nr];			            //index to link subwatersheds to LMSs
int l_end [nr];			              //index to link subwatersheds to LMSs
int d_start [nr];			            //index to link dischargers to LMSs
int d_end [nr];			              //index to link dischargers to LMSs
matrix[nr,9] pr;                  //precipitation matrix from Jan-Aug
}

transformed data{


}

parameters {  

real<lower =0> Be_a;        //Agriculture EC
real<lower =0> Be_d_pre;        //pre-1980 developed EC 
real<lower =0> Be_d_post;        //post-1980 developed EC
real<lower =0> Be_w;       //Undeveloped EC
real<lower =0, upper = 1> Be_ch;       //Chickens EC
real<lower =0, upper = 1> Be_h;       //Hogs (swine) EC
real<lower =0, upper = 1> Be_cw;       //Cows EC
real <lower =0, upper = 1> Sn;          // Stream retention
real <lower =1, upper = 60> Sn2;        //Reservoir retention
real <lower =0, upper = 0.40> PIC_q;      //PIC for retention
vector<lower =0, upper = 10> [7] pic_p;   //PIC for land classes
real <lower=0> Be_dch;      			//Point source DC 
real<lower=0, upper = 1000000> sigma_res;		//Model residual SD
real <lower=0, upper = 1000000> sigma_w;		//Random effect SD
vector [wshed_size] alpha;			//# of watersheds
real<lower = 0, upper = 2> sigma_B1;		//PIC SD
real<lower = 0, upper = 3> Bp_mean;		//PIC mean
vector [nr] ly;					// unknown true loads
real<lower = 1 , upper=9> Be_psi; //precipitation weighting coeff.
}

transformed parameters {


}

model {
vector [nr] tot; 			// Sum of all loadings from all sources
vector [nr] sigma;		//SD for watershed random effects
vector [nr] y_hat;		//total loadings plus random effects minus waterbody losses
vector [nr] A;			//To compile Agriculture with PIC
vector [nr] Dpr;			//To compile pre-1980 urban with PIC
vector [nr] Dpt;			//To compile post-1980 urban with PIC
vector [nr] W;			//To compile undeveloped urban with PIC
vector [nr] Dch; 		 	//To compile dischargers with stream/reservoir losses
vector [nr] alpha_vals;		// Watershed indicator
vector [nr] A_lc;			//Agriculture vector
vector [nr] D_lc_pre;		//pre-1980  vector		
vector [nr] D_lc_post;		//post-1980  vector
vector [nr] W_lc;		//Undeveloped  vector
vector [nr] Disch;		//point source dischargers
vector [nr] C;   			//chickens for adding PIC
vector [nr] H;			//hogs for adding PIC
vector [nr] Cw;			//cows for adding PIC
vector [nr] C_r;   		//chickens for aggregating subwatersheds
vector [nr] H_r;			//swine for aggregating subwatersheds
vector [nr] Cw_r;		//cows  for aggregating subwatersheds
vector [nr] y;				//loading
vector [8] psi;       //precipitation weight
vector [nr] prec;     //scaled precipitation
vector [nr] av_prec;  //normalized precipitation
int w;
//Loop to get precipitation weighting coefficient
for (m in 1:8){
if (m<=(Be_psi-1)){
psi[m]=0;
} else if ((Be_psi-1)<m && m<=Be_psi){
psi[m]=m+1-Be_psi;
} else if (m>=Be_psi){
psi[m]=1;
}
}
//weighted precipitation
prec=pr*psi;

av_prec=(prec-mean(prec))/sd(prec);
prec=prec/mean(prec);

// Loop to determine export for each watershed-year
for(i in 1:nr){
//Looping to aggregate subwatersheds and corresponding stream and reservoir retention.
A_lc[i]= sum((ag[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
D_lc_pre[i]= sum((devpre[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
D_lc_post[i]= sum((devpost[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
W_lc[i]= sum((wild[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Disch[i]=sum((d_vals[d_start[i]:d_end[i]]) .* exp((-Sn/(1+PIC_q*av_prec[i])) * d_loss_str[d_start[i]:d_end[i]]) .* exp((-Sn2/(1+PIC_q*av_prec[i])) ./ d_loss_res[d_start[i]:d_end[i]]));
C_r[i]= sum((chick[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
H_r[i]= sum((hog[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Cw_r[i]= sum((cow[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));

//Adding precipitation impact coefficient to land export 

A[i]=  Be_a * pow(prec[i],pic_p[1]) .* A_lc[i];
Dpr[i] = Be_d_pre * pow(prec[i],pic_p[2]) .* D_lc_pre[i];
Dpt[i] = Be_d_post * pow(prec[i],pic_p[3]) .* D_lc_post[i];
D[i]   = Dpr[i]+Dpt[i];
W[i] = Be_w * pow(prec[i],pic_p[4]) .* W_lc[i];

C[i] = Be_ch * pow(prec[i],pic_p[5]) .* C_r[i];
H[i] = Be_h * pow(prec[i],pic_p[6]) .* H_r[i];
Cw[i] = Be_cw * pow(prec[i],pic_p[7]) .* Cw_r[i];

Dch[i] = Be_dch * Disch[i];
}


for (i in 1:nr){
//Loop to determine random effect for each watershed
w= wsd[i];
sigma[i] = sqrt(pow(SD[i],2)+pow(sigma_res,2));   //with regular sigma values
alpha_vals[i] = alpha[w];

}
//Sum loadings from all sources
tot =  A + D + W + C + H + Cw + Dch;
//Add random effects to source loadings and subtract losses from upstream loads
y_hat = tot + trans * alpha_vals - up_t_load1 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load1) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load1))) - up_t_load2 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load2) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load2))) - up_t_load3 .* (1-(exp((-Sn ./(1+PIC_q*av_prec)) .* str_loss_load3) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load3)));


//priors
Be_a ~ normal(100,65);  //Prior for agriculture
Be_d_pre ~ normal(100,90);  //Prior for pre-1980 development
Be_d_post ~ normal(100,90);  //Prior forpost-1980 development
Be_w ~ normal(15,5);   //Prior for undeveloped

Be_ch ~ normal(0.005,0.0025);  //Prior for chickens
Be_h ~ normal(0.02,0.01);  //Prior for hogs (swine)
Be_cw ~ normal(0,5);  //Uninformed Prior for cows

Be_dch ~ normal(1,.03);   //Prior for point source delivery
sigma_res ~ normal(0,1000000); //st error of the model
sigma_w ~ normal(0,1000000);     //st. deviation of random effect hyperdistribution
alpha ~ normal(0,sigma_w);    //watershed random effects
sigma_B1 ~ normal(0,.5);  //st. deviation of precipitation coefficient
Bp_mean ~ normal(1,.5);    //mean PIC for hyperdistribution
pic_p[1] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for ag
pic_p[2] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for pre-1980 deve
pic_p[3] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for post-1980 dev
pic_p[4] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for undev
pic_p[5] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for chicken
pic_p[6] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for swine
pic_p[7] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for cow
Sn ~ normal(.2,.08);    //Prior for stream retention rate
Sn2 ~ normal(30,8.5);   //Prior for reservoir retention rate
PIC_q ~ normal(0,1);    //prior for PIC for retention
Be_psi~ uniform (1,8);  //prior for precipitation weighting coef.

ly ~ normal(y_hat,sigma_res);        //parameter that calibrates ly_hat with ly () 
load ~ normal(ly,SD);                // load = WRTDS estimate, SD = WRTDS sd

}

generated quantities {

}
'
Annual_TP <- readRDS("./Annual_TP.rds") #load annual input
#run the annual model (data is the list of data sets. For other parameters, return to the function description)
model_annual = stan(model_code=stanmodelcode_annual, data=Annual_TP, iter=iter, 
             warmup=warmup, thin=thin, chains=3,cores=3,
             control = list(adapt_delta =adapt_delta ,max_treedepth =max_treedepth ))

Summer_TP <- readRDS("./Summer_TP.rds") #load summer input
#run the summer model
model_summer = stan(model_code=stanmodelcode_summer, data=Summer_TP, iter=iter, 
             warmup=warmup, thin=thin, chains=3,cores=3,
             control = list(adapt_delta =adapt_delta ,max_treedepth =max_treedepth ))
