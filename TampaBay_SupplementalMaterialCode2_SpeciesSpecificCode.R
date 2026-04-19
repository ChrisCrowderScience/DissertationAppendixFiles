#Predator Fish Negative Binomial, Logistic, and Zero-Inflated Script
#2/8/25
#Prepared for Supplemental Material for Manuscript by Crowder et al 2025.
#Defining the MCMC Tampa settings######################################
rstan_options(auto_write = TRUE)
options(mc.cores=parallel::detectCores())
ni<-6000
nt<-2 
nb<-3000
nc<-3
N<-14223 #For TB
#Tampa Bay Prey Species Specific Zero-Inflated##################################
#Prey Fish Models 2C and 2G. 
StanTBData <- list(
  Gear1 = FinalTB3Data$GearType1,
  Gear2 = FinalTB3Data$GearType2,
  Gear3 = FinalTB3Data$GearType3,
  YearR = FinalTB3Data$YearR,
  SGCODEDV2 = FinalTB3Data$SGCODEDV2,
  SGCODEDV1 = FinalTB3Data$SGCODEDV1,
  Model2CN = FinalTB3Data$Model2CN,
  SGCODEDV1V2 = FinalTB3Data$SGCODEDV1V2,
  SGCODEDV2V2 = FinalTB3Data$SGCODEDV2V2,
  SGCODEDV3V2 = FinalTB3Data$SGCODEDV3V2,
  SGCODEDV4V2 = FinalTB3Data$SGCODEDV4V2,
  Dim1TB = FinalTB3Data$Dim1TB,
  Dim2TB = FinalTB3Data$Dim2TB,
  SilverPerchNumber = FinalTB3Data$SilverPerchNumber,
  PinfishNumber = FinalTB3Data$PinfishNumber,
  MulletNumber = FinalTB3Data$MulletNumber,
  MojarraNumber = FinalTB3Data$MojarraNumber,
  HerringNumber = FinalTB3Data$HerringNumber,
  AnchovyNumber = FinalTB3Data$AnchovyNumber,
  SilverSidesNumber = FinalTB3Data$SilverSidesNumber,
  PigfishNumber = FinalTB3Data$PigfishNumber
)
#Tampa Bay ZI Herring####
#Model 2C Herring TB
Herring2CTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int HerringNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real b1z;
  real b2z;
  real c2z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2 ~ normal(0,100);
  c2z ~ normal(0.5,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (HerringNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(HerringNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(HerringNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (HerringNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(HerringNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(HerringNumber[i]|lambda[i]));
}
"
Herring2CTBZIStanFile<-write_stan_file(Herring2CTBZI)
print(Herring2CTBZIStanFile)
Herring2CTBZIStanFileResults<-stan(Herring2CTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Model 2G Herring TB
Herring2GTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int HerringNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real c3;
  real c4;
  real b1z;
  real b2z;
  real c2z;
  real c3z;
  real c4z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2z ~ normal(.5,100);
  c3z ~ normal(.5,100);
  c4z ~ normal(.5,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  c2 ~ normal(0,100);
  c3 ~ normal(0,100);
  c4 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (HerringNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(HerringNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(HerringNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (HerringNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(HerringNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(HerringNumber[i]|lambda[i]));
}
"
Herring2GTBZIStanFile<-write_stan_file(Herring2GTBZI)
print(Herring2GTBZIStanFile)
Herring2GTBZIStanFileResults<-stan(Herring2GTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Tampa Bay ZI Anchovy####
#Model 2C Anchovy TB
Anchovy2CTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int AnchovyNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real b1z;
  real b2z;
  real c2z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2 ~ normal(0,100);
  c2z ~ normal(0.5,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (AnchovyNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(AnchovyNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(AnchovyNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (AnchovyNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(AnchovyNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(AnchovyNumber[i]|lambda[i]));
}
"
Anchovy2CTBZIStanFile<-write_stan_file(Anchovy2CTBZI)
print(Anchovy2CTBZIStanFile)
Anchovy2CTBZIStanFileResults<-stan(Anchovy2CTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Model 2G Anchovy TB
Anchovy2GTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int AnchovyNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real c3;
  real c4;
  real b1z;
  real b2z;
  real c2z;
  real c3z;
  real c4z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2z ~ normal(.5,100);
  c3z ~ normal(.5,100);
  c4z ~ normal(.5,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  c2 ~ normal(0,100);
  c3 ~ normal(0,100);
  c4 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (AnchovyNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(AnchovyNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(AnchovyNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (AnchovyNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(AnchovyNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(AnchovyNumber[i]|lambda[i]));
}
"
Anchovy2GTBZIStanFile<-write_stan_file(Anchovy2GTBZI)
print(Anchovy2GTBZIStanFile)
Anchovy2GTBZIStanFileResults<-stan(Anchovy2GTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Tampa Bay ZI Silver Perch####
#Model 2C Silver Perch TB
SilverPerch2CTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int SilverPerchNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real b1z;
  real b2z;
  real c2z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2 ~ normal(0,100);
  c2z ~ normal(0.5,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (SilverPerchNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverPerchNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverPerchNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (SilverPerchNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverPerchNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverPerchNumber[i]|lambda[i]));
}
"
SilverPerch2CTBZIStanFile<-write_stan_file(SilverPerch2CTBZI)
print(SilverPerch2CTBZIStanFile)
SilverPerch2CTBZIStanFileResults<-stan(SilverPerch2CTBZIStanFile,data=StanTBData, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1,init_r = 0.2)
#Model 2G Silver Perch TB
SilverPerch2GTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int SilverPerchNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real c3;
  real c4;
  real b1z;
  real b2z;
  real c2z;
  real c3z;
  real c4z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2z ~ normal(.5,100);
  c3z ~ normal(.5,100);
  c4z ~ normal(.5,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  c2 ~ normal(0,100);
  c3 ~ normal(0,100);
  c4 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (SilverPerchNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverPerchNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverPerchNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (SilverPerchNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverPerchNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverPerchNumber[i]|lambda[i]));
}
"
SilverPerch2GTBZIStanFile<-write_stan_file(SilverPerch2GTBZI)
print(SilverPerch2GTBZIStanFile)
SilverPerch2GTBZIStanFileResults<-stan(SilverPerch2GTBZIStanFile,data=StanTBData, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1,init_r = 0.2)
#Tampa Bay ZI Pinfish####
#Model 2C Pinfish TB
Pinfish2CTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int PinfishNumber;
  array[14223] int MulletNumber;
  array[14223] int MojarraNumber;
  array[14223] int HerringNumber;
  array[14223] int AnchovyNumber;
  array[14223] int SilverSidesNumber;
  array[14223] int PigfishNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real b1z;
  real b2z;
  real c2z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2 ~ normal(0,100);
  c2z ~ normal(0.5,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (PinfishNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PinfishNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PinfishNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (PinfishNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PinfishNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PinfishNumber[i]|lambda[i]));
}
"
Pinfish2CTBZIStanFile<-write_stan_file(Pinfish2CTBZI)
print(Pinfish2CTBZIStanFile)
Pinfish2CTBZIStanFileResults<-stan(Pinfish2CTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Model 2G Pinfish TB
Pinfish2GTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int PinfishNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real c3;
  real c4;
  real b1z;
  real b2z;
  real c2z;
  real c3z;
  real c4z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2z ~ normal(.5,100);
  c3z ~ normal(.5,100);
  c4z ~ normal(.5,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  c2 ~ normal(0,100);
  c3 ~ normal(0,100);
  c4 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (PinfishNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PinfishNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PinfishNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (PinfishNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PinfishNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PinfishNumber[i]|lambda[i]));
}
"
Pinfish2GTBZIStanFile<-write_stan_file(Pinfish2GTBZI)
print(Pinfish2GTBZIStanFile)
Pinfish2GTBZIStanFileResults<-stan(Pinfish2GTBZIStanFile,data=StanTBData,init_r=0.2, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Tampa Bay ZI Mullet####
#Model 2C Mullet TB
Mullet2CTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int PinfishNumber;
  array[14223] int MulletNumber;
  array[14223] int MojarraNumber;
  array[14223] int HerringNumber;
  array[14223] int AnchovyNumber;
  array[14223] int SilverSidesNumber;
  array[14223] int PigfishNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real b1z;
  real b2z;
  real c2z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2 ~ normal(0,100);
  c2z ~ normal(0.5,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  g2 ~ normal(10,10);
  g3 ~ normal(10,10);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (MulletNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MulletNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MulletNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (MulletNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MulletNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MulletNumber[i]|lambda[i]));
}
"
Mullet2CTBZIStanFile<-write_stan_file(Mullet2CTBZI)
print(Mullet2CTBZIStanFile)
Mullet2CTBZIStanFileResults<-stan(Mullet2CTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Model 2G Mullet TB
Mullet2GTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int MulletNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real c3;
  real c4;
  real b1z;
  real b2z;
  real c2z;
  real c3z;
  real c4z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2z ~ normal(.5,100);
  c3z ~ normal(.5,100);
  c4z ~ normal(.5,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  c2 ~ normal(0,100);
  c3 ~ normal(0,100);
  c4 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (MulletNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MulletNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MulletNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (MulletNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MulletNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MulletNumber[i]|lambda[i]));
}
"
Mullet2GTBZIStanFile<-write_stan_file(Mullet2GTBZI)
print(Mullet2GTBZIStanFile)
Mullet2GTBZIStanFileResults<-stan(Mullet2GTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Tampa Bay ZI Mojarra####
#Model 2C Mojarra TB
Mojarra2CTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int PinfishNumber;
  array[14223] int MulletNumber;
  array[14223] int MojarraNumber;
  array[14223] int HerringNumber;
  array[14223] int AnchovyNumber;
  array[14223] int SilverSidesNumber;
  array[14223] int PigfishNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real b1z;
  real b2z;
  real c2z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2 ~ normal(0,100);
  c2z ~ normal(0.5,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (MojarraNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MojarraNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MojarraNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (MojarraNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MojarraNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MojarraNumber[i]|lambda[i]));
}
"
Mojarra2CTBZIStanFile<-write_stan_file(Mojarra2CTBZI)
print(Mojarra2CTBZIStanFile)
Mojarra2CTBZIStanFileResults<-stan(Mojarra2CTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Model 2G Mojarra TB
Mojarra2GTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int MojarraNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real c3;
  real c4;
  real b1z;
  real b2z;
  real c2z;
  real c3z;
  real c4z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2z ~ normal(.5,100);
  c3z ~ normal(.5,100);
  c4z ~ normal(.5,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  c2 ~ normal(0,100);
  c3 ~ normal(0,100);
  c4 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (MojarraNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MojarraNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MojarraNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (MojarraNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MojarraNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(MojarraNumber[i]|lambda[i]));
}
"
Mojarra2GTBZIStanFile<-write_stan_file(Mojarra2GTBZI)
print(Mojarra2GTBZIStanFile)
Mojarra2GTBZIStanFileResults<-stan(Mojarra2GTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Tampa Bay ZI SilverSides####
#Model 2C SilverSides TB
SilverSides2CPFTB<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int PinfishNumber;
  array[14223] int MulletNumber;
  array[14223] int MojarraNumber;
  array[14223] int HerringNumber;
  array[14223] int AnchovyNumber;
  array[14223] int SilverSidesNumber;
  array[14223] int PigfishNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real b1z;
  real b2z;
  real c2z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2 ~ normal(0,100);
  c2z ~ normal(0.5,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (SilverSidesNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverSidesNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverSidesNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (SilverSidesNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverSidesNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverSidesNumber[i]|lambda[i]));
}
"
SilverSides2CPFTBStanFile<-write_stan_file(SilverSides2CPFTB)
print(SilverSides2CPFTBStanFile)
SilverSides2CPFTBStanFileResults<-stan(SilverSides2CPFTBStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Model 2G SilverSides TB
SilverSides2GTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int SilverSidesNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real c3;
  real c4;
  real b1z;
  real b2z;
  real c2z;
  real c3z;
  real c4z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2z ~ normal(.5,100);
  c3z ~ normal(.5,100);
  c4z ~ normal(.5,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  c2 ~ normal(0,100);
  c3 ~ normal(0,100);
  c4 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (SilverSidesNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverSidesNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverSidesNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (SilverSidesNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverSidesNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(SilverSidesNumber[i]|lambda[i]));
}
"
SilverSides2GTBZIStanFile<-write_stan_file(SilverSides2GTBZI)
print(SilverSides2GTBZIStanFile)
SilverSides2GTBZIStanFileResults<-stan(SilverSides2GTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Tampa Bay ZI Pigfish####
#Model 2C Pigfish TB
Pigfish2CTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int PinfishNumber;
  array[14223] int MulletNumber;
  array[14223] int MojarraNumber;
  array[14223] int HerringNumber;
  array[14223] int AnchovyNumber;
  array[14223] int SilverSidesNumber;
  array[14223] int PigfishNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real b1z;
  real b2z;
  real c2z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2 ~ normal(0,100);
  c2z ~ normal(0.5,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (PigfishNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PigfishNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PigfishNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (PigfishNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PigfishNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PigfishNumber[i]|lambda[i]));
}
"
Pigfish2CTBZIStanFile<-write_stan_file(Pigfish2CTBZI)
print(Pigfish2CTBZIStanFile)
Pigfish2CTBZIStanFileResults<-stan(Pigfish2CTBZIStanFile,data=StanTBData, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Model 2G Pigfish TB
Pigfish2GTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int PigfishNumber;
  array[14223] int Gear3;
  array[14223] int Gear2;
  array[14223] int Model2CN;
  array[14223] int YearR;
  array[14223] int SGCODEDV4V2;
  array[14223] int SGCODEDV3V2;
  array[14223] int SGCODEDV2V2;
}
parameters{
  real a;
  real b1;
  real b2;
  real c2;
  real c3;
  real c4;
  real b1z;
  real b2z;
  real c2z;
  real c3z;
  real c4z;
  real z1;
  real g2;
  real g2z;
  real g3;
  real g3z;
  vector[9] HC;
  real<lower=0> sigmaHC;
  vector[12] Yv;
  real<lower=0> sigmaYear;
  vector[9] HCz;
  real<lower=0> sigmaHCz;
  vector[12] Yvz;
  real<lower=0> sigmaYearz;
}
model{
  vector[14223] lambda;
  vector[14223] pbar;
  sigmaYear ~ cauchy( 0 , 50 );
  Yv ~ normal( 0 , sigmaYear );
  sigmaHC ~ cauchy( 0 , 50 );
  HC ~ normal( 0 , sigmaHC );
  sigmaYearz ~ cauchy( 0 , 50 );
  Yvz ~ normal( 0 , sigmaYearz );
  sigmaHCz ~ cauchy( 0 , 50 );
  HCz ~ normal( 0 , sigmaHCz );
  c2z ~ normal(.5,100);
  c3z ~ normal(.5,100);
  c4z ~ normal(.5,100);
  b1z ~ normal(0,10);
  b2z ~ normal(0,10);
  c2 ~ normal(0,100);
  c3 ~ normal(0,100);
  c4 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  g2 ~ normal(10,100);
  g3 ~ normal(10,100);
  g2z ~ normal(0.5,100);
  g3z ~ normal(0.5,100);
  a ~ normal( 5 , 100 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (PigfishNumber[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PigfishNumber[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PigfishNumber[i]|lambda[i]));
}
generated quantities{
  vector[14223] lambda;
  vector[14223] pbar;
  real dev;
  dev = 0;
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (PigfishNumber[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PigfishNumber[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PigfishNumber[i]|lambda[i]));
}
"
Pigfish2GTBZIStanFile<-write_stan_file(Pigfish2GTBZI)
print(Pigfish2GTBZIStanFile)
Pigfish2GTBZIStanFileResults<-stan(Pigfish2GTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)

#Tampa Bay Prey Species Specific Logistic######################################
FinalTB3Data$HerringPresAbs<-ifelse(FinalTB3Data$HerringNumber==0,0,1)
FinalTB3Data$AnchovyPresAbs<-ifelse(FinalTB3Data$AnchovyNumber==0,0,1)
FinalTB3Data$MulletPresAbs<-ifelse(FinalTB3Data$MulletNumber==0,0,1)
FinalTB3Data$MojarraPresAbs<-ifelse(FinalTB3Data$MojarraNumber==0,0,1)
FinalTB3Data$SilverPerchPresAbs<-ifelse(FinalTB3Data$SilverPerchNumber==0,0,1)
FinalTB3Data$PinfishPresAbs<-ifelse(FinalTB3Data$PinfishNumber==0,0,1)
FinalTB3Data$PigfishPresAbs<-ifelse(FinalTB3Data$PigfishNumber==0,0,1)
FinalTB3Data$SilverSidesPresAbs<-ifelse(FinalTB3Data$SilverSidesNumber==0,0,1)
LogStanTBData <- list(
  Gear1 = FinalTB3Data$Gear1,
  Gear2 = FinalTB3Data$Gear2,
  Gear3 = FinalTB3Data$Gear3,
  YearR = FinalTB3Data$YearR,
  SGCODEDV2 = FinalTB3Data$SGCODEDV2,
  SGCODEDV1 = FinalTB3Data$SGCODEDV1,
  Model2CN = FinalTB3Data$Model2CN,
  SGCODEDV1V2 = FinalTB3Data$SGCODEDV1V2,
  SGCODEDV2V2 = FinalTB3Data$SGCODEDV2V2,
  SGCODEDV3V2 = FinalTB3Data$SGCODEDV3V2,
  SGCODEDV4V2 = FinalTB3Data$SGCODEDV4V2,
  Dim1TB = FinalTB3Data$Dim1TB,
  Dim2TB = FinalTB3Data$Dim2TB,
  SilverPerchPresAbs = FinalTB3Data$SilverPerchPresAbs,
  PinfishPresAbs = FinalTB3Data$PinfishPresAbs,
  MulletPresAbs = FinalTB3Data$MulletPresAbs,
  MojarraPresAbs = FinalTB3Data$MojarraPresAbs,
  HerringPresAbs = FinalTB3Data$HerringPresAbs,
  AnchovyPresAbs = FinalTB3Data$AnchovyPresAbs,
  SilverSidesPresAbs = FinalTB3Data$SilverSidesPresAbs,
  PigfishPresAbs= FinalTB3Data$PigfishPresAbs
)


#Herring TB Log 2C####
Herring2CTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10, 100 );
    g2 ~ normal( 10, 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    HerringPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( HerringPresAbs[i] | 1 , p[i] );
}"

Herring2CTBLogStanFile<-write_stan_file(Herring2CTBLog)
print(Herring2CTBLogStanFile)
Herring2CTBLogStanFileResults<-stan(Herring2CTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Herring TB Log 2G####
Herring2GTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    HerringPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( HerringPresAbs[i] | 1 , p[i] );
}"

Herring2GTBLogStanFile<-write_stan_file(Herring2GTBLog)
print(Herring2GTBLogStanFile)
Herring2GTBLogStanFileResults<-stan(Herring2GTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Anchovy TB Log 2C####
Anchovy2CTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10, 100 );
    g2 ~ normal( 10, 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    AnchovyPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( AnchovyPresAbs[i] | 1 , p[i] );
}"

Anchovy2CTBLogStanFile<-write_stan_file(Anchovy2CTBLog)
print(Anchovy2CTBLogStanFile)
Anchovy2CTBLogStanFileResults<-stan(Anchovy2CTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Anchovy TB Log 2G####
Anchovy2GTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    AnchovyPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( AnchovyPresAbs[i] | 1 , p[i] );
}"

Anchovy2GTBLogStanFile<-write_stan_file(Anchovy2GTBLog)
print(Anchovy2GTBLogStanFile)
Anchovy2GTBLogStanFileResults<-stan(Anchovy2GTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Mullet TB Log 2C####
Mullet2CTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10, 100 );
    g2 ~ normal( 10, 100 );
    c2 ~ normal( 0 , 10 );
    b2 ~ normal( 0 , 10 );
    b1 ~ normal( 0 , 10 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    MulletPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( MulletPresAbs[i] | 1 , p[i] );
}"

Mullet2CTBLogStanFile<-write_stan_file(Mullet2CTBLog)
print(Mullet2CTBLogStanFile)
Mullet2CTBLogStanFileResults<-stan(Mullet2CTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Mullet TB Log 2G####
Mullet2GTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10, 100 );
    c4 ~ normal( 0 , 10 );
    c3 ~ normal( 0 , 10 );
    c2 ~ normal( 0 , 10 );
    b2 ~ normal( 0 , 10 );
    b1 ~ normal( 0 , 10 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    MulletPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( MulletPresAbs[i] | 1 , p[i] );
}"

Mullet2GTBLogStanFile<-write_stan_file(Mullet2GTBLog)
print(Mullet2GTBLogStanFile)
Mullet2GTBLogStanFileResults<-stan(Mullet2GTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Mojarra TB Log 2C####
Mojarra2CTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10, 100 );
    g2 ~ normal( 10, 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    MojarraPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( MojarraPresAbs[i] | 1 , p[i] );
}"

Mojarra2CTBLogStanFile<-write_stan_file(Mojarra2CTBLog)
print(Mojarra2CTBLogStanFile)
Mojarra2CTBLogStanFileResults<-stan(Mojarra2CTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Mojarra TB Log 2G####
Mojarra2GTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    MojarraPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( MojarraPresAbs[i] | 1 , p[i] );
}"

Mojarra2GTBLogStanFile<-write_stan_file(Mojarra2GTBLog)
print(Mojarra2GTBLogStanFile)
Mojarra2GTBLogStanFileResults<-stan(Mojarra2GTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#SilverPerch TB Log 2C####
SilverPerch2CTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10, 100 );
    g2 ~ normal( 10, 100 );
    c2 ~ normal( 0 , 10 );
    b2 ~ normal( 0 , 10 );
    b1 ~ normal( 0 , 10 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    SilverPerchPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( SilverPerchPresAbs[i] | 1 , p[i] );
}"

SilverPerch2CTBLogStanFile<-write_stan_file(SilverPerch2CTBLog)
print(SilverPerch2CTBLogStanFile)
SilverPerch2CTBLogStanFileResults<-stan(SilverPerch2CTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#SilverPerch TB Log 2G####
SilverPerch2GTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    SilverPerchPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( SilverPerchPresAbs[i] | 1 , p[i] );
}"

SilverPerch2GTBLogStanFile<-write_stan_file(SilverPerch2GTBLog)
print(SilverPerch2GTBLogStanFile)
SilverPerch2GTBLogStanFileResults<-stan(SilverPerch2GTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Pinfish TB Log 2C####
Pinfish2CTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10, 100 );
    g2 ~ normal( 10, 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    PinfishPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( PinfishPresAbs[i] | 1 , p[i] );
}"

Pinfish2CTBLogStanFile<-write_stan_file(Pinfish2CTBLog)
print(Pinfish2CTBLogStanFile)
Pinfish2CTBLogStanFileResults<-stan(Pinfish2CTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Pinfish TB Log 2G####
Pinfish2GTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    PinfishPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( PinfishPresAbs[i] | 1 , p[i] );
}"

Pinfish2GTBLogStanFile<-write_stan_file(Pinfish2GTBLog)
print(Pinfish2GTBLogStanFile)
Pinfish2GTBLogStanFileResults<-stan(Pinfish2GTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Pigfish TB Log 2C####
Pigfish2CTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10, 100 );
    g2 ~ normal( 10, 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    PigfishPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( PigfishPresAbs[i] | 1 , p[i] );
}"

Pigfish2CTBLogStanFile<-write_stan_file(Pigfish2CTBLog)
print(Pigfish2CTBLogStanFile)
Pigfish2CTBLogStanFileResults<-stan(Pigfish2CTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Pigfish TB Log 2G####
Pigfish2GTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    PigfishPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( PigfishPresAbs[i] | 1 , p[i] );
}"

Pigfish2GTBLogStanFile<-write_stan_file(Pigfish2GTBLog)
print(Pigfish2GTBLogStanFile)
Pigfish2GTBLogStanFileResults<-stan(Pigfish2GTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#SilverSides TB Log 2C####
SilverSides2CTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10, 100 );
    g2 ~ normal( 10, 100 );
    c2 ~ normal( 0 , 10 );
    b2 ~ normal( 0 , 10 );
    b1 ~ normal( 0 , 10 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    SilverSidesPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( SilverSidesPresAbs[i] | 1 , p[i] );
}"

SilverSides2CTBLogStanFile<-write_stan_file(SilverSides2CTBLog)
print(SilverSides2CTBLogStanFile)
SilverSides2CTBLogStanFileResults<-stan(SilverSides2CTBLogStanFile,data=LogStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#SilverSides TB Log 2G####
SilverSides2GTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SilverSidesPresAbs;
    array[14223] int PigfishPresAbs;
    array[14223] int PinfishPresAbs;
    array[14223] int SilverPerchPresAbs;
    array[14223] int MojarraPresAbs;
    array[14223] int AnchovyPresAbs;
    array[14223] int MulletPresAbs;
    array[14223] int HerringPresAbs;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int Gear3;
    array[14223] int Gear2;
    array[14223] int Model2CN;
    array[14223] int YearR;
    array[14223] int SGCODEDV2;
     vector[14223] Dim2TB;
     vector[14223] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
}
model{
     vector[14223] p;
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 10 );
    c3 ~ normal( 0 , 10 );
    c2 ~ normal( 0 , 10 );
    b2 ~ normal( 0 , 10 );
    b1 ~ normal( 0 , 10 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    SilverSidesPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( SilverSidesPresAbs[i] | 1 , p[i] );
}"

SilverSides2GTBLogStanFile<-write_stan_file(SilverSides2GTBLog)
print(SilverSides2GTBLogStanFile)
SilverSides2GTBLogStanFileResults<-stan(SilverSides2GTBLogStanFile,data=LogStanTBData, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1,init_r = 0.2)

#Tampa Bay Prey Species Specific Neg Binom#####################################
NonZeroHerringDf<-subset(FinalTB3Data,HerringNumber>0)
HerringStanTBData <- list(
  Gear1 = NonZeroHerringDf$Gear1,
  Gear2 = NonZeroHerringDf$Gear2,
  Gear3 = NonZeroHerringDf$Gear3,
  YearR = NonZeroHerringDf$YearR,
  SGCODEDV2 = NonZeroHerringDf$SGCODEDV2,
  SGCODEDV1 = NonZeroHerringDf$SGCODEDV1,
  Model2CN = NonZeroHerringDf$Model2CN,
  SGCODEDV1V2 = NonZeroHerringDf$SGCODEDV1V2,
  SGCODEDV2V2 = NonZeroHerringDf$SGCODEDV2V2,
  SGCODEDV3V2 = NonZeroHerringDf$SGCODEDV3V2,
  SGCODEDV4V2 = NonZeroHerringDf$SGCODEDV4V2,
  Dim1TB = NonZeroHerringDf$Dim1TB,
  Dim2TB = NonZeroHerringDf$Dim2TB,
  HerringNumber = NonZeroHerringDf$HerringNumber
)


NonZeroAnchovyDf<-subset(FinalTB3Data,AnchovyNumber>0)
AnchovyStanTBData <- list(
  Gear1 = NonZeroAnchovyDf$Gear1,
  Gear2 = NonZeroAnchovyDf$Gear2,
  Gear3 = NonZeroAnchovyDf$Gear3,
  YearR = NonZeroAnchovyDf$YearR,
  SGCODEDV2 = NonZeroAnchovyDf$SGCODEDV2,
  SGCODEDV1 = NonZeroAnchovyDf$SGCODEDV1,
  Model2CN = NonZeroAnchovyDf$Model2CN,
  SGCODEDV1V2 = NonZeroAnchovyDf$SGCODEDV1V2,
  SGCODEDV2V2 = NonZeroAnchovyDf$SGCODEDV2V2,
  SGCODEDV3V2 = NonZeroAnchovyDf$SGCODEDV3V2,
  SGCODEDV4V2 = NonZeroAnchovyDf$SGCODEDV4V2,
  Dim1TB = NonZeroAnchovyDf$Dim1TB,
  Dim2TB = NonZeroAnchovyDf$Dim2TB,
  AnchovyNumber = NonZeroAnchovyDf$AnchovyNumber
)
NonZeroMulletDf<-subset(FinalTB3Data,MulletNumber>0)
MulletStanTBData <- list(
  Gear1 = NonZeroMulletDf$Gear1,
  Gear2 = NonZeroMulletDf$Gear2,
  Gear3 = NonZeroMulletDf$Gear3,
  YearR = NonZeroMulletDf$YearR,
  SGCODEDV2 = NonZeroMulletDf$SGCODEDV2,
  SGCODEDV1 = NonZeroMulletDf$SGCODEDV1,
  Model2CN = NonZeroMulletDf$Model2CN,
  SGCODEDV1V2 = NonZeroMulletDf$SGCODEDV1V2,
  SGCODEDV2V2 = NonZeroMulletDf$SGCODEDV2V2,
  SGCODEDV3V2 = NonZeroMulletDf$SGCODEDV3V2,
  SGCODEDV4V2 = NonZeroMulletDf$SGCODEDV4V2,
  Dim1TB = NonZeroMulletDf$Dim1TB,
  Dim2TB = NonZeroMulletDf$Dim2TB,
  MulletNumber = NonZeroMulletDf$MulletNumber
)
NonZeroMojarraDf<-subset(FinalTB3Data,MojarraNumber>0)
MojarraStanTBData <- list(
  Gear1 = NonZeroMojarraDf$Gear1,
  Gear2 = NonZeroMojarraDf$Gear2,
  Gear3 = NonZeroMojarraDf$Gear3,
  YearR = NonZeroMojarraDf$YearR,
  SGCODEDV2 = NonZeroMojarraDf$SGCODEDV2,
  SGCODEDV1 = NonZeroMojarraDf$SGCODEDV1,
  Model2CN = NonZeroMojarraDf$Model2CN,
  SGCODEDV1V2 = NonZeroMojarraDf$SGCODEDV1V2,
  SGCODEDV2V2 = NonZeroMojarraDf$SGCODEDV2V2,
  SGCODEDV3V2 = NonZeroMojarraDf$SGCODEDV3V2,
  SGCODEDV4V2 = NonZeroMojarraDf$SGCODEDV4V2,
  Dim1TB = NonZeroMojarraDf$Dim1TB,
  Dim2TB = NonZeroMojarraDf$Dim2TB,
  MojarraNumber = NonZeroMojarraDf$MojarraNumber
)
NonZeroSilverPerchDf<-subset(FinalTB3Data,SilverPerchNumber>0)
SilverPerchStanTBData <- list(
  Gear1 = NonZeroSilverPerchDf$Gear1,
  Gear2 = NonZeroSilverPerchDf$Gear2,
  Gear3 = NonZeroSilverPerchDf$Gear3,
  YearR = NonZeroSilverPerchDf$YearR,
  SGCODEDV2 = NonZeroSilverPerchDf$SGCODEDV2,
  SGCODEDV1 = NonZeroSilverPerchDf$SGCODEDV1,
  Model2CN = NonZeroSilverPerchDf$Model2CN,
  SGCODEDV1V2 = NonZeroSilverPerchDf$SGCODEDV1V2,
  SGCODEDV2V2 = NonZeroSilverPerchDf$SGCODEDV2V2,
  SGCODEDV3V2 = NonZeroSilverPerchDf$SGCODEDV3V2,
  SGCODEDV4V2 = NonZeroSilverPerchDf$SGCODEDV4V2,
  Dim1TB = NonZeroSilverPerchDf$Dim1TB,
  Dim2TB = NonZeroSilverPerchDf$Dim2TB,
  SilverPerchNumber = NonZeroSilverPerchDf$SilverPerchNumber
)
NonZeroPinfishDf<-subset(FinalTB3Data,PinfishNumber>0)
PinfishStanTBData <- list(
  Gear1 = NonZeroPinfishDf$Gear1,
  Gear2 = NonZeroPinfishDf$Gear2,
  Gear3 = NonZeroPinfishDf$Gear3,
  YearR = NonZeroPinfishDf$YearR,
  SGCODEDV2 = NonZeroPinfishDf$SGCODEDV2,
  SGCODEDV1 = NonZeroPinfishDf$SGCODEDV1,
  Model2CN = NonZeroPinfishDf$Model2CN,
  SGCODEDV1V2 = NonZeroPinfishDf$SGCODEDV1V2,
  SGCODEDV2V2 = NonZeroPinfishDf$SGCODEDV2V2,
  SGCODEDV3V2 = NonZeroPinfishDf$SGCODEDV3V2,
  SGCODEDV4V2 = NonZeroPinfishDf$SGCODEDV4V2,
  Dim1TB = NonZeroPinfishDf$Dim1TB,
  Dim2TB = NonZeroPinfishDf$Dim2TB,
  PinfishNumber = NonZeroPinfishDf$PinfishNumber
)
NonZeroPigfishDf<-subset(FinalTB3Data,PigfishNumber>0)
PigfishStanTBData <- list(
  Gear1 = NonZeroPigfishDf$Gear1,
  Gear2 = NonZeroPigfishDf$Gear2,
  Gear3 = NonZeroPigfishDf$Gear3,
  YearR = NonZeroPigfishDf$YearR,
  SGCODEDV2 = NonZeroPigfishDf$SGCODEDV2,
  SGCODEDV1 = NonZeroPigfishDf$SGCODEDV1,
  Model2CN = NonZeroPigfishDf$Model2CN,
  SGCODEDV1V2 = NonZeroPigfishDf$SGCODEDV1V2,
  SGCODEDV2V2 = NonZeroPigfishDf$SGCODEDV2V2,
  SGCODEDV3V2 = NonZeroPigfishDf$SGCODEDV3V2,
  SGCODEDV4V2 = NonZeroPigfishDf$SGCODEDV4V2,
  Dim1TB = NonZeroPigfishDf$Dim1TB,
  Dim2TB = NonZeroPigfishDf$Dim2TB,
  PigfishNumber = NonZeroPigfishDf$PigfishNumber
)
NonZeroSilverSidesDf<-subset(FinalTB3Data,SilverSidesNumber>0)
SilverSidesStanTBData <- list(
  Gear1 = NonZeroSilverSidesDf$Gear1,
  Gear2 = NonZeroSilverSidesDf$Gear2,
  Gear3 = NonZeroSilverSidesDf$Gear3,
  YearR = NonZeroSilverSidesDf$YearR,
  SGCODEDV2 = NonZeroSilverSidesDf$SGCODEDV2,
  SGCODEDV1 = NonZeroSilverSidesDf$SGCODEDV1,
  Model2CN = NonZeroSilverSidesDf$Model2CN,
  SGCODEDV1V2 = NonZeroSilverSidesDf$SGCODEDV1V2,
  SGCODEDV2V2 = NonZeroSilverSidesDf$SGCODEDV2V2,
  SGCODEDV3V2 = NonZeroSilverSidesDf$SGCODEDV3V2,
  SGCODEDV4V2 = NonZeroSilverSidesDf$SGCODEDV4V2,
  Dim1TB = NonZeroSilverSidesDf$Dim1TB,
  Dim2TB = NonZeroSilverSidesDf$Dim2TB,
  SilverSidesNumber = NonZeroSilverSidesDf$SilverSidesNumber
)
#Herring TB NB 2C####
N<-399
Herring2CTBNB<-"data{
    array[399] int Gear1;
    array[399] int SGCODEDV1;
    array[399] int HerringNumber;
    array[399] int Gear3;
    array[399] int Gear2;
    array[399] int Model2CN;
    array[399] int YearR;
    array[399] int SGCODEDV2;
    array[399] int SGCODEDV4V2;
    array[399] int SGCODEDV1V2;
    array[399] int SGCODEDV2V2;
    array[399] int SGCODEDV3V2;
     vector[399] Dim2TB;
     vector[399] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[399] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:399 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    HerringNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[399] log_lik;
     vector[399] pbar;
    for ( i in 1:399 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:399 ) log_lik[i] = neg_binomial_2_lpmf( HerringNumber[i] | pbar[i] , scale );
}"

Herring2CTBNBStanFile<-write_stan_file(Herring2CTBNB)
print(Herring2CTBNBStanFile)
Herring2CTBNBStanFileResults<-stan(Herring2CTBNBStanFile,data=HerringStanData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Herring TB NB 2G####
Herring2GTBNB<-"data{
    array[399] int HerringNumber;
    array[399] int Gear1;
    array[399] int SGCODEDV2;
    array[399] int SGCODEDV1;
    array[399] int Gear3;
    array[399] int Gear2;
    array[399] int Model2CN;
    array[399] int YearR;
    array[399] int SGCODEDV4V2;
    array[399] int SGCODEDV3V2;
    array[399] int SGCODEDV2V2;
     vector[399] Dim2TB;
     vector[399] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[399] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:399 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    HerringNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[399] log_lik;
     vector[399] pbar;
    for ( i in 1:399 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:399 ) log_lik[i] = neg_binomial_2_lpmf( HerringNumber[i] | pbar[i] , scale );
}"

Herring2GTBNBStanFile<-write_stan_file(Herring2GTBNB)
print(Herring2GTBNBStanFile)
Herring2GTBNBStanFileResults<-stan(Herring2GTBNBStanFile,data=HerringStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Anchovy TB NB 2C####
N<-2441
Anchovy2CTBNB<-"data{
     array[2441]  int AnchovyNumber;
    array[2441] int Gear1;
    array[2441] int SGCODEDV2;
    array[2441] int SGCODEDV1;
    array[2441] int SGCODEDV1V2;
    array[2441] int Gear3;
    array[2441] int Gear2;
    array[2441] int Model2CN;
    array[2441] int YearR;
    array[2441] int SGCODEDV4V2;
    array[2441] int SGCODEDV3V2;
    array[2441] int SGCODEDV2V2;
     vector[2441] Dim2TB;
     vector[2441] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[2441] pbar;
    scale ~ cauchy( 0 , 10 );
    sigmaYear ~ cauchy( 0 , 10 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 10 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 50 );
    for ( i in 1:2441 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    AnchovyNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[2441] log_lik;
     vector[2441] pbar;
    for ( i in 1:2441 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:2441 ) log_lik[i] = neg_binomial_2_lpmf( AnchovyNumber[i] | pbar[i] , scale );
}"

Anchovy2CTBNBStanFile<-write_stan_file(Anchovy2CTBNB)
print(Anchovy2CTBNBStanFile)
Anchovy2CTBNBStanFileResults<-stan(Anchovy2CTBNBStanFile,data=AnchovyStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Anchovy TB NB 2G####
Anchovy2GTBNB<-"data{
    array[2441] int AnchovyNumber;
    array[2441] int Gear1;
    array[2441] int SGCODEDV2;
    array[2441] int SGCODEDV1;
    array[2441] int SGCODEDV1V2;
    array[2441] int Gear3;
    array[2441] int Gear2;
    array[2441] int Model2CN;
    array[2441] int YearR;
    array[2441] int SGCODEDV4V2;
    array[2441] int SGCODEDV3V2;
    array[2441] int SGCODEDV2V2;
     vector[2441] Dim2TB;
     vector[2441] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[2441] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:2441 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    AnchovyNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[2441] log_lik;
     vector[2441] pbar;
    for ( i in 1:2441 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:2441 ) log_lik[i] = neg_binomial_2_lpmf( AnchovyNumber[i] | pbar[i] , scale );
}"

Anchovy2GTBNBStanFile<-write_stan_file(Anchovy2GTBNB)
print(Anchovy2GTBNBStanFile)
Anchovy2GTBNBStanFileResults<-stan(Anchovy2GTBNBStanFile,data=AnchovyStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Mullet TB NB 2C####
N<-1974
Mullet2CTBNB<-"data{
    array[1974] int MulletNumber;
    array[1974] int Gear1;
    array[1974] int SGCODEDV2;
    array[1974] int SGCODEDV1;
    array[1974] int SGCODEDV1V2;
    array[1974] int Gear3;
    array[1974] int Gear2;
    array[1974] int Model2CN;
    array[1974] int YearR;
    array[1974] int SGCODEDV4V2;
    array[1974] int SGCODEDV3V2;
    array[1974] int SGCODEDV2V2;
     vector[1974] Dim2TB;
     vector[1974] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[1974] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:1974 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    MulletNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[1974] log_lik;
     vector[1974] pbar;
    for ( i in 1:1974 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:1974 ) log_lik[i] = neg_binomial_2_lpmf( MulletNumber[i] | pbar[i] , scale );
}"

Mullet2CTBNBStanFile<-write_stan_file(Mullet2CTBNB)
print(Mullet2CTBNBStanFile)
Mullet2CTBNBStanFileResults<-stan(Mullet2CTBNBStanFile,data=MulletStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Mullet TB NB 2G####
Mullet2GTBNB<-"data{
    array[1974] int MulletNumber;
    array[1974] int Gear1;
    array[1974] int SGCODEDV2;
    array[1974] int SGCODEDV1;
    array[1974] int SGCODEDV1V2;
    array[1974] int Gear3;
    array[1974] int Gear2;
    array[1974] int Model2CN;
    array[1974] int YearR;
    array[1974] int SGCODEDV4V2;
    array[1974] int SGCODEDV3V2;
    array[1974] int SGCODEDV2V2;
     vector[1974] Dim2TB;
     vector[1974] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[1974] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:1974 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    MulletNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[1974] log_lik;
     vector[1974] pbar;
    for ( i in 1:1974 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:1974 ) log_lik[i] = neg_binomial_2_lpmf( MulletNumber[i] | pbar[i] , scale );
}"

Mullet2GTBNBStanFile<-write_stan_file(Mullet2GTBNB)
print(Mullet2GTBNBStanFile)
Mullet2GTBNBStanFileResults<-stan(Mullet2GTBNBStanFile,data=MulletStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Mojarra TB NB 2C####
N<-6006
Mojarra2CTBNB<-"data{
    array[6006] int MojarraNumber;
    array[6006] int Gear1;
    array[6006] int SGCODEDV2;
    array[6006] int SGCODEDV1;
    array[6006] int SGCODEDV1V2;
    array[6006] int Gear3;
    array[6006] int Gear2;
    array[6006] int Model2CN;
    array[6006] int YearR;
    array[6006] int SGCODEDV4V2;
    array[6006] int SGCODEDV3V2;
    array[6006] int SGCODEDV2V2;
     vector[6006] Dim2TB;
     vector[6006] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[6006] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:6006 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    MojarraNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[6006] log_lik;
     vector[6006] pbar;
    for ( i in 1:6006 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:6006 ) log_lik[i] = neg_binomial_2_lpmf( MojarraNumber[i] | pbar[i] , scale );
}"

Mojarra2CTBNBStanFile<-write_stan_file(Mojarra2CTBNB)
print(Mojarra2CTBNBStanFile)
Mojarra2CTBNBStanFileResults<-stan(Mojarra2CTBNBStanFile,data=MojarraStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Mojarra TB NB 2G####
Mojarra2GTBNB<-"data{
    array[6006] int MojarraNumber;
    array[6006] int Gear1;
    array[6006] int SGCODEDV2;
    array[6006] int SGCODEDV1;
    array[6006] int SGCODEDV1V2;
    array[6006] int Gear3;
    array[6006] int Gear2;
    array[6006] int Model2CN;
    array[6006] int YearR;
    array[6006] int SGCODEDV4V2;
    array[6006] int SGCODEDV3V2;
    array[6006] int SGCODEDV2V2;
     vector[6006] Dim2TB;
     vector[6006] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[6006] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:6006 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    MojarraNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[6006] log_lik;
     vector[6006] pbar;
    for ( i in 1:6006 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:6006 ) log_lik[i] = neg_binomial_2_lpmf( MojarraNumber[i] | pbar[i] , scale );
}"

Mojarra2GTBNBStanFile<-write_stan_file(Mojarra2GTBNB)
print(Mojarra2GTBNBStanFile)
Mojarra2GTBNBStanFileResults<-stan(Mojarra2GTBNBStanFile,data=MojarraStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#SilverPerch TB NB 2C####
N<-2168
SilverPerch2CTBNB<-"data{
    array[2168] int SilverPerchNumber;
    array[2168] int Gear1;
    array[2168] int SGCODEDV2;
    array[2168] int SGCODEDV1;
    array[2168] int SGCODEDV1V2;
    array[2168] int Gear3;
    array[2168] int Gear2;
    array[2168] int Model2CN;
    array[2168] int YearR;
    array[2168] int SGCODEDV4V2;
    array[2168] int SGCODEDV3V2;
    array[2168] int SGCODEDV2V2;
     vector[2168] Dim2TB;
     vector[2168] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[2168] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:2168 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    SilverPerchNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[2168] log_lik;
     vector[2168] pbar;
    for ( i in 1:2168 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:2168 ) log_lik[i] = neg_binomial_2_lpmf( SilverPerchNumber[i] | pbar[i] , scale );
}"

SilverPerch2CTBNBStanFile<-write_stan_file(SilverPerch2CTBNB)
print(SilverPerch2CTBNBStanFile)
SilverPerch2CTBNBStanFileResults<-stan(SilverPerch2CTBNBStanFile,data=SilverPerchStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#SilverPerch TB NB 2G####
SilverPerch2GTBNB<-"data{
    array[2168] int SilverPerchNumber;
    array[2168] int Gear1;
    array[2168] int SGCODEDV2;
    array[2168] int SGCODEDV1;
    array[2168] int SGCODEDV1V2;
    array[2168] int Gear3;
    array[2168] int Gear2;
    array[2168] int Model2CN;
    array[2168] int YearR;
    array[2168] int SGCODEDV4V2;
    array[2168] int SGCODEDV3V2;
    array[2168] int SGCODEDV2V2;
     vector[2168] Dim2TB;
     vector[2168] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[2168] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:2168 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    SilverPerchNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[2168] log_lik;
     vector[2168] pbar;
    for ( i in 1:2168 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:2168 ) log_lik[i] = neg_binomial_2_lpmf( SilverPerchNumber[i] | pbar[i] , scale );
}"

SilverPerch2GTBNBStanFile<-write_stan_file(SilverPerch2GTBNB)
print(SilverPerch2GTBNBStanFile)
SilverPerch2GTBNBStanFileResults<-stan(SilverPerch2GTBNBStanFile,data=SilverPerchStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Pinfish TB NB 2C####
N<-6538
Pinfish2CTBNB<-"data{
    array[6538] int PinfishNumber;
    array[6538] int Gear1;
    array[6538] int SGCODEDV2;
    array[6538] int SGCODEDV1;
    array[6538] int SGCODEDV1V2;
    array[6538] int Gear3;
    array[6538] int Gear2;
    array[6538] int Model2CN;
    array[6538] int YearR;
    array[6538] int SGCODEDV4V2;
    array[6538] int SGCODEDV3V2;
    array[6538] int SGCODEDV2V2;
     vector[6538] Dim2TB;
     vector[6538] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[6538] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:6538 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    PinfishNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[6538] log_lik;
     vector[6538] pbar;
    for ( i in 1:6538 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:6538 ) log_lik[i] = neg_binomial_2_lpmf( PinfishNumber[i] | pbar[i] , scale );
}"

Pinfish2CTBNBStanFile<-write_stan_file(Pinfish2CTBNB)
print(Pinfish2CTBNBStanFile)
Pinfish2CTBNBStanFileResults<-stan(Pinfish2CTBNBStanFile,data=PinfishStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Pinfish TB NB 2G####
Pinfish2GTBNB<-"data{
    array[6538] int PinfishNumber;
    array[6538] int Gear1;
    array[6538] int SGCODEDV2;
    array[6538] int SGCODEDV1;
    array[6538] int SGCODEDV1V2;
    array[6538] int Gear3;
    array[6538] int Gear2;
    array[6538] int Model2CN;
    array[6538] int YearR;
    array[6538] int SGCODEDV4V2;
    array[6538] int SGCODEDV3V2;
    array[6538] int SGCODEDV2V2;
     vector[6538] Dim2TB;
     vector[6538] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[6538] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:6538 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    PinfishNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[6538] log_lik;
     vector[6538] pbar;
    for ( i in 1:6538 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:6538 ) log_lik[i] = neg_binomial_2_lpmf( PinfishNumber[i] | pbar[i] , scale );
}"

Pinfish2GTBNBStanFile<-write_stan_file(Pinfish2GTBNB)
print(Pinfish2GTBNBStanFile)
Pinfish2GTBNBStanFileResults<-stan(Pinfish2GTBNBStanFile,data=PinfishStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Pigfish TB NB 2C####
N<-2705
Pigfish2CTBNB<-"data{
    array[2705] int PigfishNumber;
    array[2705] int Gear1;
    array[2705] int SGCODEDV2;
    array[2705] int SGCODEDV1;
    array[2705] int SGCODEDV1V2;
    array[2705] int Gear3;
    array[2705] int Gear2;
    array[2705] int Model2CN;
    array[2705] int YearR;
    array[2705] int SGCODEDV4V2;
    array[2705] int SGCODEDV3V2;
    array[2705] int SGCODEDV2V2;
     vector[2705] Dim2TB;
     vector[2705] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[2705] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:2705 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    PigfishNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[2705] log_lik;
     vector[2705] pbar;
    for ( i in 1:2705 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:2705 ) log_lik[i] = neg_binomial_2_lpmf( PigfishNumber[i] | pbar[i] , scale );
}"

Pigfish2CTBNBStanFile<-write_stan_file(Pigfish2CTBNB)
print(Pigfish2CTBNBStanFile)
Pigfish2CTBNBStanFileResults<-stan(Pigfish2CTBNBStanFile,data=PigfishStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#Pigfish TB NB 2G####
Pigfish2GTBNB<-"data{
    array[2705] int PigfishNumber;
    array[2705] int Gear1;
    array[2705] int SGCODEDV2;
    array[2705] int SGCODEDV1;
    array[2705] int SGCODEDV1V2;
    array[2705] int Gear3;
    array[2705] int Gear2;
    array[2705] int Model2CN;
    array[2705] int YearR;
    array[2705] int SGCODEDV4V2;
    array[2705] int SGCODEDV3V2;
    array[2705] int SGCODEDV2V2;
     vector[2705] Dim2TB;
     vector[2705] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[2705] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:2705 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    PigfishNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[2705] log_lik;
     vector[2705] pbar;
    for ( i in 1:2705 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:2705 ) log_lik[i] = neg_binomial_2_lpmf( PigfishNumber[i] | pbar[i] , scale );
}"

Pigfish2GTBNBStanFile<-write_stan_file(Pigfish2GTBNB)
print(Pigfish2GTBNBStanFile)
Pigfish2GTBNBStanFileResults<-stan(Pigfish2GTBNBStanFile,data=PigfishStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#SilverSides TB NB 2C####
N<-759
SilverSides2CTBNB<-"data{
    array[759] int SilverSidesNumber;
    array[759] int Gear1;
    array[759] int SGCODEDV2;
    array[759] int SGCODEDV1;
    array[759] int SGCODEDV1V2;
    array[759] int Gear3;
    array[759] int Gear2;
    array[759] int Model2CN;
    array[759] int YearR;
    array[759] int SGCODEDV4V2;
    array[759] int SGCODEDV3V2;
    array[759] int SGCODEDV2V2;
     vector[759] Dim2TB;
     vector[759] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[759] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:759 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    SilverSidesNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[759] log_lik;
     vector[759] pbar;
    for ( i in 1:759 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:759 ) log_lik[i] = neg_binomial_2_lpmf( SilverSidesNumber[i] | pbar[i] , scale );
}"

SilverSides2CTBNBStanFile<-write_stan_file(SilverSides2CTBNB)
print(SilverSides2CTBNBStanFile)
SilverSides2CTBNBStanFileResults<-stan(SilverSides2CTBNBStanFile,data=SilverSidesStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)
#SilverSides TB NB 2G####
SilverSides2GTBNB<-"data{
    array[759] int SilverSidesNumber;
    array[759] int Gear1;
    array[759] int SGCODEDV2;
    array[759] int SGCODEDV1;
    array[759] int SGCODEDV1V2;
    array[759] int Gear3;
    array[759] int Gear2;
    array[759] int Model2CN;
    array[759] int YearR;
    array[759] int SGCODEDV4V2;
    array[759] int SGCODEDV3V2;
    array[759] int SGCODEDV2V2;
     vector[759] Dim2TB;
     vector[759] Dim1TB;
}
parameters{
     real a;
     real b1;
     real b2;
     real c2;
     real c3;
     real c4;
     real g2;
     real g3;
     vector[9] HC;
     real<lower=0> sigmaHC;
     vector[12] Yv;
     real<lower=0> sigmaYear;
     real scale;
}
model{
     vector[759] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 10 , 100 );
    g2 ~ normal( 10 , 100 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:759 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    SilverSidesNumber ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[759] log_lik;
     vector[759] pbar;
    for ( i in 1:759 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:759 ) log_lik[i] = neg_binomial_2_lpmf( SilverSidesNumber[i] | pbar[i] , scale );
}"

SilverSides2GTBNBStanFile<-write_stan_file(SilverSides2GTBNB)
print(SilverSides2GTBNBStanFile)
SilverSides2GTBNBStanFileResults<-stan(SilverSides2GTBNBStanFile,data=SilverSidesStanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)