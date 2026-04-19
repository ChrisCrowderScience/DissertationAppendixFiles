#Grouped Prey Fish Negative Binomial, Logistic, and Zero-Inflated Script
#2/10/25
#Prepared for Supplemental Material for Manuscript by Crowder et al 2025.
#Tampa Bay Grouped Zero-inflated###############################################################
#Prey Fish Models 2C and 2G. Predator Fish Models 2I and 2E. 
StanTBData <- list(
  Gear1 = FinalTB3Data$GearType1,
  Gear2 = FinalTB3Data$GearType2,
  Gear3 = FinalTB3Data$GearType3,
  YearR = FinalTB3Data$YearR,
  SGCODEDV2 = FinalTB3Data$SGCODEDV2,
  SGCODEDV1 = FinalTB3Data$SGCODEDV1,
  Model6CN = FinalTB3Data$Model6CN,
  Model5CN = FinalTB3Data$Model5CN,
  Model4CN = FinalTB3Data$Model4CN,
  Model3CN = FinalTB3Data$Model3CN,
  Model2CN = FinalTB3Data$Model2CN,
  SGCODEDV1V2 = FinalTB3Data$SGCODEDV1V2,
  SGCODEDV2V2 = FinalTB3Data$SGCODEDV2V2,
  SGCODEDV3V2 = FinalTB3Data$SGCODEDV3V2,
  SGCODEDV4V2 = FinalTB3Data$SGCODEDV4V2,
  SportFishCount = FinalTB3Data$SportFishCount,
  LogPreyFishCount = FinalTB3Data$LogPreyFishCount,
  PreyFishCount = FinalTB3Data$PreyFishCount,
  Dim1TB = FinalTB3Data$Dim1TB,
  Dim2TB = FinalTB3Data$Dim2TB,
  SilverPerchNumber = FinalTB3Data$SilverPerchNumber,
  PinfishNumber = FinalTB3Data$PinfishNumber,
  TroutNumber = FinalTB3Data$TroutNumber,
  SnookNumber = FinalTB3Data$SnookNumber,
  SnapperNumber = FinalTB3Data$SnapperNumber,
  MulletNumber = FinalTB3Data$MulletNumber,
  MojarraNumber = FinalTB3Data$MojarraNumber,
  JackNumber = FinalTB3Data$JackNumber,
  HerringNumber = FinalTB3Data$HerringNumber,
  AnchovyNumber = FinalTB3Data$AnchovyNumber,
  SilverSidesNumber = FinalTB3Data$SilverSidesNumber,
  PigfishNumber = FinalTB3Data$PigfishNumber,
  RedDrumNumber = FinalTB3Data$RedDrumNumber
)
#Defining the MCMC settings
rstan_options(auto_write = TRUE)
options(mc.cores=parallel::detectCores())
ni<-6000
nt<-2
nb<-3000
nc<-3
N<-14223
#Tampa Bay Prey 2C ZI####
Model2CPFTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int Model4CN;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  array[14223] int Model6CN;
  array[14223] int Model5CN;
  array[14223] int Model3CN;
  array[14223] int SportFishCount;
  array[14223] int PreyFishCount;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
  array[14223] int SilverPerchNumber;
  array[14223] int PinfishNumber;
  array[14223] int TroutNumber;
  array[14223] int SnookNumber;
  array[14223] int SnapperNumber;
  array[14223] int MulletNumber;
  array[14223] int MojarraNumber;
  array[14223] int JackNumber;
  array[14223] int HerringNumber;
  array[14223] int AnchovyNumber;
  array[14223] int SilverSidesNumber;
  array[14223] int PigfishNumber;
  array[14223] int RedDrumNumber;
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
  c2z ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  b1z ~ normal(0,100);
  b2z ~ normal(0,100);
  g2 ~ normal(1.5,10);
  g3 ~ normal(1.5,10);
  g2z ~ normal(.5,10);
  g3z ~ normal(.5,10);
  a ~ normal( 1.98 , 98.34 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+c2*SGCODEDV2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1+ b1z*Dim1TB[i]+b2z*Dim2TB[i]+c2z*SGCODEDV2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (PreyFishCount[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PreyFishCount[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PreyFishCount[i]|lambda[i]));
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
    if (PreyFishCount[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PreyFishCount[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PreyFishCount[i]|lambda[i]));
}
"
Model2CPFTBZIStanFile<-write_stan_file(Model2CPFTBZI)
print(Model2CPFTBZIStanFile)
Model2CPFTBZIStanFileResults<-stan(Model2CPFTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)

#Tampa Bay Prey 2G ZI####
Model2GPFTBZI<-"data{
  array[14223] int Gear1;
  array[14223] int SGCODEDV1V2;
  array[14223] int SGCODEDV2;
  array[14223] int SGCODEDV1;
  array[14223] int PreyFishCount;
  vector[14223] Dim2TB;
  vector[14223] Dim1TB;
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
  real g2;
  real g2z;
  real g3;
  real g3z;
  real z1;
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
  c2z ~ normal(0,100);
  c3z ~ normal(0,100);
  c4z ~ normal(0,100);
  b1z ~ normal(0,100);
  b2z ~ normal(0,100);
  c2 ~ normal(0,100);
  c3 ~ normal(0,100);
  c4 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  g2 ~ normal(1.5,10);
  g3 ~ normal(1.5,10);
  g2z ~ normal(.5,10);
  g3z ~ normal(.5,10);
  a ~ normal( 1.98 , 98.34 );
  for ( i in 1:14223 ) {
    lambda[i] = a + b1*Dim1TB[i]+b2*Dim2TB[i]+ c2*SGCODEDV2[i]+ c3*SGCODEDV3V2[i]+ c4*SGCODEDV4V2[i]+Yv[YearR[i]]+HC[Model2CN[i]]+g2*Gear2[i]+g3*Gear3[i];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:14223 ) {
    pbar[i] = z1 + b1z*Dim1TB[i]+b2z*Dim2TB[i]+ c2z*SGCODEDV2[i]+ c3z*SGCODEDV3V2[i]+ c4z*SGCODEDV4V2[i]+Yvz[YearR[i]]+HCz[Model2CN[i]]+g2z*Gear2[i]+g3z*Gear3[i];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:14223 )
    if (PreyFishCount[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PreyFishCount[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PreyFishCount[i]|lambda[i]));
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
    if (PreyFishCount[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PreyFishCount[i]|lambda[i])));
  else
    dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(PreyFishCount[i]|lambda[i]));
}
"
Model2GPFTBZIStanFile<-write_stan_file(Model2GPFTBZI)
print(Model2GPFTBZIStanFile)
Model2GPFTBZIStanFileResults<-stan(Model2GPFTBZIStanFile,data=StanTBData,init_r=0.5, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1)


#Tampa Bay Grouped Logistic Regressions####
FinalTB3Data$SportFishPresAbs<-ifelse(FinalTB3Data$SportFishCount==0,0,1)
FinalTB3Data$PreyFishPresAbs<-ifelse(FinalTB3Data$PreyFishCount==0,0,1)
LogStanTBData <- list(
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
  SportPresAbs = FinalTB3Data$SportFishPresAbs,
  PreyPresAbs = FinalTB3Data$PreyFishPresAbs,
  LogPreyFishCount = FinalTB3Data$LogPreyFishCount
)
#Tampa Bay Prey 2C Log####
Model2CPFTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int PreyPresAbs;
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
    g3 ~ normal( 1.5 , 10 );
    g2 ~ normal( 1.5 , 10 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 0.5 , 100 );
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    PreyPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( PreyPresAbs[i] | 1 , p[i] );
}"
Model2CPFTBLogStanFile<-write_stan_file(Model2CPFTBLog)
print(Model2CPFTBLogStanFile)
Model2CPFTBLogStanFileResults<-stan(Model2CPFTBLogStanFile,data=LogStanTBData, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1,init_r = 0.5)
#Tampa Bay Prey 2G Log####
Model2GPFTBLog<-"data{
    array[14223] int Gear1;
    array[14223] int SGCODEDV1;
    array[14223] int SGCODEDV4V2;
    array[14223] int SGCODEDV3V2;
    array[14223] int SGCODEDV2V2;
    array[14223] int SGCODEDV1V2;
    array[14223] int PreyPresAbs;
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
    g3 ~ normal( 1.5 , 10 );
    g2 ~ normal( 1.5 , 10 );
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
    PreyPresAbs ~ binomial( 1 , p );
}
generated quantities{
    vector[14223] log_lik;
     vector[14223] p;
    for ( i in 1:14223 ) {
        p[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:14223 ) log_lik[i] = binomial_lpmf( PreyPresAbs[i] | 1 , p[i] );
}
"
Model2GPFTBLogStanFile<-write_stan_file(Model2GPFTBLog)
print(Model2GPFTBLogStanFile)
Model2GPFTBLogStanFileResults<-stan(Model2GPFTBLogStanFile,data=LogStanTBData, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1,init_r = 0.5)

#Tampa Bay Grouped Negative Binomial#####
NonZeroPreyFishDf<-subset(FinalTB3Data,PreyFishCount>0)
TBPFStanData <- list(
  Gear1 = NonZeroPreyFishDf$GearType1,
  Gear2 = NonZeroPreyFishDf$GearType2,
  Gear3 = NonZeroPreyFishDf$GearType3,
  YearR = NonZeroPreyFishDf$YearR,
  SGCODEDV2 = NonZeroPreyFishDf$SGCODEDV2,
  SGCODEDV1 = NonZeroPreyFishDf$SGCODEDV1,
  SGCODEDV2V2 = NonZeroPreyFishDf$SGCODEDV2V2,
  SGCODEDV1V2 = NonZeroPreyFishDf$SGCODEDV1V2,
  SGCODEDV3V2 = NonZeroPreyFishDf$SGCODEDV3V2,
  SGCODEDV4V2 = NonZeroPreyFishDf$SGCODEDV4V2,
  Model2CN = NonZeroPreyFishDf$Model2CN,
  Dim1TB = NonZeroPreyFishDf$Dim1TB,
  Dim2TB = NonZeroPreyFishDf$Dim2TB,
  PreyFishCount = NonZeroPreyFishDf$PreyFishCount
)
NonZeroSportFishDf<-subset(FinalTB3Data,SportFishCount>0)
TBSFStanData <- list(
  Gear1 = NonZeroSportFishDf$GearType1,
  Gear2 = NonZeroSportFishDf$GearType2,
  Gear3 = NonZeroSportFishDf$GearType3,
  YearR = NonZeroSportFishDf$YearR,
  SGCODEDV2 = NonZeroSportFishDf$SGCODEDV2,
  SGCODEDV1 = NonZeroSportFishDf$SGCODEDV1,
  SGCODEDV2V2 = NonZeroSportFishDf$SGCODEDV2V2,
  SGCODEDV1V2 = NonZeroSportFishDf$SGCODEDV1V2,
  SGCODEDV3V2 = NonZeroSportFishDf$SGCODEDV3V2,
  SGCODEDV4V2 = NonZeroSportFishDf$SGCODEDV4V2,
  Model2CN = NonZeroSportFishDf$Model2CN,
  Dim1TB = NonZeroSportFishDf$Dim1TB,
  Dim2TB = NonZeroSportFishDf$Dim2TB,
  LogPreyFishCount = NonZeroSportFishDf$LogPreyFishCount,
  SportFishCount = NonZeroSportFishDf$SportFishCount
)
#Tampa Bay Prey 2C NegBin####
N<-10441
Model2CPFTBNegBin<-"data{
    array[10441] int Gear1;
    array[10441] int SGCODEDV1;
    array[10441] int SGCODEDV4V2;
    array[10441] int SGCODEDV3V2;
    array[10441] int SGCODEDV2V2;
    array[10441] int SGCODEDV1V2;
    array[10441] int PreyFishCount;
    array[10441] int Gear3;
    array[10441] int Gear2;
    array[10441] int Model2CN;
    array[10441] int YearR;
    array[10441] int SGCODEDV2;
     vector[10441] Dim2TB;
     vector[10441] Dim1TB;
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
     vector[10441] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 1.5 , 10 );
    g2 ~ normal( 1.5 , 10 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:10441 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    PreyFishCount ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[10441] log_lik;
     vector[10441] pbar;
    for ( i in 1:10441 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:10441 ) log_lik[i] = neg_binomial_2_lpmf( PreyFishCount[i] | pbar[i] , scale );
}"
Model2CPFTBNegBinStanFile<-write_stan_file(Model2CPFTBNegBin)
print(Model2CPFTBNegBinStanFile)
Model2CPFTBNegBinStanFileResults<-stan(Model2CPFTBNegBinStanFile,data=TBPFStanData, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1,init_r = 0.5)
#Tampa Bay Prey 2G NegBin####
Model2GPFTBNegBin<-"data{
    array[10441] int Gear1;
    array[10441] int SGCODEDV1;
    array[10441] int SGCODEDV4V2;
    array[10441] int SGCODEDV3V2;
    array[10441] int SGCODEDV2V2;
    array[10441] int SGCODEDV1V2;
    array[10441] int PreyFishCount;
    array[10441] int Gear3;
    array[10441] int Gear2;
    array[10441] int Model2CN;
    array[10441] int YearR;
    array[10441] int SGCODEDV2;
     vector[10441] Dim2TB;
     vector[10441] Dim1TB;
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
     vector[10441] pbar;
    scale ~ cauchy( 0 , 50 );
    sigmaYear ~ cauchy( 0 , 50 );
    Yv ~ normal( 0 , sigmaYear );
    sigmaHC ~ cauchy( 0 , 50 );
    HC ~ normal( 0 , sigmaHC );
    g3 ~ normal( 1.5 , 10 );
    g2 ~ normal( 1.5 , 10 );
    c4 ~ normal( 0 , 100 );
    c3 ~ normal( 0 , 100 );
    c2 ~ normal( 0 , 100 );
    b2 ~ normal( 0 , 100 );
    b1 ~ normal( 0 , 100 );
    a ~ normal( 5 , 100 );
    for ( i in 1:10441 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    PreyFishCount ~ neg_binomial_2( pbar , scale );
}
generated quantities{
    vector[10441] log_lik;
     vector[10441] pbar;
    for ( i in 1:10441 ) {
        pbar[i] = a + b1 * Dim1TB[i] + b2 * Dim2TB[i] + c2 * SGCODEDV2V2[i] + c3 * SGCODEDV3V2[i] + c4 * SGCODEDV4V2[i] + Yv[YearR[i]] + HC[Model2CN[i]] + g2 * Gear2[i] + g3 * Gear3[i];
        pbar[i] = exp(pbar[i]);
    }
    for ( i in 1:10441 ) log_lik[i] = neg_binomial_2_lpmf( PreyFishCount[i] | pbar[i] , scale );
}"
Model2GPFTBNegBinStanFile<-write_stan_file(Model2GPFTBNegBin)
print(Model2GPFTBNegBinStanFile)
Model2GPFTBNegBinStanFileResults<-stan(Model2GPFTBNegBinStanFile,data=TBPFStanData, chains=nc,iter=ni,warmup=nb,thin=nt,seed=1,init_r = 0.5)
