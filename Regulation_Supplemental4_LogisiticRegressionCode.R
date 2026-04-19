#Exploring Proposed Regulations Reveal Differences Between Freshwater and Saltwater Florida Recreational Anglers
##By: Christopher Crowder - CrowderFishScience[at]gmail.com
###Produced 11/6/25
####Supplemental Material 5: Logistic Regression Code


####About: This script contains the code used to produce the logistic regressions results to determine which characteristics best describe if a participant believes a regulation is appropriate or not. The following describe names given after the model title and the codes for variables within the model. 
##Model Title Descriptions:
###Type of fisher: classification of participant as primarily saltwater or freshwater angler based on criteria noted in manuscript. 
###Survey mode: the mode in which the participant completed the survey. 
###Years of Experience: The number of years the participant noted to have been fishing in Florida. This is a dummary variable with 0-2 years as the intercept, 3-5 as Ex2/EC2, 
####More Restrictive Question: This is the participant's response to "do you believe a more restrictive regulations would work in the conservation of your favorite fish?". 


#Set Up#####################################################################
library(rethinking)
#Model Analysis#############################################################
#Model 1 - Type of Fisher (Freshwater VS Saltwater)
M1<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+ff2*SWF,
    a ~ dnorm(.5,10),
    ff2~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 2 - Survey Mode
M2<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+s2*SM2+s3*SM3,
    a ~ dnorm(.5,10),
    s2~dnorm(0,10),
    s3~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 3 - Years of Experience
M3<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+Ex2*EC2+Ex3*EC3+Ex4*EC4+Ex5*EC5+Ex6*EC6,
    a ~ dnorm(.5,10),
    Ex2~dnorm(0,10),
    Ex3~dnorm(0,10),
    Ex4~dnorm(0,10),
    Ex5~dnorm(0,10),
    Ex6~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 4 - More Restrictive Question
M4<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+m2*MR2+m3*MR3+m4*MR4+m5*MR5,
    a ~ dnorm(.5,10),
    m2~dnorm(0,10),
    m3~dnorm(0,10),
    m4~dnorm(0,10),
    m5~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 5 - Type of Fisher + Survey Mode
M5<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+ff2*SWF+s2*SM2+s3*SM3,
    a ~ dnorm(.5,10),
    ff2~dnorm(0,10),
    s2~dnorm(0,10),
    s3~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 6 - Type of Fisher + Years of Experience
M6<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+ff2*SWF+Ex2*EC2+Ex3*EC3+Ex4*EC4+Ex5*EC5+Ex6*EC6,
    a ~ dnorm(.5,10),
    ff2~dnorm(0,10),
    Ex2~dnorm(0,10),
    Ex3~dnorm(0,10),
    Ex4~dnorm(0,10),
    Ex5~dnorm(0,10),
    Ex6~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 7 - Type of Fisher + More Restrictive
M7<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+ff2*SWF+m2*MR2+m3*MR3+m4*MR4+m5*MR5,
    a ~ dnorm(.5,10),
    ff2~dnorm(0,10),
    m2~dnorm(0,10),
    m3~dnorm(0,10),
    m4~dnorm(0,10),
    m5~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 8 - Survey Mode + Years of Experience
M8<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+Ex2*EC2+Ex3*EC3+Ex4*EC4+Ex5*EC5+Ex6*EC6+s2*SM2+s3*SM3,
    a ~ dnorm(.5,10),
    Ex2~dnorm(0,10),
    Ex3~dnorm(0,10),
    Ex4~dnorm(0,10),
    Ex5~dnorm(0,10),
    Ex6~dnorm(0,10),
    s2~dnorm(0,10),
    s3~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 9 - Survey Mode + More Restrictive
M9<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+s2*SM2+s3*SM3+m2*MR2+m3*MR3+m4*MR4+m5*MR5,
    a ~ dnorm(.5,10),
    s2~dnorm(0,10),
    s3~dnorm(0,10),
    m2~dnorm(0,10),
    m3~dnorm(0,10),
    m4~dnorm(0,10),
    m5~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 10 - Years of Experience + More Restrictive
M10<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+Ex2*EC2+Ex3*EC3+Ex4*EC4+Ex5*EC5+Ex6*EC6+m2*MR2+m3*MR3+m4*MR4+m5*MR5,
    a ~ dnorm(.5,10),
    Ex2~dnorm(0,10),
    Ex3~dnorm(0,10),
    Ex4~dnorm(0,10),
    Ex5~dnorm(0,10),
    Ex6~dnorm(0,10),
    m2~dnorm(0,10),
    m3~dnorm(0,10),
    m4~dnorm(0,10),
    m5~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 11 - Type of Fisher + Survey Mode + Years of Experience 
M11<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+ff2*SWF+Ex2*EC2+Ex3*EC3+Ex4*EC4+Ex5*EC5+Ex6*EC6+s2*SM2+s3*SM3,
    a ~ dnorm(.5,10),
    ff2~dnorm(0,10),
    Ex2~dnorm(0,10),
    Ex3~dnorm(0,10),
    Ex4~dnorm(0,10),
    Ex5~dnorm(0,10),
    Ex6~dnorm(0,10),
    s2~dnorm(0,10),
    s3~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 12 - Type of Fisher + Survey Mode + More Restrictive Question
M12<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+ff2*SWF+s2*SM2+s3*SM3+m2*MR2+m3*MR3+m4*MR4+m5*MR5,
    a ~ dnorm(.5,10),
    ff2~dnorm(0,10),
    s2~dnorm(0,10),
    s3~dnorm(0,10),
    m2~dnorm(0,10),
    m3~dnorm(0,10),
    m4~dnorm(0,10),
    m5~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 13 - Type of Fisher + Years of Experience + More Restrictive Question
M13<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+ff2*SWF+Ex2*EC2+Ex3*EC3+Ex4*EC4+Ex5*EC5+Ex6*EC6+m2*MR2+m3*MR3+m4*MR4+m5*MR5,
    a ~ dnorm(.5,10),
    ff2~dnorm(0,10),
    Ex2~dnorm(0,10),
    Ex3~dnorm(0,10),
    Ex4~dnorm(0,10),
    Ex5~dnorm(0,10),
    Ex6~dnorm(0,10),
    m2~dnorm(0,10),
    m3~dnorm(0,10),
    m4~dnorm(0,10),
    m5~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 14 - Survey mode + Years of Experience + Conservation Question
M14<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+Ex2*EC2+Ex3*EC3+Ex4*EC4+Ex5*EC5+Ex6*EC6+s2*SM2+s3*SM3+m2*MR2+m3*MR3+m4*MR4+m5*MR5,
    a ~ dnorm(.5,10),
    Ex2~dnorm(0,10),
    Ex3~dnorm(0,10),
    Ex4~dnorm(0,10),
    Ex5~dnorm(0,10),
    Ex6~dnorm(0,10),
    s2~dnorm(0,10),
    s3~dnorm(0,10),
    m2~dnorm(0,10),
    m3~dnorm(0,10),
    m4~dnorm(0,10),
    m5~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)
#Model 15 - Type of fisher + Survey mode + Years of Experience + Conservation Question
M15<- rethinking::ulam(
  alist(
    BinomRegResponse ~ dbinom(1,p),
    logit(p)  <-a+ff2*SWF+Ex2*EC2+Ex3*EC3+Ex4*EC4+Ex5*EC5+Ex6*EC6+s2*SM2+s3*SM3+m2*MR2+m3*MR3+m4*MR4+m5*MR5,
    a ~ dnorm(.5,10),
    ff2~dnorm(0,10),
    Ex2~dnorm(0,10),
    Ex3~dnorm(0,10),
    Ex4~dnorm(0,10),
    Ex5~dnorm(0,10),
    Ex6~dnorm(0,10),
    s2~dnorm(0,10),
    s3~dnorm(0,10),
    m2~dnorm(0,10),
    m3~dnorm(0,10),
    m4~dnorm(0,10),
    m5~dnorm(0,10)
  ),
  data = ModelBinomialReg,
  iter=6000,warmup=3000,chains =3, cores = 6, log_lik = TRUE, sample=TRUE
)

compare(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,M15)