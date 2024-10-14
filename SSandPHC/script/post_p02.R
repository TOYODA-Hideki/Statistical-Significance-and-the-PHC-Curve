########################################################################
(n_wd<-getwd())                # Confirmation of working directory
source('myfunc/myfunc.R')      # Loading self-made functions
library(cmdstanr)              # load and attach add-on package cmdstanr
library(posterior)             # load and attach add-on package posterior

##################  Chapter 5  ######################

# Example of Bayesian inference calculation
# In the case where words C, D, and E are in the email
(C<-(0.85*0.66)/((0.85*0.66)+0.11*(1-0.66)))
(D<-(0.83*C   )/((0.83*C   )+0.09*(1-C)))
(E<-(0.88*D   )/((0.88*D   )+0.08*(1-D)))
# In the case where words C, D, and E are not in the email
(C<-((1-0.85)*0.25)/(((1-0.85)*0.25)+(1-0.11)*(1-0.25)))
(D<-((1-0.83)*C   )/(((1-0.83)*C   )+(1-0.09)*(1-C   )))
(E<-((1-0.88)*D   )/(((1-0.88)*D   )+(1-0.08)*(1-D   )))

# Reanalysis of "Presentiment" by Bem (2011) using PHC
x<-829; n<-1560;                      # Number of successes & number of trials
Bem<-Bi01(x,n);
gqcal(Bem$theta);                                      # Table 5-3
hist(Bem$theta, breaks=100,cex.axis=2.0);              # Figure 5-1
seq01<-seq(0.50,0.60,0.01);                            # Table 5-4
PHC01(seq01,Bem$theta,0,cc="gtc",byoga="no",dedits=3)
seq02<-seq(0.50,0.55,0.005);                           # Figure 5-2
PHC01(seq02,Bem$theta,0,cex.lab=2.2, cex.axis=2.5)
lines(c(0.50,0.510,0.510),c(0.95,0.95,0.00),lty=2,lwd=2.0)
mean(abs(Bem$theta-0.5)<0.05);                         # ROPE [0.45,0.55]
seq03<-seq(0,0.1,0.005);                               # Figure/Table 5-3
PHC01(seq03,abs(Bem$theta-0.5),0,cc="ltc",cex.lab=2.2, cex.axis=2.5)
lines(c(0.05,0.05,0.00),c(0.00,0.9299,0.9299),lty=2,lwd=2.0)

# Differences in PHC curves for large and small values of n - Figure 5-4
B02<-Bi01(x=7,n=24, cha=2, war=100, ite=2100)
B03<-Bi01(x=700,n=2400, cha=2, war=100, ite=2100)
seq04<-seq(0.00,0.60,0.005)
par(mfrow=c(2,1));                                     # Draw figures in 2 rows and 1 column
PHC01(seq04,B02$theta,0,cc="ltc",cex.lab=2.2, cex.axis=1.5)
PHC01(seq04,B03$theta,0,cc="ltc",cex.lab=2.2, cex.axis=1.5)
# dev.copy2pdf(file='z0002.pdf')
par(mfrow=c(1,1))



##################  Chapter 6  ######################

# Table 6-1 Raw data on the therapeutic effect of antibiotic A
experimental_group<-c(5.5,10.0,5.5,6.0,9.0,9.5,6.5,7.0,12.5,7.0,6.5,10.5,
                      9.0,4.5,6.5,9.5,10.0,9.5,10.5,6.5,8.5,12.5,4.0,9.0)
control_group<-c(5.5,11.5,7.0,7.5,10.5,11.0,8.0,8.5,14.0,8.5,10.0,8.0,12.0,
                 10.5,6.0,8.0,11.0,11.5,11.0,12.0,8.0,10.0,14.0,5.5,10.5,9.0)

# EAP estimates, confidence intervals
outDEF<-G2Ind(control_group, experimental_group,
              EQU=0, prior=T, mL=0, mH=100, sL=0, sH=100)
gqcal(outDEF$ext[,1:6])

mean_def<-outDEF$mu1-outDEF$mu2;                           # Posterior distribution of mean difference
hist(mean_def, breaks=200,cex.axis=2.0);                   # Figure 6-1

PHC01(seq(0.0,0.5,0.1),mean_def,0,cc="gtc", byoga="no");   # Table 6-2
PHC01(seq(0.6,2.0,0.2),mean_def,0,cc="gtc", byoga="no")
PHC01(seq(2.5,4.0,0.5),mean_def,0,cc="gtc", byoga="no")

PHC01(seq(0.0,4.0,0.1),mean_def,0,cc="gtc", byoga="yes");  # Figure 6-2  # Left figure
PHC01(seq(0.0,4.0,0.1),mean_def,0,cc="rope",byoga="yes");  # Right figure

PHC01(seq(0.0,1.0,0.1),mean_def,0,cc="rope", byoga="no");  # Table 6-3
PHC01(seq(2.0,3.0,0.1),mean_def,0,cc="rope", byoga="no")

delta_con<-(outDEF$mu1-outDEF$mu2)/outDEF$sigma1;  # Posterior distribution of standardized mean difference
hist(delta_con, breaks=200,cex.axis=2.0);                  # Figure 6-3

PHC01(seq(-0.5,2.0,0.05),delta_con,0,cc="gtc",byoga="yes");# Figure 6-4
PHC01(seq(0.0,0.4,0.05), delta_con,0,cc="gtc",byoga= "no");# Table 6-4

(9.6-8.1)/2.3;                  # Effect size calculated with sample mean & sample standard deviation

U3<-pnorm(outDEF$mu1,outDEF$mu2,outDEF$sigma2);            # Posterior distribution of measure of nonoverlap
PHC01(seq(0.4,1.0,0.05),U3,0,cc="gtc", byoga="yes");       # Figure 6-6
PHC01(seq(0.5,0.6,0.01),U3,0,cc="gtc", byoga="no");        # Table 6-5


xas1<-outDEF$xaste1; xas2<-outDEF$xaste2;                  # Posterior predictive distribution
PHC01(seq(0.0,2.5,0.1),xas1,xas2,cc="gtc", byoga="yes");   # Figure 6-7
plot(xas1,xas2,pch='.',xlim=c(1,17),ylim=c(1,17));         # Figure 6-8
abline(-0.5,1.0,lwd=2.0);
PHC01(seq(0.0,2.5,0.25),xas1,xas2,cc="gtc",byoga="no");    # Table 6-6

hiritu<-outDEF$mu2/outDEF$mu1;                             # Posterior distribution of ratio of means
gqcal(hiritu,2);                                           # Summary statistics
hist(hiritu, breaks=100,cex.axis=2.0);                     # Histogram
PHC01(seq(0.9,1.0,0.01),hiritu,0,cc="ltc", byoga="no");    # Table 6-7
PHC01(seq(0.7,1.0,0.01),hiritu,0,cc="ltc", byoga="yes");   # Practical Assignment Solutions

##################  Chapter 7  ######################

# Table 7-1
x1<-c(48.5,45.6,56.9,52.6,52.6,60.5,52.2,55.3,58.8,54.0,47.3,57.4,48.2,49.9,
      47.9,54.1,55.1,49.7,48.7,53.6,52.0,42.0,58.0,49.6,54.8,50.3,49.3,51.9,
      46.0,53.2)
x2<-c(47.0,42.6,54.4,51.0,51.0,58.3,50.4,54.4,57.6,51.5,45.8,55.7,46.9,48.1,
      45.0,51.9,53.2,47.0,46.0,51.6,50.5,40.0,56.4,47.2,53.0,48.0,47.0,49.4,
      44.8,50.9)
x3<-c(148.5,147.4,164.6,158.3,158.3,169.7,154.0,162.3,170.5,164.3,150.1,165.3,
      148.0,154.1,151.0,156.8,162.0,153.8,152.3,159.8,157.4,146.9,162.4,153.7,
      161.5,154.8,149.7,161.1,148.0,155.5)

x<-cbind(x1,x2,x3); n<-nrow(x);                            # Table 7-2
(round(apply(cbind(x,x1-x2),2,mean),1));                   # Mean values
(round(apply(cbind(x,x1-x2),2,sd),2));                     # Standard deviations
(round(apply(cbind(x,x1-x2),2,median),1));                 # Median values

t.test(x1, x2,paired = T); # Paired t-test for two groups, check it's the same as below
t.test(x1-x2, mu = 0);     # One-sample t-test for mean difference of paired scores being 0


########### Stan code for the difference in obesity levels #################### Refer to the following materials
BMI_01<-'
data {
  int<lower=0> n; // Number of data points
  array[n] vector[3] x; // Data
}
parameters {
  vector[3] mu; // Mean (no range specified)
  cov_matrix[3] Sigma; // Covariance matrix
}
transformed parameters {
  
}
model {
  for (i in 1 : n) {
    x[i] ~ multi_normal(mu, Sigma);
  } // Multivariate normal distribution
}
generated quantities {
  vector[3] xaste;
  xaste = multi_normal_rng(mu, Sigma); // Predictive distribution
}
';

##### Bayesian analysis of the difference in obesity levels ########
par<-c("mu","Sigma","xaste");                        # Parameters
dat <- list(n=nrow(x), x=x);                         # Data
initi<-function(){list(mu=colMeans(x),Sigma=cov(x))};# Initial values

########### Execution with cmdstanr
modfileBMI <- write_stan_file(BMI_01)              # Write to a temporary file
modBMI <- cmdstan_model(modfileBMI)                   # Compile
csrfit_BMI <- modBMI$sample(data = dat,chains = 2,iter_sampling = 1500,
        init = initi, iter_warmup = 500,parallel_chains = 2,seed=1234)  # MCMC
      draw<-csrfit_BMI$draws(par)
      post<-as_draws_df(draw)
      (colnames(post) <- gsub("\\[|\\]|,", "", colnames(post)))
gqcal(post[,1:15]);                                  #output

### Posterior distribution of generated quantities, Predictive distribution ###
withinSD<-sqrt((post$Sigma11+post$Sigma22)/2); 
diffSD  <-sqrt(post$Sigma11+post$Sigma22-2*post$Sigma21);
mudiff<-post$mu1-post$mu2;                  
withinES<-mudiff/withinSD;                  
diffES  <-mudiff/diffSD;                    
x_diff  <-post$xaste1-post$xaste2;          
BMI01   <- post$xaste1/(post$xaste3/100)^2; 
BMI02   <- post$xaste2/(post$xaste3/100)^2; 
BMIdiff <- BMI01-BMI02;                     

# Table 7-3    ### Output of summary statistics of generated quantities ###
pro<-c(0.05, 0.5, 0.95);               # Specifying probability points
gqcal(withinSD,probs=pro,digits =2);   # SD within the group
gqcal(diffSD ,probs=pro,digits =2);    # SD of weight difference
gqcal(withinES,probs=pro,digits =2);   # Posterior distribution according to formula (3.8)
gqcal(diffES ,probs=pro,digits =2);    # Posterior distribution according to formula (3.9)
gqcal(mudiff     ,probs=pro,digits =2);# Posterior distribution of weight difference
gqcal(x_diff   ,probs=pro,digits =2);  # Predictive distribution of weight difference
gqcal(BMI01    ,probs=pro,digits =1);  # Predictive distribution of obesity level b
gqcal(BMI02    ,probs=pro,digits =1);  # Predictive distribution of obesity level a
gqcal(BMIdiff  ,probs=pro,digits =2); # Predictive distribution of the difference in obesity levels
mean(abs(mudiff)<1.0);  #ROPE [-1,1]; # PHC for ROPE of ±1kg range
hist(BMIdiff,breaks=100); # Predictive distribution histogram of obesity level difference

# Figure 7-1      # Drawing PHC curves ###########################
par(mfrow=c(3,2))
seq_a<-seq(0.2,0.6,length.out=30)
PHC01(seq_a,withinES,0,cc="gtc",xlab="a. Mean difference standardized by within-group SD")
seq_b<-seq(0.4,0.9,length.out=30)
PHC01(seq_b,diffSD,0,cc="ltc",     xlab="b. SD of weight difference (kg)")
seq_c<-seq(1.7,2.3,length.out=30)
PHC01(seq_c,mudiff,0,cc="gtc",     xlab="c. Mean weight difference (kg)")
seq_d<-seq(0,2.3,length.out=30)
PHC01(seq_d,abs(mudiff),0,cc="ltc",xlab="d. ROPE of mean weight difference (kg)")
seq_e<-seq(0.8,3.3,length.out=30)
PHC01(seq_e,x_diff,0,cc="gtc",     xlab="e. Individual weight difference (kg)")
seq_f<-seq(0.25,1.4,length.out=30)
PHC01(seq_f,BMIdiff,0,cc="gtc",    xlab="f. Difference in obesity level (BMI)")
# # dev.copy2pdf(file='z0003.pdf')
par(mfrow=c(1,1))

