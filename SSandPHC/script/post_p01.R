
##################  Chapter 1 ######################

#Figure 1-2
binom01<-function(x){dbinom(x,10,0.5)};z01<-numeric(11);    #Binomial distribution definition for n=10
for (i in 0:10){z01[i+1]<-binom01(i)}; names(z01)<-0:10
barplot(z01,ylab='',xlab='x',cex.lab=1.7,cex.axis=1.7, col=gray(0.9));#Figure1-2
round(z01*100,1);                                                  #probability indication

round(dbinom(7,10,0.5),3);                                         #Equation (1.5)
sum(round(z01*100,1)[8:11])*2;                                     #Equation (1.6)
round(dbinom(70,100,0.5),6);                                       #Equation (1.7)
pbinom(30,100,0.5)+pbinom(69,100,0.5,lower.tail =F);               #Equation (1.8)
#Note that lower.tail =T is less than, lower.tail =F is greater than

#Figure 1-4
binom01<-function(x){dbinom(x,10,0.5)};z01<-numeric(11)
for (i in 0:10){z01[i+1]<-binom01(i)}; names(z01)<-c(0:10)/10;     #Sample Ratio
barplot(z01,ylab='',xlab='',cex.lab=1.7,cex.axis=1.7, col=gray(0.9))#Figure 1-4

#Figure 1-5
binom01<-function(x){dbinom(x,100,0.5)};z01<-numeric(101);  #Binomial distribution definition for n=100
for (i in 0:100){z01[i+1]<-binom01(i)};                     #Assigning probabilities
z02<-rep(NA,101);z02[c(1:10)*10]<-c(1:10)/10;z02[1]<-0;     #Creating labels
names(z01)<-z02
barplot(z01,ylab='',xlab='',cex.lab=1.7,cex.axis=1.7, col=gray(0.9));#Figure 1-5

#Rejection range for n=100
qbinom(0.025,100,0.5);qbinom(0.975,100,0.5)

#Interval estimation with confidence intervals and the width of the confidence intervals 
#(note that this is not the width of the adoption region)
binom.test(x=7,n=10,p=0.5);        round(0.9332605-0.3475471,3) 
binom.test(x=70,n=100,p=0.5);      round(0.7875936-0.6001853,3) 

#"Feeling the future" test results Reanalysis of Bem (2011) ##############
x <- 829 ;  n <- 1560                               #Number of positive responses and data 
binom.test(x,n,p=0.5);                              #Tests using the binomial distribution
(z<-(x-n*0.5)/sqrt(n*0.5*(1-0.5)));                 #test statistic with normal approximation 
(1-pnorm(z))*2;                                     #p-value by normal approximation

#Rejection range for n=1560
n*0.5+1.96*sqrt(n*0.5*(1-0.5));                     #Upper side rounded up to 819
n*0.5-1.96*sqrt(n*0.5*(1-0.5));                     #Lower side truncated to 741
819/n; 741/n;                                       #Rejection range of sample proportions
curve(dnorm(x,780,sqrt(1560*0.5*(1-0.5))),0,1560);  #Figure 1-6

#Let's do a "telekinesis experiment."
n<-100000;    x<-53100;
n*0.5+1.96*sqrt(n*0.5*(1-0.5));                     #Upper side rounded up to 50310
n*0.5-1.96*sqrt(n*0.5*(1-0.5));                     #Lower side truncated to 49690
50310/n; 49689/n;                                   #Rejection range of sample proportions
binom.test(x,n,p=0.5);                              #Binomial distribution for accurate testing
curve(dnorm(x,50000,sqrt(100000*0.5*(1-0.5))),0,100000);  #Figure 1-7
segments(50310,0,50310,0.00125);arrows(50310,0.00125,65000,0.00125);
segments(49690,0,49690,0.00125);arrows(49690,0.00125,35000,0.00125);

##################  Chapter 2  ######################

#Table 2-1
x1<-c(53.1,51.5,45.5,55.5,49.6,50.1,59.2,54.7,53.0,48.6,
      55.3,52.6,51.7,48.6,56.4,42.9,50.3,42.4,51.2,39.1)
x2<-c(48.3,45.2,46.6,56.6,41.2,44.6,51.9,55.5,45.4,47.6,
      50.6,54.5,49.0,43.9,53.8,40.1,52.8,35.3,55.6,38.0)
xd=x1-x2;                                           #weight difference
x<-cbind(x1,x2,xd);                                 #Summary of data
(n<-nrow(x));                                       #Number of observation targets

#Figure 2-1 
# right=F is the argument to make [greater than, less than] T is the argument to make [greater than, less than]
hist(xd,breaks=seq(-5,9,2),right=F,ylim=c(0,5),col=gray(0.9));

#numerical summary 
round(mean(xd),2);                                  # mean (weight difference)
round(mean(x1),2);round(mean(x2),2);                # mean
sort(xd);                                           # sort
(me1<-median(xd));                                  # median
van<-function(x){mean((x-mean(x))^2)};              # Function of sample variance
vad<-van(xd);round(vad,2);                          # variance
round(var(xd),2);                                   # unbiased variance
sdd<-sqrt(vad);round(sdd,2);                        # standard deviation
round(sd(xd),2);                                    # SD by unbiased variance

#Test for difference of means of two corresponding groups
t.test(x1,x2,paired=T);                             # t-test (exact)
(t<-(mean(xd)/sd(xd))*sqrt(20));                    # t-value (approximate match by normal approximation)
(1-pnorm(3.2185,0,1))*2;                            # p-value (slightly off by normal approximation)

#Move n in equation (2.12)
for (n in c(20,50,100,1000,5000,100000,5)){print(round((1.96*3.81)/sqrt(n),3))}

##################  Chapter 3  ######################

#Raw data on therapeutic effect of antimicrobial A
ex_group<-c(5.5,10.0,5.5,6.0,9.0,9.5,6.5,7.0,12.5,7.0,6.5,10.5,
                 9.0,4.5,6.5,9.5,10.0,9.5,10.5,6.5,8.5,12.5,4.0,9.0)
con_group<-c(5.5,11.5,7.0,7.5,10.5,11.0,8.0,8.5,14.0,8.5,10.0,8.0,12.0,
                 10.5,6.0,8.0,11.0,11.5,11.0,12.0,8.0,10.0,14.0,5.5,10.5,9.0)

#Table 3-1
length(ex_group);       length(con_group);               # Number of observation targets
round(mean(ex_group),1);round(mean(con_group),1);        # sample mean
round(sd(ex_group),1);  round(sd(con_group),1);          # unbiased variance

#Figure 3-1
par(mfrow=c(2,1))
hist(ex_group,breaks=seq(3,15,2), col =gray(0.9),cex.lab=2.0,cex.axis=2.0)
hist(con_group,breaks=seq(3,15,2), col =gray(0.9),cex.lab=2.0,cex.axis=2.0)
par(mfrow=c(1,1))

(s_within<-sqrt((23*sd(ex_group)^2 + 25*sd(con_group)^2)/(24+26-2)));   #Equation (3.3)
(naste<-sqrt((24 * 26)/(24+26)));                                       #Equation (3.4)
(z<-((mean(con_group)-mean(ex_group))/s_within)*naste);                 #Equation (3.5)
(1-pnorm(2.159465))*2;                                        #p-value (normal approximation)
(1-pt(2.159465,df=24+26-2))*2;                                #p-value (t-distribution)
t.test(con_group,ex_group,var.equal = T);                          #t-test

# Move n in Equation (3.6)
nas<-function(n){sqrt(n^2/(2*n))};                            #Equation (3.6)
for (n in c(1000,10000,100000,5)) { print(round(1.96*2.3/nas(n),4)) }
24* 0.2016;                                          #Convert days to units of time
60*24*0.0638;                                        #Convert days to minutes
60*24*0.0202;                                        #Convert days to minutes

#Statistical Power Analysis #################################
#install.packages("pwr");                            #Run once if no pwr
library(pwr);                                        #Load and attach add-on package pwr

#Table 3-4   Determining n for Tests of Proportions
#h1<-seq(0.45,0.05,-0.05);          # To calculate from 0.45 to 0.05 (values are the same)
h1<-seq(0.55,0.95,0.05);            # Scope of calculation
h2<-ES.h(h1,0.5);                   # Effect size for ratio tests
box1<-matrix(0,3,length(h1));       # Creating a matrix to store the results
ii<-0;                              # counter initialization
for (i in h2){
 ii<-ii+1
 box1[1,ii]<-ceiling(
       pwr.p.test(h=i,sig.level=0.05,power=0.7,alternative="two.sided")$n)
 box1[2,ii]<-ceiling(
       pwr.p.test(h=i,sig.level=0.05,power=0.8,alternative="two.sided")$n)
 box1[3,ii]<-ceiling(
       pwr.p.test(h=i,sig.level=0.05,power=0.9,alternative="two.sided")$n)
}
rowname01<-c("1-\beta=0.7","1-\beta=0.8","1-\beta=0.9")
rownames(box1)<-rowname01
colnames(box1)<-h1
print(box1)

#Figure 3-3
nor01<-function(x,y){dnorm(x,y,1)};                #Functional definition of normal probability density
par(mfrow=c(1,3))
x01<-c(-3,3.3);                                    #δ=0.3
curve(nor01(x,y=0),xlim=x01,ylim=c(0,0.4))
par(new=T);curve(nor01(x,y=0.3),xlim=x01,ylim=c(0,0.4))
x01<-c(-3,3.5);                                    #δ=0.5
curve(nor01(x,y=0),xlim=x01,ylim=c(0,0.4))
par(new=T);curve(nor01(x,y=0.5),xlim=x01,ylim=c(0,0.4))
x01<-c(-3,3.9);                                    #δ=0.8
curve(nor01(x,y=0),xlim=x01,ylim=c(0,0.4))
par(new=T);curve(nor01(x,y=0.9),xlim=x01,ylim=c(0,0.4))
par(mfrow=c(1,1))

#Table 3-6  Number of Observations Required for the Test of Mean Differences Between Two Independent Groups
h1<-seq(0.1,1.0,0.1);               #Scope of calculation
box2<-matrix(0,3,length(h1));       #Creating a matrix to store the results
ii<-0;                              #counter initialization
for (i in h1){
 ii<-ii+1
 box2[1,ii]<-ceiling(
       pwr.t.test(d=i,sig.level=0.05,power=0.7,type="two.sample")$n)
 box2[2,ii]<-ceiling(
       pwr.t.test(d=i,sig.level=0.05,power=0.8,type="two.sample")$n)
 box2[3,ii]<-ceiling(
       pwr.t.test(d=i,sig.level=0.05,power=0.9,type="two.sample")$n)
}
rownames(box2)<-rowname01
colnames(box2)<-h1
print(box2)

#Table 3-7   Number of Observations Required for the Test of Mean Differences Between Two Related Groups 
h1<-seq(0.1,1.0,0.1);               #Scope of calculation
box3<-matrix(0,3,length(h1));       #Creating a matrix to store the results
ii<-0;                              #counter initialization
for (i in h1){
 ii<-ii+1
 box3[1,ii]<-ceiling(pwr.t.test(d=i,sig.level=0.05,power=0.7,type="paired")$n)
 box3[2,ii]<-ceiling(pwr.t.test(d=i,sig.level=0.05,power=0.8,type="paired")$n)
 box3[3,ii]<-ceiling(pwr.t.test(d=i,sig.level=0.05,power=0.9,type="paired")$n)
}
rownames(box3)<-rowname01
colnames(box3)<-h1
print(box3)

##################  Chapter 4  ######################

# One-tailed probability and binomial test with intention 1
pbinom(7,24,0.5);                                # Probability of having less than 7 people in binomial distribution
binom.test(7,24,p=0.5,alternative = "two.sided");# Binomial test (intention 1 not significant)

# One-sided probability of negative binomial distribution by intention2.
# The grammar of the distribution function is pnbinom(number of cured, number of noncured, mother ratio).
# It took 24 treatments before we observed 7 nonhealers.
# The probability that the nonheal rate is lower than that is 24-7=17 or more healed patients.
# Subtract the probability (distribution function) of having less than 16 cured patients from 1.
# Since doubling the number is less than 0.05, intention 2 is a significant difference.
1-pnbinom(16,7,0.5)

#Figure 4-1　Sample Distribution of Patient Numbers by Negative Binomial Distribution
#Let the number of treatments be a function of the number of treatments as 
# pnbinom(Number of Treatments-Number of Non-Cured, Number of Non-Cured, Population Ratio)
dnbinom01<-function(x){dnbinom(x-7,7,0.5)};z01<-numeric(24);   #function definition
for (i in 7:30){z01[i-6]<-dnbinom01(i)}; names(z01)<-7:30;     #probability calculation
barplot(z01,col=gray(0.9));                                    #Figure 4-1

# Simulation of true risk rate with intention 4
# Calculate from 3 cases to 1 case each with a cure rate that is larger is more desirable according to the test and bem
# Smaller value of s =10000 will allow for shorter computation time.
set.seed(1234)
a<-365*2;                                           # Number of certifications
s<-100;                                             # Number of simulations
proba<-0.50;                                        # Population Ratio
r<-numeric(s);                                      # Vector Preparation
for (j in 1:s) {                                    # Counting of follow-up examinations
  x<-rbinom(a,1,prob=proba);                        # Data generation
  for (i in 3:a){                                   # 3 case or later certification
    no1<-sum(x[1:i]);hi<-no1/i;                     # Number of successes /success rate
    pval<-binom.test(no1,i,p=0.5)$p.value;           # Every morning examination
    if ((pval<0.05)&(hi>=0.5)) {r[j]<-1; break};     # If more than 0.5 and significant, write
  }
}
#print(r);                                          # Raw data Significant difference is 1
(rm<-mean(r));                                      # True risk rate

#The true risk rate for 2 years of hard work was 0.1994.
#The true risk ratio for 3 years of hard work was 0.2185.
#If we set the mother ratio proba=0.53 to match Bem,
# The true risk ratio for 2 years of hard work was 0.6091.
# The true risk rate for 3 years of hard work was 0.7314.
