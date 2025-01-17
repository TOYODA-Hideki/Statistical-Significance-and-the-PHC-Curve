//JeSŞkÌäĤÌvIŞ
data { 
  int<lower=0>    k;                       //JeS 
  int<lower=0> x[k];                       //½ 
}
transformed data{
}
parameters {
  simplex[k]    pi;                        //êäĤ
}
model {
  x ~ multinomial(pi);                     //½Şz
}
generated quantities{
  real<lower=0,upper=1>   U2[k,k];         //2är
  real log_lik;
  log_lik  = multinomial_lpmf(x | pi);
  for (i in 1:k){  U2[i,i] =0;}
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      U2[i,j] = pi[i]-pi[j]>0 ? 1 : 0;
      U2[j,i] = !(U2[i,j]);  }  }
}
