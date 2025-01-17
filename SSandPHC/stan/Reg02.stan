//dρAf
data { 
  int<lower=0>  n;                     //f[^ 
  int<lower=0>  p;                     //\ͺΟ   
  vector[n]     y;                     //ξΟ   
  matrix[n,p]   X;                     //\ͺΟsρ   
}
parameters {
  real        a;                       //ΨΠ
  vector[p]   b;                       //ΞρAWxNg  
  real<lower=0> sigma;                 //λ·WΞ·
}
transformed parameters {
  vector[n]  yhat;                     //yhat  
  yhat = a + X*b;
}
model {
    y ~ normal(yhat, sigma);           //³Kͺzf
}
generated quantities{
  real log_lik;                        //Ξήx
  real<lower=0> vyhat;                 //yhatΜͺU
  real<lower=0,upper=1>    r2;         //θW
  vyhat =(variance(yhat)*(n-1))/n;
  r2 = vyhat/(vyhat+sigma^2);
  log_lik  =  normal_lpdf(y | a + X*b, sigma);    
}
