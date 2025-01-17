//ÎÌ éa~bÌNX\ÉÖ·évIª
data { 
  int<lower=0> a;                          //s 
  int<lower=0> b;                          //ñ 
  int<lower=0> x[a,b];                     //½Ìsñ`® 
}
transformed data{
  int<lower=0> N;                          //v½ 
  int<lower=0> ab;                         //Zv 
  int<lower=0> xv[a*b];                    //½ÌxNg`®
  ab  = a*b;
  for (i in 1:a){
    for (j in 1:b){
      xv[(i-1)*b+j]  = x[i,j];}}
  N  = sum(xv);
}
parameters {
  simplex[ab]    pi;                        //aª1Ìêä¦ÌxNg`®
}
transformed parameters{
  real<lower=0,upper=1> pim[a,b];           //aª1Ìêä¦Ìsñ`®
  for (i in 1:a){
    for (j in 1:b){
      pim[i,j]  = pi[(i-1)*b+j];}}
}
model {
  xv ~ multinomial(pi);                     //½ªz
}
generated quantities{
  int<lower=0> xastev[ab];                  //\ªªzxNg`®
  int<lower=0> xaste[a,b];                  //\ªªzsñ`®
  real log_lik;
  real pa[a]; real pb[b];  real V;  real res[a,b];
  real<lower=0,upper=1>   Up[a,b];          //sA\c·ª{
  real<lower=0,upper=1>   Um[a,b];          //sA\c·ª[
  real<lower=0>            L[a,b];          //¯m¦ÆüÓm¦ÌÏÆÌä
    xastev    = multinomial_rng(pi,N);      //\ªªz
    V  = 0;
    for (i in 1:a){pa[i] =0;}               //üÓm¦
    for (j in 1:b){pb[j] =0;}
    for (i in 1:a){
      for (j in 1:b){
        pa[i] =pa[i]+pim[i,j];
        pb[j] =pb[j]+pim[i,j];
    }}
    V  = 0;
    for (i in 1:a){
      for (j in 1:b){
        xaste[i,j]  = xastev[(i-1)*b+j];
        res[i,j]    =(pim[i,j]-pa[i]*pb[j])/sqrt(pa[i]*pb[j]);
        L[i,j]      =pim[i,j]/(pa[i]*pb[j]);
        V  = V+ pow(res[i,j],2);
        Up[i,j] = res[i,j]>0 ? 1 : 0;
        Um[i,j] = !(Up[i,j]);
    }}
    V  = sqrt(V /(min(a,b)-1));
    log_lik  = multinomial_lpmf(xv | pi);      //Îm¦
}
