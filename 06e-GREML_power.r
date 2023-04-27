# Functions used in the functions above
var_vg_func <- function(N, var_pi=2e-5){
    return(2/(N^2*var_pi))
}

var_rg_func <- function(N1, N2, hsq1, hsq2, rg, rp, overlap=TRUE, var_pi=2e-5){
    if(overlap==T) var_rg=((1-rg*rp)^2+(rg-rp)^2)/(hsq1*hsq2*N1^2*var_pi)
    if(overlap==F) var_rg=(rg^2*(N1^2*hsq1^2+N2^2*hsq2^2)+2*hsq1*hsq2*N1*N2)/(2*hsq1^2*hsq2^2*N1^2*N2^2*var_pi)
    return(var_rg)
}

power_func <- function(ncp, alpha){
    pchisq(qchisq(alpha, df=1,lower.tail=F), ncp=ncp, df=1, lower.tail=F)
}

h2O_func <- function(ncase, ncontrol, K, h2L, var_pi=2e-5){
    n=ncase+ncontrol
    v=ncase/(ncase+ncontrol)
    z=dnorm(qnorm(K))
    c=(K*(1-K))^2/(v*(1-v)*z^2)
    h2O=h2L/c
    var_h2O=var_vg_func(n, var_pi)
    var_h2L=c^2*var_h2O
    return(list(h2L=h2L, var_h2L=var_h2L, h2O=h2O, var_h2O=var_h2O))
}

# Function for a quantitative trait
# n = sample size 
# hsq = variance explained by all SNPs
# alpha = significance level
# var_pi = variance of the off-diagonal elements of the GRM
# The output are: se (standard error), ncp (non-centrality parameter) and power
calcUniQt <- function(
    n     =1000, 
    hsq   =0.5, 
    alpha =0.05,
  var_pi=2e-5
){
    l <- list()
    var_vg <- var_vg_func(n, var_pi)
    l$se <- sqrt(var_vg)
    l$ncp <- hsq^2/var_vg;
    l$power <- power_func(l$ncp, alpha)
    return(l)
}

calcUniQt(5976, 0.15, alpha=0.05, 2e-5)
