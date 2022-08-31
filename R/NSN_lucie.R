library(clinfun)

ph2single(pu=0.1,pa=0.2,ep1=0.05,ep2=0.2,nsoln=5)



ph2simon(pu=0.15, pa=0.3, ep1=0.05, ep2=0.1, nmax=200)
ph2simon(pu=0.15, pa=0.25, ep1=0.05, ep2=0.1, nmax=200)

ph2simon(pu=0.08, pa=0.2, ep1=0.05, ep2=0.1, nmax=200)
ph2simon(pu=0.1, pa=0.2, ep1=0.05, ep2=0.1, nmax=200)

###########################################
# TOX only
###########################################

ptox <- 0.4
a0 <- 1
b0 <- 1


seuil <- NULL

for (i in 1:10000){
x <- c(rbinom(4, 5, ptox), rbinom(7, 10, ptox))
n <- c(rep(5,4), rep(10,7))

res <- pbeta(0.25, a0+cumsum(x), b0+cumsum(n)-cumsum(x) ,lower.tail = F)

seuil <- c(seuil ,sum(res>0.9))

} 
table(seuil>0)


###########################################
# TOP + TOX only
###########################################

ptox <- 0.3
prep <- 0.15
a0 <- 1
b0 <- 1


seuil_tox <- seuil_tox2 <- seuil_tox3 <- seuil_tox4 <- stop_fut <- efficacite <-  NULL
essais <- gen_patients_multinomTOP(
  n_sim = 10000,
  ana_inter = c(30, 51),
  ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
  interpatient = 6, 
  max_teff = 180,
  max_ttox = 180,
  multinom_ttt = c(.09, .06, .21, .64),
  seed = 1993
)
n <- cumsum(c(rep(5,4), rep(10,5), 11))
nr <- c(30,51)

for (i in 1:10000){
  x <- c(rbinom(4, 5, ptox), rbinom(5, 10, ptox), rbinom(1,11,ptox))
  
  y <- c(rbinom(1, 30, prep), rbinom(1, 51, prep))
  
  
  z <- cbind(rmultinom(4, 5, c(.09, .06, .21, .64)), rmultinom(5, 10, c(.09, .06, .21, .64)), rmultinom(1, 11, c(.09, .06, .21, .64)))
  z <- t(apply(z, 1, cumsum))
  tox <- z[1, ] + z[3, ]
  
  w <- rmultinom(81, 1, c(.09, .06, .21, .64))
  tox2 <- w[1, ] + w[3, ]
  tox2 <- map_dbl(c(5, 10, 15, 20, 30, 40, 50, 60, 70, 81), ~ sum(tox2[seq_len(.x)]))
  
  tox3 <- map_dbl(c(5, 10, 15, 20, 30, 40, 50, 60, 70, 81), ~ sum(essais[[i]]$tox[seq_len(.x)]))
  
  res <- pbeta(0.25, a0+cumsum(x), b0+n-cumsum(x) ,lower.tail = F)
  res2 <- pbeta(.25, a0 + tox, b0 + n - tox, lower.tail = F)
  res3 <- pbeta(.25, a0 + tox2, b0 + n - tox2, lower.tail = F)
  res4 <- pbeta(.25, a0 + tox3, b0 + n - tox3, lower.tail = F)
  
  seuil_tox <- c(seuil_tox ,sum(res>0.95))
  stop_fut <- c(stop_fut, I(y[1]<4)*1)
  efficacite <- c(efficacite, I(sum(y)>17)*1)
  seuil_tox2 <- c(seuil_tox2, sum(res2 > .95))
  seuil_tox3 <- c(seuil_tox3, sum(res3 > .95))
  seuil_tox4 <- c(seuil_tox4, sum(res4 > .95))
  
} 
table(seuil_tox>0)
table(seuil_tox2 > 0)
table(seuil_tox3 > 0)
table(seuil_tox4 > 0)
table(efficacite==1 & stop_fut==0)
table(efficacite==1 & stop_fut==0 & seuil_tox==0)

test_res <- map_dfr(
  essais,
  ~ realisation_essai_bopiva(ana_inter_eff = c(30, 51),
                             ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                             tableau = .x,
                             phitox = c(.15, .25),
                             prior_eff = .15, prior_tox = c(1, 1),
                             C_ = .92, gamm = .97,
                             p_n = c(.09, .06, .21, .64),
                             critere_tox = .95, interpatient = 6, 
                             max_teff = 180,
                             max_ttox = 180)
)

###############################################################
binom.test(5, 10, 0.5)


a<- 1
b <- 1
res <- res2 <- res3 <- NULL
for (n in 10:15){
  
  res <- cbind(res, 
               c(paste(round((a+(0:n))/(a+(0:n)+ b+n-(0:n)), 3), " (", 
                              round(qbeta(0.025, a+(0:n), b+n-(0:n)),3),
                              "-",
                              round(qbeta(0.975, a+(0:n), b+n-(0:n)),3), ")",
                              sep=""), 
                      rep("", 16-n)),
               c(round(pbeta(0.5, a+(0:n), b+n-(0:n), lower.tail=F),3), rep("", 16-n)), 
   c(round(pbeta(0.4, a+(0:n), b+n-(0:n), lower.tail=F),3), rep("", 16-n)))
  
}
write.csv2(res,"X:/METIMGAST/PROTOCOLE/betabin.csv" )

