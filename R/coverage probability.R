data("df")
#Gaussian
bell<-link(t = df$Time, mu = 2.29, k1 = 2.39, k2 = 21.67, t0 = 0.59)
mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)

k1_gau <- NULL
k2_gau <- NULL
t0_gau <- NULL
for (i in 1:1000){
  sim_gau <- rnorm(500, mean = mut, sd = 2) #Gaussian
  sim_gau[sim_gau < 0] <- 0

  result_gau <- scGTM(t=df$Time, y1 = sim_gau, marginal= "Gaussian", iter = 200)
  if(2.39 > result_gau$k1_lower && 2.39 < result_gau$k1_upper){
    k1_gau <- c(k1_gau, T)
  }else{
    k1_gau <- c(k1_gau, F)
  }
  if(21.67 > result_gau$k2_lower && 21.67 < result_gau$k2_upper){
    k2_gau <- c(k2_gau, T)
  }else{
    k2_gau <- c(k2_gau, F)
  }
  if(0.59 > result_gau$t0_lower && 0.59 < result_gau$t0_upper){
    t0_gau <- c(t0_gau, T)
  }else{
    t0_gau <- c(t0_gau, F)
  }

}

sum(k1_gau)/length(k1_gau)
sum(k2_gau)/length(k2_gau)
sum(t0_gau)/length(t0_gau)



#Poisson
bell<-link(t = df$Time, mu = 2.39, k1 = 11.81, k2 = 4.05, t0 = 0.42)
mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)

k1_pois <- NULL
k2_pois <- NULL
t0_pois <- NULL
for(i in 1:1000){
  sim_pois <- rpois(500, lambda = mut)

  result_pois <- scGTM(t=df$Time, y1 = sim_pois, marginal="Poisson", iter = 200, seed = i)
  if(11.81 > result_pois$k1_lower && 11.81 < result_pois$k1_upper){
    k1_pois <- c(k1_pois, T)
  }else{
    k1_pois <- c(k1_pois, F)
  }
  if(4.05 > result_pois$k2_lower && 4.05 < result_pois$k2_upper){
    k2_pois <- c(k2_pois, T)
  }else{
    k2_pois <- c(k2_pois, F)
  }
  if(0.42 > result_pois$t0_lower && 0.42 < result_pois$t0_upper){
    t0_pois <- c(t0_pois, T)
  }else{
    t0_pois <- c(t0_pois, F)
  }

}
sum(k1_pois)/length(k1_pois)
sum(k2_pois)/length(k2_pois)
sum(t0_pois)/length(t0_pois)


#NB
bell<-link(t = df$Time, mu = 2.23, k1 = 8.97, k2 = 3.56, t0 = 0.45)
mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)

k1_nb <- NULL
k2_nb <- NULL
t0_nb <- NULL
for(i in 1:1000){
  sim_nb <- rnbinom(500, 1.5, 1 - mut/(mut+1.5))

  result_nb <- scGTM(t=df$Time, y1 = sim_nb, marginal='NB', iter = 200, seed = i)
  if(8.97 > result_nb$k1_lower && 8.97 < result_nb$k1_upper){
    k1_nb <- c(k1_nb, T)
  }else{
    k1_nb <- c(k1_nb, F)
  }
  if(3.56 > result_nb$k2_lower && 3.56 < result_nb$k2_upper){
    k2_nb <- c(k2_nb, T)
  }else{
    k2_nb <- c(k2_nb, F)
  }
  if(0.45 > result_nb$t0_lower && 0.45 < result_nb$t0_upper){
    t0_nb <- c(t0_nb, T)
  }else{
    t0_nb <- c(t0_nb, F)
  }

}
sum(k1_nb)/length(k1_nb)
sum(k2_nb)/length(k2_nb)
sum(t0_nb)/length(t0_nb)








