#' Compute log_likelihood cost function of one gene based on given parameters
#' @description
#' The log_likelihood cost of fitting a gene's expression counts and pseudotime into a
#' self-defined marginal distribution and parameters
#' @param b A numeric vector of the parameters of the predefined distribution,
#' length equals 6 or 7 depending on distribution(take 6 for Poisson or ZIP, 7 for NB and ZINB)
#' @param y A integer vector of the input expression counts data of a given gene,
#' length equals the numbers of cells
#' @param t A numeric vector of the input normalized pseodotime data of a given gene,
#' length equals the numbers of cells
#' @param marginal A string of the distribution name. One of \code{Poisson}, \code{ZIP}, \code{NB} and \code{ZINB}.
#'
#' @return The log_likelihood cost of the gene
#'
#' @importFrom stats dnbinom dpois
#' @export pso_obj_fct
#'
#' @examples
#' b<-c(0.1,0.2,0.3,0.4,0.5,0.6)
#' y<-floor(runif(100, min = 0, max = 20))
#' t<-runif(100, min = 0, max = 1)
#' marginal<-"ZIP"
#' pso_obj_fct(b,y=y, t=t, marginal=marginal)
#'
#' @author Shiyu Ma
#'
pso_obj_fct<-function(b, y, t, marginal){
  d<- dim(b)

  cost <- 0
  if(marginal == "ZINB"){
    mu<-b[1]
    k1<-b[2]
    k2<-b[3]
    t0<-b[4]
    phi<-b[5]
    alpha<-b[6]
    beta<-b[7]
    cost <- single_gene_log_likelihood_ZINB(y, t, mu, k1, k2, t0, phi, alpha, beta)
  }
  else if(marginal == "NB"){
    mu<-b[1]
    k1<-b[2]
    k2<-b[3]
    t0<-b[4]
    phi<-b[5]
    alpha<-b[6]
    beta<-b[7]
    cost <- single_gene_log_likelihood_NB(y, t, mu, k1, k2, t0, phi)
  }
  else if(marginal == "ZIP"){
    mu<-b[1]
    k1<-b[2]
    k2<-b[3]
    t0<-b[4]
    alpha<-b[5]
    beta<-b[6]
    cost <- single_gene_log_likelihood_ZIP(y, t, mu, k1, k2, t0, alpha, beta)
  }
  else{
    mu<-b[1]
    k1<-b[2]
    k2<-b[3]
    t0<-b[4]
    alpha<-b[5]
    beta<-b[6]
    cost <- single_gene_log_likelihood_Poisson(y, t, mu, k1, k2, t0)
  }
  cost
}

link<-function(t, mu, k1, k2, t0){
  part1<-mu * exp(- abs(k1) * (t - t0) ** 2) * (sign(k1) + (k1 == 0))
  part2<-mu * exp(- abs(k2) * (t - t0) ** 2) * (sign(k2) + (k2 == 0))

  part1 * (t <= t0) + part2 * (t > t0)
}

## Poisson
single_gene_log_likelihood_Poisson<-function(y, t, mu, k1, k2, t0){
  bell<-link(t, mu, k1, k2, t0)
  mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)

  sum(-log(dpois(y, lambda=mut) + 1e-300))
}

## Zero-inflated Poisson
single_gene_log_likelihood_ZIP<-function(y, t, mu, k1, k2, t0, alpha, beta){
  bell<-link(t, mu, k1, k2, t0)
  mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)
  cache<-dpois(y, lambda=mut) + 1e-300

  ## Zero-inflation
  p <- 1 / (1 + exp(alpha * log(mut) + beta))

  sum(- log(cache * (1 - p) + p * (y == 0)))
}

## Negative Binomial
single_gene_log_likelihood_NB<-function(y, t, mu, k1, k2, t0, phi){
  bell<-link(t, mu, k1, k2, t0)
  mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)

  phi <- ifelse(phi>1, phi, 1)
  p0 <- mut / (mut + phi)
  cache <- dnbinom(y, phi, 1 - p0) + 1e-300

  sum(-log(cache))
}


## Zero-inflated Negative Binomial
single_gene_log_likelihood_ZINB<-function(y, t, mu, k1, k2, t0, phi, alpha, beta){
  bell<-link(t, mu, k1, k2, t0)
  mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)

  phi <- ifelse(phi>1, phi, 1)
  p0 <- phi / (mut + phi)
  cache <- dnbinom(y, phi, p0) + 1e-300

  ## Zero-inflation
  p <- 1 / (1 + exp(alpha * log(mut) + beta))

  sum(- log(cache * (1 - p) + p * (y == 0)))
}


#' Compute Fisher information matrix of interested parameters
#'
#' @param t A numeric vector of the input normalized pseodotime data of a given gene,
#' length equals the numbers of cells
#' @param para A numeric vector of the Estimated parameters of the predefined distribution,
#' length equals 6 or 7 depending on distribution(take 6 for Poisson or ZIP, 7 for NB and ZINB)
#' @param marginal A string of the distribution name. One of \code{Poisson}, \code{ZIP}, \code{NB} and \code{ZINB}.
#'
#' @return A 4x4 Fisher information matrix.
#' @export Fisher_info
#'
#' @examples
#' para<-c(0.1,0.2,0.3,0.4,0.5,0.6)
#' t<-runif(100, min = 0, max = 1)
#' marginal<-"ZIP"
#' Fisher_info(t=t, para=para,marginal=marginal)
#'
Fisher_info <- function(t, para, marginal){
  if (marginal == "ZIP" | marginal == "Poisson"){
    mu_fit <- para[1]
    k1_fit <- para[2]
    k2_fit <- para[3]
    t0_fit <- para[4]
    alpha_fit <- para[5]
    beta_fit <- para[6]
  }else{
    mu_fit <- para[1]
    k1_fit <- para[2]
    k2_fit <- para[3]
    t0_fit <- para[4]
    phi_fit <- para[5]
    alpha_fit <- para[6]
    beta_fit <- para[7]
  }

  log_mut_fit <- link(t, mu_fit, k1_fit, k2_fit, t0_fit)

  mut <- ifelse(exp(log_mut_fit)>0.1, exp(log_mut_fit), 0.1)

  t0_deri <- 2 * (k1_fit * (t <= t0_fit) + k2_fit * (t > t0_fit)) * (t - t0_fit) * log_mut_fit  # * (1 + 1/mut)
  k1_deri <- (t - t0_fit) ** 2 * log_mut_fit * (t <= t0_fit)  # * (1 + 1/mut)
  k2_deri <- (t - t0_fit) ** 2 * log_mut_fit * (t > t0_fit)  # * (1 + 1/mut)
  mu_deri <- log_mut_fit / mu_fit  # * (1 + 1/mut)

  if (marginal == "Poisson"){
    t0_deri <- t0_deri*(1 + 1 / mut)
    k1_deri <- k1_deri*(1 + 1 / mut)
    k2_deri <- k2_deri*(1 + 1 / mut)
    mu_deri <- mu_deri*(1 + 1 / mut)
  }else if (marginal == "ZIP"){
    t0_deri <- t0_deri*(1 + 1 / mut)
    k1_deri <- k1_deri*(1 + 1 / mut)
    k2_deri <- k2_deri*(1 + 1 / mut)
    mu_deri <- mu_deri*(1 + 1 / mut)
  }else{
    t0_deri <- t0_deri*(phi_fit * (mut + 1) / (mut + phi_fit))
    k1_deri <- k1_deri*(phi_fit * (mut + 1) / (mut + phi_fit))
    k2_deri <- k2_deri*(phi_fit * (mut + 1) / (mut + phi_fit))
    mu_deri <- mu_deri*(phi_fit * (mut + 1) / (mut + phi_fit))
  }
  cache <- rbind(t0_deri, k1_deri, k2_deri, mu_deri)
  cache %*% t(cache)
}
