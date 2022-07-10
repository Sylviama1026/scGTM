#' Estimate Parameters in Single-cell Gene Expression Generalized Trend Model
#' @description
#' The model fits a gene's expression counts and its normalized pseudotime into one of the four marginal distributions.
#' This function estimates corresponding parameters and confidence intervals in the corresponding marginal distributions.
#'
#' @param gene_index A single integer vector, indicates which gene to estimate in the model,start from 1,
#' default=100
#' @param t A numeric vector of the input normalized pseudotime data of a given gene,
#' length equals the numbers of cells
#' @param y1 A vector of integers, representing the input expression counts of a given gene,
#' length equals the numbers of cells
#' @param gene_name A single string vector, indicates the gene name used in the model,
#' default=NULL
#' @param marginal A string of the distribution name. One of \code{Poisson}, \code{ZIP}, \code{NB}, \code{ZINB}, and  \code{Gaussian}.
#' default=\code{ZIP}
#' @param iter_num A single integer vector, indicates max number of iteration used in the PSO algorithm
#' that estimates model parameters
#' @param seed A numeric variable of the random seed, affecting parametric fitting of the marginal distribution.
#' default=123
#' @param k_design A single positive integer, indicates the number of variables in Design matrix,
#' default=NULL
#' @param Design_X A numerical matrix whose number of rows equals the length of y1,
#' number of columns equals k_design, default=NULL
#'
#' @return The log_likelihood cost, estimated parameters, and their confidence intervals of a gene
#'
#' @importFrom stats dnbinom dpois optim runif cor dnorm
#' @export scGTM
#'
#' @examples
#' y1<-floor(runif(100, min = 0, max = 20))
#' t<-runif(100, min = 0, max = 1)
#' marginal<-"ZIP"
#' scGTM(gene_index=1, t=t, y1=y1, marginal=marginal)
#'
#' data("df")
#' t<-df$Time
#' marginal<-"ZIP"
#' y1<-df$Gene11
#' scGTM(gene_index=11, t=t, y1=y1, marginal=marginal)
#'
#' data("df")
#' t<-df$Time
#' marginal<-"Gaussian"
#' y1<-df$Gene1
#' scGTM(gene_index=1, t=t, y1=y1, marginal=marginal)
#'
#' @author Shiyu Ma, Lehan Zou
#'
scGTM<-function(gene_index = 100, t, y1, gene_name=NULL, marginal="ZIP", iter_num=50, seed=123,
                k_design=NULL, Design_X=NULL){
  #This function automatically determines Hill- or Valley- trend
  ## Flag calculation
  flag <- (cor(t[t<0.5], y1[t<0.5]) < 0) && (cor(t[t>0.5], y1[t>0.5]) > 0) #slope of 1st half and 2nd half
  cat("The need of transformation: " , flag)  #flag=False=hill

  ## ESTIMATION
  cat("\nWe are estimating gene", gene_index,"with marginal",marginal,".\n")

  #transformation if valley
  if(flag){
    y1 <- log(y1+1)
    y1 <- -y1+max(y1) #enforce b=0, align valley with hill estimation
    y1 <- floor(exp(y1)-1)
  }

  #Add Design Matrix
  if(!is.null(k_design)){
    est<-estimation(y1, t, marginal, seed, iter_num, k_design, Design_X)
  }else{
    est<-estimation(y1, t, marginal, seed, iter_num)
  }

  gcost <- est[[1]]
  gbest <- est[[2]]

  result<-NULL
  result['negative_log_likelihood'] <- gcost

  if (gcost > 1e2){
    paste("Best negative log-likelihood: ", round(gcost, 2))
  }else{
    paste("Algorithm fails to find reasonable estimation.")
  }

  if (marginal == "ZIP"){
    mu <- gbest[1]; k1 <- gbest[2]
    k2 <- gbest[3]; t0 <- gbest[4]; phi <- NA; sd<-NA
    alpha <- gbest[5]; beta <- gbest[6]

    cat("Best parameter estimation:\n mu , k1 , k2 , t0 , p:\n",round(gbest, 2))
  }else if(marginal == "ZINB"){ #use 7 para,
    gbest[6] <- ifelse(gbest[6]>1,gbest[6], 1) #alpha>1
    mu <- gbest[1]; k1 <- gbest[2];k2 <- gbest[3]
    t0 <- gbest[4]; phi <- gbest[5]; sd<-NA
    alpha <- gbest[6]; beta <- gbest[7]

    cat("Best parameter estimation:\n mu , k1 , k2 , t0 , phi, p:\n",round(gbest, 2))
  }else if(marginal == "Poisson"){ #use 4 para
    mu <- gbest[1]; k1 <- gbest[2];k2 <- gbest[3]
    t0 <- gbest[4]; phi <- NA; sd<-NA; alpha<-NA; beta<-NA

    cat("Best parameter estimation:\n mu , k1 , k2 , t0:\n",round(gbest[-length(gbest)], 2))
  }else if(marginal == "NB"){ #use 5 para
    gbest[6] <- ifelse(gbest[6]>1,gbest[6], 1) #alpha>1
    mu <- gbest[1]; k1 <- gbest[2];k2 <- gbest[3]
    t0 <- gbest[4]; phi <- gbest[5]; sd<-NA; alpha<-NA; beta<-NA

    cat("Best parameter estimation:\n mu , k1 , k2 , t0, phi:\n",round(gbest[-length(gbest)], 2))
  }else{ #Gaussian, use 5 para
    mu <- gbest[1]; k1 <- gbest[2];k2 <- gbest[3]
    t0 <- gbest[4]; phi<-NA; sd <- gbest[5]; alpha<-NA; beta<-NA

    cat("Best parameter estimation:\n mu , k1 , k2, t0, sd:\n",round(gbest[-length(gbest)], 2))
  }

  ## FISHER INFORMATION
  fisher<-inference(t, gbest, marginal)[[1]] #4x4 matrix
  var<-inference(t, gbest, marginal)[[2]] #4x4 matrix
  t0_lower<-inference(t, gbest, marginal)[[3]]
  t0_upper<-inference(t, gbest, marginal)[[4]]

  cat("\nThe 95% confidence interval of the activation time t0:\n" ,
      "t0 : (" , t0_lower, ", " , t0_upper , ")\n")
  if(length(dim(var))>1){
    t0_std <- sqrt(var[1,1])
    k1_lower <- round(gbest[2] - 1.96 * sqrt(var[2, 2]), 3)
    k1_upper <- round(gbest[2] + 1.96 * sqrt(var[2, 2]), 3)
    k2_lower <- round(gbest[3] - 1.96 * sqrt(var[3, 3]), 3)
    k2_upper <- round(gbest[3] + 1.96 * sqrt(var[3, 3]), 3)
    mu_lower <- round(gbest[1] - 1.96 * sqrt(var[4, 4]), 3)
    mu_upper <- round(gbest[1] + 1.96 * sqrt(var[4, 4]), 3)

    cat("\nThe 95% CIs for activation strength k1 and k2:\n" ,
        "k1 : (" , k1_lower , ", " , k1_upper , ")\n",
        "k2 : (" , k2_lower , ", " , k2_upper , ")\n")

    k1_std <- sqrt(var[2, 2])
    k2_std <- ifelse(is.na(k2_lower),NA,sqrt(var[3, 3]))
    mu_std <- sqrt(var[4, 4])
    Fisher <- 'Non-singular'
  }else{
    var <- fisher
    var[1, 1] <- 1 / (var[1, 1] + 1e-100)
    var[2, 2] <- 1 / (var[2, 2] + 1e-100)
    var[3, 3] <- 1 / (var[3, 3] + 1e-100)
    var[4, 4] <- 1 / (var[4, 4] + 1e-100)

    t0_std <- sqrt(var[1,1])
    k1_lower <- round(gbest[2] - 1.96 * sqrt(var[2, 2]), 3)
    k1_upper <- round(gbest[2] + 1.96 * sqrt(var[2, 2]), 3)
    k2_lower <- round(gbest[3] - 1.96 * sqrt(var[3, 3]), 3)
    k2_upper <- round(gbest[3] + 1.96 * sqrt(var[3, 3]), 3)
    mu_lower <- round(gbest[1] - 1.96 * sqrt(var[4, 4]), 3)
    mu_upper <- round(gbest[1] + 1.96 * sqrt(var[4, 4]), 3)

    cat("\nThe 95% CIs for activation strength k1 and k2:\n" ,
        "k1 : (" , k1_lower , ", " , k1_upper , ")\n",
        "k2 : (" , k2_lower , ", " , k2_upper , ")\n")

    k1_std <- sqrt(var[2, 2])
    k2_std <- ifelse(is.na(k2_lower),NA,sqrt(var[3, 3]))
    mu_std <- sqrt(var[4, 4])
    Fisher <- 'Singular'
  }
  Transform <- as.integer(flag)

  #Add Design Matrix
  if(!is.null(k_design)){
    if(marginal %in% c("Poisson", "ZIP") ){
      Design_para<-gbest[7:length(gbest)]
    } else{ #marginal %in% c("NB", "ZINB", "Gaussian")
      Design_para<-gbest[8:length(gbest)]
    }
  }else{
    Design_para<-NA
  }

  list(negative_log_likelihood = gcost,
       mu = mu,
       k1 = k1,
       k2 = k2,
       t0 = t0,
       phi = phi,
       sd = sd, #Gaussian
       alpha = alpha,
       beta = beta,
       t0_lower = t0_lower,
       t0_upper = t0_upper,
       t0_std = t0_std,
       k1_lower = k1_lower,
       k1_upper = k1_upper,
       k1_std = k1_std,
       k2_lower = k2_lower,
       k2_upper = k2_upper,
       k2_std = k2_std,
       mu_lower = mu_lower,
       mu_upper = mu_upper,
       mu_std = mu_std,
       Fisher = Fisher,
       Transform = Transform,
       Design_para = Design_para)
}


# I. Objective_function
##Compute log_likelihood cost function of one gene based on given parameters
pso_obj_fct<-function(b, y, t, marginal,Design_X=NULL){
  cost <- 0
  if(marginal == "ZINB"){
    mu<-b[1]
    k1<-b[2]
    k2<-b[3]
    t0<-b[4]
    phi<-b[5]
    alpha<-b[6]
    beta<-b[7]
    if(!is.null(Design_X)){
      k_design_para<-b[8:length(b)]
      cost <- single_gene_log_likelihood_ZINB(y, t, mu, k1, k2, t0, phi, alpha, beta, Design_X, k_design_para)
    }else{
      cost <- single_gene_log_likelihood_ZINB(y, t, mu, k1, k2, t0, phi, alpha, beta)
    }
  }
  else if(marginal == "NB"){
    mu<-b[1]
    k1<-b[2]
    k2<-b[3]
    t0<-b[4]
    phi<-b[5]
    alpha<-b[6]
    beta<-b[7]
    if(!is.null(Design_X)){
      k_design_para<-b[8:length(b)]
      cost <- single_gene_log_likelihood_NB(y, t, mu, k1, k2, t0, phi, Design_X, k_design_para)
    }else{
      cost <- single_gene_log_likelihood_NB(y, t, mu, k1, k2, t0, phi)
    }
  }
  else if(marginal == "ZIP"){
    mu<-b[1]
    k1<-b[2]
    k2<-b[3]
    t0<-b[4]
    alpha<-b[5]
    beta<-b[6]
    if(!is.null(Design_X)){
      k_design_para<-b[7:length(b)]
      cost <- single_gene_log_likelihood_ZIP(y, t, mu, k1, k2, t0, alpha, beta, Design_X, k_design_para)
    }else{
      cost <- single_gene_log_likelihood_ZIP(y, t, mu, k1, k2, t0, alpha, beta)
    }
  }
  else if(marginal == "Poisson"){
    mu<-b[1]
    k1<-b[2]
    k2<-b[3]
    t0<-b[4]
    alpha<-b[5]
    beta<-b[6]
    if(!is.null(Design_X)){
      k_design_para<-b[7:length(b)]
      cost <- single_gene_log_likelihood_Poisson(y, t, mu, k1, k2, t0, Design_X, k_design_para)
    }else{
      cost <- single_gene_log_likelihood_Poisson(y, t, mu, k1, k2, t0)
    }
  }
  else{ #Gaussian
    mu<-b[1]
    k1<-b[2]
    k2<-b[3]
    t0<-b[4]
    sd<-b[5]
    alpha<-b[6]
    beta<-b[7]
    if(!is.null(Design_X)){
      k_design_para<-b[8:length(b)]
      cost <- single_gene_log_likelihood_Gaussian(y, t, mu, k1, k2, t0, sd, Design_X, k_design_para)
    }else{
      cost <- single_gene_log_likelihood_Gaussian(y, t, mu, k1, k2, t0, sd)
    }
  }
  cost
}

link<-function(t, mu, k1, k2, t0, Design_X=NULL, k_design_para=NULL){
  part1<-mu * exp(- abs(k1) * (t - t0) ** 2) * (sign(k1) + (k1 == 0))
  part2<-mu * exp(- abs(k2) * (t - t0) ** 2) * (sign(k2) + (k2 == 0))

  if(!is.null(Design_X)){
    out<-part1 * (t <= t0) + part2 * (t > t0) + Design_X %*% k_design_para
  }else{
    out<-part1 * (t <= t0) + part2 * (t > t0)
  }
  out
}

## Gaussian
single_gene_log_likelihood_Gaussian<-function(y, t, mu, k1, k2, t0, sd, Design_X=NULL, k_design_para=NULL){
  bell<-link(t, mu, k1, k2, t0, Design_X, k_design_para)
  mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)

  sd <- ifelse(sd>0, sd, 1e-300)

  sum(-log(dnorm(y,  mean=mut, sd=sd) + 1e-300))
}

## Poisson
single_gene_log_likelihood_Poisson<-function(y, t, mu, k1, k2, t0, Design_X=NULL, k_design_para=NULL){
  bell<-link(t, mu, k1, k2, t0, Design_X, k_design_para)
  mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)

  sum(-log(dpois(y, lambda=mut) + 1e-300))
}

## Zero-inflated Poisson
single_gene_log_likelihood_ZIP<-function(y, t, mu, k1, k2, t0, alpha, beta, Design_X=NULL, k_design_para=NULL){
  bell<-link(t, mu, k1, k2, t0, Design_X, k_design_para)
  mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)
  cache<-dpois(y, lambda=mut) + 1e-300

  ## Zero-inflation
  p <- 1 / (1 + exp(alpha * log(mut) + beta))

  sum(- log(cache * (1 - p) + p * (y == 0)))
}

## Negative Binomial
single_gene_log_likelihood_NB<-function(y, t, mu, k1, k2, t0, phi, Design_X=NULL, k_design_para=NULL){
  bell<-link(t, mu, k1, k2, t0, Design_X, k_design_para)
  mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)

  phi <- ifelse(phi>1, phi, 1)
  p0 <- mut / (mut + phi)
  cache <- dnbinom(y, phi, 1 - p0) + 1e-300

  sum(-log(cache))
}


## Zero-inflated Negative Binomial
single_gene_log_likelihood_ZINB<-function(y, t, mu, k1, k2, t0, phi, alpha, beta, Design_X=NULL, k_design_para=NULL){
  bell<-link(t, mu, k1, k2, t0, Design_X, k_design_para)
  mut<-ifelse(exp(bell)>0.1, exp(bell), 0.1)

  phi <- ifelse(phi>1, phi, 1)
  p0 <- phi / (mut + phi)
  cache <- dnbinom(y, phi, p0) + 1e-300

  ## Zero-inflation
  p <- 1 / (1 + exp(alpha * log(mut) + beta))

  sum(- log(cache * (1 - p) + p * (y == 0)))
}




# II. Fisher information matrix
##Compute Fisher information matrix of interested parameters
Fisher_info <- function(t, para, marginal){
  if (marginal == "ZIP" | marginal == "Poisson"){
    mu_fit <- para[1]
    k1_fit <- para[2]
    k2_fit <- para[3]
    t0_fit <- para[4]
    alpha_fit <- para[5]
    beta_fit <- para[6]
  }else if(marginal == "Gaussian"){
    mu_fit <- para[1]
    k1_fit <- para[2]
    k2_fit <- para[3]
    t0_fit <- para[4]
    sd_fit <- para[5]
    alpha_fit <- para[6]
    beta_fit <- para[7]
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

  if (marginal == "ZIP" | marginal == "Poisson" | marginal == "Gaussian"){
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



# III. Parameters Estimation
##Using PSO algorithm to estimate parameters within the fitting distribution
estimation<-function(y, t, marginal, seed, iter=50, k_design=NULL, Design_X=NULL){
  n<-30
  if(marginal %in% c("Poisson", "ZIP") ){
    d<-6
    bounds<-list(c(log(min(y+1)),-99,-99,min(t),-99,-99),
                 c(log(max(y)),99,99,max(t),99,99))
  } else if(marginal %in% c("NB", "ZINB")){
    d<-7
    bounds<-list(c(log(min(y+1)),-99,-99,min(t), 1, 0, 0),
                 c(log(max(y+1)),99,99,max(t), 100, 100, 100))
  } else if(marginal=="Gaussian"){
    d<-7
    bounds<-list(c(log(min(y+1)),-99,-99,min(t), 1e-300, 0, 0),
                 c(log(max(y+1)),99,99,max(t), 100, 100, 100))
  }else{
    stop("Enter a valid marginal distribution: [NB, ZINB, Poisson, ZIP, Gaussian]!")
  }

  #Add Design matrix
  if(!is.null(k_design)){
    d_all<-d+k_design
    bounds<-list(c(bounds[[1]],rep(0,k_design)),
                 c(bounds[[2]],rep(99,k_design))) #KEY！！
  }else{
    d_all<-d
  }
  b <- runif(d_all)

  # Set-up hyperparameters and correct initial position
  options <- list('c.p'= 1.2, 'c.g'= 0.3, 'w'= 0.9 , 'maxit'=iter,'s'=n)#,'trace'=10,'trace.stats'=T)
  b[1] <- log(mean(y) + 1)
  b[d] <- 0.1 #set beta=0.1
  b[2] <- b[2]+5
  b[3] <- b[3]+5
  if (d == 7){ #phi or sd
    b[5] <- b[5]+1
  }

  #Add Design Matrix
  if(!is.null(k_design)){
    gmodel<-my_psoptim(par=b, fn = pso_obj_fct, y=y, t=t, marginal=marginal,Design_X=Design_X, seed=seed,
                       lower = bounds[[1]], upper = bounds[[2]], control=options)
  }else{
    gmodel<-my_psoptim(par=b, fn = pso_obj_fct, y=y, t=t, marginal=marginal,seed=seed,
                       lower = bounds[[1]], upper = bounds[[2]], control=options)
  }
  gbest<-gmodel$par
  gcost<-gmodel$value

  list(gcost,gbest)
}




# IV. Parameters Inference
##Compute estimated Fisher information matrix and covariance matrix for estimated parameters
inference <- function(t, para, marginal){
  fisher <- Fisher_info(t, para, marginal)
  tryCatch(
    expr = {
      cov <- solve(fisher)
      lower <- round(para[4] - 1.96 * sqrt(cov[1, 1]), 3)
      upper <- round(para[4] + 1.96 * sqrt(cov[1, 1]), 3)
    }, error = function(e){
      print("The Fisher information matrix is singular: estimated t0 is either less than 0 or greater than 1.\n",
            "Calculating the sub-Fisher information matrix of t0 instead.\n")
      cov <- 1 / fisher[1, 1]
      lower <- round(para[4] - 1.96 * sqrt(cov), 3)
      upper <- round(para[4] + 1.96 * sqrt(cov), 3)

    })
  list(fisher, cov, lower, upper)
}




###### my_psoptim function used in estimation#######
#### This function is modified from Claus Bendtsen's psoptim function in the pso package.
my_psoptim <- function (par, fn, gr = NULL, ..., lower=-1, upper=1, seed, Design_X=NULL,
                        control = list()) {
  fn1 <- function(par) fn(par, ...)/p.fnscale
  mrunif <- function(n,m,lower,upper) {
    return(matrix(runif(n*m,0,1),nrow=n,ncol=m)*(upper-lower)+lower)
  }
  norm <- function(x) sqrt(sum(x*x))
  rsphere.unif <- function(n,r) {
    temp <- runif(n)
    return((runif(1,min=0,max=r)/norm(temp))*temp)
  }
  svect <- function(a,b,n,k) {
    temp <- rep(a,n)
    temp[k] <- b
    return(temp)
  }
  mrsphere.unif <- function(n,r) {
    m <- length(r)
    temp <- matrix(runif(n*m),n,m)
    return(temp%*%diag(runif(m,min=0,max=r)/apply(temp,2,norm)))
  }
  npar <- length(par)
  lower <- as.double(rep(lower, ,npar))
  upper <- as.double(rep(upper, ,npar))
  con <- list(trace = 0, fnscale = 1, maxit = 1000L, maxf = Inf,
              abstol = -Inf, reltol = 0, REPORT = 10,
              s = NA, k = 3, p = NA, w = 1/(2*log(2)),
              c.p = .5+log(2), c.g = .5+log(2), d = NA,
              v.max = NA, rand.order = TRUE, max.restart=Inf,
              maxit.stagnate = Inf,
              vectorize=FALSE, hybrid = FALSE, hybrid.control = NULL,
              trace.stats = FALSE, type = "SPSO2007")
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  ## Argument error checks
  if (any(upper==Inf | lower==-Inf))
    stop("fixed bounds must be provided")

  p.type <- pmatch(con[["type"]],c("SPSO2007","SPSO2011"))-1
  if (is.na(p.type)) stop("type should be one of \"SPSO2007\", \"SPSO2011\"")

  p.trace <- con[["trace"]]>0L # provide output on progress?
  p.fnscale <- con[["fnscale"]] # scale funcion by 1/fnscale
  p.maxit <- con[["maxit"]] # maximal number of iterations
  p.maxf <- con[["maxf"]] # maximal number of function evaluations
  p.abstol <- con[["abstol"]] # absolute tolerance for convergence
  p.reltol <- con[["reltol"]] # relative minimal tolerance for restarting
  p.report <- as.integer(con[["REPORT"]]) # output every REPORT iterations
  p.s <- ifelse(is.na(con[["s"]]),ifelse(p.type==0,floor(10+2*sqrt(npar)),40),
                con[["s"]]) # swarm size
  p.p <- ifelse(is.na(con[["p"]]),1-(1-1/p.s)^con[["k"]],con[["p"]]) # average % of informants
  p.w0 <- con[["w"]] # exploitation constant
  if (length(p.w0)>1) {
    p.w1 <- p.w0[2]
    p.w0 <- p.w0[1]
  } else {
    p.w1 <- p.w0
  }
  p.c.p <- con[["c.p"]] # local exploration constant
  p.c.g <- con[["c.g"]] # global exploration constant
  p.d <- ifelse(is.na(con[["d"]]),norm(upper-lower),con[["d"]]) # domain diameter
  p.vmax <- con[["v.max"]]*p.d # maximal velocity
  p.randorder <- as.logical(con[["rand.order"]]) # process particles in random order?
  p.maxrestart <- con[["max.restart"]] # maximal number of restarts
  p.maxstagnate <- con[["maxit.stagnate"]] # maximal number of iterations without improvement
  p.vectorize <- as.logical(con[["vectorize"]]) # vectorize?
  if (is.character(con[["hybrid"]])) {
    p.hybrid <- pmatch(con[["hybrid"]],c("off","on","improved"))-1
    if (is.na(p.hybrid)) stop("hybrid should be one of \"off\", \"on\", \"improved\"")
  } else {
    p.hybrid <- as.integer(as.logical(con[["hybrid"]])) # use local BFGS search
  }
  p.hcontrol <- con[["hybrid.control"]] # control parameters for hybrid optim
  if ("fnscale" %in% names(p.hcontrol))
    p.hcontrol["fnscale"] <- p.hcontrol["fnscale"]*p.fnscale
  else
    p.hcontrol["fnscale"] <- p.fnscale
  p.trace.stats <- as.logical(con[["trace.stats"]]) # collect detailed stats?

  if (p.trace) {
    message("S=",p.s,", K=",con[["k"]],", p=",signif(p.p,4),", w0=",
            signif(p.w0,4),", w1=",
            signif(p.w1,4),", c.p=",signif(p.c.p,4),
            ", c.g=",signif(p.c.g,4))
    message("v.max=",signif(con[["v.max"]],4),
            ", d=",signif(p.d,4),", vectorize=",p.vectorize,
            ", hybrid=",c("off","on","improved")[p.hybrid+1])
    if (p.trace.stats) {
      stats.trace.it <- c()
      stats.trace.error <- c()
      stats.trace.f <- NULL
      stats.trace.x <- NULL
    }
  }
  ## Initialization
  if (p.reltol!=0) p.reltol <- p.reltol*p.d
  if (p.vectorize) {
    lowerM <- matrix(lower,nrow=npar,ncol=p.s)
    upperM <- matrix(upper,nrow=npar,ncol=p.s)
  }
  #################################################MODIFY!!!!!!!!!!!!!!!!
  set.seed(seed)
  if(!is.null(Design_X)){
    d_npar<-npar-ncol(Design_X)
  }else{
    d_npar<-npar
  }
  X <- mrunif(npar,p.s,0,1)
  if (!any(is.na(par)) && all(par>=lower) && all(par<=upper)){ #specify initial pos
    X[1,] <- par[1]
    X[d_npar,] <- 0.1
    X[2,] <- X[2,]+5
    X[3,] <- X[3,]+5
    if(par[4]>0.5){
      X[4,] <- runif(1, min = 0.5, max = 1)
    }
    else{
      X[4,] <- runif(1, min = 0, max = 0.5)
    }
    if (d_npar == 7){
      X[5,] <- X[5,]+1
    }
  }
  #################################################
  if (p.type==0) {
    V <- (mrunif(npar,p.s,lower,upper)-X)/2
  } else { ## p.type==1
    V <- matrix(runif(npar*p.s,min=as.vector(lower-X),max=as.vector(upper-X)),npar,p.s)
    p.c.p2 <- p.c.p/2 # precompute constants
    p.c.p3 <- p.c.p/3
    p.c.g3 <- p.c.g/3
    p.c.pg3 <- p.c.p3+p.c.g3
  }
  if (!is.na(p.vmax)) { # scale to maximal velocity
    temp <- apply(V,2,norm)
    temp <- pmin.int(temp,p.vmax)/temp
    V <- V%*%diag(temp)
  }
  f.x <- apply(X,2,fn1) # first evaluations
  stats.feval <- p.s
  P <- X
  f.p <- f.x
  P.improved <- rep(FALSE,p.s)
  i.best <- which.min(f.p)
  error <- f.p[i.best]
  init.links <- TRUE
  if (p.trace && p.report==1) {
    message("It 1: fitness=",signif(error,4))
    if (p.trace.stats) {
      stats.trace.it <- c(stats.trace.it,1)
      stats.trace.error <- c(stats.trace.error,error)
      stats.trace.f <- c(stats.trace.f,list(f.x))
      stats.trace.x <- c(stats.trace.x,list(X))
    }
  }
  ## Iterations
  stats.iter <- 1
  stats.restart <- 0
  stats.stagnate <- 0
  while (stats.iter<p.maxit && stats.feval<p.maxf && error>p.abstol &&
         stats.restart<p.maxrestart && stats.stagnate<p.maxstagnate) {
    stats.iter <- stats.iter+1
    if (p.p!=1 && init.links) {
      links <- matrix(runif(p.s*p.s,0,1)<=p.p,p.s,p.s)
      diag(links) <- TRUE
    }
    ## The swarm moves
    if (!p.vectorize) {
      if (p.randorder) {
        index <- sample(p.s)
      } else {
        index <- 1:p.s
      }
      for (i in index) {
        if (p.p==1)
          j <- i.best
        else
          j <- which(links[,i])[which.min(f.p[links[,i]])] # best informant
        temp <- (p.w0+(p.w1-p.w0)*max(stats.iter/p.maxit,stats.feval/p.maxf))
        V[,i] <- temp*V[,i] # exploration tendency
        if (p.type==0) {
          V[,i] <- V[,i]+runif(npar,0,p.c.p)*(P[,i]-X[,i]) # exploitation
          if (i!=j) V[,i] <- V[,i]+runif(npar,0,p.c.g)*(P[,j]-X[,i])
        } else { # SPSO 2011
          if (i!=j)
            temp <- p.c.p3*P[,i]+p.c.g3*P[,j]-p.c.pg3*X[,i] # Gi-Xi
          else
            temp <- p.c.p2*P[,i]-p.c.p2*X[,i] # Gi-Xi for local=best
          V[,i] <- V[,i]+temp+rsphere.unif(npar,norm(temp))
        }
        if (!is.na(p.vmax)) {
          temp <- norm(V[,i])
          if (temp>p.vmax) V[,i] <- (p.vmax/temp)*V[,i]
        }
        X[,i] <- X[,i]+V[,i]
        ## Check bounds
        temp <- X[,i]<lower
        if (any(temp)) {
          X[temp,i] <- lower[temp]
          V[temp,i] <- 0
        }
        temp <- X[,i]>upper
        if (any(temp)) {
          X[temp,i] <- upper[temp]
          V[temp,i] <- 0
        }
        ## Evaluate function
        if (p.hybrid==1) {
          temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower,
                        upper=upper,control=p.hcontrol)
          V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
          X[,i] <- temp$par
          f.x[i] <- temp$value
          stats.feval <- stats.feval+as.integer(temp$counts[1])
        } else {
          f.x[i] <- fn1(X[,i])
          stats.feval <- stats.feval+1
        }
        if (f.x[i]<f.p[i]) { # improvement
          P[,i] <- X[,i]
          f.p[i] <- f.x[i]
          if (f.p[i]<f.p[i.best]) {
            i.best <- i
            if (p.hybrid==2) {
              temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower,
                            upper=upper,control=p.hcontrol)
              V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
              X[,i] <- temp$par
              P[,i] <- temp$par
              f.x[i] <- temp$value
              f.p[i] <- temp$value
              stats.feval <- stats.feval+as.integer(temp$counts[1])
            }
          }
        }
        if (stats.feval>=p.maxf) break
      }
    } else {
      if (p.p==1)
        j <- rep(i.best,p.s)
      else # best informant
        j <- sapply(1:p.s,function(i)
          which(links[,i])[which.min(f.p[links[,i]])])
      temp <- (p.w0+(p.w1-p.w0)*max(stats.iter/p.maxit,stats.feval/p.maxf))
      V <- temp*V # exploration tendency
      if (p.type==0) {
        V <- V+mrunif(npar,p.s,0,p.c.p)*(P-X) # exploitation
        temp <- j!=(1:p.s)
        V[,temp] <- V[,temp]+mrunif(npar,sum(temp),0,p.c.p)*(P[,j[temp]]-X[,temp])
      } else { # SPSO 2011
        temp <- j==(1:p.s)
        temp <- P%*%diag(svect(p.c.p3,p.c.p2,p.s,temp))+
          P[,j]%*%diag(svect(p.c.g3,0,p.s,temp))-
          X%*%diag(svect(p.c.pg3,p.c.p2,p.s,temp)) # G-X
        V <- V+temp+mrsphere.unif(npar,apply(temp,2,norm))
      }
      if (!is.na(p.vmax)) {
        temp <- apply(V,2,norm)
        temp <- pmin.int(temp,p.vmax)/temp
        V <- V%*%diag(temp)
      }
      X <- X+V
      ## Check bounds
      temp <- X<lowerM
      if (any(temp)) {
        X[temp] <- lowerM[temp]
        V[temp] <- 0
      }
      temp <- X>upperM
      if (any(temp)) {
        X[temp] <- upperM[temp]
        V[temp] <- 0
      }
      ## Evaluate function
      if (p.hybrid==1) { # not really vectorizing
        for (i in 1:p.s) {
          temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower,
                        upper=upper,control=p.hcontrol)
          V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
          X[,i] <- temp$par
          f.x[i] <- temp$value
          stats.feval <- stats.feval+as.integer(temp$counts[1])
        }
      } else {
        f.x <- apply(X,2,fn1)
        stats.feval <- stats.feval+p.s
      }
      temp <- f.x<f.p
      if (any(temp)) { # improvement
        P[,temp] <- X[,temp]
        f.p[temp] <- f.x[temp]
        i.best <- which.min(f.p)
        if (temp[i.best] && p.hybrid==2) { # overall improvement
          temp <- optim(X[,i.best],fn,gr,...,method="L-BFGS-B",lower=lower,
                        upper=upper,control=p.hcontrol)
          V[,i.best] <- V[,i.best]+temp$par-X[,i.best] # disregards any v.max imposed
          X[,i.best] <- temp$par
          P[,i.best] <- temp$par
          f.x[i.best] <- temp$value
          f.p[i.best] <- temp$value
          stats.feval <- stats.feval+as.integer(temp$counts[1])
        }
      }
      if (stats.feval>=p.maxf) break
    }
    if (p.reltol!=0) {
      d <- X-P[,i.best]
      d <- sqrt(max(colSums(d*d)))
      if (d<p.reltol) {
        X <- mrunif(npar,p.s,lower,upper)
        V <- (mrunif(npar,p.s,lower,upper)-X)/2
        if (!is.na(p.vmax)) {
          temp <- apply(V,2,norm)
          temp <- pmin.int(temp,p.vmax)/temp
          V <- V%*%diag(temp)
        }
        stats.restart <- stats.restart+1
        if (p.trace) message("It ",stats.iter,": restarting")
      }
    }
    init.links <- f.p[i.best]==error # if no overall improvement
    stats.stagnate <- ifelse(init.links,stats.stagnate+1,0)
    error <- f.p[i.best]
    if (p.trace && stats.iter%%p.report==0) {
      if (p.reltol!=0)
        message("It ",stats.iter,": fitness=",signif(error,4),
                ", swarm diam.=",signif(d,4))
      else
        message("It ",stats.iter,": fitness=",signif(error,4))
      if (p.trace.stats) {
        stats.trace.it <- c(stats.trace.it,stats.iter)
        stats.trace.error <- c(stats.trace.error,error)
        stats.trace.f <- c(stats.trace.f,list(f.x))
        stats.trace.x <- c(stats.trace.x,list(X))
      }
    }
  }
  if (error<=p.abstol) {
    msg <- "Converged"
    msgcode <- 0
  } else if (stats.feval>=p.maxf) {
    msg <- "Maximal number of function evaluations reached"
    msgcode <- 1
  } else if (stats.iter>=p.maxit) {
    msg <- "Maximal number of iterations reached"
    msgcode <- 2
  } else if (stats.restart>=p.maxrestart) {
    msg <- "Maximal number of restarts reached"
    msgcode <- 3
  } else {
    msg <- "Maximal number of iterations without improvement reached"
    msgcode <- 4
  }
  if (p.trace) message(msg)
  o <- list(par=P[,i.best],value=f.p[i.best],
            counts=c("function"=stats.feval,"iteration"=stats.iter,
                     "restarts"=stats.restart),
            convergence=msgcode,message=msg)
  if (p.trace && p.trace.stats) o <- c(o,list(stats=list(it=stats.trace.it,
                                                         error=stats.trace.error,
                                                         f=stats.trace.f,
                                                         x=stats.trace.x)))
  return(o)
}

