#' Estimate Parameters in Single-cell Gene Expression Generalized Trend Model on a List of Genes
#'
#' @param t A numeric vector of the input normalized pseudotime data of a given gene,
#' length equals the numbers of cells
#' (If sce is not null, t is a string of gene names to use in the model)
#' @param y A tibble, representing the input expression counts of corresponding lists of genes,
#' number of rows equals the numbers of cells, number of columns equals the numbers of targeted genes,
#' (If sce is not null, y is a SingleCellExperiment object with counts data)
#' @param sce A character vector, indicates the assay name applied to the SingleCellExperiment object
#' default=NULL
#' @param marginal A string of the distribution name. One of \code{Poisson}, \code{ZIP}, \code{NB}, \code{ZINB} and \code{Gaussian}.
#' default=\code{ZIP}
#' @param iter_num A single integer vector, indicates max number of iteration used in the PSO algorithm
#' that estimates model parameters
#' @param seed A numeric variable of the random seed, affecting parametric fitting of the marginal distribution.
#' default=123
#' @param mc.cores Number of cores used for computing.
#'
#' @return A tibble of summary results of genes
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom tibble as_tibble
#' @importFrom tidyr unnest
#' @importFrom SummarizedExperiment assays
#' @export runscGTM
#'
#'
#' @author Shiyu Ma
#'
#' @examples
#' data("df")
#' res <- scGTM::runscGTM(t=df$Time, y=df[,3:5])
#'
#' data("sce")
#' t_sce<-rownames(sce)
#' res_sce <- runscGTM(t=t_sce, y=sce, sce="logcounts", marginal="Gaussian")
#'
runscGTM<-function(t, #gene_names to run or pseudotime
                   y, #object with data or matrix with data
                   sce = NULL, #assay name
                   marginal = "ZIP",
                   iter_num = 50,
                   seed = 123,
                   mc.cores = 2){
  set.seed(seed)

  BPPARAM <- BiocParallel::bpparam()
  BPPARAM$workers <- mc.cores

  if(!is.null(sce)){
    if(class(t)!="numeric"){
      gene_name<-t #all genes' name in y(sce object)
      t<-y$pseudotime #time of all cells
      values<-SummarizedExperiment::assays(y)[[sce]][gene_name,] #value of selected genes and assays
      y<-tibble::as_tibble(t(values))
    }else{
      stop("If sce is not NULL, t must be a string to indicate gene names, not numeric")
    }

  }else{
    gene_name<-colnames(y)
  }

  res <- BiocParallel::bplapply(gene_name, function(x, ...) {
    cur_res <- tryCatch(expr = scGTM(t=t,
                                     y1=y[,x][[1]],
                                     gene_name=x,
                                     marginal=marginal,
                                     iter_num=iter_num,
                                     seed=seed),
                        error = function(e) {
                          list(negative_log_likelihood = NA,
                               mu = NA,
                               k1 = NA,
                               k2 = NA,
                               t0 = NA,
                               phi = NA,
                               sd = NA, #Gaussian
                               alpha = NA,
                               beta = NA,
                               t0_lower = NA,
                               t0_upper = NA,
                               t0_std = NA,
                               k1_lower = NA,
                               k1_upper = NA,
                               k1_std = NA,
                               k2_lower = NA,
                               k2_upper = NA,
                               k2_std = NA,
                               mu_lower = NA,
                               mu_upper = NA,
                               mu_std = NA,
                               Fisher = NA,
                               Transform = NA,
                               Design_para = NA)
                        })
    cur_res
  },
  BPPARAM = BPPARAM)

  res <- simplify2array(res)
  res <- t(res)
  rownames(res) <- gene_name
  res <- tibble::as_tibble(res,rownames="gene") #,.name_repair = 'unique'
  res <- tidyr::unnest(res,cols=colnames(res)[-1])

  res
}
