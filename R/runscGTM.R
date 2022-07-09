#' Estimate Parameters in Single-cell Gene Expression Generalized Trend Model on a List of Genes
#'
#' @param gene.vec A vector of integers, indicating targeted genes' index in the data.
#' @param t A numeric vector of the input normalized pseudotime data of a given gene,
#' length equals the numbers of cells
#' @param y1 A vector of integers, representing the input expression counts of corresponding lists of genes,
#' number of rows equals the numbers of cells, number of columns equals the numbers of targeted genes
#' @param gene_name A vector of strings, indicates the genes' name used in the model
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
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom tidyr unnest
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment counts
#' @export runscGTM
#'
#'
#' @author Shiyu Ma
#'
#' @examples
#' data("df")
#' res <- scGTM::runscGTM(gene.vec=1:5, t=df$Time, y1=df[,3:7],gene_name=colnames(df[,3:7]))
#'
#' data("sce")
#' t_sce<-SummarizedExperiment::colData(sce)$pseudotime
#' d_sce<-t(SingleCellExperiment::counts(sce)[1:5,])
#' sce_sort<-cbind('Time'= t_sce,d_sce)
#' sce_sort <- tibble::as_tibble(sce_sort)
#' sce_sort <- sce_sort[order(sce_sort$Time),]
#' sce_sort <- tibble::rownames_to_column(sce_sort, "Index")
#' if(class(sce)[1]=="SingleCellExperiment"){t<-sce_sort$Time;y1<-sce_sort[,3:7]}
#' name<-rownames(SingleCellExperiment::counts(sce)[1:5,])
#' res1 <- scGTM::runscGTM(gene.vec=1:5, t=t, y1=y1, gene_name=name)
#'
runscGTM<-function(gene.vec,
                   t,
                   y1,
                   gene_name=NULL,
                   marginal="ZIP",
                   iter_num=50,
                   seed=123,
                   mc.cores = 2){
  set.seed(seed)

  BPPARAM <- BiocParallel::bpparam()
  BPPARAM$workers <- mc.cores

  res <- BiocParallel::bplapply(gene.vec,
                                function(x, ...) {
                                  cur_res <- tryCatch(expr = scGTM::scGTM(gene_index=x,
                                                                  t=t,
                                                                  y1=y1[,x][[1]],
                                                                  gene_name=gene_name[x],
                                                                  marginal=marginal,
                                                                  iter_num=iter_num),
                                                      error = function(e) {
                                                        list(negative_log_likelihood = NA,
                                                             mu = NA,
                                                             k1 = NA,
                                                             k2 = NA,
                                                             t0 = NA,
                                                             phi = NA,
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
                                                             Transform = NA)
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
