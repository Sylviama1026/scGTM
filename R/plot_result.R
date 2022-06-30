#' Plot the Fitted Model
#'
#' @param para A numeric vector of the Estimated parameters of the predefined distribution,
#' length equals 6 or 7 depending on distribution(take 6 for Poisson or ZIP, 7 for NB and ZINB)
#' @param t A numeric vector of the input normalized pseudotime data of a given gene,
#' length equals the numbers of cells
#' @param color A string vector of length 4 to define plot color, default=\code{c('red', 'blue', 'orange', 'darkgreen')}
#' @param marginal A string of the distribution name. One of \code{Poisson}, \code{ZIP}, \code{NB} and \code{ZINB}.
#' @param flag A boolean variable, flag=T indicates Valley shape, flag=F indicates Hill shape
#' @param y1 A vector of integers, representing the input expression counts of a given gene after transformation,
#' length equals the numbers of cells, y1=raw if data hs Hill shape(flag=F)
#' @param raw A vector of integers, representing the input expression counts of a given gene,
#' length equals the numbers of cells
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @import cowplot
#' @export plot_result
#'
#' @examples
#' data("df")
#' y1<-df$Gene1
#' t<-df$Time
#' flag<-FALSE
#' para<-c(2.29,3.27,11.79,0.58,30.4,60.82)
#' marginal<-"ZIP"
#' color<-c('red', 'blue', 'orange', 'darkgreen')
#' plot_result(para, t, color, marginal, flag, y1)
#'
#' @author Shiyu Ma, Lehan Zou
plot_result <- function(para, t, color, marginal, flag, y1, raw){
  mu_fit <- para[1]
  k1_fit <- para[2]
  k2_fit <- para[3]
  t0_fit <- para[4]
  log_mut_fit <- link(sort(t), mu_fit, k1_fit, k2_fit, t0_fit)

  if (flag){
    log_mut_fit = -log_mut_fit + log(max(y1) + 1)
  }

  p_fit <- 1 / (1 + exp(para[length(para)] + para[length(para)-1]* exp(log_mut_fit)))
  ZIlog_mut_fit <- ifelse(log_mut_fit + log(1 - p_fit) > 0, log_mut_fit + log(1 - p_fit), 0)

  y_1<-ZIlog_mut_fit - 0.1
  y_2<-p_fit
  coff<-y_2/y_1

  data<-as.data.frame(cbind(t,log_mut_fit,log(y1+1)))

  if (flag){
    y1<-raw
  }

  p <-ggplot(data)+
    geom_point(aes(x = t, y = log(y1+1),color = log(y1+1)))+
    scale_color_gradient(low="blue", high="red",name = "Counts")+
    ylim(min(log(y1+1))-1, max(log(y1+1))+1)+
    geom_line(aes(x= sort(t), y = log_mut_fit),colour=color[1], size=1.2)+
    geom_text(aes(x = max(t)-0.025,
                  y = log(y1+1)[length(log(y1+1))]+0.25,
                  label="Fitted"), size = 3) +
    xlab("Pseudotime") +
    ylab("Expression log(count +1)") +
    theme_bw()

  if(t0_fit <= 1 & t0_fit >= 0){
    p<-p+
      geom_vline(xintercept = t0_fit, linetype="dashed",
                 color = color[3], size=1.5)
    p
  }

  if (marginal == 'ZIP'|marginal == 'ZINB'){
    p1<-p+
      geom_line(aes( x= sort(t), y = ZIlog_mut_fit - 0.1),color=color[2], size=0.8)+
      geom_text(aes(x = max(t)-0.1,
                    y = ZIlog_mut_fit[length(t)]-0.1,
                    label="W/dropout"), size = 3)
    p2<-ggplot(data)+
      geom_line(aes( x= sort(t), y = p_fit),color=color[4], size=1)+
      xlab("Pseudotime") +
      ylab("Dropout Rate") +
      theme_bw()
    cowplot::plot_grid(p1,
                       p2,
                       nrow = 2,
                       align = "hv",
                       axis = "tblr",
                       rel_heights=c(2.5,1))
  }
}
