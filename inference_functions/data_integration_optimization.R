# loads the pwm occurrence (prior TFBM matrix)
load("rdata/pwm_occurrences_N_response_varala.rdata")

# fills with 1/2 as a nuetral prior when TFBM is missing
pwm_imputed <- pwm_occurrence
pwm_imputed[is.na(pwm_imputed)] <- 0.5


# parallel sapply to parallelise the computing of optimal alphas
mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}


#' Draw gene effective data integration (EDI)
#' 
#' Plots the rank of the importance of TFs that have a specific prior value, here one, (supported by a TFBM)
#' as alpha is increased
#'
#' @param gene Gene to plot
#' @param mats List of matrices giving the TF-target genes importances values of an inference model
#' @param prior TFs thtat have this prior value will be shown
#' @param type Weather to plot the average importance directly ("imp") or the average rank of the importance
#' @param return Weather to return the average effective data integration of a gene
#' as a function of alpha.
#'
#' @return A ggplot object, or a data frame of effective data integration as a function of alpha
#' @export
#'
draw_gene_effective_integration <- function(gene, mats, prior=1, type = "rank", return = F){
  tfs_with_motif <- names(which(pwm_imputed[gene,]== prior))
  if(type == "rank")
    data <- data.frame(lapply(mats, function(mat){mean(rank(mat[,gene])[tfs_with_motif])}))
  if(type == "imp")
    data <- data.frame(lapply(mats, function(mat){mean(mat[,gene][tfs_with_motif])}))
  data <- data %>%
    gather(key = "setting", value = "summed_importance") %>%
    separate(setting, into = c("method", "alpha", "dataset", "rep"), sep = "_")
  
  
  if(return) return(data %>%
                      filter(dataset == "trueData") %>%
                      group_by(alpha, dataset) %>%
                      summarise(mean_imp = mean(summed_importance),
                                sd_imp = sd(summed_importance)))
  
  plot <- data %>%
    group_by(alpha, dataset) %>%
    mutate(mean_imp = mean(summed_importance),
           sd_imp = sd(summed_importance),
           alpha = as.numeric(alpha))%>%
    ggplot(aes(x=alpha, y = summed_importance, color = dataset, fill = dataset)) +
    geom_point(alpha = 0.1) + geom_line(aes(y=mean_imp))+
    geom_ribbon(aes(ymin = mean_imp - sd_imp , 
                    ymax = mean_imp + sd_imp  ), 
                alpha = .4)  +theme_pubr(legend = "none")+
    xlab(expression(alpha))+ylab("EDI")+ ggtitle(gene)+
    scale_color_manual(values = setNames(c("grey", "#70AD47"), c("shuffled", "trueData")))+
    scale_fill_manual(values = setNames(c("grey", "#70AD47"), c("shuffled", "trueData")))
  plot + xlab(expression(alpha))
}


#' Draw the MSE of a target gene
#' 
#' Plots the MSE as alpha is increased for a given gene.
#'
#' @param gene Gene to plot.
#' @param lmses Matrix of MSE values.
#' @param title Title of the plot.
#'
#' @return a ggplot object
#' @export
draw_gene_mse <- function(gene, lmses, title = NULL){
  data <-
    lmses[gene, ] %>%
    gather() %>%
    separate(key,
             into = c("model", "alpha", "dataset",  "rep"),
             sep = "_") 
  
  data %>%
    group_by(alpha, dataset) %>%
    mutate(mean_mse = mean(value, na.rm = T),
           sd_mse = sd(value, na.rm = T)) %>%
    ggplot(aes(
      x = as.numeric(alpha),
      y = value,
      color = dataset,
      fill = dataset
    )) +ylab("MSE")+ xlab(expression(alpha))+
    geom_ribbon(aes(ymin = mean_mse - sd_mse , 
                    ymax = mean_mse + sd_mse  ), 
                alpha = .4)  +theme_pubr(legend = "top")+ggtitle(gene)+
    geom_point(alpha = 0.1) + geom_line(aes(y=mean_mse))+xlab(expression(alpha))+ 
    scale_color_manual(values = setNames(c("grey", "#70AD47"), c("shuffled", "trueData")))+
    scale_fill_manual(values = setNames(c("grey", "#70AD47"), c("shuffled", "trueData")))
}


#' Finds the optimal value of alpha for a given gene
#' 
#' @param gene Gene of interest
#' @param mats List of matrices giving the TF-target genes importance values of an inference model
#' @param lmses Matrix of MSE values.
#' @param type Weather EDI is the average importance directly ("imp") or the average rank of the importance
#' of TFBM supported regulators
#' @param dev Deviation measure used to normalise the difference of means between true and null datasets.
#' "mean" means the mean between true and shuffled, else "true" or "shuff".
#' @param return_alpha Weather to return the value of optimal alpha.
#' @param metric Type of metric : "dev" for the maximal deviation from the shuffled baseline, "min" for
#' minimal MSE regardless of permutations.
#'
#' @return A ggplot object, or a value of alpha.
#' @export
get_opt_alpha_per_gene <- function(gene, mats, lmses, type = "rank", dev = "true",
                                   return_alpha = F, metric = 'div', pval.adjust = "fdr"){
  tfs_with_motif <- names(which(pwm_occurrence[gene,]==1))
  
 
  # gene with no TFBS has an optimal alpha of 0
  if(return_alpha & length(tfs_with_motif)==0) return(0)
  # else
  if(length(tfs_with_motif)>0){
    
    if(type == "rank")
     data <- data.frame(lapply(mats, function(mat){mean(rank(mat[,gene])[tfs_with_motif])}))
    if(type == "imp")
      data <- data.frame(lapply(mats, function(mat){mean(mat[,gene][tfs_with_motif])}))
    
    # get effective data integration
    inte <- data %>%
      gather(key = "setting", value = "summed_importance") %>%
      separate(setting, into = c("method", "alpha", "dataset", "rep"), sep = "_") %>%
      group_by(alpha, dataset) %>%
      mutate(mean_imp = mean(summed_importance, na.rm=T),
             sd_imp = sd(summed_importance, na.rm=T),
             alpha = as.numeric(alpha))
    
    # joins effective data integration with MSE data
    curves <- lmses[gene, ] %>%
      gather() %>%
      separate(key,
               into = c("model", "alpha", "dataset",  "rep"),
               sep = "_") %>%
      group_by(alpha, dataset) %>%
      mutate(mean_mse = mean(value, na.rm = T),
             sd_mse = sd(value, na.rm = T)) %>%
      mutate(alpha = as.numeric(alpha))%>%
      full_join(inte, by = c("alpha", "dataset", "rep")) %>%
      group_by(alpha, mean_imp, dataset) %>%
      summarise(mean_mse = mean(value, na.rm = T),
                sd_mse = sd(value, na.rm = T)) 
    
    # Approximates shuffled mean MSE and sd MSE for all true values of data integration
    # ( as the two curves are not aligned in terms of effectiva data integration)
    curves <- curves %>%
      group_by(dataset) %>%
      arrange(dataset, mean_imp)%>%
      mutate(imps=curves[curves$dataset=="trueData", ]$mean_imp) %>%
      mutate(approx_mse = approx(mean_imp,mean_mse,curves[curves$dataset=="trueData", ]$mean_imp, rule=2)$y,
             approx_sd = approx(mean_imp,sd_mse,curves[curves$dataset=="trueData", ]$mean_imp, rule=2)$y) 
    
    
    true <- curves[curves$dataset=="trueData",]
    shuff <- curves[curves$dataset!="trueData",]
    
    # true$div <- (shuff$approx_mse - shuff$approx_sd) - (true$mean_mse)
    
    # to get mse difference divided by permutations sd :
    if(dev == "shuff")
      true$div <- ifelse((shuff$approx_mse - true$mean_mse)/shuff$approx_sd>1, 
                         (shuff$approx_mse - true$mean_mse)/shuff$approx_sd, 
                         0)
    
    # to get mse difference divided by true sd :
    if(dev == "true")
      true$div <- ifelse((shuff$approx_mse - true$mean_mse)/true$sd_mse>1, 
                         (shuff$approx_mse - true$mean_mse)/true$sd_mse, 
                         0)
    
    # to get mse difference divided by the mean of true and shuffled sd :
    if(dev=="mean")
      true$div <- ifelse(2*(shuff$approx_mse - true$mean_mse)/(true$sd_mse+shuff$approx_sd)>1, 
                         2*(shuff$approx_mse - true$mean_mse)/(true$sd_mse+shuff$approx_sd), 
                         0)
    
    if(dev=="student"){
      # number of observations in the compared samples
      N = length(mats)/nrow(true)/2
      # quantile to be exceeded to be significant at the 5% threshold
      # critical_value = qt(p=.05/nrow(true), df=N-2, lower.tail=F)
      # student statistic is maximized
      pvals <- p.adjust(pt((shuff$approx_mse - true$mean_mse)/sqrt((true$sd_mse^2+shuff$approx_sd^2)/N), 
                  df = N-2, lower.tail = F), method = pval.adjust)
      true$div <- ifelse(pvals < 0.05, 
                         (shuff$approx_mse - true$mean_mse)/sqrt((true$sd_mse^2+shuff$approx_sd^2)/N), 
                         0)
    }
      
    
    # gets the value of alpha optimizing the chosen criterion
    if(metric == "div"){
      if(max(true$div)>0) alpha_opt <- true[true$div == max(true$div),]$alpha
      else alpha_opt <- 0
    }
    
    # If the criterion is instead the minimal MSE (regardless of shuffled data)
    if(metric == "min"){
      alpha_opt <- true[true$mean_mse == min(true$mean_mse),]$alpha
    }
    
    # value of effective data integration corresponding to the optimal alpha value
    eff_opt = true[true$alpha==alpha_opt,]$mean_imp
    
    # if asked, returns the optimal value of alpha for this gene
    if(return_alpha) return(alpha_opt)
    
    annotation = c(expression(alpha[t][",opt"]))
    
    # else, draws the figure :
    padding = ifelse(type=="rank", ifelse(alpha_opt<0.2, 12, -12), 0)
    curves%>%
      ggplot(aes(y=mean_mse, x = mean_imp, color = dataset, fill = dataset))+
      # geom_ribbon(aes(ymin =mean_mse-sd_mse , 
      #                 ymax = mean_mse + sd_mse),  alpha = .4)+ 
      geom_ribbon(aes(x=imps,ymin =approx_mse-approx_sd ,
                      ymax = approx_mse + approx_sd), alpha = .25)+
      geom_point(alpha = 0.1, size = 0.5) +
      geom_line(aes(x=imps, y = approx_mse), size=1) + 
      # geom_line(aes(x=mean_imp, y = mean_mse), size=1) + 
      # geom_line(aes(x=imps, y = approx_mse), col="black")+
      theme_pubr(legend = "none") +
      ylab("MSE") + xlab("EDI") + ggtitle(gene)+
      scale_color_manual(values = setNames(c("grey", "#70AD47"), c("shuffled", "trueData")))+
      scale_fill_manual(values = setNames(c("grey", "#70AD47"), c("shuffled", "trueData")))+
      geom_vline(xintercept = eff_opt, size = 2, col="#4670CD") +
      annotate("text", x=eff_opt+padding, y=max(shuff$mean_mse), 
               label=annotation, parse = T, 
               angle=0, col = "#4670CD", size =3.5 )
  }
  else ggplot()
}

