ungroup_tfs <- function(TFs){
  TFs_spread <- TFs[!grepl("mean_", TFs)]
  groups <- setdiff(TFs, TFs_spread)
  for (group in groups) {
    TFs_spread <- c(TFs_spread,
                    strsplit(stringr::str_split_fixed(group, "_", 2)[, 2],
                             '-')[[1]])
  }
  return(TFs_spread)
}