#' Differential methylated CpGs detection
#'
#' @param dat GenomicRatioSet object
#'
#' @return DiffMeth object
#' @export
#'
#' @importFrom minfi dropLociWithSnps getM pData getAnnotation getBeta
#' @importFrom siggenes qvalue.cal pi0.est
#' @importFrom dplyr as_tibble mutate select inner_join filter
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate
#' @importFrom limma lmFit contrasts.fit makeContrasts eBayes
#' @importFrom magrittr %>%
#' @importFrom stats setNames model.matrix na.omit
#'
#' @import data.table
#' @import methods
#' @import BiocGenerics
#'
diffMeth <- function(dat) {

  dat <- dropLociWithSnps(dat)
  M <- getM(dat)

  ## Generate contrasts
  sampleNames <- pData(dat)$Sample_Name
  sampleGroups <- pData(dat)$Sample_Group
  if (is.numeric(sampleGroups)){
    sampleGroups <- make.names(sampleGroups)
  }

  sampleGroups <- factor(as.character(sampleGroups))
  sg <- levels(sampleGroups)
  contrasts <- c()
  it_val <- 1
  for (i in 1:(length(sg) - 1)) {
    for (z in 1:(length(sg) - it_val)) {
      contrasts <- c(contrasts, paste(sg[i], sg[i+z], sep = '-'))
    }
    it_val <- it_val + 1
  }
  rm(it_val)

  ## Modeling
  design <- model.matrix(~sampleGroups)
  colnames(design) <- sg
  fit <- lmFit(M, design)  ## <- CHECK TO REDUCE EXECUTION TIME
  fit <- contrasts.fit(fit, makeContrasts(contrasts = contrasts, levels = sg))
  fit <- eBayes(fit)
  fit$p.value <- na.omit(fit$p.value)

  ## Ectracting topTables
  tt_results <- as_tibble()
  for (contrast in contrasts) {
    tt <- data.frame(pvalue = fit$p.value[,contrast]) %>%
      rownames_to_column(var = 'cpg') %>%
      as_tibble() %>%
      mutate(qval = qvalue.cal(pvalue, pi0.est(pvalue)$p0),
             comparation = contrast)
    tt_results <- rbind(tt_results, tt)
  }

  ## Calculate means
  B <- getBeta(dat) %>%
    as.data.frame()  %>%
    rownames_to_column('cpg') %>%
    as.data.table() %>%
    setNames(c('cpg', sampleNames)) %>%
    melt(id.vars = c('cpg'), measure.vars = sampleNames) %>%
    separate(variable, c('group','rep'))

  means <- B[, mean(value), by = c('cpg', 'group')] %>%
   dcast(cpg ~ group, value.var = 'V1')

  ## Calculate differentials
  contrasts <- list()
  it_val <- 1
  for (i in 1:(length(sg) - 1)) {
    for (z in 1:(length(sg) - it_val)) {
      contrasts[[paste(sg[i], sg[i+z], sep = '-')]] <- c(sg[i], sg[i+z])
    }
    it_val <- it_val + 1
  }
  rm(it_val)

  diff_results <- as_tibble()
  for (i in names(contrasts)) {
    diff <- select(means, c('cpg', contrasts[[i]]))
    diff$diff <-  diff[, .SD[,2] - .SD[,3]]
    diff <- as_tibble(diff) %>%
      select(c('cpg', 'diff')) %>%
      mutate(comparation = i)
    diff_results <- rbind(diff_results,diff)
  }

  ## Get Betas
  beta <- getBeta(dat)

  ## Get Annotation
  ann <- getAnnotation(dat) %>%
    as_tibble()

  ## stats_table
  stats_table <- inner_join(tt_results, diff_results, by = c('cpg', 'comparation')) %>%
    filter(cpg %in% rownames(beta))

  stats_table <- inner_join(stats_table, means, by = 'cpg')

  ## Returining output
  output <- list(stats = stats_table,
                 beta = beta,
                 ann = ann,
                 pData = pData(dat))
  return(output)
}
