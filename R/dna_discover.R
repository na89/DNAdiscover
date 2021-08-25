#' Run predictor for presence of the platform bias
#'
#' @import pROC
#' @import ROCR
#' @param data data.frame with each column representing specific variant
#' annotation type. Annotations that could be included:
#' group, lcr,decoy, segdup,nonpar, variant_type, allele_type, was_mixed,
#' has_star,qd,info_MQRankSum, info_SOR, info_InbreedingCoeff, info_FS,
#' info_QD, info_MQ, info_DP, rf_probability, was_split, score,
#' qual, BaseQRankSum, ClippingRankSum,
#' FS, InbreedingCoeff, MQ, MQRankSum, QD, ReadPosRankSum.
#' For specific annotation descriptions please look up gnomad.broadinstitute.org
#' @param P a threshold for significance of discordance that was observed in
#' gnomad (exomes vs genomes in NFE population). Default value is 1e-5.
#' Smaller values increase specificity, greater values increase sensitivity.
#' @param auc boolean. If the ROC AUC estimate on training data should be printed
#' @return boolean vector of predictions
#'
#'
#' @export

dna_discover <- function(data, P = 1e-5, auc = TRUE){
 stopifnot(is.data.frame(data))
  load("R/sysdata.rda")
  data <- data[, colnames(data) %in% colnames(train_data)]

  train_data$group <- rep('bad', nrow(train_data))
  train_data$group[which(train_data$P > P)] <- 'good'
  good_features <- c('group', 'lcr', 'decoy',
                     'segdup', 'nonpar', 'variant_type',
                     'allele_type', 'was_mixed', 'has_star',
                     'qd', 'info_MQRankSum', 'info_SOR',
                     'info_InbreedingCoeff', 'info_FS', 'info_QD',
                     'info_MQ', 'info_DP', 'rf_probability',
                     'was_split', 'score', 'qual',
                     'BaseQRankSum', 'ClippingRankSum', 'FS',
                     'InbreedingCoeff', 'MQ', 'MQRankSum', 'QD',
                     'ReadPosRankSum')
  train_data <- train_data[, colnames(train_data) %in% good_features]
  train_data <- train_data[, colnames(train_data) %in% c('group', colnames(data))]
  stopifnot(all(colnames(train_data) %in% c('group', colnames(data))))
  rf <- randomForest::randomForest(as.factor(group) ~ ., data = train_data )
  if(auc){
    rf_p_train <- stats::predict(rf, type="prob")[,2]
    rf_pr_train <- ROCR::prediction(rf_p_train, train_data$group)
    r_auc_train <- ROCR::performance(rf_pr_train, measure = "auc")@y.values[[1]]
  }
  res <- stats::predict(rf, data, type = 'class')
  message(paste0('Test data AUC=', r_auc_train))
  res
}
