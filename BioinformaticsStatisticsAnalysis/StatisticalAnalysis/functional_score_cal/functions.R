kernelDensityThreshold = function(score, label, weights, posteriors) {
  # use precisions to choose the threshold for both sides
  
  # input
  # score: n*1 vector, the functional score
  # label: n*1 factor vector the label of the variant effects, 
  #   assume only two classes, the levels should be specified in the order
  #   of increasing scores
  # weights: 2*1 vector weights for each class
  # posteriors: a m*1 vector for different posterior cutoffs
  
  # output
  # cutoffs: m*2 matrix corresponding to specified posteriors for both classes
  
  scoreRange = range(score)
  
  levels = levels(label)
  hist(score[label == levels[2]], prob = T)
  density1 = density(score[label == "Benign"], 
                     from = scoreRange[1], to = scoreRange[2], n = 1000)
  # lines(density1, col = "green")
  
  hist(score[label == levels[1]], prob = T)
  density2 = density(score[label == levels[1]], 
                     from = scoreRange[1], to = scoreRange[2], n = 1000)
  # lines(density2, col = "red")
  
  # calculate the likelihood ratio
  L1 = density1$y + 1e-10
  L2 = density2$y + 1e-10
  
  prior = sum(label == levels[1]) * weights[1] / 
   (sum(label == levels[1]) * weights[1] + sum(label == levels[2]) * weights[2])
  posteriorProb = prior * L2 / (prior * L2 + (1 - prior) * L1)
  LR = L2 / L1
  plot(density1$x, posteriorProb)
  
  # find cutoffs for posteriors
  cutoffs = matrix(NA, length(posteriors), 2)
  rownames(cutoffs) = posteriors
  colnames(cutoffs) = levels(label)
  uniqueValues = density1$x
  for (i in 1:length(posteriors)) {
    # for class1, choose the rightmost
    cutoffs[i, 1] = uniqueValues[max(which(posteriorProb >= posteriors[i]))]
    # for class2, choose the leftmost
    cutoffs[i, 2] = uniqueValues[min(which(1 - posteriorProb >= posteriors[i]))]
  }
  
  return(list(cutoffs = cutoffs, LR = LR, posteriorProb = posteriorProb))
}

precisionBasedThreshold = function(score, label, weights, precisions) {
  # use precisions to choose the threshold for both sides
  
  # input
  # score: n*1 vector, the functional score
  # label: n*1 factor vector the label of the variant effects, 
  #   assume only two classes, the levels should be specified in the order
  #   of increasing scores
  # weights: 2*1 vector weights for each class
  # precisions: a m*1 vector for different precision cutoffs
  
  # output
  # cutoffs: m*2 matrix corresponding to specified precisions for both classes  
  
  # sort scores
  indexSort = order(score, decreasing = F)
  score = score[indexSort]
  label = label[indexSort]
  
  # assign the precision values for each score for both classes
  uniqueValues = unique(score)
  k = length(uniqueValues)
  precision1 = numeric(length(uniqueValues)) 
  precision2 = numeric(length(uniqueValues))
  score1 = score[label == levels(label)[1]]
  score2 = score[label == levels(label)[2]]
  
  for (i in 1:k) { # not efficient, but easy to write
    # for class 1
    TP1 = sum(score1 <= uniqueValues[i]) * weights[1]
    FP1 = sum(score2 <= uniqueValues[i]) * weights[2]
    
    # for class 2
    TP2 = sum(score2 >= uniqueValues[i]) * weights[2]
    FP2 = sum(score1 >= uniqueValues[i]) * weights[1]
    
    precision1[i] = TP1 / (TP1 + FP1)
    precision2[i] = TP2 / (TP2 + FP2)
  }
  
  # find cutoffs for precisions
  cutoffs = matrix(NA, length(precisions), 2)
  rownames(cutoffs) = precisions
  colnames(cutoffs) = levels(label)
  for (i in 1:length(precisions)) {
    # for class1, choose the rightmost
    cutoffs[i, 1] = uniqueValues[max(which(precision1 >= precisions[i]))]
    # for class2, choose the leftmost
    cutoffs[i, 2] = uniqueValues[min(which(precision2 >= precisions[i]))]
  }
  
  return(cutoffs)
}

randomEffectsBasedModels = function(effect, standardError, classLabel) {
  # generate two Gaussian models, assuming random effects for each class
  
  # input 
  # effect: the estimated effects
  # standardError: the estiamted standard error of the estimated effects
  # classLabel: a factor of two levels
  
  library(meta)
  factorLevels = levels(classLabel)
  mu = numeric(2)
  sd = numeric(2)
  for (i in 1:2) {
    class = factorLevels[i]
    index = (classLabel == class)
    effectEach = effect[index]
    standardErrorEach = standardError[index]
    result = metagen(effectEach, standardErrorEach, 
                     control=list(stepadj=0.5, maxiter=10000))
    mu[i] = result$TE.random
    sd[i] = result$tau
  }
  
  names(mu) = factorLevels
  names(sd) = factorLevels
  
  return(list(mu = mu, sd = sd))
}

mixtureModelEM = function(effect, classLabel) {
  
  library(mixtools)
  
  factorLevels = levels(classLabel)
  mu = numeric(2)
  sd = numeric(2)
  
  for (i in 1:2) {
    class = factorLevels[i]
    index = (classLabel == class)
    effectEach = effect[index]
    mu[i] = mean(effectEach)
    sd[i] = sd(effectEach)
  }
  
  x = normalmixEM(effect, 
                  lambda = 0.5, k = 2,
                  # fix the man and sd
                  mean.constr = mu,
                  sd.constr = sd,
                  mu = mu, 
                  sigma = sd)
  
  return(x)
}

posteriorProbabilityFromGaussianModels = function(effect, standardError,
                                                  meanAndSD, prior) {
  # calculate the posterior probabilities
  
  # effect: the estimated effect
  # standardError: the estimated standard error of the estimated effect
  # meanAndSD: a list with mu and sd for each component
  # prior: the prior vector for each variant
  
  mu = meanAndSD$mu
  sd = meanAndSD$sd
  # calculate the likelihood 
  logLikelihood1 = dnorm(effect,  mu[1], sqrt(sd[1]^2 + standardError^2), 
                         log = T)
  logLikelihood2 = dnorm(effect,  mu[2], sqrt(sd[2]^2 + standardError^2), 
                         log = T)
  
  logLikelihoodRatio = logLikelihood1 - logLikelihood2
  likelihoodRatio = exp(logLikelihoodRatio)
  
  posteriorP = (likelihoodRatio * prior) / ((likelihoodRatio - 1) * prior + 1)
  
  return(list(posteriorP = posteriorP, likelihoodRatio = likelihoodRatio))
}


ROCOut = function(ROCResult, truth) {
  levelName = levels(truth)                 
  confusion = c(sum(ROCResult[[2]] == levelName[1] & truth == levelName[1]),
                sum(ROCResult[[2]] == levelName[2] & truth == levelName[1]),
                sum(ROCResult[[2]] == levelName[1] & truth == levelName[2]),
                sum(ROCResult[[2]] == levelName[2] & truth == levelName[2]))
  confusion = paste(confusion, collapse = " ")
  out = data.frame(ROCResult[[1]]$auc, ROCResult[[3]], ROCResult[[4]], 
                   ROCResult[[5]], confusion)
  colnames(out) = c("AUC", "cutoff", "Sensitivity", "Specificity", 
                    "ConfusionMatrix")
  return(out)
}

ROCBasedPrediction = function(score, truth) {
  # calculate AUC and the predicted labels based on the optimal cut
  # assume a binary classification
  
  # input
  # score: the score
  # truth: the class label, a factor with lower level corresponds to lower score
  
  roc_score = roc(predictor = score, response = truth,
                  levels = rev(levels(truth)))
  out = data.frame(functional.score = roc_score$"threshold", 
                   sens =  roc_score$"sensitivities", 
                   specs = roc_score$"specificities")
  
  out$avg = out$sens+out$specs
  index = which(out$avg == max(out$avg))
  if (length(index) > 1) {
    out$F1 = 1 / (1 / out$sens + 1 / out$specs)
    F1Scores = out$F1[index]
    index = index[which.max(F1Scores)] # only keep the first with the max F1
    
    if (length(index) > 1) {
      sensitivity = out$sens[index]
      index = index[which.max(sensitivity)]
    }
  }
  
  cutoff0 = out$functional.score[index]
  sen0 = out$sens[index]
  spec0 = out$specs[index]
  
  labels = levels(truth)
  predicted = ifelse(score < cutoff0, labels[1], labels[2])
  predicted = factor(predicted, levels = levels(truth))
  confusion = table(truth, predicted)
  out = list(roc_score, predicted, cutoff0, sen0, spec0, confusion)
  names(out) = c("roc", "predicted", "cutoff", "sensitivity", "specificity",
                 "confusion_matrix")
  return(out)
}

medianIRQ = function(data, referenceType) {
  # normalize the data using the minus median and divided by interquartile
  
  # input
  # data: a data frame
  # referenceType: a logical index selecting the rows that are used as the 
  #  reference points
  
  median = apply(X = data[referenceType, , drop = F], MARGIN = 2, FUN = median)
  quartiles = apply(X = data[referenceType, , drop = F], MARGIN = 2, 
                    FUN = quantile, probs = c(0.25, 0.75))
  IRQ = quartiles[2, ] - quartiles[1, ]
  normalized = data
  for (i in 1:dim(data)[2]) {
    normalized[, i] = (data[, i] - median[i]) / IRQ[i]
  }
  
  return(normalized)
}

median2Normalize = function(data, referenceType1, referenceType2) {
  # normalize the data using the median of two reference types
  
  # input
  # data: a data frame
  # referenceType1, referenceType2: 
  #    a logical index selecting the rows that are used as the 
  #    reference points
  
  D14libColumns = grepl("D14lib", colnames(data))
  median1 = apply(X = data[referenceType1, D14libColumns, drop = F], MARGIN = 2, 
                  FUN = median)
  median2 = apply(X = data[referenceType2, D14libColumns, drop = F], MARGIN = 2, 
                  FUN = median)
  medianMean1 = mean(median1)
  medianMean2 = mean(median2)
  
  nReplicate = length(median1)
  a = numeric(nReplicate)
  b = numeric(nReplicate)
  
  for (i in 1:nReplicate) {
    a[i] = (medianMean2 - medianMean1) / (median2[i] - median1[i])
    b[i] = medianMean2 - a[i] * median2[i]
  }
  
  normalized = data
  
  replicateID = sub("_.*$", "", colnames(data))
  replicateNumber = as.integer(sub("R", "", replicateID))
  
  for (i in 1:dim(data)[2]) {
    normalized[, i] = a[replicateNumber[i]] * data[, i] + b[replicateNumber[i]]
  }
  
  return(normalized)
}

plotAcrossReplicates = function(data, nReplicate, featureID, levelsOfEvent) {
  dataForPlot = NULL
  for (i in 1:nReplicate) {
    names = paste0("R", i, "_", featureID)
    dataPerReplicate = data[, c(names, "EventType")]
    colnames(dataPerReplicate)[1] = "log2Ratio"
    replicate = paste0("R", i)
    dataPerReplicate = cbind(dataPerReplicate, replicate)
    dataForPlot = rbind(dataForPlot, dataPerReplicate)
  }
  g1 = ggplot(dataForPlot, aes(replicate, log2Ratio, 
                    fill = factor(EventType, levels = levelsOfEvent))) + 
    geom_violin()
  return(g1)
}

functionalClass = function(score, pathogenicT, benignT) {
  FSClass = character(length(score))
  FSClass[score <= pathogenicT] = "abnormal"
  FSClass[score >= benignT] = "normal"
  FSClass[score < benignT & score > pathogenicT] = 
    "I/U"
  
  FSClass = factor(FSClass, levels = c("normal", "I/U", "abnormal"))
  return(FSClass)
}

BayesDelNoAFClass = function(score) {
  FSClass = character(length(score))
  FSClass[] = "NA"

  FSClass[score >= 0.5] = "P_Strong"
  FSClass[score < 0.5 & score >= 0.27] = 
    "P_Moderate"
  FSClass[score < 0.27 & score >= 0.13] = 
    "P_Supporting"
  FSClass[score <= -0.18 & score > -0.36] = 
    "B_Supporting"
  FSClass[score <= -0.36] = 
    "B_Moderate"

  FSClass = factor(FSClass, levels = 
    c("P_Strong", "P_Moderate", "P_Supporting", "B_Supporting", "B_Moderate"))
  return(FSClass)
}
