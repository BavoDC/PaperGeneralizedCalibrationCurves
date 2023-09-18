#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Calibration prediction models #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### 0. Settings ####
library(magrittr)
library(MASS)
library(mgcv)
library(gbm)


Op = par(no.readonly = T)

###### 0.1 Custom functions #####
.GCVgbm <- function(f, Argz, nfolds, Df, OptFun, seedGCV = 1) {
  set.seed(seedGCV)
  N       = nrow(Df)
  FoldID  = sample(rep(seq(nfolds), length = N))

  ResFold = lapply(seq(nfolds), function(i) {
    which = FoldID == i
    DfTr  = Df[!which, , drop = F]
    Argz$data = DfTr
    set.seed(i)
    ModelFit  = do.call("f", Argz)
    DfTe  = Df[which, , drop = F]
    Pred  = predict(ModelFit, newdata = DfTe, type = "response", n.trees = Argz$n.trees)
    y     = DfTe[, all.vars(Argz$formula)[1]]
    Res =
      if(length(OptFun) == 1) {
        OptFun[[1]](y = y, mu = Pred)
      } else {
        do.call("cbind", lapply(OptFun, function(f) f(y = y, mu = Pred)))
      }
    return(Res)
  })
  ResFold = do.call("rbind", ResFold)
  PerfGCV = apply(ResFold, 2, mean)
  return(list(Performance = PerfGCV, Results = ResFold))
}

GCVgbm <- function(f, Argz, Grid, nfolds = 10, Df, OptFun, seed = 1,
                   Parallel = F, NrCores = parallel::detectCores() - 1,
                   ProgressBar = T, importPkgs = NULL, importObjects = NULL) {
  N       = nrow(Df)
  if (N / nfolds < 3)
    stop("< 3 observations per fold!")
  if(!all(colnames(Grid) %in% names(formals(f))))
    warning("Not all column names of Grid found in the arguments of f.")
  if(!is.list(OptFun))
    stop("OptFun has to be of type list.")
  if(!all(sapply(OptFun, is.function)))
    stop("Not all entries in OptFun are functions.")
  if(Parallel) {
    if(!"doSNOW" %in% names(sessionInfo()$otherPkgs))
      library(doSNOW)
    Cl = makeCluster(NrCores, type = "SOCK")
    registerDoSNOW(Cl)

    if(!missing(importPkgs)) {
      cat("\n\nLoading packages into workers.\n\n")
      clusterExport(Cl, "importPkgs", envir = environment())
      clusterEvalQ(Cl, lapply(importPkgs, library, character.only = T))
    }
    clusterExport(Cl, ".GCVgbm")
    if(!missing(importObjects))
      clusterExport(Cl,  importObjects, envir = environment())
    on.exit(stopCluster(Cl))
    if(ProgressBar) {
      pb       = txtProgressBar(max = nrow(Grid), style = 3)
      progress = function(n) setTxtProgressBar(pb, n)
      opts     = list(progress = progress)
    } else {
      opts = NULL
    }
    AllRes = foreach(i = seq_len(nrow(Grid)), .options.snow = opts,
                     .combine = "rbind") %dopar%
      {
        Param = Grid[i, , drop = T]
        for(j in names(Param))
          Argz[[j]] = Param[[j]]
        Res   = .GCVgbm(f, Argz, nfolds, Df, OptFun, seedGCV = seed)$Performance
        return(Res)
      }
  } else {
    if(ProgressBar)
      pb = txtProgressBar(min = 0, max = nrow(Grid), style = 3)
    AllRes = do.call("rbind",
                     lapply(seq_len(nrow(Grid)), function(i) {
                       Param = Grid[i, , drop = T]
                       for(j in names(Param))
                         Argz[[j]] = Param[[j]]
                       Res   = .GCVgbm(f, Argz, nfolds, Df, OptFun, seedGCV = seed)$Performance
                       if(ProgressBar)
                         setTxtProgressBar(pb, i)
                       return(Res)
                     }))
  }
  AllRes = cbind.data.frame(Grid, AllRes, row.names = NULL)
  return(AllRes)
}

poisson.dev <- function(y, mu, wt = rep(1, length(y))) {
  sum(poisson()$dev.resids(y = y, mu = mu, wt = wt))
}

genCalibration <- function(y, yHat, family, plot = FALSE, Smooth = FALSE, GLMCal = TRUE, lwdIdeal = 2, colIdeal = "gray", ltyIdeal = 1,
                           lwdSmooth = 1, colSmooth = "blue", ltySmooth = 1, argzSmooth = alist(degree = 2),
                           lwdGLMCal = 1, colGLMCal = "red", ltyGLMCal = 1,
                           AddStats = TRUE, Digits = 3, posStats = NULL, cexStats = 1, lwdLeg = 1.5, Legend = TRUE, legendPos = "bottomright",
                           xLim = NULL, yLim = NULL,
                           confLimitsSmooth = c("none", "bootstrap", "pointwise"), confLevel = 0.95,
                           Title = "Calibration plot",
                           xlab = "Predicted", ylab = "Observed",
                           EmpiricalDistribution = TRUE, length.seg = 1, ...) {
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  confLimitsSmooth = match.arg(confLimitsSmooth)
  a = 1 - confLevel

  Eta    = family$linkfun(yHat)
  ClInt = tryCatch(
    glm(
      y ~ offset(Eta),
      family = family,
      control = glm.control(maxit = 1e2)
    ),
    error = function(e)
      T,
    warning = function(w)
      T
  )
  if(is.logical(ClInt)) {
    # https://stackoverflow.com/questions/8212063/glm-starting-values-not-accepted-log-link
    ClInt = glm(I(y + .Machine$double.eps) ~ offset(Eta), family = family, control = glm.control(maxit = 1e2))
  }
  ClSl  =
    tryCatch(
      glm(y ~ Eta, family = family, control = glm.control(maxit = 1e2)),
      error = function(e)
        T,
      warning = function(w)
        T
    )
  if(is.logical(ClSl)) {
    lmFit = lm(y ~ Eta)
    ClSl  = glm(y ~ Eta, family = family, control = glm.control(maxit = 1e2), start = coef(lmFit))
  }
  ClSl2 = tryCatch(
    glm(y ~ Eta - 1, family = family),
    error = function(e)
      T,
    warning = function(w)
      T
  )
  if(is.logical(ClSl2)) {
    lmFit = lm(y ~ Eta - 1)
    ClSl2  = glm(y ~ Eta - 1, family = family, control = glm.control(maxit = 1e2), start = coef(lmFit))
  }
  CalibrStats = c("Calibration intercept" = unname(coef(ClInt)), "Calibration slope" = unname(coef(ClSl)[2]))


  if(plot) {
    y    = y[order(yHat)]
    Eta  = Eta[order(yHat)]
    yHat = sort(yHat)
    if(GLMCal) {
      glmFit = glm(y ~ Eta, family = family)
      rangeY = range(glmFit$fitted)
    }
    if(Smooth) {
      argzSmooth$formula = y ~ yHat
      SmFit <- Sm <- do.call("loess", argzSmooth)
      Sm     = data.frame(Sm$x, Sm$fitted)
      rangeY = if(GLMCal) c(min(rangeY, SmFit$fitted), max(rangeY, SmFit$fitted)) else range(SmFit$fitted)
    }
    xLim = if(is.null(xLim)) range(yHat) else xLim
    yLim = if(is.null(yLim)) c(min(c(xLim, rangeY)), max(c(xLim, rangeY))) else yLim
    yLim[1] =
      if(yLim[1] <= 0.5) {
        0 - 0.1 * diff(range(yLim))
      } else {
        yLim[1] * 0.9
      }
    plot(mean(xLim), mean(yLim), col = "white", pch = 1, xlab = xlab, ylab = ylab,
         xlim = xLim, ylim = yLim, main = Title, ...)
    clip(min(c(xLim, yHat)), max(c(xLim, yHat)), min(c(yLim, rangeY)), max(c(yLim, rangeY)))

    labLeg = "Ideal"
    colLeg = colIdeal
    ltyLeg = ltyIdeal
    lwdLeg = lwdIdeal

    if(Smooth) {
      lines(Sm, lty = ltySmooth, lwd = lwdSmooth, col = colSmooth)
      if(confLimitsSmooth != "none") {
        if(confLimitsSmooth == "bootstrap") {
          yHatGrid = seq(min(yHat), max(yHat), length = 200)
          resBoot  = replicate(2000, bootSamples(y, yHat, yHatGrid))
          clBoot   = apply(resBoot, 1, quantile, c(0.025, 0.975))
          dfCL     = data.frame(x = yHatGrid, ymin = clBoot[1, ], ymax = clBoot[2, ])
          rownames(dfCL) = NULL
        } else {
          cl.loess = predict(SmFit, type = "fitted", se = TRUE)
          dfCL     = data.frame(x = yHat, ymin = with(cl.loess, fit - qnorm(1 - a / 2) * se.fit),
                                ymax = with(cl.loess, fit + qnorm(1 - a / 2) * se.fit))
        }
        with(dfCL,
             polygon(
               x = c(x, rev(x)),
               y = c(ymax,
                     rev(ymin)),
               col = rgb(177, 177, 177, 177, maxColorValue = 255),
               border = NA
             )
        )
      }
      labLeg = c(labLeg, "Flexible calibration")
      colLeg = c(colLeg, colSmooth)
      ltyLeg = c(ltyLeg, ltySmooth)
      lwdLeg = c(lwdLeg, lwdSmooth)
    }
    if(GLMCal) {
      lines(yHat, fitted(glmFit), lty = ltyGLMCal, lwd = lwdGLMCal, col = colGLMCal)
      labLeg = c(labLeg, "GLM calibration")
      colLeg = c(colLeg, colGLMCal)
      ltyLeg = c(ltyLeg, ltyGLMCal)
      lwdLeg = c(lwdLeg, lwdGLMCal)
    }
    abline(0, 1, col = colIdeal, lwd = lwdIdeal, lty = ltyIdeal)
    do.call("clip", as.list(par()$usr))
    if(EmpiricalDistribution) {
      x     <- yHat
      bins  <- seq(min(x), max(x), length = 101)
      f0	  <- table(cut(x, bins))
      bins  <- (bins[-101])
      maxf  <- max(f0)
      f0	  <- (0.1 * f0) / maxf

      segments(bins, yLim[1], bins, yLim[1] + length.seg * f0)
      lines(c(min(bins) - 0.01, max(bins) + 0.01), c(yLim[1], yLim[1]))
    }
    if(AddStats) {
      StatsPlot = paste0('Calibration\n',
                         '...intercept: ',
                         sprintf(paste0("%.", Digits, "f"), CalibrStats[1]), '\n',
                         '...slope: ',
                         sprintf(paste0("%.", Digits, "f"), CalibrStats[2]), '\n')
      if(is.null(posStats))
        text(xLim[1], xLim[2] * 0.85, StatsPlot, pos = 4, cex = cexStats)
      else
        text(posStats[1], posStats[2], StatsPlot, pos = 4, cex = cexStats)
    }
    if(Legend)
      if(is.character(legendPos))
        legend(legendPos, legend = labLeg, col = colLeg, lty = ltyLeg, bty = "n", lwd = lwdLeg)
    else
      legend(legendPos[1], legendPos[2], legend = labLeg, col = colLeg, lty = ltyLeg, bty = "n", lwd = lwdLeg)
  }
  return(CalibrStats)
}

HistGGplot <- function(var, data, Digits = options()$digits,
                       freq = T, ylim, xlim, box = T, Bins = 30, BinAdj = F,
                       Summary = c("none", "topleft", "topright"), txt.x, txt.y, xticks = NULL, yticks = NULL,
                       main = "You forgot the text", xlab = "Number thingies",
                       ylab = if(freq) "Frequency" else "Probability",
                       axisCex = 12, axisTitleCex = 14, plotTitleCex = 14,
                       SummaryCex = 5, Col = "white", alpha = 1){
  Argz = as.list(match.call())[-1]
  Summary = match.arg(Summary)

  if(! "ggplot2" %in% names(sessionInfo()$otherPkgs))
    library(ggplot2)
  if(missing(data)) {
    x = var
  } else {
    x = eval(Argz$var, data)
  }

  if(!is.numeric(x))
    stop("Variable has to be of type numeric.")

  if(anyNA(x)) {
    warning("There are missing values and these will be removed.")
    x = x[!is.na(x)]
  }

  Avg = round(mean(x), Digits)
  SD  = round(sd(x), Digits)
  Med = round(median(x), Digits)
  Range = round(range(x), Digits)


  # ggplot
  if(!freq){
    df = as.data.frame(x)
    p1 = ggplot(df, aes(x = x)) + geom_histogram(aes(y = ..density..),
                                                 bins = if(BinAdj) length(unique(x)) else Bins,
                                                 fill = I(Col),
                                                 col = I("black"),
                                                 alpha = alpha)
    p1 = p1 + xlab(xlab) + ylab(ylab) + ggtitle(main)
  }else{
    p1 = qplot(x, geom = "histogram", alpha = I(alpha), bins = if(BinAdj) length(unique(x)) else Bins,
               fill = I(Col), col = I("black"),
               main = main, xlab = xlab, ylab = ylab)
  }

  if(!missing(ylim))
    p1 = p1 + scale_y_continuous(limits = ylim, breaks = if(!is.null(yticks)) yticks else NULL)
  if(!missing(xlim))
    p1 = p1 + scale_x_continuous(limits = xlim, breaks = if(!is.null(xticks)) xticks else NULL)

  p3 = p1 +theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  if(box)
    p3 = p3 + theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5))


  p4 =
    if(Summary != "none") {
      xTmp = ggplot_build(p3)$layout$panel_ranges[[1]]$x.range
      yTmp = ggplot_build(p3)$layout$panel_ranges[[1]]$y.range
      xCoord = if(Summary == "topleft") min(xTmp) else 0.8 * max(xTmp)
      yCoord = if(Summary == "topleft") max(yTmp) else 0.8 * max(yTmp)
      p3  + annotate("text", label = paste0("Mean (SD) = ", Avg," (", SD,")",
                                            "\n","Median = ", Med,
                                            "\n","Range = ", Range[1]," to ", Range[2]),
                     x = xCoord, y = yCoord, size = SummaryCex, hjust = 0)
    } else {
      p3
    }

  FinalPlot <- p4 + theme(axis.text = element_text(size = axisCex),
                          axis.title = element_text(size = axisTitleCex),
                          plot.title = element_text(size = plotTitleCex, face = "bold"))
  FinalPlot
}

ScaleRange <- function(x, xmin = -1, xmax = 1) {
  xRange = range(x)
  (x - xRange[1]) / diff(xRange) * (xmax - xmin) + xmin
}

#### 1. Poisson regression ####
KULbg <- "#116E8A"

##### 1.1 Examples paper ####

###### 1.1.1 Generate population data ####
set.seed(144)
p    = 5
N    = 1e6
n    = 5e3
nOOS = 1e3
S    = matrix(NA, 5, 5)
rho  = c(0.025, 0, 0, 0.05, 0.075, 0, 0, 0.025, 0, 0)
S[upper.tri(S)] = rho
S[lower.tri(S)] = t(S)[lower.tri(S)]
diag(S) = 1
Matrix::isSymmetric(S)


X  = mvrnorm(N, rep(0, p), Sigma = S, empirical = T)
cor(X)
X  = apply(X, 2, ScaleRange)
apply(X, 2, range)
B  = c(-2.3, 1.5, 2, -1, -2, -1.5)
mu = poisson()$linkinv(cbind(1, X) %*% B)
Y  = rpois(N, mu)
table(Y)

Df = data.frame(Y, X)
colnames(Df)[-1] %<>% tolower()
glm(Y ~ ., data = Df, family = poisson)
(yDistr = HistGGplot(Y, Df, xlab = expression(Y[i]), main = "",
                     Col = KULbg, alpha = 0.5, freq = T) + theme(legend.position = "none"))

###### 1.1.2 Illustration with GLMs ####
ModelFits = list()

####### 1.1.2.1 Near perfect #######
set.seed(2)
DfS    = Df[sample(1:nrow(Df), n, F), ]
glmFit = glm(Y ~ ., data = DfS, family = poisson)
Bhat   = coef(glmFit)

## Validation dataset ##
DfOOS = Df[sample(1:nrow(Df), nOOS, F), ]
yOOS  = DfOOS$Y
XOOS  = model.matrix(glmFit, data = DfOOS)
yHat  = as.vector(exp(XOOS %*% Bhat))
muOOS = as.vector(exp(XOOS %*% B))

ModelFits[["Theoretical"]] = list(y = yOOS, yHat = muOOS)
ModelFits[["NearPerfect"]] = list(y = yOOS, yHat = yHat)


####### 1.1.2.2 Overfit ######
gamFit = gam(Y ~ x1 + x3 + x1:x3 + s(x5), data = DfS, family = poisson)

## Validation dataset ##
yHat  = predict(gamFit, newdata = DfOOS, type = "response")

ModelFits[["Overfit"]] = list(y = yOOS, yHat = yHat)


####### 1.1.2.3 Underfit ######
glmFit = glm(Y ~ x2, data = DfS, family = poisson)
yHat   = predict(glmFit, newdata = DfOOS, type = "response")

ModelFits[["Underfit"]] = list(y = yOOS, yHat = yHat)


###### 1.1.2.4 All models #######
par(mfrow = c(2, 2), cex.axis = 1.25, cex.lab = 1.5, cex.main = 1.75, oma=c(2,2,0,0))
lapply(seq_along(ModelFits), function(i) {
  x = ModelFits[[i]]
  list2env(x, envir = environment())
  ArgzPlot = alist(
    y = y,
    yHat = yHat,
    family = poisson,
    plot = T,
    xlab = "",
    ylab = "",
    colGLMCal = "black", ltyGLMCal = 5, Smooth = TRUE, argzSmooth = alist(degree = 1, span = 0.75),
    colSmooth = "black", colIdeal = "black", ltyIdeal = 3,
    lwdLeg = 2, Legend = FALSE, Title = paste0("(", letters[i], ")"), legendPos = c(0.6, 0.25),
    confLimitsSmooth = "pointwise", las = 1,
    cexStat = 1.5, cex.axis = 1.5,
    xLim = if(i < 3) NULL else if(i == 3) c(0, 0.6) else c(0, 0.5)
  )
  do.call("genCalibration", ArgzPlot)
})
mtext("Predicted value", line=0, side=1, outer=TRUE, cex=1.85)
mtext("Empirical average", line=0, side=2, outer=TRUE, cex=1.85)



###### 1.1.3 GBM ######
gbmFits = list()

###### 1.1.3.1 Optimal fit #####
gridGBM =
  expand.grid(
    interaction.depth = 1:10,
    n.trees = seq(1e2, 5e3, by = 1e2)
  )
p    = ncol(DfS)
nDfS = nrow(DfS)
argzGBM =
  alist(
    formula = Y ~ .,                           # formula
    data = DfS,                                # data set
    var.monotone = rep(0, p - 1),              # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
    distribution = "poisson",                  # see the help for other choices
    n.trees = 3500,                            # number of trees
    shrinkage = 0.01,                          # shrinkage or learning rate
    interaction.depth = 1,                     # 1: additive model, 2: two-way interactions, etc.
    bag.fraction = 0.75,                       # subsampling fraction
    train.fraction = 1,                        # fraction of data for training
    n.minobsinnode =
      ceiling(nDfS * 0.01 * 0.75),             # minimum total weight needed in each node
    cv.folds = 0,                              # no k-fold cross-validation
    keep.data = TRUE,                          # keep a copy of the dataset with the object
    verbose = FALSE,                           # don't print out progress
    n.cores = 1
  )


resGCV   = GCVgbm(gbm, argzGBM, gridGBM, Df = DfS, OptFun = list(poisson.dev), Parallel = TRUE, NrCores = 10,
                  importObjects = c("DfS", "p", "nDfS"), seed = 1234)
optParam = resGCV[which(resGCV$AllRes == min(resGCV$AllRes)), 1:2, drop = TRUE]
for(j in names(optParam))
  argzGBM[[j]] = optParam[[j]]

set.seed(1234)
optFit   = do.call("gbm", argzGBM)

yHat = predict(optFit, DfOOS, n.trees = optFit$n.trees, type = 'response')
gbmFits[["NearPerfect"]] = list(y = yOOS, yHat = yHat)



###### 1.1.3.2 Overfit ######
set.seed(54321)
gbmFit <- gbm(
  formula = Y ~ .,
  data = DfS,
  distribution = 'poisson',
  n.trees = 5000, # T in Table 3
  interaction.depth = 5, # d in Table 3
  shrinkage = 0.01, # lambda in Table 1
  bag.fraction = 0.75, # delta in Table 1
  n.minobsinnode = 0.01 * 0.75 * nrow(DfS), # kappa * delta in Table 1
  verbose = FALSE
)

yHat = predict(gbmFit, DfOOS, n.trees = gbmFit$n.trees, type = 'response')
gbmFits[["Overfit"]] = list(y = yOOS, yHat = yHat)


####### 1.1.3.3 Underfit ######
set.seed(54321)
gbmFit <- gbm(
  formula = Y ~ .,
  data = DfS,
  distribution = 'poisson',
  n.trees = 200, # T in Table 3
  interaction.depth = 1, # d in Table 3
  shrinkage = 0.01, # lambda in Table 1
  bag.fraction = 0.75, # delta in Table 1
  n.minobsinnode = 0.01 * 0.75 * nrow(DfS), # kappa * delta in Table 1
  verbose = FALSE
)

yHat = predict(gbmFit, DfOOS, n.trees = gbmFit$n.trees, type = 'response')
gbmFits[["Underfit"]] = list(y = yOOS, yHat = yHat)


###### 1.1.3.4 Plot ######
m =
  rbind(
    c(1, 1, 2, 2),
    c(4, 3, 3, 5)
  )

par(cex.axis = 1.25, cex.lab = 1.5, cex.main = 1.75)
layout(m)
lapply(seq_along(gbmFits), function(i) {
  x = gbmFits[[i]]
  list2env(x, envir = environment())
  ArgzPlot = alist(
    y = y,
    yHat = yHat,
    family = poisson,
    plot = T,
    xlab = "",
    ylab = "",
    colGLMCal = "black", ltyGLMCal = 5, Smooth = TRUE, argzSmooth = alist(degree = 1, span = 0.75),
    colSmooth = "black", colIdeal = "black", ltyIdeal = 3,
    lwdLeg = 2, Legend = FALSE, Title = paste0("(", letters[i], ")"),
    posStats = if(i < 3) NULL else c(0.625, 0),
    cexStat = 1.5, cex.axis = 1.5,
    confLimitsSmooth = "pointwise", las = 1,
    xLim = if(i == 1) c(0, 1.2) else c(0, 1),
    yLim = if(i == 1) c(0, 1.2) else c(0, 1)
  )
  do.call("genCalibration", ArgzPlot, envir = environment())
})
mtext("Predicted value", line=0, side=1, outer=TRUE, cex=1.85)
mtext("Empirical average", line=0, side=2, outer=TRUE, cex=1.85)
layout(1)
