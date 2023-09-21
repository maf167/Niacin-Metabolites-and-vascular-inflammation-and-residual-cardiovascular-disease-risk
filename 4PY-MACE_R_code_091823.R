#
# Custom R code used for the article: "Terminal metabolites of niacin, a 
# component of cereal fortification, are associated with vascular inflammation 
# and residual cardiovascular disease risk"
#  
# Note custom code was used to generate the following display items:
#   * Figure 2a-b
#   * Figure 3a-h
#   * Extended Data Figure 6c
#   * Extended Data Figure 7a-b
#   * Extended Data Figure 8e     
#   * Supplemental  Figure 5a-b
#   * Supplemental  Table  2
#   * Supplemental  Table  3
#   * Supplemental  Table  4
# 
# Dependencies: dplyr, rms, survival, reshape2, topGO, org.Hs.eg.db, ggplot2
#
# Note d is a data frame containing the fields referenced below.
#
# Author: Marc Ferrell

library(dplyr)
library(ggplot2)
library(rms)
library(survival)
library(reshape2)

library(topGO)
library(org.Hs.eg.db)

drawKMQuartiles <- function(t, v, x, ymin){
  # t is Time
  # v is vital status (event = 1)
  # x is quartile
  
  require(survival)
  
  # Remove NAs
  comp <- complete.cases(data.frame(a = t, b = v, c = x))
  t <- t [comp]; v <- v[comp]; x <- x[comp]
  
  ## Make a survfit object for each quantile
  
  surv.fit <- survfit(Surv(t,v) ~ x)
  # plot(surv.fit)
  plot(surv.fit, col = c("black", "darkgreen", "blue", "red"), 
       ylim = c(ymin, 1))
  
}

my.coxph <- function(my.formula){
  
  t <- summary( robcov( coxph(my.formula) ) )
  
  out <- as.data.frame(signif(t$conf.int, 3))
  
  out$P <- signif(as.vector(t$coef[,"Pr(>|z|)", drop = TRUE]), 3)
  
  (out)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 2a: Kaplan-Meir estimates for MACE ranked by quartiles of a plasma 
#            analyte with m/z = 153.0656 Da. 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Transform Unknown # 221 concentrations to quartiles

u221   <- d$HILIC_posESI..221

q   <- quantile(u221, c(0.25, 0.50, 0.75), na.rm=TRUE)
df1 <- which(!is.na(u221) & u221 <= q[1])
df2 <- which(!is.na(u221) & u221 <= q[2] & u221 > q[1])
df3 <- which(!is.na(u221) & u221 <  q[3] & u221 > q[2])
df4 <- which(!is.na(u221) & u221 >= q[3])

a      <- rep(NA, length(u221))
a[df1] <- 1; a[df2] <- 2; a[df3] <- 3; a[df4] <- 4
d$Unknown.221.Quartile  <- a


## Kaplan-Meier plot

drawKMQuartiles(d$DTDMS3_YU, d$DMS3, d$Unknown.221.Quartile, 0.8)

## Log-rank test
diff <- survdiff(Surv(d$DTDMS3_YU,d$DMS3) ~ d$Unknown.221.Quartile)
signif(
  pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE),
  digits = 2 )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 2b: Cox proportional hazards for MACE ranked by quartiles of a plasma 
#            analyte with m/z = 153.0656 Da. 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Matrix of Covariants for Cox models

## Ethanol consumption (Yes/No to >= 1 drink per day)

d$EtohDrinksPerDayme[d$EtohDrinksPerDayme == "(Unknown)"] <- NA
d$EtohDrinksPerDayme[d$EtohDrinksPerDayme == "<1, occasiol or social"] <- "0"

d$EtohDrinksPerDayme <- as.numeric(d$EtohDrinksPerDayme)
d$EtohDrinksPerDayme <- ifelse(d$EtohDrinksPerDayme <= 1, 0, 1)

## Sex (Male is 1, Female is 0)

d$Gender <- ifelse(d$Gender == "M", 1, 0)

my.covar <- as.matrix(
  d[,c("AgeAtBloodDraw", "Gender", "BPSystolic",
       "LDL.Priority", "HDL.Priority", "TG.Priority",
       "DIABETICS","CRP16", "CurrentSmoker", "EtohDrinksPerDayme",
       "Education.Level")]
)

### Unadjusted

u221.unadjusted <- rbind(
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$Unknown.221.Quartile, `1` = 0, `2` = 1, .default = NaN)),
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$Unknown.221.Quartile, `1` = 0, `3` = 1, .default = NaN)),
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$Unknown.221.Quartile, `1` = 0, `4` = 1, .default = NaN))
)

rownames(u221.unadjusted) <- paste("Q", 2:4, sep = "")

print(u221.unadjusted) 


### Adjusted

u221.adjusted <- rbind(
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$Unknown.221.Quartile, `1` = 0, `2` = 1, .default = NaN) + my.covar)[1,],
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$Unknown.221.Quartile, `1` = 0, `3` = 1, .default = NaN) + my.covar)[1,],
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$Unknown.221.Quartile, `1` = 0, `4` = 1, .default = NaN) + my.covar)[1,]
)

rownames(u221.adjusted) <- paste("Q", 2:4, sep = "")

print(u221.adjusted)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 3a: Kaplan-Meir estimates for MACE ranked by quartiles of plasma 4PY 
#            in the US Validation Cohort
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



### Transform 4PY concentrations to quartiles

py4   <- d$X4PY

q   <- quantile(py4, c(0.25, 0.50, 0.75), na.rm=TRUE)
df1 <- which(!is.na(py4) & py4 <= q[1])
df2 <- which(!is.na(py4) & py4 <= q[2] & py4 > q[1])
df3 <- which(!is.na(py4) & py4 <  q[3] & py4 > q[2])
df4 <- which(!is.na(py4) & py4 >= q[3])

a      <- rep(NA, length(py4))
a[df1] <- 1; a[df2] <- 2; a[df3] <- 3; a[df4] <- 4
d$PY4.Quartile  <- a

## Kaplan-Meier plot

drawKMQuartiles(d$DTDMS3_YU, d$DMS3, d$PY4.Quartile, 0.75)

## Log-rank test
diff <- survdiff(Surv(d$DTDMS3_YU,d$DMS3) ~ d$PY4.Quartile)
signif(
  pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE),
  digits = 2 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 3b: Kaplan-Meir estimates for MACE ranked by quartiles of plasma 2PY 
#            in the US Validation Cohort 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Transform 2PY concentrations to quartiles

py2   <- d$X2PY

q   <- quantile(py2, c(0.25, 0.50, 0.75), na.rm=TRUE)
df1 <- which(!is.na(py2) & py2 <= q[1])
df2 <- which(!is.na(py2) & py2 <= q[2] & py2 > q[1])
df3 <- which(!is.na(py2) & py2 <  q[3] & py2 > q[2])
df4 <- which(!is.na(py2) & py2 >= q[3])

a      <- rep(NA, length(py2))
a[df1] <- 1; a[df2] <- 2; a[df3] <- 3; a[df4] <- 4
d$PY2.Quartile  <- a

## Kaplan-Meier plot

drawKMQuartiles(d$DTDMS3_YU, d$DMS3, d$PY2.Quartile, 0.75)

## Log-rank test
diff <- survdiff(Surv(d$DTDMS3_YU,d$DMS3) ~ d$PY2.Quartile)
signif(
  pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE),
  digits = 2 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 3c: Cox proportional hazards for MACE ranked by quartiles of 4PY 
#            in the US Validation Cohort
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Matrix of Covariants for Cox models

## Ethanol consumption (Yes/No to >= 1 drink per day)

d$EtohDrinksPerDayme[d$EtohDrinksPerDayme == "(Unknown)"] <- NA
d$EtohDrinksPerDayme[d$EtohDrinksPerDayme == "<1, occasiol or social"] <- "0"

d$EtohDrinksPerDayme <- as.numeric(d$EtohDrinksPerDayme)
d$EtohDrinksPerDayme <- ifelse(d$EtohDrinksPerDayme <= 1, 0, 1)

## Sex (Male is 1, Female is 0)

d$Gender <- ifelse(d$Gender == "M", 1, 0)

my.covar.gb <- as.matrix(
  d[,c("AgeAtBloodDraw", "Gender", "BPSystolic",
       "LDL.Priority", "HDL.Priority", "TG.Priority",
       "DIABETICS","CRP16", "CurrentSmoker", "EtohDrinksPerDayme",
       "Education.Level")]
)

### Unadjusted

py4.unadjusted <- rbind(
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY4.Quartile, `1` = 0, `2` = 1, .default = NaN)),
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY4.Quartile, `1` = 0, `3` = 1, .default = NaN)),
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY4.Quartile, `1` = 0, `4` = 1, .default = NaN))
)

rownames(py4.unadjusted) <- paste("Q", 2:4, sep = "")

print(py4.unadjusted) 


### Adjusted

py4.adjusted <- rbind(
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY4.Quartile, `1` = 0, `2` = 1, .default = NaN) + my.covar.gb)[1,],
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY4.Quartile, `1` = 0, `3` = 1, .default = NaN) + my.covar.gb)[1,],
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY4.Quartile, `1` = 0, `4` = 1, .default = NaN) + my.covar.gb)[1,]
)

rownames(py4.adjusted) <- paste("Q", 2:4, sep = "")

print(py4.adjusted)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 3d: Cox proportional hazards for MACE ranked by quartiles of 2PY  
#            in the US Validation Cohort
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Unadjusted

py2.unadjusted <- rbind(
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY2.Quartile, `1` = 0, `2` = 1, .default = NaN)),
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY2.Quartile, `1` = 0, `3` = 1, .default = NaN)),
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY2.Quartile, `1` = 0, `4` = 1, .default = NaN))
)

rownames(py2.unadjusted) <- paste("Q", 2:4, sep = "")

print(py2.unadjusted) 


### Adjusted

py2.adjusted <- rbind(
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY2.Quartile, `1` = 0, `2` = 1, .default = NaN) + my.covar.gb)[1,],
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY2.Quartile, `1` = 0, `3` = 1, .default = NaN) + my.covar.gb)[1,],
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$PY2.Quartile, `1` = 0, `4` = 1, .default = NaN) + my.covar.gb)[1,]
)

rownames(py2.adjusted) <- paste("Q", 2:4, sep = "")

print(py2.adjusted)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 3e: Kaplan-Meir estimates for MACE ranked by quartiles of plasma 4PY 
#            in the European Validation Cohort 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Transform 4PY concentrations to quartiles

py4   <- d$X4PY

q   <- quantile(py4, c(0.25, 0.50, 0.75), na.rm=TRUE)
df1 <- which(!is.na(py4) & py4 <= q[1])
df2 <- which(!is.na(py4) & py4 <= q[2] & py4 > q[1])
df3 <- which(!is.na(py4) & py4 <  q[3] & py4 > q[2])
df4 <- which(!is.na(py4) & py4 >= q[3])

a      <- rep(NA, length(py4))
a[df1] <- 1; a[df2] <- 2; a[df3] <- 3; a[df4] <- 4
d$PY4.Quartile  <- a

## Kaplan-Meier plot

drawKMQuartiles(d$D_MI_Stroke, d$D_MI_StrokeYN, d$PY4.Quartile, 0.55)

## Log-rank test
diff <- survdiff(Surv(d$D_MI_Stroke,d$D_MI_StrokeYN) ~ d$PY4.Quartile)
signif(
  pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE),
  digits = 2 )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 3f: Kaplan-Meir estimates for MACE ranked by quartiles of plasma 2PY 
#            in the European Validation Cohort 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### Transform 2PY concentrations to quartiles

py2   <- d$X2PY

q   <- quantile(py2, c(0.25, 0.50, 0.75), na.rm=TRUE)
df1 <- which(!is.na(py2) & py2 <= q[1])
df2 <- which(!is.na(py2) & py2 <= q[2] & py2 > q[1])
df3 <- which(!is.na(py2) & py2 <  q[3] & py2 > q[2])
df4 <- which(!is.na(py2) & py2 >= q[3])

a      <- rep(NA, length(py2))
a[df1] <- 1; a[df2] <- 2; a[df3] <- 3; a[df4] <- 4
d$PY2.Quartile  <- a

## Kaplan-Meier plot

drawKMQuartiles(d$D_MI_Stroke, d$D_MI_StrokeYN, d$PY2.Quartile, 0.55)

## Log-rank test
diff <- survdiff(Surv(d$D_MI_Stroke,d$D_MI_StrokeYN) ~ d$PY2.Quartile)
signif(
  pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE),
  digits = 2 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 3g: Cox proportional hazards for MACE ranked by quartiles of 4PY  
#            in the European Validation Cohort
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Matrix of Covariants for Cox models

my.covar.ber <- d[,c("Age", "Sex", "Diabetes",
                     "Smoker", "Hypertension", "hsCRP", "ETOH")]

##### Dyslipidemia variable

# Sex must be 1 for male, 2 for female
my.dyslipid <- function(ldl, hdl, tg, sex){
  if(all(is.na(c(ldl, hdl, tg)))){
    (NA)
  }
  else{
    dysl <- c(FALSE, FALSE, FALSE)
    if(!is.na(ldl) & ldl >= 100) dysl[1] <- TRUE
    if(!is.na(hdl) & ((hdl <= 40 & sex == 1) | (hdl <= 50 & sex == 2))) dysl[2] <- TRUE
    if(!is.na(tg)  & tg >= 150) dysl[3] <- TRUE
    
    (as.numeric(any(dysl)))
  }
}

my.covar.ber$Dyslipid <- sapply(1:nrow(d), function(x){
  my.dyslipid(d$LDL..mg.dl.[x], d$HDL..mg.dl.[x],
              d$Triglyceride..mg.dl.[x], d$Sex..M.1..F.2.[x])
})

my.covar.ber <- as.matrix(my.covar.ber)

### Unadjusted

py4.unadjusted <- rbind(
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY4.Quartile, `1` = 0, `2` = 1, .default = NaN)),
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY4.Quartile, `1` = 0, `3` = 1, .default = NaN)),
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY4.Quartile, `1` = 0, `4` = 1, .default = NaN))
)

rownames(py4.unadjusted) <- paste("Q", 2:4, sep = "")

print(py4.unadjusted) 


### Adjusted

py4.adjusted <- rbind(
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY4.Quartile, `1` = 0, `2` = 1, .default = NaN) + my.covar.ber)[1,],
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY4.Quartile, `1` = 0, `3` = 1, .default = NaN) + my.covar.ber)[1,],
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY4.Quartile, `1` = 0, `4` = 1, .default = NaN) + my.covar.ber)[1,]
)

rownames(py4.adjusted) <- paste("Q", 2:4, sep = "")

print(py4.adjusted)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Figure 3h: Cox proportional hazards for MACE ranked by quartiles of 2PY  
#            in the European Validation Cohort
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Unadjusted

py2.unadjusted <- rbind(
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY2.Quartile, `1` = 0, `2` = 1, .default = NaN)),
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY2.Quartile, `1` = 0, `3` = 1, .default = NaN)),
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY2.Quartile, `1` = 0, `4` = 1, .default = NaN))
)

rownames(py2.unadjusted) <- paste("Q", 2:4, sep = "")

print(py2.unadjusted) 


### Adjusted

py2.adjusted <- rbind(
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY2.Quartile, `1` = 0, `2` = 1, .default = NaN) + my.covar.ber)[1,],
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY2.Quartile, `1` = 0, `3` = 1, .default = NaN) + my.covar.ber)[1,],
  my.coxph(Surv(d$D_MI_Stroke, d$D_MI_StrokeYN) ~ 
             recode(d$PY2.Quartile, `1` = 0, `4` = 1, .default = NaN) + my.covar.ber)[1,]
)

rownames(py2.adjusted) <- paste("Q", 2:4, sep = "")

print(py2.adjusted)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extended Data Figure 6c: Spearman correlations among 2PY, 4PY, 
#                          and traditional risk factors for MACE.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


cormat <- cor(merge.dat[,c("X4PY", "X2PY", "Age", "LDL", "HDL", "TG","CRP", 
                           "EGFR")], method = "spearman", 
              use="pairwise.complete.obs")

melted_cormat <- melt(cormat)
melted_cormat$P <- sapply(1:nrow(melted_cormat), function(i){
  cor.test(merge.dat[,as.character(melted_cormat[i, "Var1"])] , 
           merge.dat[,as.character(melted_cormat[i, "Var2"])] ,
           method = "spearman")$p.value
})

melted_cormat$P_adjusted <- p.adjust(melted_cormat$P, method = "fdr")

print(melted_cormat)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extended Data Figure 7a: Sensitivity analysis of 4PY association 
#                          with MACE in the merged cohort
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Support functions

CoxInt <- function(x, t, e, c, thres = 0){
  
  # Establish a numeric threshold
  if(is.character(thres)){
    g = c >= median(c, na.rm = TRUE)
  }
  else if(thres == 0){
    g = c != 0
  }
  else{
    g = c >= thres
  }
  
  # Convert to quartiles within blocks
  z <- ifelse(x <= quantile(x, 0.25, na.rm=TRUE), 0,
              ifelse(x >= quantile(x, 0.75, na.rm=TRUE),1,
                     NA))
  z.hi <- ifelse(x[g] <= quantile(x[g], 0.25, na.rm=TRUE), 0,
                 ifelse(x[g] >= quantile(x[g], 0.75, na.rm=TRUE),1,
                        NA))
  z.lo <- ifelse(x[!g] <= quantile(x[!g], 0.25, na.rm=TRUE), 0,
                 ifelse(x[!g] >= quantile(x[!g], 0.75, na.rm=TRUE),1,
                        NA))
  
  pint <- coef(summary(robcov(coxph(Surv(t, e) ~ z * c))))[3,5]
  
  m.hi <- robcov(coxph(Surv(t[g], e[g]) ~ z.hi))
  m.lo <- robcov(coxph(Surv(t[!g], e[!g]) ~ z.lo))
  
  summ = data.frame(
    coef = c(exp(coef(m.hi)), exp(coef(m.lo))), 
    lo = c(exp(confint(m.hi))[1], exp(confint(m.lo))[1]), 
    hi = c(exp(confint(m.hi))[2], exp(confint(m.lo))[2]),
    N = rev(as.numeric(table(g)))
  )
  rownames(summ) = c("High", "Low")
  
  ( list(c(pint = pint), summ = summ ) )
  
}

my.sens <- function(x, t, e, c, thres, data){
  pint <- c()
  sens.tab <- data.frame(coef = numeric(), lo = numeric(), hi = numeric(), 
                         N = numeric())
  
  for(i in 1:length(c)){
    the.thres <- if(is.na(as.numeric(thres[i]))) thres[i] else as.numeric(thres[i])
    tmp <- CoxInt(data[, x], data[, t], data[, e], data[, c[i]], the.thres)
    
    pint <- c(pint, tmp[[1]])
    sens.tab <- rbind(sens.tab, tmp$summ)
    
  }
  
  names(pint)  <- c
  sens.tab$Var <- rep(c, each = 2)
  
  (list( pint, sens.tab ))
  
  
}

py4.int <- my.sens("X4PY", "DTDMS3", "DMS3", 
                   c=c("Age", "Sex", "DM", "HTN", "Smoke", "LDL", "HDL", "TG",
                       "VCAM", "BMI", "EGFR", "CVDYN", "CAD", "PAD"),
                   thres=c("med", 0, 0, 140, 0, 100, 40, 150, "med", 27, 60, 0,
                           0, 0),
                   data = merge.dat)

### Hazard ratios

print(py4.int[[2]])

### Interaction P values

print(p.adjust(py4.int[[1]], method = "fdr"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extended Data Figure 7b: Sensitivity analysis of 2PY association 
#                          with MACE in the merged cohort
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

py2.int <- my.sens("X2PY", "DTDMS3", "DMS3", 
                   c=c("Age", "Sex", "DM", "HTN", "Smoke", "LDL", "HDL", "TG",
                       "VCAM", "BMI", "EGFR", "CVDYN", "CAD", "PAD"),
                   thres=c("med", 0, 0, 140, 0, 100, 40, 150, "med", 27, 60, 0,
                           0, 0),
                   data = merge.dat)

### Hazard ratios

print(py2.int[[2]])

### Interaction P values

print(p.adjust(py2.int[[1]], method = "fdr"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extended Data Figure 8e: Pathway analysis of transcript levels in human 
#                          endothelial cells exposed to 2PY, 4PY, or Vehicle 
#                          control, as determined by RNA sequencing.
#
#  * Note custom R code was used to perform pathway analysis and and generate 
#    the heatmap shown in Extended Data Figure 8e. The heatmap shows
#    normalized transcript levels for genes anotated with the Gene Ontology term
#    "cellular response to tumor necrosis factor" (GO:0071356).
#
#  * Note V_4PY_dge is a data frame containing the differentially expressed
#    genes between Vehicle and 4PY.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#########################
#
# Pathway Enrichment Analaysis in TopGO
#
#########################

library(topGO)
library(org.Hs.eg.db)

## org.Hs.egGO is an R object with GO annotations of each gene

mapped_genes <- mappedkeys(org.Hs.egGO)

# Convert to a list
my.genes2GO <- as.list(org.Hs.egGO[mapped_genes])

my.genes2GO <- lapply(my.genes2GO, names)


genes.entrez <- read.table("Human_entrez_genes_tested_biomart.txt", sep=",",
                           header = TRUE)

gene.entrez <- unique(genes.entrez[,5:6]) ## All genes tested

geneNames <- unique(gene.entrez$NCBI.gene..formerly.Entrezgene..ID)

myInterestingGenes <- V_4PY_dge$entrez
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

GOdata <- new("topGOdata", allGenes = geneList, ontology = "BP", 
              annot = annFUN.gene2GO, gene2GO = my.genes2GO)

resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

FisTab <- GenTable(GOdata, classic = resultFis, numChar = 200,
                   orderBy = "weight", ranksOf = "classic", topNodes = 15620)

FisTab$classic <- as.numeric(FisTab$classic)

FisTab.sig <- FisTab[FisTab$classic < 0.05,]

## Label "tip" terms, ie. terms without children  

my.tips <- sapply(FisTab.sig$GO.ID, function(x){
  if(is.na(x)) {(FALSE)}
  else {!any(get(x, GOBPOFFSPRING) %in% FisTab.sig$GO.ID) }
)

FisTab.sig.tips <- FisTab.sig[my.tips,]

FisTab.sig.tips$Enrich <- FisTab.sig.tips$Significant / FisTab.sig.tips$Expected

rownames(FisTab.sig.tips) <- 1:nrow(FisTab.sig.tips)

FisTab.sig.tips[grep("kappaB", FisTab.sig.tips$Term),]
FisTab.sig.tips[grep("necrosis", FisTab.sig.tips$Term),]

GO.for.plt <- FisTab.sig.tips[order(FisTab.sig.tips$classic)[1:15],]

GO.for.plt$Term <- factor(GO.for.plt$Term, 
                          levels = rev(GO.for.plt$Term))

#########################
#
# Extended Figure 8e TNF response heatmap
#
#########################

## GO2genes

GO2genes <- as.list(org.Hs.egGO2EG)

## DGEs in TNF response
tnfgenes <- unique(
  c(V_4PY_dge[V_4PY_dge$entrez %in% GO2genes[["GO:0071356"]],"symbol"])
) 

tpm_tnf <- tpm[tpm$Gene %in% tnfgenes,]

rownames(tpm_tnf) <- tpm_tnf$Gene

my.heatmap <- function(my.tpm){
  ## Scale TPM
  my.tpm.scl <- t(scale(t(my.tpm)))
  
  # clustering
  row.order <- hclust(dist(my.tpm.scl))$order
  
  my.tpm.scl <- my.tpm.scl[row.order,]
  
  my.tpm.scl_l <- pivot_longer(as.data.frame(my.tpm.scl), cols = 1:18,
                               names_to = "Sample", values_to = "Val")
  
  my.tpm.scl_l$Gene <- rep(rownames(my.tpm.scl), each =18)
  
  my.tpm.scl_l$Sample <- factor(my.tpm.scl_l$Sample, 
                                levels = colnames(my.tpm.scl)[1:18])
  
  my.tpm.scl_l$Gene <- factor(my.tpm.scl_l$Gene, 
                              levels = rownames(my.tpm.scl))
  
  (my.tpm.scl_l)
}

hm_tnf <- my.heatmap(tpm_tnf[,c(2:13, 18:23)])

hm.plt <- ggplot(hm_tnf, aes(Sample, Gene)) +
  geom_tile(aes(fill = Val)) +
  scale_fill_gradient2(low = "blue", high = "red", limits = c(-4,4))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Supplemental Figure 5a: Kaplan-Meir estimates for MACE ranked by quartiles 
#                         of plasma sVCAM-1 in the US Validation Cohort 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Transform sVCAM-1 concentrations to quartiles

vcam   <- d$VCAM.ng.mL

q   <- quantile(vcam, c(0.25, 0.50, 0.75), na.rm=TRUE)
df1 <- which(!is.na(vcam) & vcam <= q[1])
df2 <- which(!is.na(vcam) & vcam <= q[2] & vcam > q[1])
df3 <- which(!is.na(vcam) & vcam <  q[3] & vcam > q[2])
df4 <- which(!is.na(vcam) & vcam >= q[3])

a      <- rep(NA, length(vcam))
a[df1] <- 1; a[df2] <- 2; a[df3] <- 3; a[df4] <- 4
d$VCAM.Quartile  <- a

## Kaplan-Meier plot

drawKMQuartiles(d$DTDMS3_YU, d$DMS3, d$VCAM.Quartile, 0.70)

## Log-rank test
diff <- survdiff(Surv(d$DTDMS3_YU,d$DMS3) ~ d$VCAM.Quartile)
signif(
  pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE),
  digits = 2 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Supplemental Figure 5b: Cox proportional hazards for MACE ranked by quartiles 
#                         of plasma sVCAM-1 in the US Validation Cohort 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Unadjusted

vcam.unadjusted <- rbind(
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$VCAM.Quartile, `1` = 0, `2` = 1, .default = NaN)),
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$VCAM.Quartile, `1` = 0, `3` = 1, .default = NaN)),
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$VCAM.Quartile, `1` = 0, `4` = 1, .default = NaN))
)

rownames(vcam.unadjusted) <- paste("Q", 2:4, sep = "")

print(vcam.unadjusted) 


### Adjusted

vcam.adjusted <- rbind(
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$VCAM.Quartile, `1` = 0, `2` = 1, .default = NaN) + my.covar.gb)[1,],
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$VCAM.Quartile, `1` = 0, `3` = 1, .default = NaN) + my.covar.gb)[1,],
  my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ 
             recode(d$VCAM.Quartile, `1` = 0, `4` = 1, .default = NaN) + my.covar.gb)[1,]
)

rownames(vcam.adjusted) <- paste("Q", 2:4, sep = "")

print(vcam.adjusted)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Supplemental Table 2: Baseline characteristics of patients in clinical 
#                       cohorts.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

############## d is the Discovery Cohort

### Age (median, (IQR))
signif(quantile(d$AgeAtBloodDraw, c(0.25, 0.5, 0.75)), 3)

### Sex (% male)
100*signif(table(d$Gender)["M"] / 1162, 3)

### Race (% White)
100*signif(table(d$p_White)["1"] / 1162, 3)

### Alcohol Use (% >1 drink per day)
100 * signif(table(d$EtohDrinksPerDayme > 1)["TRUE"] / 1162, 3)

### Education Level (% 4 year degree)
100 * signif(table(d$Education.Level > 3)["TRUE"] / 1162, 3)

### Diabetes Status (%)
100 * signif(table(d$DIABETICS)["1"] / 1162, 3)

### Current Smoking Status (%)
100 * signif(table(d$CurrentSmoker)["1"] / 1162, 3)

### History of CAD (%)
100 * signif(table(d$ALLCADEVER)["1"] / 1162, 3)

### History of Hypertension (%)
100 * signif(table(d$HxHtn)["1"] / 1162, 3)

### History of Hyperlipidemia (%)
100 * signif(table(d$HxHyperlipidemia)["1"] / 1162, 3)

### LDL (median, (IQR))
signif(quantile(d$LDL.Priority, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### HDL (median, (IQR))
signif(quantile(d$HDL.Priority, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### TG (median, (IQR))
signif(quantile(d$TG.Priority, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### hsCRP (median, (IQR))
signif(quantile(d$CRP16, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### eGFR (median, (IQR))
signif(quantile(d$EGFR_2021, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### Baseline aspirin (%)
100 * signif(table(d$Aspirin)["1"] / 1162, 3)

### Baseline ACE inhibitor (%)
100 * signif(table(d$CV.Ace.I.A2.Meds)["1"] / 1162, 3)

### Baseline beta blocker (%)
100 * signif(table(d$CV.Beta.Blockers)["1"] / 1162, 3)

### Baseline statin (%)
100 * signif(table(d$CV.Chol.Lowering...Statin)["1"] / 1162, 3)

### 3-year Myocardial Infarction (%)
100 * signif(table(d$MI3)["1"] / 1162, 3)

### 3-year Stroke (%)
100 * signif(table(d$STROKE3)["1"] / 1162, 3)

### 3-year Death (%)
100 * signif(table(d$DTH3)["1"] / 1162, 3)

############## gb is the US Validation Cohort

### Age (median, (IQR))
signif(quantile(gb$AgeAtBloodDraw, c(0.25, 0.5, 0.75)), 3)

### Sex (% male)
100*signif(table(gb$Gender)["M"] / 2331, 3)

### Race (% White)
100*signif(table(gb$p_White)["1"] / 2331, 3)

### Alcohol Use (% >1 drink per day)
100 * signif(table(gb$EtohDrinksPerDayme > 1)["TRUE"] / 2331, 3)

### Education Level (% 4 year degree)
100 * signif(table(gb$Education.Level > 3)["TRUE"] / 2331, 3)

### Diabetes Status (%)
100 * signif(table(gb$DIABETICS)["1"] / 2331, 3)

### Current Smoking Status (%)
100 * signif(table(gb$CurrentSmoker)["1"] / 2331, 3)

### History of CAD (%)
100 * signif(table(gb$ALLCADEVER)["1"] / 2331, 3)

### History of Hypertension (%)
100 * signif(table(gb$HxHtn)["1"] / 2331, 3)

### History of Hyperlipidemia (%)
100 * signif(table(gb$HxHyperlipidemia)["1"] / 2331, 3)

### LDL (median, (IQR))
signif(quantile(gb$LDL.Priority, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### HDL (median, (IQR))
signif(quantile(gb$HDL.Priority, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### TG (median, (IQR))
signif(quantile(gb$TG.Priority, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### hsCRP (median, (IQR))
signif(quantile(gb$CRP16, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### eGFR (median, (IQR))
signif(quantile(gb$EGFR_2021, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### Baseline aspirin (%)
100 * signif(table(gb$Aspirin)["1"] / 2331, 3)

### Baseline ACE inhibitor (%)
100 * signif(table(gb$CV.Ace.I.A2.Meds)["1"] / 2331, 3)

### Baseline beta blocker (%)
100 * signif(table(gb$CV.Beta.Blockers)["1"] / 2331, 3)

### Baseline statin (%)
100 * signif(table(gb$CV.Chol.Lowering...Statin)["1"] / 2331, 3)

### 3-year Myocardial Infarction (%)
100 * signif(table(gb$MI3)["1"] / 2331, 3)

### 3-year Stroke (%)
100 * signif(table(gb$STROKE3)["1"] / 2331, 3)

### 3-year Death (%)
100 * signif(table(gb$DTH3)["1"] / 2331, 3)

############## ber is the European Validation Cohort

### Age (median, (IQR))
signif(quantile(ber$Age, c(0.25, 0.5, 0.75)), 3)

### Sex (% male)
100*signif(table(ber$Sex)["1"] / 832, 3)

### Race (% White)
# 100*signif(table(ber$)["1"] / 832, 3)

### Alcohol Use (% >1 drink per day)
100 * signif(table(ber$ETOH)["1"] / 832, 3)

### Education Level (% 4 year degree)
100 * signif(table(ber$Education_23)["1"] / 434, 3)

### Diabetes Status (%)
100 * signif(table(ber$Diabetes)["1"] / 832, 3)

### Current Smoking Status (%)
100 * signif(table(ber$Smoker)["1"] / 832, 3)

### History of CAD (%)
100 * signif(table(ber$CAD01)["1"] / 832, 3)

### History of Hypertension (%)
100 * signif(table(ber$Hypertension)["1"] / 832, 3)

### History of Hyperlipidemia (%)
100 * signif(table(ber$Hypercholesterinaemia)["1"] / 832, 3)

### LDL (median, (IQR))
signif(quantile(ber$LDL..mg.dl., c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### HDL (median, (IQR))
signif(quantile(ber$HDL..mg.dl., c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### TG (median, (IQR))
signif(quantile(ber$Triglyceride..mg.dl., c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### hsCRP (median, (IQR))
signif(quantile(ber$hsCRP, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### eGFR (median, (IQR))
signif(quantile(ber$EGFR_2021, c(0.25, 0.5, 0.75), na.rm = TRUE), 3)

### Baseline aspirin (%)
100 * signif(table(ber$Aspirin)["1"] / 832, 3)

### Baseline ACE inhibitor (%)
100 * signif(table(ber$ACE_Inhi_AT1_R_Bloc)["1"] / 832, 3)

### Baseline beta blocker (%)
100 * signif(table(ber$Beta.Blocker)["1"] / 832, 3)

### Baseline statin (%)
100 * signif(table(ber$Statin)["1"] / 832, 3)

### 3-year Myocardial Infarction (%)
100 * signif(table(ber$MIYN)["1"] / 832, 3)

### 3-year Stroke (%)
100 * signif(table(ber$StrokeYN)["1"] / 832, 3)

### 3-year Death (%)
100 * signif(table(ber$DeathYN)["1"] / 832, 3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Supplemental Table 3/4: Previously reported plasma metabolites with known 
#                         chemical structures and those associated with 
#                         prospective residual CVD risk from untargeted 
#                         metabolomics studies of the Discovery Cohort.
#
#  * All metabolite levels were analyzed with the code below. Knowns are
#    reported in Supplemental Table 3, while Unknowns are reported in 
#    Supplemental Table 4.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my.stat <- data.frame(Name = character(), mz = numeric(), hr = numeric(), 
                      hi = numeric(),lo = numeric(), p = numeric(), 
                      model = character())

for( i in mzs[,1]){
  
  if(!(i %in% c("HILIC_posESI..385"))) # Number 385 has too little variance 
    # Compte Quartiles
    q <- quantile(d[,"i"], c(0.25, 0.75), na.rm = TRUE)
  i.quar <- ifelse(d[,"i"] <= q[1], 0, ifelse(d[,"i"] >= q[2], 1, NA))
  
  # Fit Cox Model
  my.hr <- my.coxph(Surv(d$DTDMS3_YU, d$DMS3) ~ i.quar + my.covar.discov)[1,]
  
  # Add to table
  my.stat <- rbind(my.stat, 
                   data.frame(Name = i, mz = mzs[mzs$P20ID == i,"mzP20"], 
                              hr = my.hr[,"exp(coef)"], 
                              hi = my.hr[,"upper .95"], 
                              lo = my.hr[,"lower .95"], p = my.hr[,"P"], 
                              model = "Adjusted"
                   )
  )
}

my.stat$HR.pretty <- apply(my.stat, 1, function(x){
  paste(signif(as.numeric(x[3]), 2), " (", 
        signif(as.numeric(x[5]), 2), " - ", 
        signif(as.numeric(x[4]), 2), ")", 
        sep = "")
})

my.stat$p.fdr <- p.adjust(my.stat$p, method = "fdr")

my.stat$p.round <- signif(my.stat$p, 1)

my.sig <- my.stat[my.stat$model == "Adjusted" & my.stat$p < 0.05,"Name"]


