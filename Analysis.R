setwd()
library(car)
library(reporttools)
library(xtable)
library(corrplot)
library(MASS)

## Data

# Read in the data set
cancer <- read.csv("cancer_reg.csv", header = TRUE)

# Histogram of cancer mortality rates
hist(cancer$TARGET_deathRate, 
     xlab = "Cancer Mortality Rates per 100,000 People",
     main = "",
     col = "gray")

# Make a table of descriptive statistics
descstats <- summary(cancer)
contcancer <- cancer[, !(names(cancer) %in% c("binnedInc", "Geography", "avgAnnCount", "avgDeathsPerYear", "PctSomeCol18_24", "PctEmployed16_Over", "PctPrivateCoverageAlone"))]
tableContinuous(
  vars = contcancer,
  stats = c("mean", "median", "min", "max"),
  cap = "Descriptive Statistics of Socioeconomic Factors Affecting Cancer Mortality",
  lab = "tab:descriptivestats"
)

## Exploratory Data Analysis Grouped by Related Predictor Variables

# Correlation Matrix Heatmap
corrplot(cor(contcancer), 
         method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(100), 
         tl.col = "black", 
         tl.srt = 60, 
         cl.pos = "r",
         cl.cex = 1, 
         cl.ratio = 0.2)

# Removed: Some College, Employment Rate, and Private Coverage Alone are missing values

# Initial Model
g1 = lm(TARGET_deathRate ~ incidenceRate + medIncome + popEst2015 + povertyPercent
        + studyPerCap + MedianAge + MedianAgeMale + MedianAgeFemale + AvgHouseholdSize
        + PercentMarried + PctNoHS18_24 + PctHS18_24 + PctBachDeg18_24 + PctHS25_Over 
        + PctBachDeg25_Over + PctUnemployed16_Over + PctPrivateCoverage 
        + PctEmpPrivCoverage + PctPublicCoverage + PctPublicCoverageAlone + PctWhite 
        + PctBlack + PctAsian + PctOtherRace + PctMarriedHouseholds + BirthRate, 
        data = cancer)
summary(g1)

# ANOVA Test of Initial Model with Reduced Model
g0 = lm(TARGET_deathRate~1, data = cancer)
anova(g1, g0)

# Export table of initial model to Latex 
digits_matrix <- matrix(3, nrow = nrow(summary(g1)$coefficients), ncol = ncol(summary(g1)$coefficients)+1)
digits_matrix[, 5] <- -1
g1tab <- xtable(g1, digits = digits_matrix)
print(g1tab)

## Model Diagnostics of Initial Model

# Generate Residual Plot
plot(summary(g1)$residuals,
     main = "", 
     xlab = "",
     ylab = "Residuals", 
     pch = 20, 
     ylim = c(-150,150))
abline(h = 0, lty = 2)

# Generate Studentized Residual Plot
plot(studres(g1), 
     main = "",
     xlab = "",
     ylab = "Studentized Residuals",
     pch = 20, 
     ylim = c(-7, 7))
abline(h = 0, lty = 2)

# Plot the Residuals vs. the Fitted Values
plot(fitted.values(g1), 
     residuals(g1),
     main = "",
     xlab = expression(hat(y)), 
     ylab = "Residuals", 
     pch = 20)

# Plot the Normal Plot
qqnorm(residuals(g1), 
       pch = 20,
       main = "", 
       xlab = "Normal Scores", 
       ylab = "Ordered Residuals")

# Remove MedianAge due unusual values and insignificant contribution to the model
gt = lm(TARGET_deathRate ~ incidenceRate + medIncome + popEst2015 + povertyPercent
        + studyPerCap + MedianAgeMale + MedianAgeFemale + AvgHouseholdSize
        + PercentMarried + PctNoHS18_24 + PctHS18_24 + PctBachDeg18_24 + PctHS25_Over 
        + PctBachDeg25_Over + PctUnemployed16_Over + PctPrivateCoverage 
        + PctEmpPrivCoverage + PctPublicCoverage + PctPublicCoverageAlone + PctWhite 
        + PctBlack + PctAsian + PctOtherRace + PctMarriedHouseholds + BirthRate, 
        data = cancer)
summary(gt)

# Export table of model with MedianAge removed to Latex 
digits_matrix <- matrix(3, nrow = nrow(summary(gt)$coefficients), ncol = ncol(summary(gt)$coefficients)+1)
digits_matrix[, 5] <- -1
gttab <- xtable(gt, digits = digits_matrix)
print(gttab)

# ANOVA Test to Check Significance of Removing MedianAge
anova(gt, g1)

# Check for multiocollinearity and adjust the model
vif(gt)

# Remove PctPublicCoverage due to high VIF and insignificant contribution to the model
gt2 = lm(TARGET_deathRate ~ incidenceRate + medIncome + popEst2015 + povertyPercent
         + studyPerCap + MedianAgeMale + MedianAgeFemale + AvgHouseholdSize
         + PercentMarried + PctNoHS18_24 + PctHS18_24 + PctBachDeg18_24 + PctHS25_Over 
         + PctBachDeg25_Over + PctUnemployed16_Over + PctPrivateCoverage 
         + PctEmpPrivCoverage + PctPublicCoverageAlone + PctWhite 
         + PctBlack + PctAsian + PctOtherRace + PctMarriedHouseholds + BirthRate, 
         data = cancer)
summary(gt2)

# Export table of model with PctPublicCoverage removed to Latex 
digits_matrix <- matrix(3, nrow = nrow(summary(gt2)$coefficients), ncol = ncol(summary(gt2)$coefficients)+1)
digits_matrix[, 5] <- -1
gt2tab <- xtable(gt2, digits = digits_matrix)
print(gt2tab)

# ANOVA Test to Check Significance of Removing PctPublicCoverage
anova(gt2, gt)

vif(gt2)

## Detect and Remove Outliers from the Multicollinearity-adjusted Model

# Find the influential observations (leverage values)
gt2inf = influence(gt2)
sum(gt2inf$hat)
infindt2 = seq(1, length(gt2inf$hat))[gt2inf$hat > (2*(length(coef(gt2)))/length(gt2inf$hat))]

# Find the estimate of the standard deviation of the errors
summary(gt2)$sig

# Find the outliers using studentized residuals
studres = residuals(gt2)/(summary(gt2)$sig*sqrt(1-gt2inf$hat))
sroutt2 = seq(1, length(studres))[abs(studres) > 3]

# Remove data points that are influential observations and outliers
b = intersect(infindt2, sroutt2)
cancertest_clean <- cancer[-b, ]

# Model with Outliers Removed
gtc = lm(TARGET_deathRate ~ incidenceRate + medIncome + popEst2015 + povertyPercent
         + studyPerCap + MedianAgeMale + MedianAgeFemale + AvgHouseholdSize
         + PercentMarried + PctNoHS18_24 + PctHS18_24 + PctBachDeg18_24 + PctHS25_Over 
         + PctBachDeg25_Over + PctUnemployed16_Over + PctPrivateCoverage 
         + PctEmpPrivCoverage + PctPublicCoverageAlone + PctWhite 
         + PctBlack + PctAsian + PctOtherRace + PctMarriedHouseholds + BirthRate, 
         data = cancertest_clean)
summary(gtc)

# Export table of model with outliers removed to Latex 
digits_matrix <- matrix(3, nrow = nrow(summary(gtc)$coefficients), ncol = ncol(summary(gtc)$coefficients)+1)
digits_matrix[, 5] <- -1
gtctab <- xtable(gtc, digits = digits_matrix)
print(gtctab)

# ANOVA Test of Model with Outliers Removed and Reduced Model
gt0 = lm(TARGET_deathRate~1, data = cancertest_clean)
anova(gtc, gt0)

## Model Diagnostics of Final Model

# Generate Residual Plot
plot(summary(gtc)$residuals,
     main = "", 
     xlab = "",
     ylab = "Residuals", 
     pch = 20, 
     ylim = c(-150,150))
abline(h = 0, lty = 2)

# Generate Studentized Residual Plot
plot(studres(gtc), 
     main = "",
     xlab = "",
     ylab = "Studentized Residuals",
     pch = 20, 
     ylim = c(-7, 7))
abline(h = 0, lty = 2)

# Plot the Residuals vs. the Fitted Values
plot(fitted.values(gtc), 
     residuals(gtc),
     main = "",
     xlab = expression(hat(y)), 
     ylab = "Residuals", 
     pch = 20)

# Plot the Normal Plot
qqnorm(residuals(gtc), 
       pch = 20,
       main = "", 
       xlab = "Normal Scores", 
       ylab = "Ordered Residuals")

## Building Confidence Intervals for the True Values of Regression Coefficients
confint(gtc)