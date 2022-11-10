# Read functions
source("functions.R")

# Set sample size
n <- 10000
# Set the number of the explanatory and the response variables
p <- 3
k <- 2
# Set the values of the coefficients used in data generation
coefX <- c(1, 1, 1)
coefY <- c(-20)

# Generate clean data based on the classifier chain model - synthetic data 1 
# as described in section 4.2.1. The coefficients used in this generation are
# specified by coefX and coefY.
cleandata <- generatesd1(n, p, k, coefX, coefY)
# Assign the values of X and Y
X <- cleandata$X
Y <- cleandata$Y

# Add outliers to the data. The function addoutlierssd1 takes as arguments the 
# clean data, percentage of outliers to be added to the data, and the values of
# the coefficients used in generation. The indexes of the observations that 
# will become outliers are chosen at random and their labels are altered
# as described in section 4.2.1. 
outlierspercent <- 0.01
outliersinfo <- addoutlierssd1(X, Y, outlierspercent, coefX, coefY)
# Get the indexes of the outliers and the new Y
outliers <- outliersinfo$outliers
Y <- outliersinfo$Ynew

# Fit the classifier chain model to the data with outliers. The underlying 
# model is ridge logistic regression. 
cc <- classchain(X, Y)
# Calculate conditional outlier scores based on the cc model
outlierscores <- -log(estprobcc(cc, X, Y))

# Plot conditional outlier scores in decreasing order
(figure <- ggplot() + 
    geom_point(aes(x = 1:length(outlierscores), 
                   y = sort(outlierscores, decreasing = T))) +
    theme_bw() + 
    labs(x = "Observation", y = "Score"))

# Specify the grid of the alert rates
alertrates <- seq(0.005, 0.045, 0.001)
# Calculate the precision and recall for the grid of alert rates
precisionrecall <- precisionrecall(outliers, alertrates, outlierscores)
# Plot PAR, RAR, and PR curves
plotprecisionrecall(precisionrecall, alertrates)
