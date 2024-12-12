# Data Location: https://www.humaneborders.org/migrant-death-mapping
 
# Install the required packages
Packages <- c(
  "ggplot2", "viridis", "corrplot", "party", "FactoMineR", "factoextra", 
  "MuMIn", "AICcmodavg", "lme4", "drat", "xgboost", 
  "caret", "tidyverse", "dplyr", "plyr", "stringr", 
  "gsubfn", "lubridate", "xts", "fpp2", "sf", 
  "mapview", "raster", "RColorBrewer", "terra", "nnet", "tinytex")

for (pkg in Packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  } else {
    library(pkg, character.only = TRUE)  }
  print("Package Available")
}

if (!requireNamespace("tinytex", quietly = TRUE)) install.packages("tinytex")
tinytex::install_tinytex() #does not always install with others above.

setwd("C:/Users/Justin.White/OneDrive - afacademy.af.edu/Desktop/DownloadTime/HarvardStats_2024/HarvardStats_2024/Final_Capstone_Files_AZ_Border")
AZ_Data <- read.csv("./Capstone_Mortality_Dataset.csv", header=TRUE)
head(AZ_Data)
  

#Read character variables as factors
AZ_Data$PointType <- factor(AZ_Data$PointType)


# Filter out blank and 'Pending' values from Cause_of_Death
filtered_data <- AZ_Data %>%
  filter(!is.na(Cause_of_Death) & Cause_of_Death != "" & Cause_of_Death != "Pending" & Cause_of_Death != "Skeletal Remains" & Cause_of_Death != "Undetermined" & Cause_of_Death != "Not Reported")

# Create a summary of the counts for each Cause_of_Death
cause_counts <- filtered_data %>%
  dplyr::count(Cause_of_Death) %>%
  dplyr::arrange(desc(n))


# Create a bar plot using plasma color gradient
ggplot(cause_counts, aes(x = reorder(Cause_of_Death, n), y = n, fill = n)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  coord_flip() +  # Flip the bar plot horizontally
  labs(title = "Cause of Death Among \n Undocumented Immigrants in \n Arizona, USA",
       x = "Cause of Death",
       y = "Frequency" ) +
  theme_minimal(base_size = 12) +   
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(size = 9)) +
  scale_fill_viridis(option = "plasma")  


#Map the data to visualize spatial arrangement of mortality sites.
loc_AZ_Data <- st_as_sf(AZ_Data, coords=c("Longitude","Latitude"), crs=4326)
mapview(loc_AZ_Data,
        zcol="PointType",
        cex=4,
        burst=TRUE,
        alpha = 0,
        alpha.regions=0.2)
 

 
#Investigate patterns in the data

#Basic correlation plot with continuous variables to inform upcoming regression models and avoid multicollinearity.
AZ_Data_Cor <- AZ_Data[6:17]
head(AZ_Data_Cor)

CorTest <- cor(AZ_Data_Cor, use = "complete.obs")

corrplot(CorTest, method = 'circle', order = 'FPC', type = 'lower', diag = FALSE)
 



##### Generalized Linear Model Tests. Models assembled based on Corrplot above #####

#Use AIC models to examine the relative quality of different statistical models by comparing their goodness-of-fit while penalizing for model complexity, thereby helping to identify the most parsimonious model that adequately explains the data.

AZ_Data_AIC <- AZ_Data[5:18]

#Split the data for GLM_training purposes. This yields an 85/15% split of the data
GLM_train <- createDataPartition(AZ_Data_AIC$PointType, p = 0.85, list = FALSE, times = 1)  
GLM_training <- AZ_Data_AIC[GLM_train,]
GLM_testing <-AZ_Data_AIC[ -GLM_train,]  

AZ_Mod.1 <-glm(PointType ~ Rough + I(Rough ^2), data = GLM_training, family = binomial) #include quadratic pattern of variable
AZ_Mod.2 <-glm(PointType ~ DiffuseRad + I(DiffuseRad ^2) + MSAVI2, data = GLM_training, family = binomial)
AZ_Mod.3 <-glm(PointType ~ Rough + I(Rough ^2) + TotalSolar + I(TotalSolar ^2), data = GLM_training, family = binomial)
AZ_Mod.4 <-glm(PointType ~ RelCont9am + I(RelCont9am ^2), data = GLM_training, family = binomial)
AZ_Mod.5 <-glm(PointType ~ RelCont9am + I(RelCont9am ^2) + Rough, data = GLM_training, family = binomial)
AZ_Mod.6 <-glm(PointType ~ RelCont6pm + I(RelCont6pm ^2) + Slope, data = GLM_training, family = binomial)
AZ_Mod.7 <-glm(PointType ~ TotalSolar + I(TotalSolar ^2), data = GLM_training, family = binomial)
AZ_Mod.8 <-glm(PointType ~ DirecDurat + I(DirecDurat ^2), data = GLM_training, family = binomial)
AZ_Mod.9 <-glm(PointType ~ DirecDurat + I(DirecDurat ^2) + Rough , data = GLM_training, family = binomial)
AZ_Mod.11 <-glm(PointType ~ TotalSolar, data = GLM_training, family = binomial)
AZ_Mod.10 <-glm(PointType ~ TotalSolar + TotAveArid, data = GLM_training, family = binomial)
AZ_Mod.12 <-glm(PointType ~ MSAVI2 + I(MSAVI2 ^2),  data = GLM_training, family = binomial)
AZ_Mod.13 <-glm(PointType ~ TotAveArid, data = GLM_training, family = binomial)
AZ_Mod.14 <-glm(PointType ~ -1, data = GLM_training, family = binomial)

AZ_Cands <- list(AZ_Mod.1, AZ_Mod.2, AZ_Mod.3, AZ_Mod.4, AZ_Mod.5, AZ_Mod.6, AZ_Mod.7, AZ_Mod.8, AZ_Mod.9, AZ_Mod.10, AZ_Mod.11, AZ_Mod.12, AZ_Mod.13, AZ_Mod.14)

AZ_Mods <- c("AZ_Mod.1", "AZ_Mod.2", "AZ_Mod.3", "AZ_Mod.4", "AZ_Mod.5", "AZ_Mod.6",  "AZ_Mod.7", "AZ_Mod.8", "AZ_Mod.9", "AZ_Mod.10", "AZ_Mod.11", "AZ_Mod.12", "AZ_Mod.13", "AZ_Mod.14")

AZ_aictable <- aictab(AZ_Cands, AZ_Mods, sort=T)

AZ_aictable

#VIFs
car::vif(AZ_Mod.9)

#Results:
summary(AZ_Mod.9)

# Ensure test labels are read as factors
GLM_test_labels <- as.factor(GLM_testing$PointType) 

# Predicted probabilities
GLM_pred_probs <- predict(AZ_Mod.9, newdata = GLM_testing, type = "response")

GLM_pred_labels <- ifelse(GLM_pred_probs > 0.1, "BodyLocation", "RandomLocation")

# Convert predictions to factor with consistent levels
GLM_pred_labels <- factor(GLM_pred_labels, levels = c("BodyLocation", "RandomLocation"))

GLM_conf_matrix <- confusionMatrix(GLM_pred_labels, GLM_test_labels)
print(GLM_conf_matrix)

 






################################################
################################################


#XGBoost Machine Learning Application - useful as it handles large datasets efficiently and supports regularization which helps avoid overfitting.

AZ_Data_Boost <- AZ_Data_AIC
head(AZ_Data_Boost)
 
#Split the data for Boost_training purposes. This yields a 70/30% split of the data
Boost_train <- createDataPartition(AZ_Data_Boost$PointType, p = 0.8, list = FALSE, times = 1) #Set list to false to return matrix or array.
Boost_training <- AZ_Data_Boost[Boost_train,]
dim(Boost_training)
Boost_testing <-AZ_Data_Boost[ -Boost_train,] #Rename the data
dim(Boost_testing)

#Create matrix
Boost_train_labels <- as.numeric(as.factor(Boost_training$PointType)) - 1
Boost_train_matrix <- xgb.DMatrix(data = as.matrix(Boost_training[, -which(names(Boost_training) == "PointType")]), label = Boost_train_labels)

Boost_test_labels <- as.numeric(as.factor(Boost_testing$PointType)) - 1
Boost_test_matrix <- xgb.DMatrix(data = as.matrix(Boost_testing[, -which(names(Boost_testing) == "PointType")]), label = Boost_test_labels)


#parameters
nc <- length(unique(Boost_train_labels))#get two binary reponse outputs
xgb_params <- list("objective" = "multi:softprob", "eval_metric" = "mlogloss", "num_class" = nc)#change eval_metric based on regression/model type.

watchlist <- list(Boost_train = Boost_train_matrix, Boost_test = Boost_test_matrix)#gives us updates on what's happening when the model is running
bst_model <-xgb.train(params = xgb_params,
                      data = Boost_train_matrix,
                      nrounds = 200,#iterations
                      watchlist = watchlist,
                      eta = 0.005,
                      max.depth = 7, 
                      subsample = 0.85,
                      verbose = 0)# Use 'eta' learning rate function (which should have values from 0-1) with lower values reducing overfitting. Depth is tree level, 7 resulted in lowest Boost_testing error. Gamma default is 0 (values are 0 to inf): larger values are more conservative and help avoid overfitting. Subsample (values 0-1): lower values also aim to reduce overfitting. The idea is that it controls the number of observations supplied to a tree. We can also determine which variables are allowed to interact, see https://xgboost.readthedocs.io/en/stable/tutorials/feature_interaction_consBoost_traint.html).

bst_model#We can use the specs from the model to create a graph of the Boost_training and Boost_testing errors
#Boost_training error with Boost_test error as lines
e <- data.frame(bst_model$evaluation_log)
plot(e$iter, e$Boost_train_mlogloss, col='blue', ylim = c(0.0, 1))#Set y-axis limit
lines(e$iter, e$Boost_test_mlogloss, col='red')


#Show which iteration had the lowest Boost_testing error.
min(e$Boost_test_mlogloss)

#Feature Importance 
importance <- xgb.importance(colnames(Boost_train_matrix), model = bst_model)
importance

#Plot based on Gain: 
# Sort the data by Gain (descending) for better visualization
importance <- importance %>%
  arrange(desc(Gain))

# Plot using ggplot2, ChatGPT was used to help build this plot.
# Create the bar plot with plasma gradient
ggplot(importance, aes(x = reorder(Feature, Gain), y = Gain, fill = Gain)) +
  geom_bar(stat = "identity", color = "black", show.legend = TRUE) +
  scale_fill_viridis_c(option = "plasma") +  # Use the plasma color gradient
  coord_flip() +  # Flip the coordinates for a horizontal bar plot
  labs(
    title = "Feature Importance by Gain",
    x = "Features",
    y = "Gain"
  ) +
  theme_minimal(base_size = 15) +  # Apply a minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 15), # Center and bold title
    axis.title.x = element_text(),
    axis.title.y = element_text(),
    axis.text = element_text(size = 12))
# Make predictions on the Boost_test set
pred_probs <- predict(bst_model, Boost_test_matrix)
pred_labels <- max.col(matrix(pred_probs, ncol = length(unique(Boost_train_labels)))) - 1

# Assess performance
Boost_conf_matrix <- confusionMatrix(factor(pred_labels), factor(Boost_test_labels))
print(Boost_conf_matrix)


####################################################
# Let's try with a slightly more complicated technique such as neural networks in 'nnet'

# Assume AZ_Data_Boost is already loaded

# Set seed for reproducibility
set.seed(123)

# Partition the data (80% NN_training, 20% NN_testing)
NN_train_index <- createDataPartition(AZ_Data_Boost$PointType, p = 0.8, list = FALSE)
NN_train_data <- AZ_Data_Boost[NN_train_index, ]
NN_test_data <- AZ_Data_Boost[-NN_train_index, ]


# Normalize the continuous variables
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

NN_train_data[, -which(names(NN_train_data) == "PointType")] <- as.data.frame(lapply(
  NN_train_data[, -which(names(NN_train_data) == "PointType")], normalize
))
NN_test_data[, -which(names(NN_test_data) == "PointType")] <- as.data.frame(lapply(
  NN_test_data[, -which(names(NN_test_data) == "PointType")], normalize
))

# Convert PointType to binary numeric (0, 1)
NN_train_data$PointType <- ifelse(NN_train_data$PointType == "BodyLocation", 1, 0)
NN_test_data$PointType <- ifelse(NN_test_data$PointType == "BodyLocation", 1, 0)

# NN_train the Neural Network
nn_model <- nnet(
  PointType ~ ., 
  data = NN_train_data, 
  size = 10,  # Number of hidden nodes
  decay = 0.01,  # Regularization
  maxit = 200  # Maximum iterations
)

# Make predictions
NN_test_preds <- predict(nn_model, NN_test_data, type = "raw")
pred_labels <- ifelse(NN_test_preds > 0.5, 1, 0)

# Evaluate the model
NN_conf_matrix <- confusionMatrix(factor(pred_labels), factor(NN_test_data$PointType))
print(NN_conf_matrix)