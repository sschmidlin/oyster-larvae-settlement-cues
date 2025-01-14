install.packages('lme4')
install.packages('Matrix')
install.packages('ggeffects')
install.packages('DHARMa')
install.packages("stringr")
install.packages("dplyr")
install.packages("ggeffects")
require(Matrix)
require(lme4)
require(ggeffects)
require(DHARMa)
require(stringr)
require(dplyr)

data <- read.csv2(file="Data_pilot.csv", sep=",")

#deleting 10hr, 20hr data 
data <- select(data, -settled_10hrs, -unattached_10hrs, -settled_20hrs, -unattached_20hrs)

#deleting 10hrs, 30hr data
data20 <- select(data, -settled_10hrs, -unattached_10hrs, -settled_30hrs, -unattached_30hrs)
#deleting 20hrs, 30hr data
data10 <- select(data, -settled_30hrs, -unattached_30hrs, -settled_20hrs, -unattached_20hrs)

colnames(data)[5] <- 'Settled'
colnames(data)[6] <- 'Unattached'

mult_s <- rep(1:nrow(data), data[, 'Settled'])
mult_u <- rep(1:nrow(data), data[, 'Unattached'])
data_s <- data[mult_s,]
data_s[, 'Settled'] <- 1
data_s <- data_s[, !(names(data_s) %in% c('Unattached'))]
data_u <- data[mult_u,]
data_u[, 'Unattached'] <- 0
data_u <- data_u[, !(names(data_u) %in% c('Settled'))]
colnames(data_u)[5] <- 'Settled'
data <- rbind(data_s, data_u)
data['conspecific_cue'] <- data['Cue']
data['predator_cue'] <- data['Cue']

data$conspecific_cue <- ifelse(grepl("conspecific cue", data$conspecific_cue), TRUE, FALSE)
data$predator_cue <- ifelse(grepl("predator cue", data$predator_cue), TRUE, FALSE)

data$predator_cue <-sub("TRUE", "Present_pred", data$predator_cue)
data$predator_cue <-sub("FALSE", "Absent_pred", data$predator_cue)
data$conspecific_cue <-sub("TRUE", "Present_con", data$conspecific_cue)
data$conspecific_cue <-sub("FALSE", "Absent_con", data$conspecific_cue)




###Forward selection, compare model ACI values adding all possible main effects## 


# Define the full formula with all potential main effects and interactions
full_formula <- Settled ~ conspecific_cue + predator_cue + conspecific_cue:predator_cue + Larvae.age + (1|Larvae.batch)

# Define the null model with only random effects
null_formula <- Settled ~ 1 + (1|Larvae.batch) 

# Function to fit a model
fit_model <- function(formula, data_combined) {
  model <- glmer(formula, data = data_combined, family = binomial)
  return(model)
}

# Function to calculate AIC
calculate_aic <- function(model) {
  return(AIC(model))
}

# Initialize variables
aic_values <- data.frame(Model = character(), AIC = numeric())
current_formula <- null_formula
current_model <- fit_model(current_formula, data)
best_aic <- calculate_aic(current_model)

# Predictors to consider adding
predictors_to_add <- c("conspecific_cue", "predator_cue", "conspecific_cue:predator_cue", "Larvae.age")



# Iterate through predictors to add
for (predictor in predictors_to_add) {
  # Add the predictor
  new_formula <- update(current_formula, as.formula(paste(". ~ . + ", predictor)))
  new_model <- fit_model(new_formula, data)
  new_aic <- calculate_aic(new_model)
  
  # Store AIC values in the dataframe
  aic_values <- rbind(aic_values, data.frame(Model = as.character(new_formula), AIC = new_aic))
  
  # Compare AIC values
  if (new_aic < best_aic) {
    current_formula <- new_formula
    current_model <- new_model
    best_aic <- new_aic
  }
}


# Print dataframe with AIC values
print(aic_values)


######final model, chosen for the lowest ACI value######### 

model <-glmer(Settled ~ conspecific_cue + predator_cue + conspecific_cue:predator_cue + Larvae.age + (1 | Larvae.batch), data = data, family = binomial)
summary(model)

#post hoc
emm_interaction1 <- emmeans(model, ~ conspecific_cue + predator_cue,  adjust = "tukey", type = "response")
# Perform pairwise comparisons between the levels of the interaction
pairs(emm_interaction1)



#visuals 
m <- ggpredict(model, terms = c("conspecific_cue", "predator_cue"))

plot(m) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("dodgerblue3", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))





