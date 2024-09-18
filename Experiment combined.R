#install.packages('lme4')
#install.packages('Matrix')
#install.packages('ggeffects')
#install.packages('DHARMa')
#install.packages("stringr")
install.packages("dplyr")
#install.packages("ggeffects")
#install.packages("emmeans")
require(Matrix)
require(lme4)
require(ggeffects)
require(DHARMa)
require(stringr)
require(dplyr)
require(emmeans)

#Combining data from MAY and AUG 


#loading data from May

data <- read.csv2(file="Cue interactions 05-2022.csv", sep=",")
data <- select(data, -settled_10hr, -unattached_10hr, -settled_20hr, -unattached_20hr)
colnames(data)[5] <- 'Settled'
colnames(data)[6] <- 'Unattached'
data$Larvae.batch = as.factor(data$Larvae.batch)

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

data['Shell'] <- data['Cue']
data['conspecific_cue'] <- data['Cue']
data['predator_cue'] <- data['Cue']
data$Shell <- ifelse(grepl("conspecific shell", data$Shell), TRUE, FALSE)
data$conspecific_cue <- ifelse(grepl("conspecific cue", data$conspecific_cue), TRUE, FALSE)
data$predator_cue <- ifelse(grepl("predator cue", data$predator_cue), TRUE, FALSE)
data$Shell <-sub("TRUE", "Untreated", data$Shell)
data$Shell <-sub("FALSE", "Sterilized", data$Shell)
data$predator_cue <-sub("TRUE", "Present", data$predator_cue)
data$predator_cue <-sub("FALSE", "Absent", data$predator_cue)
data$conspecific_cue <-sub("TRUE", "Present", data$conspecific_cue)
data$conspecific_cue <-sub("FALSE", "Absent", data$conspecific_cue)


#adding data from august
data_aug <- read.csv2(file="Cue interactions 08-2022.csv", sep=",")
colnames(data_aug)[4] <- 'Well_Number'
colnames(data_aug)[2] <- 'Age'
colnames(data_aug)[3] <- 'Batch'
colnames(data_aug)[5] <- 'Cue'
colnames(data_aug)[6] <- 'Settled'
colnames(data_aug)[7] <- 'Unattached'
data_aug[, 'Cue'] <- as.factor(data_aug[, 'Cue'])


amult_s <- rep(1:nrow(data_aug), data_aug[, 'Settled'])
amult_u <- rep(1:nrow(data_aug), data_aug[, 'Unattached'])
adata_s <- data_aug[amult_s,]
adata_s[, 'Settled'] <- 1
adata_s <- adata_s[, !(names(adata_s) %in% c('Unattached'))]
adata_u <- data_aug[amult_u,]
adata_u[, 'Unattached'] <- 0
adata_u <- adata_u[, !(names(adata_u) %in% c('Settled'))]
colnames(adata_u)[6] <- 'Settled'
data_aug <- rbind(adata_s, adata_u)

data_aug['Shell'] <- data_aug['Cue']
data_aug['conspecific_cue'] <- data_aug['Cue']
data_aug['predator_cue'] <- data_aug['Cue']
data_aug['biofilm'] <- data_aug['Cue']
data_aug$Shell <- ifelse(grepl("conspecific shell", data_aug$Shell), TRUE, FALSE)
data_aug$conspecific_cue <- ifelse(grepl("conspecific cue", data_aug$conspecific_cue), TRUE, FALSE)
data_aug$predator_cue <- ifelse(grepl("predator cue", data_aug$predator_cue), TRUE, FALSE)
data_aug$biofilm <- ifelse(grepl("biofilm", data_aug$biofilm), TRUE, FALSE)
data_aug$biofilm <-sub("TRUE", "Present", data_aug$biofilm)
data_aug$biofilm <-sub("FALSE", "Absent", data_aug$biofilm)
data_aug$Shell <-sub("TRUE", "Untreated", data_aug$Shell)
data_aug$Shell <-sub("FALSE", "Sterilized", data_aug$Shell)
data_aug$predator_cue <-sub("TRUE", "Present", data_aug$predator_cue)
data_aug$predator_cue <-sub("FALSE", "Absent", data_aug$predator_cue)
data_aug$conspecific_cue <-sub("TRUE", "Present", data_aug$conspecific_cue)
data_aug$conspecific_cue <-sub("FALSE", "Absent", data_aug$conspecific_cue)


#make levels the same name 
data$Cue <- as.factor(data$Cue) 
levels(data$Cue)[levels(data$Cue) == "conspecific shell"] <- "conspecific shell_FSW"
levels(data$Cue)[levels(data$Cue) == "conspecific shell_conspecific cue"] <- "conspecific cue_conspecific shell"
levels(data$Cue)[levels(data$Cue) == "conspecific shell_predator cue"] <- "predator cue_conspecific shell"
levels(data$Cue)[levels(data$Cue) == "steralized shell"] <- "sterilized shell_FSW"
levels(data$Cue)[levels(data$Cue) == "steralized shell_conspecific cue"] <- "conspecific cue_sterlized shell"
levels(data$Cue)[levels(data$Cue) == "steralized shell_conspecific cue_predator cue"] <- "conspecific cue_predator cue_sterilized shell"
levels(data$Cue)[levels(data$Cue) == "steralized shell_predator cue"] <- "predator cue_sterlized shell"


#combining data 
data_aug_new = data_aug[,-1] #remove "Date.started" variable
colnames(data_aug_new) = c("Larvae.age", "Larvae.batch", "Tray.Number", "Cue", "Settled", "Shell", "conspecific_cue", "predator_cue", "biofilm")
colnames(data_aug_new)
data_aug_new$Larvae.batch = as.factor(data_aug_new$Larvae.batch)
exp = rep("aug", nrow(data_aug_new))
data_aug_new = cbind(data_aug_new, exp)

biofilm = rep("Absent", nrow(data))
data_may_new = cbind(data, biofilm)
exp = rep("may", nrow(data_may_new))
data_may_new = cbind(data_may_new, exp)
levels(data_aug_new$Larvae.batch) = c("3", "4")
data_combined = rbind(data_may_new, data_aug_new)
data_combined$larvae.batch <- as.factor(data_combined$Larvae.batch)


#############################
##  creating glmer models  ##
#############################

###Forward selection, compare model ACI values adding all possible main effects## 


# Define the full formula with all potential main effects and interactions
full_formula <- Settled ~ conspecific_cue + predator_cue + Shell + 
  conspecific_cue:predator_cue + predator_cue:Shell + 
  conspecific_cue:Shell + Larvae.age + (1|Larvae.batch) + (1|biofilm)

# Define the null model with only random effects
null_formula <- Settled ~ 1 + conspecific_cue + predator_cue + Shell +(1|Larvae.batch) + (1|biofilm)

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
current_model <- fit_model(current_formula, data_combined)
best_aic <- calculate_aic(current_model)

# Predictors to consider adding
predictors_to_add <- c("conspecific_cue:predator_cue", "Shell:conspecific_cue", "Shell:predator_cue", "Larvae.age")



# Iterate through predictors to add
for (predictor in predictors_to_add) {
  # Add the predictor
  new_formula <- update(current_formula, as.formula(paste(". ~ . + ", predictor)))
  new_model <- fit_model(new_formula, data_combined)
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

final_model <-glmer(Settled ~ conspecific_cue + predator_cue + Shell + conspecific_cue:predator_cue + Larvae.age + (1 | Larvae.batch) + (1|biofilm), data = data_combined, family = binomial)
summary(final_model)

#Post Hoc
emm_interaction1 <- emmeans(final_model, ~ conspecific_cue + predator_cue,  adjust = "tukey", type = "response")
emm_interaction2 <- emmeans(final_model, ~ Shell + predator_cue, adjust = "tukey", type = "response")
emm_interaction3 <- emmeans(final_model, ~ Shell + conspecific_cue, adjust = "tukey", type = "response")


# Perform pairwise comparisons between the levels of the interaction
pairs(emm_interaction1)
pairs(emm_interaction2)
pairs(emm_interaction3)


#visuals
m<-ggpredict(final_model, terms = c('conspecific_cue', 'predator_cue'))
plot(m)

m1 <-ggpredict(final_model, terms = c('Shell','predator_cue'))
plot(m1)

m2<-ggpredict(final_model, terms = c('Shell', 'conspecific_cue'))
plot(m2)


#visuals 
plot(m) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("dodgerblue3", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))


plot(m1) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("Khaki2", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))


plot(m2) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Conspecific Cue (waterbourne)")) + 
  scale_color_manual(values = c("Khaki2","dodgerblue3"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))


#####################################################
#########CREATING IMAGE FOR PUBLICATION##############
#####################################################

#Adding experiment from May so that plots can be combined##################


data1 <- read.csv2(file="Cue interactions 05-2022.csv", sep=",")
data1 <- select(data1, -settled_10hr, -unattached_10hr, -settled_20hr, -unattached_20hr)
colnames(data1)[5] <- 'Settled'
colnames(data1)[6] <- 'Unattached'
data1$Larvae.batch = as.factor(data1$Larvae.batch)

mult_s1 <- rep(1:nrow(data1), data1[, 'Settled'])
mult_u1 <- rep(1:nrow(data1), data1[, 'Unattached'])
data_s1 <- data1[mult_s1,]
data_s1[, 'Settled'] <- 1
data_s1 <- data_s1[, !(names(data_s1) %in% c('Unattached'))]
data_u1 <- data1[mult_u1,]
data_u1[, 'Unattached'] <- 0
data_u1 <- data_u1[, !(names(data_u1) %in% c('Settled'))]
colnames(data_u1)[5] <- 'Settled'
data1 <- rbind(data_s1, data_u1)

data1['Shell'] <- data1['Cue']
data1['conspecific_cue'] <- data1['Cue']
data1['predator_cue'] <- data1['Cue']
data1$Shell <- ifelse(grepl("conspecific shell", data1$Shell), TRUE, FALSE)
data1$conspecific_cue <- ifelse(grepl("conspecific cue", data1$conspecific_cue), TRUE, FALSE)
data1$predator_cue <- ifelse(grepl("predator cue", data1$predator_cue), TRUE, FALSE)
data1$Shell <-sub("TRUE", "Untreated", data1$Shell)
data1$Shell <-sub("FALSE", "Sterilized", data1$Shell)

data1$predator_cue <-sub("TRUE", "Present", data1$predator_cue)
data1$predator_cue <-sub("FALSE", "Absent", data1$predator_cue)
data1$conspecific_cue <-sub("TRUE", "Present", data1$conspecific_cue)
data1$conspecific_cue <-sub("FALSE", "Absent", data1$conspecific_cue)


model1 <-glmer(Settled ~ conspecific_cue + predator_cue + Shell + conspecific_cue:predator_cue + conspecific_cue:Shell + predator_cue:Shell + Larvae.age + (1|Larvae.batch), data=data1, family=binomial)
#visuals
m3<-ggpredict(model1, terms = c('conspecific_cue', 'predator_cue'))

m4<-ggpredict(model1, terms = c('Shell','predator_cue'))

m5<-ggpredict(model1, terms = c('Shell', 'conspecific_cue'))



####adding august experiment so plots can be combined##########################

data_aug <- read.csv2(file="Cue interactions 08-2022.csv", sep=",")
colnames(data_aug)[4] <- 'Well_Number'
colnames(data_aug)[2] <- 'Age'
colnames(data_aug)[3] <- 'Batch'
colnames(data_aug)[5] <- 'Cue'
colnames(data_aug)[6] <- 'Settled'
colnames(data_aug)[7] <- 'Unattached'
data_aug[, 'Cue'] <- as.factor(data_aug[, 'Cue'])


amult_s <- rep(1:nrow(data_aug), data_aug[, 'Settled'])
amult_u <- rep(1:nrow(data_aug), data_aug[, 'Unattached'])
adata_s <- data_aug[amult_s,]
adata_s[, 'Settled'] <- 1
adata_s <- adata_s[, !(names(adata_s) %in% c('Unattached'))]
adata_u <- data_aug[amult_u,]
adata_u[, 'Unattached'] <- 0
adata_u <- adata_u[, !(names(adata_u) %in% c('Settled'))]
colnames(adata_u)[6] <- 'Settled'
data_aug <- rbind(adata_s, adata_u)

data_aug['Shell'] <- data_aug['Cue']
data_aug['conspecific_cue'] <- data_aug['Cue']
data_aug['predator_cue'] <- data_aug['Cue']
data_aug['biofilm'] <- data_aug['Cue']
data_aug$Shell <- ifelse(grepl("conspecific shell", data_aug$Shell), TRUE, FALSE)
data_aug$conspecific_cue <- ifelse(grepl("conspecific cue", data_aug$conspecific_cue), TRUE, FALSE)
data_aug$predator_cue <- ifelse(grepl("predator cue", data_aug$predator_cue), TRUE, FALSE)
data_aug$biofilm <- ifelse(grepl("biofilm", data_aug$biofilm), TRUE, FALSE)
data_aug$biofilm <-sub("TRUE", "Present", data_aug$biofilm)
data_aug$biofilm <-sub("FALSE", "Absent", data_aug$biofilm)
data_aug$Shell <-sub("TRUE", "Untreated", data_aug$Shell)
data_aug$Shell <-sub("FALSE", "Sterilized", data_aug$Shell)
data_aug$predator_cue <-sub("TRUE", "Present", data_aug$predator_cue)
data_aug$predator_cue <-sub("FALSE", "Absent", data_aug$predator_cue)
data_aug$conspecific_cue <-sub("TRUE", "Present", data_aug$conspecific_cue)
data_aug$conspecific_cue <-sub("FALSE", "Absent", data_aug$conspecific_cue)

model_aug <-glmer(Settled ~ conspecific_cue + Shell + biofilm + predator_cue + conspecific_cue:predator_cue + conspecific_cue:biofilm + Shell:biofilm + Age + (1|Batch), data=data_aug, family=binomial)

m6<-ggpredict(model_aug, terms = c('conspecific_cue', 'predator_cue'))


m7<-ggpredict(model_aug, terms = c('conspecific_cue', 'biofilm'))


m8<-ggpredict(model_aug, terms = c('Shell', 'conspecific_cue'))


m9<-ggpredict(model_aug, terms = c('biofilm','predator_cue'))


m10<-ggpredict(model_aug, terms = c('Shell', 'predator_cue'))


m11<-ggpredict(model_aug, terms = c('Shell', 'biofilm'))




### making one image ############################
install.packages("grid")
install.packages("gridExtra")
library(grid)
library(gridExtra)

# Convert your ggeffects objects into ggplot objects
plot_m <- plot(m) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("dodgerblue3", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot_m1 <- plot(m1) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("Khaki2", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot_m2 <- plot(m2) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Conspecific Cue (waterbourne)")) + 
  scale_color_manual(values = c("Khaki2","dodgerblue3"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot_m3 <- plot(m3) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= '',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("dodgerblue3", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot_m4 <- plot(m4) + 
  labs(x = 'Conspecific Shell', 
       y= '',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("Khaki2", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot_m5 <- plot(m5) + 
  labs(x = 'Conspecific Shell', 
       y= '',
       title = "") +
  guides(color = guide_legend(title = "Conspecific Cue (waterbourne)")) + 
  scale_color_manual(values = c("Khaki2","dodgerblue3"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot_m6 <- plot(m6) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= '',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("dodgerblue3", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))
  

plot_m7 <- plot(m7) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Biofilm")) + 
  scale_color_manual(values = c("dodgerblue3", "chartreuse4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot_m8 <- plot(m8) + 
  labs(x = 'Conspecific Shell', 
       y= '',
       title = "") +
  guides(color = guide_legend(title = "Conspecific Cue (waterbourne)")) + 
  scale_color_manual(values = c("Khaki2","dodgerblue3"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot_m9 <- plot(m9) + 
  labs(x = 'Biofilm', 
       y= '',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("chartreuse4", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot_m10 <- plot(m10) + 
  labs(x = 'Conspecific Shell', 
       y= '',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("Khaki2","orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot_m11 <- plot(m11) + 
  labs(x = 'Conspecific Shell', 
       y= '',
       title = "") +
  guides(color = guide_legend(title = "Biofilm")) + 
  scale_color_manual(values = c("Khaki2","chartreuse4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

labels <- c("1A", "1B", "1C", "2A", "2B", "2C", "3A", "3B", "3C", "4A", "4B", "4C")

# List of ggplot objects
plot_list <- list(
  plot_m, plot_m3, plot_m6,
  plot_m1, plot_m4, plot_m10,
  plot_m2, plot_m5, plot_m8,
  plot_m7, plot_m11, plot_m9
)

labeled_plots <- mapply(function(p, label) {
  arrangeGrob(p, bottom = textGrob(label, gp = gpar(fontsize = 10, fontface = "bold")))
}, plot_list, labels, SIMPLIFY = FALSE)

# Arrange the plots in a 4x3 grid
grid.arrange(
  grobs = labeled_plots,
  nrow = 4,
  ncol = 3
)






########OLD CODE DO NOT RUN###################################################################################################################################
###############################################################################################################################################################



###Forward selection, compare model ACI values adding all possible main effects## 


# Define the full formula with all potential main effects and interactions
full_formula <- Settled ~ conspecific_cue + predator_cue + Shell + 
  conspecific_cue:predator_cue + predator_cue:Shell + 
  conspecific_cue:Shell + biofilm + Larvae.age + (1|Larvae.batch)

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
current_model <- fit_model(current_formula, data_combined)
best_aic <- calculate_aic(current_model)

# Predictors to consider adding
predictors_to_add <- c("conspecific_cue", "Shell", "predator_cue", "biofilm", "conspecific_cue:predator_cue", "Shell:conspecific_cue", "Shell:predator_cue", "Larvae.age")



# Iterate through predictors to add
for (predictor in predictors_to_add) {
  # Add the predictor
  new_formula <- update(current_formula, as.formula(paste(". ~ . + ", predictor)))
  new_model <- fit_model(new_formula, data_combined)
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

final_model <-glmer(Settled ~ conspecific_cue + predator_cue + Shell + conspecific_cue:predator_cue + Larvae.age + biofilm + (1 | Larvae.batch), data = data_combined, family = binomial)
summary(final_model)



#Post Hoc
emm_interaction1 <- emmeans(final_model, ~ conspecific_cue + predator_cue,  adjust = "tukey", type = "response")
emm_interaction2 <- emmeans(final_model, ~ Shell + predator_cue, adjust = "tukey", type = "response")
emm_interaction3 <- emmeans(final_model, ~ Shell + conspecific_cue, adjust = "tukey", type = "response")


# Perform pairwise comparisons between the levels of the interaction
pairs(emm_interaction1)
pairs(emm_interaction2)
pairs(emm_interaction3)


#visuals
m<-ggpredict(final_model, terms = c('conspecific_cue', 'predator_cue'))
plot(m)

m1 <-ggpredict(final_model, terms = c('Shell','predator_cue'))
plot(m1)

m2<-ggpredict(final_model, terms = c('Shell', 'conspecific_cue'))
plot(m2)


#visuals 
plot(m) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("dodgerblue3", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))


plot(m1) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("Khaki2", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))


plot(m2) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Conspecific Cue (waterbourne)")) + 
  scale_color_manual(values = c("Khaki2","dodgerblue3"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))































##########################OLD CODE DONT RUN###############################################

# combined model
model <-glmer(Settled ~ conspecific_cue + predator_cue + conspecific_cue:predator_cue + (1|Shell) + (1 | biofilm) + (1 | Larvae.batch), data = data_combined, family = binomial)
model1 <-glmer(Settled ~ predator_cue + Shell + predator_cue:Shell + (1|conspecific_cue) + (1 | biofilm) + (1 | Larvae.batch), data = data_combined, family = binomial)
model2 <-glmer(Settled ~ conspecific_cue + Shell + conspecific_cue:Shell + (1|predator_cue) + (1 | biofilm) + (1 | Larvae.batch), data = data_combined, family = binomial)

summary(model)
summary(model1)
summary(model2)

m<-ggpredict(model, terms = c('conspecific_cue', 'predator_cue'))
plot(m)

m1 <-ggpredict(model1, terms = c('Shell','predator_cue'))
plot(m1)

m2<-ggpredict(model2, terms = c('Shell', 'conspecific_cue'))
plot(m2)

#Post Hoc
emm_interaction1 <- emmeans(model, ~ conspecific_cue + predator_cue,  adjust = "tukey", type = "response")
emm_interaction2 <- emmeans(model1, ~ Shell + predator_cue, adjust = "tukey", type = "response")
emm_interaction3 <- emmeans(model2, ~ Shell + conspecific_cue, adjust = "tukey", type = "response")


# Perform pairwise comparisons between the levels of the interaction
pairs(emm_interaction1)
pairs(emm_interaction2)
pairs(emm_interaction3)

#visuals 
plot(m) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("dodgerblue3", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))


plot(m1) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("Khaki2", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))


plot(m2) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Conspecific Cue (waterbourne)")) + 
  scale_color_manual(values = c("Khaki2","dodgerblue3"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))


#figure out percent settled per treatment 

new_df <- data_combined[, c("Cue", "Settled")]
levels(new_df$Cue)
subset_data <- new_df[new_df$Cue == "predator cue_sterlized shell", ]
subset_data$Settled <- as.numeric(subset_data$Settled)
sum(subset_data$Settled)
level <- "0"
level_count <- sum(subset_data$Settled == level)
print(level_count)
(sum(subset_data$Settled))/(print(level_count))

