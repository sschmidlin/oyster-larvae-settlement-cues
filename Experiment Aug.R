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

#data from august
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

#############################
##  creating glmer models  ##
#############################

###Forward selection, compare model ACI values adding all possible main effects## 

# Define the full formula with all potential main effects and interactions
full_formula <- Settled ~ conspecific_cue + Shell + biofilm + predator_cue +
  conspecific_cue:predator_cue + predator_cue:Shell + 
  conspecific_cue:Shell + biofilm:conspecific_cue + biofilm:predator_cue + biofilm:Shell + Age (1|Batch)

# Define the null model with only random effects
null_formula <- Settled ~ 1 + conspecific_cue + predator_cue + Shell + biofilm + (1|Batch) 

# Function to fit a model
fit_model <- function(formula, data_aug) {
  model <- glmer(formula, data = data_aug, family = binomial)
  return(model)
}

# Function to calculate AIC
calculate_aic <- function(model) {
  return(AIC(model))
}

# Initialize variables
aic_values <- data.frame(Model = character(), AIC = numeric())
current_formula <- null_formula
current_model <- fit_model(current_formula, data_aug)
best_aic <- calculate_aic(current_model)

# Predictors to consider adding
predictors_to_add <- c("Age", "conspecific_cue:predator_cue", "Shell:conspecific_cue", "Shell:predator_cue", "biofilm:conspecific_cue", " biofilm:predator_cue", "biofilm:Shell")

# Iterate through predictors to add
for (predictor in predictors_to_add) {
  # Add the predictor
  new_formula <- update(current_formula, as.formula(paste(". ~ . + ", predictor)))
  new_model <- fit_model(new_formula, data_aug)
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

model <-glmer(Settled ~ conspecific_cue + Shell + biofilm + predator_cue + conspecific_cue:predator_cue + conspecific_cue:biofilm + Shell:biofilm + Age + (1|Batch), data=data_aug, family=binomial)
summary(model)



#Post Hoc
emm_interaction <- emmeans(model, ~ conspecific_cue + predator_cue, adjust = "tukey", type = "response")
emm_interaction1 <- emmeans(model, ~ conspecific_cue + biofilm,  adjust = "tukey", type = "response")
emm_interaction2 <- emmeans(model, ~ conspecific_cue + Shell, adjust = "tukey", type = "response")
emm_interaction3 <- emmeans(model, ~ predator_cue + biofilm, adjust = "tukey", type = "response")
emm_interaction4 <- emmeans(model, ~ predator_cue + Shell, adjust = "tukey", type = "response")
emm_interaction5 <- emmeans(model, ~ Shell + biofilm, adjust = "tukey", type = "response")


# Perform pairwise comparisons between the levels of the interaction
pairs(emm_interaction)
pairs(emm_interaction1)
pairs(emm_interaction2)
pairs(emm_interaction3)
pairs(emm_interaction4)
pairs(emm_interaction5)


m<-ggpredict(model, terms = c('conspecific_cue', 'predator_cue'))
plot(m)

m1 <-ggpredict(model, terms = c('conspecific_cue', 'biofilm'))
plot(m1)

m2<-ggpredict(model, terms = c('Shell', 'conspecific_cue'))
plot(m2)

m3<-ggpredict(model, terms = c('biofilm','predator_cue'))
plot(m3)

m4<-ggpredict(model, terms = c('Shell', 'predator_cue'))
plot(m4)

m5<-ggpredict(model, terms = c('Shell', 'biofilm'))
plot(m5)

#visuals 
plot(m) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("dodgerblue3", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))+
  geom_errorbar(data = pairwise_emm, aes(x = conspecific_cue, ymin = lower.CL, ymax = upper.CL, color = predator_cue), width = 0.2, position = position_dodge(width = 0.2))


plot(m1) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Biofilm")) + 
  scale_color_manual(values = c("dodgerblue3", "chartreuse4"))+
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

plot(m3) + 
  labs(x = 'Biofilm', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("chartreuse4", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot(m4) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("Khaki2","orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot(m5) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Biofilm")) + 
  scale_color_manual(values = c("Khaki2","chartreuse4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

























##########OLD CODE DO NOT RUN##############################

#model with con and predator cue interaction 
model <-glmer(Settled ~ conspecific_cue + predator_cue + conspecific_cue:predator_cue + (1|biofilm) + (1|Shell) + (1|Batch), data=data_aug, family=binomial)
model1 <-glmer(Settled ~ conspecific_cue + biofilm + conspecific_cue:biofilm + (1|Shell) + (1|predator_cue) + (1|Batch), data=data_aug, family=binomial)
model2 <-glmer(Settled ~ conspecific_cue + Shell + conspecific_cue:Shell + (1|predator_cue) + (1|biofilm) + (1|Batch), data=data_aug, family=binomial)
model3 <-glmer(Settled ~ predator_cue + biofilm + predator_cue:biofilm + (1|conspecific_cue) + (1|Shell) + (1|Batch), data=data_aug, family=binomial)
model4 <-glmer(Settled ~ predator_cue + Shell + predator_cue:Shell + (1|conspecific_cue) + (1|biofilm) + (1|Batch), data=data_aug, family=binomial)
model5 <-glmer(Settled ~ Shell + biofilm + Shell:biofilm + (1|conspecific_cue) + (1|predator_cue) + (1|Batch), data=data_aug, family=binomial)


summary(model)
summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)

m<-ggpredict(model, terms = c('conspecific_cue', 'predator_cue'))
plot(m)

m1 <-ggpredict(model1, terms = c('conspecific_cue', 'biofilm'))
plot(m1)

m2<-ggpredict(model2, terms = c('Shell', 'conspecific_cue'))
plot(m2)

m3<-ggpredict(model3, terms = c('biofilm','predator_cue'))
plot(m3)

m4<-ggpredict(model4, terms = c('Shell', 'predator_cue'))
plot(m4)

m5<-ggpredict(model5, terms = c('Shell', 'biofilm'))
plot(m5)

#Post Hoc
emm_interaction <- emmeans(model, ~ conspecific_cue + predator_cue, adjust = "tukey", type = "response")
emm_interaction1 <- emmeans(model1, ~ conspecific_cue + biofilm,  adjust = "tukey", type = "response")
emm_interaction2 <- emmeans(model2, ~ conspecific_cue + Shell, adjust = "tukey", type = "response")
emm_interaction3 <- emmeans(model3, ~ predator_cue + biofilm, adjust = "tukey", type = "response")
emm_interaction4 <- emmeans(model4, ~ predator_cue + Shell, adjust = "tukey", type = "response")
emm_interaction5 <- emmeans(model5, ~ Shell + biofilm, adjust = "tukey", type = "response")


# Perform pairwise comparisons between the levels of the interaction
pairs(emm_interaction)
pairs(emm_interaction1)
pairs(emm_interaction2)
pairs(emm_interaction3)
pairs(emm_interaction4)
pairs(emm_interaction5)

#visuals 
plot(m) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("dodgerblue3", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))+
 geom_errorbar(data = pairwise_emm, aes(x = conspecific_cue, ymin = lower.CL, ymax = upper.CL, color = predator_cue), width = 0.2, position = position_dodge(width = 0.2))


plot(m1) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Biofilm")) + 
  scale_color_manual(values = c("dodgerblue3", "chartreuse4"))+
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

plot(m3) + 
  labs(x = 'Biofilm', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("chartreuse4", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot(m4) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("Khaki2","orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

plot(m5) + 
  labs(x = 'Conspecific Shell', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Biofilm")) + 
  scale_color_manual(values = c("Khaki2","chartreuse4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))
