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
data_aug$biofilm <-sub("TRUE", "biofilm_present", data_aug$biofilm)
data_aug$biofilm <-sub("FALSE", "biofilm_absent", data_aug$biofilm)
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

biofilm = rep("absent", nrow(data))
data_may_new = cbind(data, biofilm)
exp = rep("may", nrow(data_may_new))
data_may_new = cbind(data_may_new, exp)
levels(data_aug_new$Larvae.batch) = c("3", "4")
data_combined = rbind(data_may_new, data_aug_new)

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

