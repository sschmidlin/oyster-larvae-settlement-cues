install.packages('lme4')
install.packages('Matrix')
install.packages('ggeffects')
install.packages('DHARMa')
install.packages("stringr")
require(lme4)
require(ggeffects)
require(DHARMa)
require(stringr)


#loading data from May
setwd("~/GitHub/Conspecific-Predator-Settlement-Cue")

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

#models
model <-glmer(Settled ~ conspecific_cue + predator_cue + conspecific_cue:predator_cue + (1|Shell) + (1|Larvae.batch), data=data, family=binomial)
model1 <-glmer(Settled ~ predator_cue + Shell + predator_cue:Shell + (1|conspecific_cue) + (1|Larvae.batch), data=data, family=binomial)
model2 <- glmer(Settled ~ conspecific_cue + Shell + conspecific_cue:Shell + (1|predator_cue) + (1|Larvae.batch), data=data, family=binomial)

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

