install.packages('lme4')
install.packages('Matrix')
install.packages('ggeffects')
install.packages('DHARMa')
install.packages("stringr")
require(lme4)
require(ggeffects)
require(DHARMa)
require(stringr)

data <- read.csv2(file="Data_pilot.csv", sep=",")

#deleting 10hr, 20hr data 
data <- select(data, -settled_10hrs, -unattached_10hrs, -settled_20hrs, -unattached_20hrs)
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

data$predator_cue <-sub("TRUE", "Present", data$predator_cue)
data$predator_cue <-sub("FALSE", "Absent", data$predator_cue)
data$conspecific_cue <-sub("TRUE", "Present", data$conspecific_cue)
data$conspecific_cue <-sub("FALSE", "Absent", data$conspecific_cue)


# making the models
model <- glmer(Settled ~ predator_cue + conspecific_cue + conspecific_cue:predator_cue + (1|Larvae.batch), data = data, family = binomial)
summary(model)
m <- ggpredict(model, terms = c("conspecific_cue", "predator_cue"))
plot(m)

#post hoc
emm_interaction1 <- emmeans(model, ~ conspecific_cue + predator_cue,  adjust = "tukey", type = "response")
# Perform pairwise comparisons between the levels of the interaction
pairs(emm_interaction1)

#visuals 
plot(m) + 
  labs(x = 'Conspecific Cue (waterbourne)', 
       y= 'Larvae Settled (%)',
       title = "") +
  guides(color = guide_legend(title = "Predator Cue")) + 
  scale_color_manual(values = c("dodgerblue3", "orangered4"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  scale_y_continuous(labels= function(x) paste0(x*100), limits = c(0,1))

