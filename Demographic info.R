#Demographics table data
#Note: completion rate and bounce rate calculated by hand to combine across multiple HIT calls: see OneNote page for scratch pad

summary(as.factor(participantdata$Q3.3)) #gender

age = 2020-participantdata$Q3.2 #age
sum(age<25)
sum(age>=25 & age<32)
sum(age>=32 & age<39)
sum(age>=39 & age<46)
sum(age>=46 & age<50)
mean(age)

summary(as.factor(participantdata$Q152)) #country of residence

summary(as.factor(participantdata$Q3.6)) #first language

summary(as.factor(participantdata$Q3.8)) #education
#all 4 "other"s are "some college"

summary(as.factor(participantdata$Q3.9)) #employment

summary(as.factor(participantdata$Q3.15)) #autism and nominal variants
summary(as.factor(participantdata$Q3.16_1_TEXT)) #other diagnoses
sum(is.na(participantdata$Q3.15) & is.na(participantdata$Q3.16_1_TEXT)) #has no diagnoses in either list

summary(participantdata$AQ)
summary(participantdata$SATQ)
summary(participantdata$SCC)
summary(participantdata$avAccuracy)
