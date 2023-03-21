# Create a model using 2 factors to predict RFFT
library(oibiostat)
data("prevend")
boxplot(prevend$RFFT, horizontal = T)
hist(prevend$RFFT)
sd(prevend$RFFT)

hist(prevend$Age)
prevend$log.Age = log(prevend$Age)
hist(prevend$log.Age)
# data exploration
prevend$Gender = factor(prevend$Gender, levels = c(0,1), labels = c('Males', 'Females'))
barplot(table(prevend$Gender))

prevend$Education = factor(prevend$Education, levels = c(0,1,2,3), labels = c('Primary School', 'Lower Secondary Education', 'Higher Secondary School', 'University'))
barplot(table(prevend$Education))

prevend$DM = factor(prevend$DM, levels = c(0,1), labels = c('Absent', 'Present'))
barplot(table(prevend$DM))

prevend$Statin = factor(prevend$Statin, levels = c(0,1), labels = c('Non-users', 'Users'))
barplot(table(prevend$Statin))

prevend$Smoking = factor(prevend$Smoking, levels = c(0,1), labels = c('Non-smokers', 'Smokers'))
barplot(table(prevend$Smoking))

hist(prevend$BMI)
prevend$log.BMI = log(prevend$BMI)
hist(prevend$log.BMI)

hist(prevend$FRS)


# check association
pairs(~ RFFT + log.Age + Gender + Education + DM + Statin + Smoking + log.BMI + FRS, data=prevend, pch = 21, cex = 0.7)

RFFT.subset = subset(prevend, select = c(RFFT, log.Age, log.BMI, FRS))
cor(RFFT.subset)
# test models
model0 = lm(RFFT ~ log.Age + Education, data=prevend)
summary(model0)

model1 = lm(RFFT ~ log.Age + log.BMI, data=prevend)
summary(model1)

model2 = lm(RFFT ~ log.Age + FRS, data=prevend)
summary(model2)

model3 = lm(RFFT ~ Education + log.BMI, data=prevend)
summary(model3)

model4 = lm(RFFT~ Education + FRS, data=prevend)
summary(model4)

model5 = lm(RFFT ~ log.Age * Education, data=prevend)
summary(model5)

# declare best model and use it to predict
final.model = model5
qqnorm(resid(final.model))
qqline(resid(final.model))

summary(final.model)

log.Age = c(log(55))
Education = c('University')
new.frame = data.frame(log.Age, Education)

predict.lm(final.model, new.frame)

# Investigate the role female bird ornamentation in natural selection

# predict nestling fate from female ornamentation characteristics
load('rubythroats.Rdata')
nestling.fate.model = glm(nestling.fate ~ carotenoid.chroma + bib.area + total.brightness + weight + wing.length + tarsus.length, data=rubythroats, family = binomial(link='logit'))
summary(nestling.fate.model)

# Check for association between nestling fate and if a bird lays a second clutch
fledged = subset(rubythroats, nestling.fate == 'Fledged')
predated = subset(rubythroats, nestling.fate == 'Predated')
fledged.yes = subset(fledged, second.clutch == 'Yes')
fledged.no = subset(fledged, second.clutch == 'No')
predated.yes = subset(predated, second.clutch == 'Yes')
predated.no = subset(predated, second.clutch == 'No')
second.clutch.vs.nestling.fate.table = matrix(c(nrow(fledged.yes), nrow(fledged.no), nrow(predated.yes), nrow(predated.no)), nrow = 2, ncol=2, byrow = T)
dimnames(second.clutch.vs.nestling.fate.table)= list('Outcome' = c('Fledged', 'Predated'), 'Second Clutch' = c('Yes', 'No'))
chisq.test(second.clutch.vs.nestling.fate.table)

# Identify the two most significant predictors to a bird laying a second clutch
summary(glm(second.clutch ~ nestling.fate + carotenoid.chroma + bib.area + total.brightness, data=rubythroats, family = binomial(link='logit')))

# Predict if a bird lays a second clutch from the top two predictors
summary(glm(second.clutch ~ total.brightness * nestling.fate, data=rubythroats, family = binomial(link='logit')))

# Investigate factors associated with a female returning to the nestling site the subsequent year
summary(glm(survival ~ carotenoid.chroma + bib.area + total.brightness + weight + wing.length + tarsus.length + first.clutch.size + second.clutch, data=rubythroats, family = binomial(link = 'logit')))
# create a model using only the top features
more.parsimonious = glm(survival ~ carotenoid.chroma + total.brightness + wing.length + first.clutch.size, data=rubythroats, family = binomial(link = 'logit'))
summary(more.parsimonious)

# compare the odds of survival for a bird that lays 3 eggs vs 5 eggs in the first clutch
carotenoid.chroma = c(0.991766)
total.brightness = c(17.25410)
wing.length = c(53.0)
first.clutch.size = c(5)
bird1 = data.frame(carotenoid.chroma, total.brightness, wing.length, first.clutch.size)
first.clutch.size = c(3)
bird2 = data.frame(carotenoid.chroma, total.brightness, wing.length, first.clutch.size)

log.odds.bird1 = predict(more.parsimonious, bird1)
log.odds.bird2 = predict(more.parsimonious, bird2)

exp(log.odds.bird1)/exp(log.odds.bird2)
