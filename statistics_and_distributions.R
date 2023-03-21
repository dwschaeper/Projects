# Investigate if T1R1-T1R3 play a role in hummingbird taste behavior
load('hummingbirds.Rdata')
t.test(hummingbirds$asp.vs.sucr.2, hummingbirds$asp.vs.sucr.1, alternative = 'greater')
t.test(hummingbirds$ery.vs.sucr.2, hummingbirds$ery.vs.sucr.1, alternative = 'greater')
t.test(hummingbirds$wat.vs.sucr.2, hummingbirds$wat.vs.sucr.1, alternative = 'greater')
t.test(hummingbirds$sucr.vs.sucr.2, hummingbirds$sucr.vs.sucr.1, alternative = 'greater')

# Determine if coffee consumption is associated with physical activity
load('coffee_exercise.Rdata')
boxplot(coffee.exercise$met ~ coffee.exercise$coffee.consumption)

qqnorm(coffee.exercise$met[coffee.exercise$coffee.consumption =='A'])
qqline(coffee.exercise$met[coffee.exercise$coffee.consumption =='A'])
qqnorm(coffee.exercise$met[coffee.exercise$coffee.consumption =='B'])
qqline(coffee.exercise$met[coffee.exercise$coffee.consumption =='B'])
qqnorm(coffee.exercise$met[coffee.exercise$coffee.consumption =='C'])
qqline(coffee.exercise$met[coffee.exercise$coffee.consumption =='C'])
qqnorm(coffee.exercise$met[coffee.exercise$coffee.consumption =='D'])
qqline(coffee.exercise$met[coffee.exercise$coffee.consumption =='D'])
qqnorm(coffee.exercise$met[coffee.exercise$coffee.consumption =='E'])
qqline(coffee.exercise$met[coffee.exercise$coffee.consumption =='E'])

summary(aov(coffee.exercise$met ~ coffee.exercise$coffee.consumption))
pairwise.t.test(coffee.exercise$met, coffee.exercise$coffee.consumption, p.adj='bonf')

# Investigate the relationship between caries and flouride content in water

# fails underlying model assumptions
load('water.Rdata')
fit = lm(water$caries ~ water$fluoride)
plot(water$fluoride, water$caries)
abline(fit)
res = resid(fit)
plot(fitted(fit), res)
abline(0,0)
# passes
load('water_new.Rdata')
fit = lm(water.new$caries ~ water.new$fluoride)
plot(water.new$fluoride, water.new$caries)
abline(fit)
res = resid(fit)
plot(fitted(fit), res)
abline(0,0)
cor(water.new$fluoride, water.new$caries)

# use binomial distribution to calculate probabilites where 37.1\% of adults and
# 57.9\% of children received a vaccine. 

# exactly 20 adults of 50
dbinom(20, 50, prob=.371)
# at most 10 children
pbinom(10, 20, .579, lower.tail = T)
# at least 11
pbinom(10, 20, .579, lower.tail = F)


# use a poisson distribution with lambda = 2
#P(X=2)
dpois(2, 2)
#P(X<=2)
ppois(2,2, lower.tail = T)
#P(X>=3)
ppois(3, 5, lower.tail = F)


# Osteosarcoma is diagnosed in approximately 8 per 1,000,000 individuals per year
# in an age group. In New York City (including all five boroughs), the number of
# young adults in this age range is approximately 1,400,000.
lambda = (8/1000000)*1400000
# P(15 or more cases)
ppois(15, lambda, lower.tail = F)
# P(10 or more in Brooklyn)
lambda = ((8/1000000)*450000)
ppois(10, lambda, lower.tail = F)


# The cholesterol levels for women aged 20 to 34 years follow an approximately 
# normal distribution with mean 185 milligrams per deciliter (mg/dl) and standard
# deviation 39 mg/dl.

# percent of women with cholesterol above 240 mg/dl
pnorm(240, mean = 185, sd = 39, lower.tail = F)
# percent with cholesterol between 240 and 200 mg/dl
pnorm(240, mean = 185, sd = 39) - pnorm(200, mean = 185, sd = 39)
# probability that no more than 5 women have levels that demand medical attention
pbinom(5, 150, .2710292, lower.tail = T)


# Hemophilia affects 1 in 5,000 male births. In the United States, there are 
# approximately 4,000,000 births per year. Assume that there are equal numbers of
# males and females born each year.

# probality at most 390 newborns have hemophilia
lambda = (1/5000)*(4000000/2)
ppois(390, lambda, lower.tail = T)
# probability 425 or more have hemophilia
ppois(425, lambda, lower.tail = F)
# calculate the standard deviation
lambda = (1/5000)*(2000000/2)
sd = lambda**(1/2)








