######### Autor: Alicja Wolny-Dominiak
######### Uniwersytet Ekonomiczny w Katowicach, mailto: woali@ue.katowice.pl
######### Modele mieszane w wybranych zagadnieniach ubezpieczeniwych
#########Model częsości szkód
#########Model wartości szkody
#########Składka czysta
#########Taryfikacja a posterioru 

library(maptools)
library(spdep)
library(sp)
library(taRifx)
library(pscl)
library(tweedie)
library(hglm)
library(glmmML)
library(gamlss)
library(sfsmisc) 
library(lattice)  
library(cplm)
library(dglm)
options(OutDec=",")

auto1 = read.csv2(file="dane Polska 1.csv") 
panel= read.csv2(file="panel1.csv")
claims.nz=auto1$CLAIM_AMOUNT>0 ## niezerowe wartości szkód
claims.nzp=panel$CLAIM_AMOUNT>0 ## niezerowe wartości szkód

#######Regresory
#> names(auto1)
# [1] "OBS"           "EXPOSURE"      "PREMIUM_SPLIT" "SEX"           "VOIVODESHIP"   "CLIENT_AGE"    "CAR_MAKE"      "ENGINE"        "POWER"         "CAPACITY"      "CAR_AGE"       "CLAIM_COUNT"  
#[13] "CLAIM_AMOUNT"

### Mapa szkodowości
woj<-readShapePoly("POL1.shp")
plot(woj)
title(main="Mapa województw w Polsce")
colors<-cm.colors(16,alpha=1)
plot(woj, col=colors[1:16])
woj.nb<- poly2nb(as(woj,"SpatialPolygons")) # macierz sąsiedztwa według kryterium wspólnej granicy
woj.listw<-nb2listw(woj.nb, style="W") # macierz wag pierwszego rzędu standaryzowana wierszami
coord=coordinates(woj)
plot(woj.nb,coord, add=TRUE)
summary(woj.nb)
summary(woj.listw)
is.symmetric.nb(woj.nb)

m=tapply(auto1$CLAIM_AMOUNT[claims.nz],auto1$VOIVODESHIP[claims.nz], mean)

brks<-fivenum(auto1$CLAIM_AMOUNT[claims.nz])
cols<-c("lightblue", "red", "yellow", "green")
plot(woj, col=cols[findInterval(auto1$CLAIM_AMOUNT[claims.nz],brks)])
legend("bottomleft",legend=c(paste("poniżej ","Q1=",round(brks[2])), paste("Q1=",round(brks[2]),"-","Q2=",round(brks[3])), paste("Q2=",round(brks[3]),"-","Q3=",round(brks[4])), paste("powyżej ","Q3=",round(brks[4]))) ,fill=cols)
pointLabel(coord,as.character(woj$NAME_1), cex=0.7)

### Model wartości szkody GLM-IG
Y=glm(CLAIM_AMOUNT[claims.nz]~as.factor(SEX)[claims.nz]+as.factor(CLIENT_AGE)[claims.nz]+as.factor(POWER)[claims.nz]+as.factor(CAPACITY)[claims.nz], 
weights=CLAIM_COUNT[claims.nz], family=inverse.gaussian(log), data=auto1) 
summary(Y)

Y.d=dglm(CLAIM_AMOUNT[claims.nz]~as.factor(SEX)[claims.nz]+as.factor(CLIENT_AGE)[claims.nz]+as.factor(POWER)[claims.nz]+as.factor(CAPACITY)[claims.nz], 
weights=CLAIM_COUNT[claims.nz], family=inverse.gaussian(log), data=auto1) 
summary(Y.d)

### Model częstości szkód GLM-POIS
N=glm(CLAIM_COUNT[claims.nz]~as.factor(SEX)[claims.nz]+as.factor(CLIENT_AGE)[claims.nz]+as.factor(POWER)[claims.nz]+as.factor(CAPACITY)[claims.nz], 
offset=EXPOSURE[claims.nz], family=poisson(log), data=auto1) 
summary(N)

N.d=glm(CLAIM_COUNT[claims.nz]~+as.factor(SEX)[claims.nz]+as.factor(CLIENT_AGE)[claims.nz]+as.factor(POWER)[claims.nz]+as.factor(CAPACITY)[claims.nz], 
offset=EXPOSURE[claims.nz], family=poisson(log), data=auto1) 
summary(N.d)

##### Model częstosci szkód NLM-ZIP
tau=zeroinfl(CLAIM_COUNT~1, data=auto1)
logit.tau=tau$coefficients$zero
tau1=exp(logit.tau)/(1+exp(logit.tau)) ### prawdopodobieństwo nie wystąpienia szkody

zip=zeroinfl(CLAIM_COUNT~as.factor(SEX)+as.factor(CLIENT_AGE)+as.factor(POWER)+as.factor(CAPACITY)|as.factor(CLIENT_AGE)+as.factor(POWER)+as.factor(CAPACITY), 
offset=log(auto1$EXPOSURE), data=auto1)
summary(zip)

##### Model mieszany wartości szkody HGLM-Tweedie (p=3, rozkład IG)
Y.V=hglm(fixed=CLAIM_AMOUNT[claims.nz]~as.factor(SEX)[claims.nz]+as.factor(CLIENT_AGE)[claims.nz]+as.factor(POWER)[claims.nz]+as.factor(CAPACITY)[claims.nz],
random=~1|VOIVODESHIP[claims.nz],
weights=CLAIM_COUNT[claims.nz], family=inverse.gaussian(log), data=auto1, vcovmat = TRUE)
summary(Y.V)
plot(Y.V)
summary(Y.V)$FixCoefMa
summary(Y.V)$RandCoefMa
summary(Y.V)$varFix
summary(Y.V)$varRanef
summary(Y.V)$SummVC1
summary(Y.V)$SummVC2
summary(Y.V)

plot(quantile(auto1$CLAIM_AMOUNT[claims.nz], probs=seq(0.4,0.7,0.0001)), quantile(Y$fitted.values, probs=seq(0.4,0.7,0.0001)),xlab="Obserwacje wartości szkód dla pojedynczej polisy", 
ylab="Wartości oszacowane", main="")
abline(0,1,lwd=3,col="gray")

plot(Y)

##### Model mieszany wartości szkody HGLM-Tweedie (p=1, rozkład POIS)
N.V=hglm(fixed=CLAIM_COUNT[claims.nz]~as.factor(SEX)[claims.nz]+as.factor(CLIENT_AGE)[claims.nz]+as.factor(POWER)[claims.nz]+as.factor(CAPACITY)[claims.nz],
random=~1|VOIVODESHIP[claims.nz],
offset=EXPOSURE[claims.nz], family=poisson(log), data=auto1, vcovmat = TRUE)
summary(N.V)
plot(N.V)
summary(N.V)$FixCoefMa
summary(N.V)$RandCoefMa
summary(N.V)$varFix
summary(N.V)$varRanef
summary(N.V)$SummVC1
summary(N.V)$SummVC2
summary(N.V)

##### Model łącznej wartości szkód $GLM-CPOIS$
NC=cpglm(CLAIM_AMOUNT[claims.nz]~VOIVODESHIP[claims.nz]+as.factor(SEX)[claims.nz]+as.factor(CLIENT_AGE)[claims.nz]+as.factor(POWER)[claims.nz]+as.factor(CAPACITY)[claims.nz], 
offset=log(auto1$EXPOSURE)[claims.nz], link="log", data=auto1)
summary(NC)

NC.zi=zcpglm(CLAIM_AMOUNT~VOIVODESHIP+as.factor(SEX)+as.factor(CLIENT_AGE)+as.factor(POWER)+as.factor(CAPACITY)||as.factor(SEX)+as.factor(CLIENT_AGE)+as.factor(POWER)+as.factor(CAPACITY), 
offset=log(auto1$EXPOSURE), link="log", data=auto1)
summary(NC.zi)

logit.tau=NC.zi$coefficients$zero
tau1=exp(logit.tau)/(1+exp(logit.tau)) ### prawdopodobieństwo nie wystąpienia szkody

##### Model łącznej wartości szkód $HGLM-CPOIS$
NC.V=cpglmm(CLAIM_AMOUNT[claims.nz]~as.factor(SEX)[claims.nz]+as.factor(CLIENT_AGE)[claims.nz]+as.factor(POWER)[claims.nz]+as.factor(CAPACITY)[claims.nz]+(1|auto1$VOIVODESHIP[claims.nz]), 
offset=log(auto1$EXPOSURE)[claims.nz], link="log", data=auto1)
summary(NC.V)

##### Model HGLM-POIS
N.p=hglm(fixed=CLAIM_COUNT[claims.nzp]~as.factor(SEX)[claims.nzp]+as.factor(ENGINE)[claims.nzp]+as.factor(POWER_RANGE)[claims.nzp]+as.factor(CAR_AGE_RANGE)[claims.nzp], 
random=~1|CLIENT_NUMBER[claims.nzp],rand.family=Gamma(log),
offset=log(panel$EXPOSURE[claims.nzp]), family=poisson(log), data=panel)
summary(N.p)

####### Model HGLM-POIS modelowanie dyspersji
N.p.d=hglm(fixed=CLAIM_COUNT[claims.nzp]~as.factor(SEX)[claims.nzp]+as.factor(ENGINE)[claims.nzp]+as.factor(POWER_RANGE)[claims.nzp]+as.factor(CAR_AGE_RANGE)[claims.nzp], 
random=~1|CLIENT_NUMBER[claims.nzp],rand.family=Gamma(log), disp=~as.factor(SEX)[claims.nzp],
offset=log(panel$EXPOSURE[claims.nzp]), family=poisson(log), data=panel)
summary(N.p.d)

##### Model częstosci szkód NLM-ZIP
tau=zeroinfl(CLAIM_COUNT~1, data=auto1)
logit.tau=tau$coefficients$zero
tau1=exp(logit.tau)/(1+exp(logit.tau)) ### prawdopodobieństwo nie wystąpienia szkody

zip=zeroinfl(CLAIM_COUNT~as.factor(SEX)+as.factor(POWER)+as.factor(CAPACITY), 
offset=log(auto1$EXPOSURE), data=auto1)
summary(zip)

##### Model wartości szkód  GLM-ZA
ind <- auto1[,7]!=0 ## indykator wartości niezerowych
ind2 <- list(ind2=!(ind)) ## indykator wartości zerowych
auto1_logit<-cbind(auto1, ind2)

tau.binomial=glm(as.numeric(ind2)~1,
family=binomial(link="logit"), data=auto1_logit)
summary(tau.binomial)

tau.binomial1=glm(as.numeric(ind2)~as.factor(Wiek.kier)+as.factor(Klasa.MC)+as.factor(Wiek.poj),
family=binomial(link="logit"), data=auto1_logit)
summary(tau.binomial1)







