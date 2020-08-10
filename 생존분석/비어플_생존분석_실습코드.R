#생명표
#install.packages("KMsurv")
library(KMsurv)
tis <- 0:10 #그룹 수
ninit <- 146 #총 n
ncen <- c(3,10,10,3,3,11,5,8,1,6) #censored 수
nevent <- c(27,18,21,9,1,2,3,1,2,2) #관측된 사건 발생 수
tab<-lifetab(tis,ninit,ncen,nevent)
tab # nrisk : the estimated number of individuals at risk of experiencing the event.

plot(tab$surv,type="s",xlab="time (year)",ylab="survival prob")
plot(tab$hazard,type="s",xlab="time (year)",ylab="hazard rate")


#p.92 백혈병 처치그룹과 위약그룹에 대한 Kaplan-Meier 생존함수 및 그래프
library(survival)
#status : 0 = right censored, 1 = 관측된 값 , 2 = left censored, 3 = interval censored
time1<-c(6,6,6,6,7,9,10,10,11,13,16,17,19,20,22,23,25,25,32,32,34)
status1<-c(1,1,1,0,1,0,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0)
time2<-c(1,1,2,2,3,4,4,5,5,8,8,8,8,11,11,12,12,15,17,22,23)
status2<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

Surv(time1,status1)
fit1<-survfit(Surv(time1,status1)~1)
summary(fit1)
fit2<-survfit(Surv(time2,status2)~1)
summary(fit2)

par(mfrow=c(1,2))
plot(fit1,main="KM Curve for Group1 (treatment)",ylab="Survival probability",xlab="Time (weeks)")  
plot(fit2,main="KM Curve for Group2 (placebo)",ylab="Survival probability",xlab="Time (weeks)")

#로그-순위 검정 (그룹간에 생존함수에 차이가 있는가)
time<-c(time1,time2)
status<-c(status1,status2)
group<-c(rep(1,21),rep(2,21))
survdiff(Surv(time,status)~group) #귀무가설 기각 -> group 간에 생존함수에 차이가 있다.

#COX비례위험모형 적합 및 동점처리
data(kidney)
str(kidney)
kidney$sex<-as.factor(kidney$sex)
table(with(kidney,Surv(time,status))) # 동점자 존재

fit2<-coxph(Surv(time,status)~age+sex+disease+frail,data=kidney,ties="efron") #동점자 처리 "breslow", "exact"
summary(fit2) #frail이 한 단위 증가할 때 위험률이 6배정도 증가한다, diseasePKD가 other 일때보다 위험률이 0.114배로 감소
table(kidney$disease)

#https://stats.stackexchange.com/questions/46532/cox-baseline-hazard predicted value 계산과정
predict(fit2,type="risk") # h(t | age,sex,disease,frail)=h0(t)*exp(beta*X)

#비모수적 Kaplan-Meier 생존함수 , 다양한 변수를 사용할 경우 cox 비례위험모형을 사용함 
fit3<-survfit(Surv(time,status)~sex,data=kidney)
summary(fit3)
plot(fit3,conf.int=T,lty=c(1,2),col=c("red","blue"),main="Kaplan-Meier survival function")
legend("topright",legend=c("sex = 1","sex = 2"),lty=c(1,2),col=c("red","blue"))

#단계형 모형선택
library(MASS)
fit.step<-stepAIC(fit2)
fit.step

#잔차와 잔차그림 , status==2 는 status가 2인 경우 관측된 값, 나머지는 censoring으로 인지
data(pbc)
str(pbc)
with(pbc,Surv(time,status==2))
pbc$time ; pbc$status

fit4<-coxph(Surv(time,status==2)~age,data=pbc)
summary(fit4)

r1<-residuals(fit4,type="schoenfeld",data=pbc)
plot(pbc$age[pbc$status==2],r1,ylab="schoenfeld",xlab="age")
lines(lowess(pbc$age[pbc$status==2],r1,iter=0),lty=2,col="red")
#cox.zph 함수를 사용하면 Schoenfeld 잔차를 이용한 비례위험률 가정에 대한 검정을 실시한다. Ho: 가정만족 X
par(mfrow=c(1,2))
a<-cox.zph(fit4)
a
plot(a)

fit5<-coxph(Surv(time,status==2)~strata(sex),data=pbc)
# strata는 모든 조합을 계산해줌 (sex는 factor라 sex=1인 경우가 intercept로 들어가기 때문에 그래프에서 안나타남)
#a <- factor(rep(1:3,4), labels=c("low", "medium", "high"))
#b <- factor(rep(1:4,3))
#levels(strata(b))
#levels(strata(a,b,shortlabel=TRUE))

plot(survfit(fit5),fun="cumhaz",col=c(1,2),lty=c(1,2),xlab="time",ylab="log cumulative hazard") #그래프가 겹쳐지므로 비례위험성 가정을 만족하지않는다
legend("topleft",legend=c("sex=1","sex=2"),col=c(1,2),lty=c(1,2))






