bjs_data <- read.csv("../data/binomial-juvenile-survival-data.csv")

names(bjs_data)

RH<-unique(bjs_data$rh)
l<-length(RH)

par(mfrow=c(2,3), bty="n")
for(i in 1:l){
  w<-which(bjs_data$rh==RH[i])
  plot(Trait/n~Temp, data=bjs_data[w,],
       main = paste("RH = ", RH[i], sep=""),
       xlim=c(10,45), xlab="Temp (deg C)",
       ylim=c(0,1), ylab = "prop emerging",
       col = i)
}

library("nimble")
library("bayesTPC")
reset_models()


my_name<- "LinExp_binom"
LinExp_formula <- expression(-exp( -a*(Topt-Temp) ) + a*(Temp-Topt) + cc)
LinExp_priors <- c( cc = "dexp(0.1)",
                    a = "dnorm(0, 1)",
                    Topt  = "dnorm(25, 1/1000)"
)
LinExp_binom <- specify_binomial_model( my_name,
                                        parameters = LinExp_priors,
                                        formula = LinExp_formula
)


sim.linexp.fit.moz<-sim.quad.fit.moz<-NULL
myWAIC<-data.frame(matrix(NA, nrow=l, ncol=3))
names(myWAIC) <-c("RH", "LINEX", "QUAD")

for(i in 1:l){
  w<-which(bjs_data$rh==RH[i])
  sim.linexp.fit.moz[[i]]<-b_TPC(data = bjs_data[w,], 
                          model = "LinExp_binom",
                          inits = list(a=0.01, Topt=28),
                          burn = 15000,
                          niter = 25000)
  
  sim.quad.fit.moz[[i]] <- bayesTPC::b_TPC(data = bjs_data[w,],
                          "binomial_glm_quad", 
                          samplerType = "AF_slice",
                          priors = list(B0 = "dunif(-100,100)", 
                                        B1 = "dunif(-10,10)", 
                                        B2 = "dunif(-1,1)"),
                          inits = list(B0 = 10, B1 = 1.2, 
                                       B2 = -0.02),
                          burn = 5000, niter = 15000)
  
  par(mfrow=c(1, 2), bty="n")
  plot(sim.linexp.fit.moz[[i]], ylim=c(0,1), legend=F, 
       main = paste("LINEX, RH = ", RH[i], sep=""))
  lines(sim.quad.fit.moz[[i]], ylim=c(0,1), legend=F, 
       main=paste("Quadratic, RH = ", RH[i], sep=""))
  
  myWAIC[i,] <- c(RH[i], 
                  get_WAIC(sim.linexp.fit.moz[[i]])[1],
                  get_WAIC(sim.quad.fit.moz[[i]])[1])
                          
}


myWAIC

#########################################
## data on competence, from marta's 2020 aggregation

marta_data <- read.csv("../data/Shocket_etal_2020.csv")#, stringsAsFactors = TRUE)

## pick out non-empty columns
w<-which(!is.na(marta_data[1,]))
marta_data <- marta_data[,w]
head(marta_data)

## choose 1 mozzie species
ct_dat<-droplevels(subset(marta_data, interactor1=="Culex tarsalis"))
summary(ct_dat)

## choose 1 virus
ct_WNV<-droplevels(subset(ct_dat, interactor2=="WNV"))

Ts<-unique(ct_WNV$interactor1temp)
l <-length(Ts)

mycols<-c(4,8,18)
ct_WNV_max<-NULL

for(i in 1:l){
  tempdat<-ww<-NULL
  tempdat <- droplevels(subset(ct_WNV, interactor1temp == Ts[i]))
  print(tempdat[,mycols])
  ww <- which.max(tempdat$totalobstimevalue)
  #print(ww)
  if(is.null(ct_WNV_max)){
    ct_WNV_max <- tempdat[ww,]
  }else ct_WNV_max <- rbind(ct_WNV_max, tempdat[ww,])
  #print(ct_WNV_max[,mycols])
}

ct_WNV_max[,mycols]

plot(originaltraitvalue ~ interactor1temp, data = ct_WNV_max)

## not sure about the correct number of mozzies
ns<-c(20,40, 20, 40, 20, 20)
success <- round(ct_WNV_max$originaltraitvalue*ns)


ct.WNV.dat<-list(Trait = success,
                 Temp = ct_WNV_max$interactor1temp, 
                 n = ns)

ct.WNV.linexp.fit<-ct.WNV.quad.fit<-NULL
myWAIC<-data.frame(matrix(NA, nrow=1, ncol=2))
names(myWAIC) <-c("LINEX", "QUAD")

ct.WNV.linexp.fit<-b_TPC(data = ct.WNV.dat, 
                                 model = "LinExp_binom",
                                 inits = list(a=0.01, Topt=28),
                                 burn = 15000,
                                 niter = 25000)
  
ct.WNV.quad.fit<- bayesTPC::b_TPC(data = ct.WNV.dat,
                                  "binomial_glm_quad",
                                  samplerType = "AF_slice",
                                  priors = list(B0 = "dunif(-100,100)",
                                                B1 = "dunif(-10,10)",
                                                B2 = "dunif(-1,1)"),
                                  inits = list(B0 = 10, B1 = 1.2, 
                                               B2 = -0.02),
                                  burn = 5000, niter = 15000)
  
par(mfrow=c(1, 2), bty="n")
  plot(ct.WNV.linexp.fit, 
       ylim=c(0,1), 
       legend=F, 
       main ="LINEX")
  plot(ct.WNV.quad.fit, 
       ylim=c(0,1), 
       legend=F, 
       main="Quadratic")
  
  myWAIC<- c(LINEX = get_WAIC(ct.WNV.linexp.fit)[1],
             QUAD = get_WAIC(ct.WNV.quad.fit)[1])
  


myWAIC

