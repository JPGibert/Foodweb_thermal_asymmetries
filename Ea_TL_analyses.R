setwd("")

library('latex2exp')

data <- read.csv("Dell_data.csv")

error.bar <- function(x, y, upper, lower=upper, length=0.05,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length	(lower) != length(upper))
	stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

dat_omn <- data[which(data$ConTrophic=="omnivore"),]
dat_herb <- data[which(data$ConTrophic=="herbivore"),]
dat_prod <- data[which(data$ConTrophic=="producer"),]
head(data)
data$ResPhylum
#---------------------------------------------------------------------------------------------------------------------------
# 1) Finding thermal sensitivities for all OMNIVORES
## Find the number of unique datasets
num_omn <- unique(dat_omn$DataSeriesID)
store_omn <- rep(0,length(num_omn))
mass_omn <- rep(0,length(num_omn))
for(i in seq(1,length(num_omn))){
	# find all datapoints for a given ID
	datID <- dat_omn[which(dat_omn$DataSeriesID==num_omn[i]),]
	## Calculate averages among rows to find max values to do the regression-
	res <- aggregate(datID[,c(5,6)], by = list(datID$AmbientTemp),mean, na.rm = TRUE)
	# if no temperature data, then move on
	if(dim(res)[1]!=0){ 
		## Find the temp at which the response is maximal so we can use that to calculate Ea
		max_T <- res$AmbientTemp[which.max(res$TraitValueSI)]
		res <- datID[which(datID$AmbientTemp <= max_T),c(5,6)]
		if(dim(res)[1]>=2){
			# Calculate Ea as the slope of log(rate) vs 1/kT
			logTrait <- log(res$TraitValueSI+0.0000001)             
			oneoverKT <- 1/((8.62*10^-5)*(273+res$AmbientTemp))             
			to_lm <- cbind(logTrait,oneoverKT)
			store_omn[i] <- lm(logTrait~oneoverKT)$coefficients[2]
			mass_omn[i] <- mean(datID$ConMassValueSI, na.rm=TRUE)
		}else{
			store_omn[i] <- NA
			mass_omn[i] <- NA
		}
	}else{
		store_omn[i] <- NA
		mass_omn[i] <- NA
	}
}

mean_omn <- -mean(store_omn,na.rm=TRUE) # we just want the magnitude
sd_omn <- sd(store_omn,na.rm=TRUE)
dat_length <- length(store_omn[!is.na(store_omn)])

c(mean_omn - sd_omn/sqrt(dat_length),mean_omn,mean_omn + sd_omn/sqrt(dat_length))

#---------------------------------------------------------------------------------------------------------------------------
# 2) Finding thermal sensitivities for all HERBIVORES
## Find the number of unique datasets
num_herb <- unique(dat_herb$DataSeriesID)
store_herb <- rep(0,length(num_herb))
mass_herb <- rep(0,length(num_herb))
for(i in seq(1,length(num_herb))){
	# find all datapoints for a given ID
	datID <- dat_herb[which(dat_herb$DataSeriesID==num_herb[i]),]
	## Calculate averages among rows to find max values to do the regression-
	res <- aggregate(datID[,c(5,6)], by = list(datID$AmbientTemp),mean, na.rm = TRUE)
	# if no temperature data, then move on
	if(dim(res)[1]!=0){ 
		## Find the temp at which the response is maximal so we can use that to calculate Ea
		max_T <- res$AmbientTemp[which.max(res$TraitValueSI)]
		res <- datID[which(datID$AmbientTemp <= max_T),c(5,6)]
		if(dim(res)[1]>=2){
			# Calculate Ea as the slope of log(rate) vs 1/kT
			logTrait <- log(res$TraitValueSI+0.0000001)             
			oneoverKT <- 1/((8.62*10^-5)*(273+res$AmbientTemp))             
			to_lm <- cbind(logTrait,oneoverKT)
			store_herb[i] <- lm(logTrait~oneoverKT)$coefficients[2] 
			mass_herb[i] <- mean(datID$ConMassValueSI, na.rm=TRUE)
		}else{
			store_herb[i] <- NA
			mass_herb[i] <- NA
		}
	}else{
		store_herb[i] <- NA
		mass_herb[i] <- NA
	}
}

mean_herb <- -mean(store_herb,na.rm=TRUE) # we just want the magnitude
sd_herb <- sd(store_herb,na.rm=TRUE)
dat_length <- length(store_herb[!is.na(store_herb)])

c(mean_herb - sd_herb/sqrt(dat_length),mean_herb,mean_herb + sd_herb/sqrt(dat_length))

#---------------------------------------------------------------------------------------------------------------------------
# 3) Finding thermal sensitivities for all PRODUCERS
## Find the number of unique datasets
num_prod <- unique(dat_prod$DataSeriesID)
store_prod <- rep(0,length(num_prod))
mass_prod <- rep(0,length(num_prod))
for(i in seq(1,length(num_prod))){
	# find all datapoints for a given ID
	datID <- dat_prod[which(dat_prod$DataSeriesID==num_prod[i]),]
	## Calculate averages among rows to find max values to do the regression-
	res <- aggregate(datID[,c(5,6)], by = list(datID$AmbientTemp),mean, na.rm = TRUE)
	# if no temperature data, then move on
	if(dim(res)[1]!=0){ 
		## Find the temp at which the response is maximal so we can use that to calculate Ea
		max_T <- res$AmbientTemp[which.max(res$TraitValueSI)]
		res <- datID[which(datID$AmbientTemp <= max_T),c(5,6)]
		if(dim(res)[1]>=2){
			# Calculate Ea as the slope of log(rate) vs 1/kT
			logTrait <- log(res$TraitValueSI+0.0000001)             
			oneoverKT <- 1/((8.62*10^-5)*(273+res$AmbientTemp))             
			to_lm <- cbind(logTrait,oneoverKT)
			store_prod[i] <- lm(logTrait~oneoverKT)$coefficients[2] 
			mass_prod[i] <- mean(datID$ConMassValueSI, na.rm=TRUE)
		}else{
			store_prod[i] <- NA
			mass_prod[i] <- NA
		}
	}else{
		store_prod[i] <- NA
		mass_prod[i] <- NA
	}
}

mean_prod <- -mean(store_prod,na.rm=TRUE) # we just want the magnitude
sd_prod <- sd(store_prod,na.rm=TRUE)
dat_length <- length(store_prod[!is.na(store_prod)])

c(mean_prod - sd_prod/sqrt(dat_length),mean_prod,mean_prod + sd_prod/sqrt(dat_length))
    
dim(across_guilds_df)    

#---------------------------------------------------------------------------------------------------------------------------
# Differences Across Trophic levels, assuming that Producers have TL=1, Herb have TL=2, and Omn have TL=3
summary(mod <- lm(log_Ea~TL,data=across_guilds_df ))$r.squared

df <- data.frame(log_Ea=seq(min(across_guilds_df$log_Ea, na.rm=TRUE),max(across_guilds_df$log_Ea, na.rm=TRUE), length.out=1000),
TL=seq(min(across_guilds_df$TL, na.rm=TRUE),max(across_guilds_df$TL, na.rm=TRUE), length.out=1000))
predict_df <- predict(mod, df,type = "response", se.fit=TRUE)
high_95 <- predict_df$fit + 1.96*predict_df$se.fit
low_95 <-predict_df$fit - 1.96*predict_df$se.fit

plot(log_Ea~jitter(TL),data=across_guilds_df, pch=16, cex=0.8, col="gray", axes=FALSE, ylab="", xlab="")
polygon(c(df$TL,rev(df$TL)),c(low_95,rev(high_95)), col=rgb(0.8,0.8,0.8,0.5), border=NA)
lines(predict_df$fit~df$TL, lwd=2)
box(lwd=2, bty='l')
axis(1,at=seq(1,3),labels=c("1","2","3"), tck=0.015, cex.axis=1.25, lwd.ticks=2, mgp=c(3, 0.5, 0))
axis(2,at=seq(-6,2), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
mtext(TeX('log Activation Energy $(E_{a})'),2, line=2,cex=1.5)
mtext(TeX('Trophic Level'),1, line=2.2,cex=1.5)


# Allowing omnivorous species to occupy random TLs from 2 to 3

sdf <- across_guilds_df
sdf_length <- length(sdf$log_Ea[which(sdf$Guild=="Omnivores")])
total <- 1000
slope <- rep(0,total)
slope <- rep(0,total)
intercept <- rep(0,total)
R_2 <- rep(0,total)
par(new = FALSE)
for(i in 1:total){
	# The first time, make the first plot
	if(i == 1){
		# Generate random TLs
		rand <- runif(sdf_length,2,3)
		sdf$TL[which(sdf$Guild=="Omnivores")] <- rand
		# Run regressions
		mod_run <- lm(log_Ea~TL,sdf)
		R_2[i] <- summary(mod_run)$r.squared
		intercept[i] <- mod_run$coefficients[1]
		slope[i] <- mod_run$coefficients[2]
		# Predict
		df_run <- data.frame(log_Ea=seq(min(sdf$log_Ea, na.rm=TRUE),max(sdf$log_Ea, na.rm=TRUE), length.out=1000),
		TL=seq(min(sdf$TL, na.rm=TRUE),max(sdf$TL, na.rm=TRUE), length.out=1000))
		predict_df_run <- predict(mod_run, df_run,type = "response", se.fit=TRUE)
		
		#plot(log_Ea~jitter(TL),data=sdf, pch=16, cex=0.5, axes=FALSE, ylab="", xlab ="", ylim=c(-0.2,0.8))
		plot(predict_df_run$fit~df_run$TL, lwd=0.5, col=rgb(0.75, 0.75, 0.75, alpha = 0.75), type="l",axes=FALSE, ylab="", xlab ="", ylim=c(-0.3,0.6))
		box(lwd=2, bty='l')
		axis(1,at=seq(1,3),labels=c("1","2","3"), tck=0.015, cex.axis=1.25, lwd.ticks=2, mgp=c(3, 0.5, 0))
		axis(2,at=seq(-0.4,0.8,0.2), tck=0.015, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))
		mtext(TeX('log Activation Energy $(E_{a})'),2, line=2,cex=1.5)
		mtext(TeX('Trophic Level'),1, line=2.2,cex=1.5)
				
	}else{
		# Generate random TLs
		rand <- runif(sdf_length,2,3)
		sdf$TL[which(sdf$Guild=="Omnivores")] <- rand
		# Run regressions
		mod_run <- lm(log_Ea~TL,sdf)
		R_2[i] <- summary(mod_run)$r.squared
		intercept[i] <- mod_run$coefficients[1]
		slope[i] <- mod_run$coefficients[2]
		# Predict
		df_run <- data.frame(log_Ea=seq(min(sdf$log_Ea, na.rm=TRUE),max(sdf$log_Ea, na.rm=TRUE), length.out=1000),
		TL=seq(min(sdf$TL, na.rm=TRUE),max(sdf$TL, na.rm=TRUE), length.out=1000))
		predict_df_run <- predict(mod_run, df_run,type = "response", se.fit=TRUE)
		lines(predict_df_run$fit~df_run$TL,lwd=0.5, col=rgb(0.75, 0.75, 0.75, alpha = 0.75))
	}	
}
abline(mean(intercept),mean(slope), lwd=2, lty="dashed")

par(new = TRUE)
par(fig = c(0.47, 0.97, 0.5, 1))
hist(slope, freq=TRUE, xlim=c(-0.35,0.05), breaks=100, las=TRUE, axes=FALSE, ylab="",xlab="", col="black", main="")
	box(bty="l", lwd=2)
	axis(1,tick=c(-3,-2,-1,0,1,2,3), tck=0.015, cex.axis=0.8,lwd.ticks=2,mgp=c(3, .2, 0))
	axis(2,tick=c(0,10,20,30,40), tck=0.015,las=TRUE, cex.axis=0.8,lwd.ticks=2,mgp=c(3, .2, 0))
	abline(v=0, col="grey", lwd=2, lty=2)
	mtext(TeX('Frequency'),2, line=1,cex=1.25)
	mtext(TeX('Fitted slope'),1, line=1.6,cex=1.25)

mean(R_2)
length(slope[slope>0])/length(slope)





## THE END