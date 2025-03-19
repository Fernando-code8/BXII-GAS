# Reference: Ray-GAS
# Created by Miguel Peña-Ramírez (miguepe@gmail.com), march/2023

####################
#### R packages ####
####################
library(forecast)
library(foreign)
library(moments)

###########################
#### R source archives ####
###########################
source("Ray-Gas_fit.R")
source("Ray-Gas.R")
source("rarma.r")
source("predict_rarma.R")
source("predict_Ray-GAS.R")

##################################
#### Organizing the data sets ####
##################################
#tempo <- proc.time()
start_time <- Sys.time()
data<-read.octave("carabas.mat")
imagem<-as.matrix(data$img)
name_patch<-c("Forest","Lake")
selectedPatch="Lake"
if (selectedPatch==name_patch[1]){
  #forest
  patch_to_fit<-imagem[101:400,][,201:325]
} else if (selectedPatch==name_patch[2]){
  #lake
  patch_to_fit<-imagem[2301:2600,][,1851:1975]
} 
rows=dim(patch_to_fit)[1]
cols=dim(patch_to_fit)[2]*0.8

############################
#### initial quantities ####
############################
mape_RayGAS<-mape_Rarma<-rmse_RayGAS<-rmse_Rarma<-
  aic_RayGAS<-aic_Rarma<-bic_RayGAS<-bic_Rarma<-
  ljBox_Rarma<-ljBox_RayGAS<-shap_RayGAS<-
  shap_Rarma<-prevRayGAS_rmse<-prevRarma_rmse<-
  prevRayGAS_MAPE<-prevRarma_MAPE<-rep(NA,rows)
p_rmse_RayGAS<-p_mape_RayGAS<-p_aic_RayGAS<-
  p_bic_RayGAS<-p_ljBox_RayGAS<-p_ljBox_Rarma<-
  p_shap_RayGAS<-p_shap_Rarma<-pprevRayGAS_rmse<-
  pprevRayGAS_mape<-0
Coeff_rarma<-Coeff<-matrix(nrow =rows ,ncol=3)

for (i in (1:rows)){  
  ############### fitting the models ############
  y<-ts(patch_to_fit[i,0:cols],start=c(1,1),frequency=1)
  rarma1<-rarma(y,ar=1,ma=1,h=1,diag=0,resid=1) 
  ray_GAS<-GASRAY.fit(y,1,1,h1=1) 
  Coeff[i,]=ray_GAS$coeff
  Coeff_rarma[i,]=rarma1$coeff
  tam_amostra=length(y)
  ############### one-step-ahead forecast ############
  # Ray-GAS
  y_out<-patch_to_fit[i,(cols+1):dim(patch_to_fit)[2]]
  y_prev_GAS<-predict_rayGAS(ray_GAS,y_out)
  prevRayGAS_rmse[i]<-forecast::accuracy(y_prev_GAS, y_out)[2]
  prevRayGAS_MAPE[i]<-forecast::accuracy(y_prev_GAS, y_out)[5]
  # RARMA
  y_prev_RARMA<-predict_rarma(rarma1,y_out)
  prevRarma_rmse[i]<- forecast::accuracy(y_prev_RARMA, y_out)[2]
  prevRarma_MAPE[i]<- forecast::accuracy(y_prev_RARMA, y_out)[5]
  
  ############### measures in the training set ############
  rmse_RayGAS[i]=forecast::accuracy(ray_GAS$fitted,y)[2];  
  mape_RayGAS[i]=forecast::accuracy(ray_GAS$fitted,y)[5]
  rmse_Rarma[i]=forecast::accuracy(rarma1$fitted,y)[2]
  mape_Rarma[i]=forecast::accuracy(rarma1$fitted,y)[5]
  
  aic_RayGAS[i]=ray_GAS$aic;  aic_Rarma[i]=rarma1$aic
  bic_RayGAS[i]=ray_GAS$bic;  bic_Rarma[i]=rarma1$bic
  shap_RayGAS[i]= as.numeric(shapiro.test(ray_GAS$residuals)[2])
  shap_Rarma[i]=  as.numeric(shapiro.test(rarma1$resid1)[2])
  ljBox_Rarma[i]= as.numeric(Box.test(rarma1$resid1, lag = 20,type = "Ljung", fitdf = 2)[3])
  ljBox_RayGAS[i]= as.numeric(Box.test(ray_GAS$residuals, lag = 20,type = "Ljung", fitdf = 2)[3])
}

# computing results
pprevRayGAS_rmse<-sum(prevRayGAS_rmse<prevRarma_rmse)
pprevRayGAS_mape<-sum(prevRayGAS_MAPE<prevRarma_MAPE)
p_mape_RayGAS<-sum(mape_RayGAS<mape_Rarma)
p_rmse_RayGAS<-sum(rmse_RayGAS<rmse_Rarma)
p_aic_RayGAS<-sum(aic_RayGAS<aic_Rarma)
p_bic_RayGAS<-sum(bic_RayGAS<bic_Rarma)
p_ljBox_both<-sum(ljBox_RayGAS > 0.05 & ljBox_Rarma> 0.05)
p_ljBox_RayGAS<-sum(ljBox_RayGAS > 0.05 & ljBox_Rarma< 0.05)
p_ljBox_Rarma<-sum(ljBox_RayGAS < 0.05 & ljBox_Rarma> 0.05)
p_shap_both<-sum(shap_RayGAS > 0.05 & shap_Rarma> 0.05)
p_shap_RayGAS<-sum(shap_RayGAS > 0.05 & shap_Rarma< 0.05)
p_shap_Rarma<-sum(shap_RayGAS < 0.05 & shap_Rarma> 0.05)

# TABLE III - NUMBER OF TIMES THE NULL HYPOTHESIS IS NOT REJECTED 
data.frame(LB=c(p_ljBox_Rarma,p_ljBox_RayGAS,p_ljBox_both),
           LB_perc=c(p_ljBox_Rarma/3,p_ljBox_RayGAS/3,p_ljBox_both/3),
           SW=c(p_shap_Rarma,p_shap_RayGAS,p_shap_both),
           SW_perc=c(p_shap_Rarma/3,p_shap_RayGAS/3,p_shap_both/3)
           )

# TABLE IV - NUMBER OF TIMES THE RAY-GAS MODEL PRESENTED THE BEST FIT
data.frame(
  results=c(p_aic_RayGAS,p_bic_RayGAS,
            p_rmse_RayGAS,p_mape_RayGAS,
            pprevRayGAS_rmse, pprevRayGAS_mape),
  results_perc=c(p_aic_RayGAS,p_bic_RayGAS,
            p_rmse_RayGAS,p_mape_RayGAS,
            pprevRayGAS_rmse, pprevRayGAS_mape)/3
)

# TABLE V - ACCURACY MEASURES FOR THE RAY-GAS AND RARMA MODELS
data.frame(
  RayGAS=c(aic_RayGAS[175],bic_RayGAS[175],
            rmse_RayGAS[175],mape_RayGAS[175],
            prevRayGAS_rmse[175], prevRayGAS_MAPE[175]),
  RARMA=c(aic_Rarma[175],bic_Rarma[175],
            rmse_Rarma[175],mape_Rarma[175],
            prevRarma_rmse[175], prevRarma_MAPE[175])
)

# Fig. 4. Observed and predicted values
y<-ts(patch_to_fit[175,0:cols],start=c(1,1),frequency=1)
y_out<-patch_to_fit[175,(cols+1):dim(patch_to_fit)[2]]
rarma1<-rarma(y,ar=1,ma=1,h=1,diag=0,resid=1) 
ray_GAS<-GASRAY.fit(y,1,1,h1=1) 
y_prev_GAS<-predict_rayGAS(ray_GAS,y_out)
y_prev_RARMA<-predict_rarma(rarma1,y_out)
par(cex.lab = 1.5, cex.axis = 1.5)
plot(y, main="",cex.lab=4, cex.axis=2, ylab="Amplitude", xlab="Signal",xlim=c(1,125))
abline(v=100,lty=2,lwd=2)
lines(101:125,y_out)
lines(ray_GAS$fitted,col="red")
lines(101:125,y_prev_GAS,col="red")
lines(rarma1$fitted,col="blue")
lines(101:125,y_prev_RARMA,col="blue")
legend("topright", #box.col = "black",box.lwd = 2,
       legend=c("Serie", "Ray-GAS", "RARMA"),
       col=c("black","red", "blue"), bty="n", 
       lty=1:2, cex=1.5)

# Fig 5 or 6 Autocorrelation functions of the residuals 
acf(ray_GAS$residuals, main="Ray-GAS")
acf(rarma1$resid1[-1], main="RARMA")
pacf(ray_GAS$residuals, main="Ray-GAS")
pacf(rarma1$resid1[-1], main="RARMA")
end_time <- Sys.time()
duration<-(end_time - start_time)
print(duration)