w1<-4 # width for plots 
h11<-4 # height for plots

setwd("~/GitHub/BXII-GAS/Apply BXII-GAS")
j=i=5 # 1,4,5,11,12,13,16
resid_plot<-paste0("Plots/Acf",i,".pdf")
pdf(resid_plot,width = w1, height = h11)
# acf BXII residuals 
acf(final[[j]]$residuals,main="")
dev.off()


n <- length(final[[j]]$residuals)
t<-seq(-5,n+6,by=1)
# res_indice
res_indice<-paste0("Plots/res_indice",i,".pdf")
pdf(res_indice,width = w1, height = h11)
# resid vs index BXII residuals 
#par(mar=c(5,6,4,1)+.1)
plot(final[[j]]$residuals,xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
lines(t,rep(-3,n+12),lty=2,col=1)
lines(t,rep(3,n+12),lty=2,col=1)
lines(t,rep(-2,n+12),lty=3,col=1)
lines(t,rep(2,n+12),lty=3,col=1)
dev.off()

# Acf plot
Acf_plot<-paste0("Plots/Acf_plot",i,".pdf")
pdf(Acf_plot,width =4, height = 4)
# acf  
par(mfrow=c(1,1))
acf(y,main=nome[j]#,main=""
)
dev.off()

# 1,4,5,11,12,13,16
################################################################################
h=12
s=12
j=i=16
y_prev<-get(paste0("y_prev",j))
y_prev<-ts(y_prev,frequency = s,start=c(2023,9))
forecast<-paste0("Plots/forecast",i,".pdf")
pdf(forecast,width = w1, height = h11)
plot(y_prev,ylab="RF",ylim=c(min(y_prev,out_forecast[,1,j]),max(y_prev,out_forecast[,1,j])+0.15),type = "l")
lines(ts(out_forecast[,1,j],frequency = s,start=c(2023,9)),col=2,lty=2,lwd=2)
# lines(ts(out_forecast[,2,j],frequency = s,start=c(2023,9)),col=4,lty=4,lwd=1.6)
legend("topleft", 
       c("BXII-GAS"#,"Ray-GAS"
         ),
       col = c(2#,4
               ),
       lty= c(2#,4
              ),
       lwd = c(2#,1.6
               ), bty="n", cex = 1.5)
dev.off()

# 1,2,11,12,16
########################################################

j=5
# data<-final[[j]]$model[,c(1,2)]
# data<-final_Ray[[j]]$model[,c(1,2)]
#data<-final_GAMMA[[j]]$model[,c(1,2)]
#data<-finalBXIIARMA[[j]]$model[,c(1,2)]
#data<-final_ARMA[[j]]
#data<-final_SARMA[[j]]
data<-final_ETS[[j]]
data
data<-round(data,3)

# Função para formatar valores
format_value <- function(est, se) {
  if (est < 0) {
    sprintf("$%.3f$ (%.3f)", est, se) # Para valores negativos
  } else {
    sprintf("%.3f (%.3f)", est, se)  # Para valores positivos
  }
}

# Aplicar formatação a cada linha
formatted <- apply(data, 1, function(row) format_value(row[1], row[2]))

# Criar dataframe formatado
formatted_output <- data.frame(Parameter = rownames(data), Estimate_SE = formatted)
print(formatted_output, row.names = FALSE)
formatted_output<-formatted_output[,2]
