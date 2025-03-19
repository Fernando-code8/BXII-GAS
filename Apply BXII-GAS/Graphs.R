w1<-4 # width for plots 
h11<-4 # height for plots

setwd("C:/Users/fernando.monteiro/Downloads/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001-20250228T054140Z-001/UG-GAS e BXII-GAS 20241120-20241217T004454Z-001/UG-GAS e BXII-GAS 20241120/Scripts Atual BXII-GAS/Apply BXII-GAS")
j=i=16 # 1,4,5,11,12,13,16
resid_plot<-paste0("Plots/Acf",i,".pdf")
pdf(resid_plot,width = w1, height = h11)
# acf BXII residuals 
acf(final[[j]]$residuals,main="")
dev.off()

# 1,4,5,11,12,13,16
################################################################################
h=12
s=12
j=i=16
y_prev<-get(paste0("y_prev",j))
y_prev<-ts(y_prev,frequency = s,start=c(2023,9))
forecast<-paste0("Plots/forecast",i,".pdf")
pdf(forecast,width = w1*2, height = h11)
plot(y_prev,ylab="RF",ylim=c(min(y_prev,out_forecast[,1,j]),0.15+max(y_prev,out_forecast[,1,j])),type = "l")
lines(ts(out_forecast[,1,j],frequency = s,start=c(2023,9)),col=2,lty=2,lwd=2)
legend("topright", 
       c("BXII-GAS"),
       col = c(2),
       lty= c(2),
       lwd = c(2), bty="n", cex = 1)
dev.off()

# 1,2,11,12,16
########################################################

j=16
# data<-final[[j]]$model[,c(1,2)]
# data<-final_Ray[[j]]$model[,c(1,2)]
data<-final_GAMMA[[j]]$model[,c(1,2)]
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
