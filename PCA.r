###############################################################
# Entra en la carpeta EVENTOS y 
# guarda los eventos centrados para utilizar en rr y en este cod
# Luego, realiza pca a cada evento en C02 y C05. 
# Calcula el rankin segun la influencia de cada PCA
################################################################

library(dplyr)

##########################Scale the events and save it ####################################################
scale_data <- function(file,cont_files){
  feature = c("Fecha",	"Value (%)",	"Value (hPa)",	"Value (W/m2)",	"Value (°C)", "Value (mm)",  "Monthly sin",	"Monthly cos",	"Bi-monthly sin", "Bi-monthly cos", "Wind X",	"Wind Y","Label")
  #Read each File and make dataset
  Dataset<-read.csv(file,col.names =feature)
  
# Fecha	Value (%)	Value (hPa)	Value (W/m2)	Value (°C) Value (mm)	
#  Monthly sin	Monthly cos	Bi-monthly sin  Bi-monthly cos	
#  Wind X	Wind Y	Label

  #CENTER DATA FOR RR
  path<-paste("/home/juan/Desktop/TESIS/Codes/codesTesis/Paso3/CenteredEvents/centered",cont_files,".csv")
  
  fecha<-filter(Dataset, )[1:24,1:1]
  data<-filter(Dataset, )[1:24,2:12]

  rr_data<-scale(data,center = TRUE)
  rr_data[is.na(rr_data)]<-0

  rr_centered_data<-cbind(fecha,rr_data)
  
  colnames(rr_centered_data)<-c("Fecha",	"Value (%)",	"Value (hPa)",	"Value (W/m2)",	"Value (°C)", "Value (mm)",  "Monthly sin",	"Monthly cos",	"Bi-monthly sin", "Bi-monthly cos", "Wind X",	"Wind Y")
  write.csv(rr_centered_data,path,row.names = FALSE)
  
  #CENTER DATA FOR PCA

  parte1<-filter(Dataset, )[1:24,2:5]
  #mm<- filter(Dataset, )[1:24,6:6]
  parte2<- filter(Dataset, )[1:24,7:12]
  #labels<- filter(Dataset, )[1:24,13:13]
  
  X<- cbind(parte1,parte2)
  colnames(X)<-c("Value (%)", "Value (hPa)","Value (W/m2)", "Value (°C)", "MonthlY sin","MonthLY cos","Bi-monthly sin","Bi-monthly cos","Wind X","Wind Y")
  
  W<-scale(X,center = TRUE)

  W[is.na(W)]<-0
  which(apply(W, 2, var)==0)
  W<-W[ , which(apply(W, 2, var) != 0)]
  
  return(W)
}
##############################################################################
# si cualquirea de los tres primeros valor en max_pc_sv
# rank_feature_sum es igual al del feature[i]
# entonce cont[i]++
rank_according_to_sum<-function(rank_feature_sum,variance,cont_rank_sum){
  max_pc_sv <- sort(rank_feature_sum, decreasing = TRUE)
  for(i in 1:length(variance)){
    if( isTRUE(max_pc_sv[1] == rank_feature_sum[i]) ) {
      cont_rank_sum[i]<-cont_rank_sum[i]+1
    }
    if(isTRUE(max_pc_sv[2] == rank_feature_sum[i])){
      cont_rank_sum[i]<-cont_rank_sum[i]+1
    }
    if(isTRUE(max_pc_sv[3] == rank_feature_sum[i]) ){
      cont_rank_sum[i]<-cont_rank_sum[i]+1
    }
  }
  return(cont_rank_sum)
}

#####################################################################3
rank_according_to_max<-function(rank_feature_max,variance,cont_rank_max){
  max_pc_sv <- sort(rank_feature_max, decreasing = TRUE)
  
  for(i in 1:length(variance)){
    if( isTRUE(max_pc_sv[1] == rank_feature_max[i]) ) {
      cont_rank_max[i]<-cont_rank_max[i]+1
    }
    if(isTRUE(max_pc_sv[2] == rank_feature_max[i])){
      cont_rank_max[i]<-cont_rank_max[i]+1
    }
    if(isTRUE(max_pc_sv[3] == rank_feature_max[i]) ){
      cont_rank_max[i]<-cont_rank_max[i]+1
    }
  }
  
  return(cont_rank_max)
  
}

##############################################################################
# leemos la data
files_Rumihurco <- list.files(path="/home/juan/Desktop/TESIS/Codes/codesTesis/Paso2/Rumihurco", pattern="*.csv", full.names=TRUE, recursive=FALSE);
files_Rumipamba <- list.files(path="/home/juan/Desktop/TESIS/Codes/codesTesis/Paso2/Rumipamba", pattern="*.csv", full.names=TRUE, recursive=FALSE);


##############################################################################
#combine them in one list
files<- append(files_Rumipamba,files_Rumihurco);

##############################################################################
#make variables for the pc ranks

feature <- c("Value (%)", "Value (hPa)","Value (W/m2)", "Value (°C)","Value (mm)","Monthly sin","Monthly cos",	"Bi-monthly sin", "Bi-monthly cos","Wind X","Wind Y")

cont_rank_max<-c(0,0,0,0,0,0,0,0,0,0)
cont_rank_sum<-c(0,0,0,0,0,0,0,0,0,0)

##############################################################################

cont_files=0
for (file in files) {
  
  cont_files=cont_files+1
  feature<-c("Value (%)", "Value (hPa)","Value (W/m2)", "Value (°C)", "MonthLY sin","MonthLY cos","Bi-monthly sin","Bi-monthly cos","Wind X","Wind Y")
  
  # Scale the data to make PCA and save it to use later in relieff
  W<- scale_data(file,cont_files)
  
  #make PCA 
  pcy<-prcomp(W)
  pc<-pcy$rotation
  #print("PCs")
  #print(pc)
  new_matrix<-pcy$x
  
  ### important results
  importance<-summary(pcy)$importance
  variance <- (importance[1,])^2
  Proportion_of_Variance<-importance[2,]*100
  Cumulative_Proportion<-importance[3,]*100
  
  
  ######################## make ranks ############################
  
  # CALCULAMOS LA RELEVANCIA
  # relevancia de variable[i] = sum(cada PC de la variable de variable "xyz" * variance)
  # relevancia de variable[i] = prod vectoria del pc[i] y la variancia
  
  relevancia<-c(0)

  for (i in 1:length(variance)) {
    relevancia[i]=abs(pc[i,])%*%variance
  }
  
  # CALCULAMOS EL RANK WRT LOS OTROS
  # para captar su rank usar: relevancia[i]/ sum(relevancia)
  
  rank_feature_sum<-c(0)
  for (i in 1:length(variance)) {
    rank_feature_sum[i]<-(relevancia[i]/sum(relevancia))
  }
  
  # para captar su rank usar: relvancia/ max(relevancia)
  rank_feature_max<-c(0)
  for (i in 1:length(variance)) {
    rank_feature_max[i]<-(relevancia[i]/max(relevancia))
  }
  
  
  #IMPRIMIMOS RESULTADOS PARA ESETE EVENTO
  rank_m<-matrix(nrow = 3, ncol = length(variance))
  
  #make a matrix to print values for each rank
  for (i in 1:length(variance)) {
    rank_m[1,i]<-feature[i]
    rank_m[2,i]<-round(rank_feature_sum[i], digits = 5)*100
    rank_m[3,i]<-round(rank_feature_max[i],digits = 5)*100
  }  
  #print("pcy")
  #print(pcy)
  #print("importance")
  #print(round(importance,digits = 5)*100)
  
  print("Rank Matrix: Analisis de relevancia con PCA")
  print(t(rank_m))
  cat(paste("\nfile: ", file), sep = "") 

  
  ######################## make counter  ############################
  
  print("Rank acording to relvancia/ sum(relevancia)")
  cont_rank_sum <- rank_according_to_sum(rank_feature_sum,variance,cont_rank_sum)
  
  # rank accroding to max ields same results
  
}

for (i in 1:length(variance)) {
  cat(paste("\nFeature ", feature[i], "count=",cont_rank_sum[i]), sep = "") 
}






