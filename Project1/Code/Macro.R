missing_row<-function(y,x){
  
  miss.all<-paste0(length(y[is.na(y)==T])," (",round((length(y[is.na(y)==T])/length(y))*100),'%',")")
  
  miss.groups<-tapply(y,x, function(x) paste0(length(x[is.na(x)==T])," (",round((length(x[is.na(x)==T])/length(x))*100),'%',")"))
  
  m<-matrix(0,ncol=length(miss.groups)+2,nrow=1)
  
  #put row together
  t<-data.frame(rbind(m,c(Variable="",All=miss.all,miss.groups)))
  
  return(t[-1,])
}

missing_table<-function(data,variables,group,col.names=T){
  
  #create table
  temp<-data.frame(do.call(rbind,lapply(data[variables],missing_row,group)),row.names=NULL)
  temp$Variable<-as.character(temp$Variable)
  temp$Variable<-Hmisc::label(data[variables])
  
  #set column names with sample sizes
  if (col.names==T){
    colnames(temp)<-column_label(group)[-length(column_label(group))]
  }
  
  return(temp)
}


easy_p<-function(out,pred){
  p<-plot(out ~ pred, data = data_final, main = paste(label(out), "vs.", label(pred)),xlab=paste(label(pred)),ylab=paste(label(out)))
  if (is.numeric(pred)==TRUE){
  y <- predict(lm(out ~ pred, data = data_final))
  x <- pred[is.na(out)==FALSE]
  l<-lines(y = y, x = x, col = "navy")
  }
}