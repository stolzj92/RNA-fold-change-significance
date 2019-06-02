library(combinat)
#Table2.csv is the normalized,floored and filtered hypoxia values merged into one table. 
# this file has already been floored using Table2[Table2<20]<-20
Table<-read.table(file="Table2.csv",sep = ",")
#The following comments show how the table was made.
#Eliminate rows with zero expression (Table is the original raw values)
#gene.sum<-apply(Table,1,sum)
#Table<-Table[gene.sum>0,]
#Normalize the data to the smallest column total
#TotalExpression<-colSums(Table)
#minExpression<-min(TotalExpression)
#Table<-sweep(Table,2,TotalExpression,'/')
#Table2<-Table*minExpression
#floor the data
#Table2[Table2<20]<-20
#filter out rows that have all 20 as expression
row.names(Table)<-Table[,1]
Table[,1]<-NULL
gene.sum<-apply(Table,1,sum)
Table2<-Table[gene.sum>160,]

times<-c(0,5,10,30,60,120,180,240)
Output<-rep(1,nrow(Table2))
for (i in 1:nrow(Table2)){ 
  f<-permn(as.numeric(Table2[i,]))
  R_2prime<-summary(lm((f[[1]])~poly(times,2)))$adj.r.squared
  total<-1
  print(i)
  for (n in 1:length(f)){
    r.squared<-summary(lm((f[[n]])~poly(times,2)))$adj.r.squared
    if ((R_2prime)>0){
      if (R_2prime<=r.squared){
        total<-total+1}
    }
    else{
      if (R_2prime>=r.squared){
        total<-total+1}
    }
  }
  Output[i]<-((total)/40320)
}
Table2$p_value<-Output
fdr<-p.adjust(Output,"fdr",n=length(Output))
Table2$adjusted_p<-fdr
Table3<-Table2[order(fdr),]
plot(times,log2(Table3[4,1:8]),ylab= "YGR131W Log2 Expression")
plot(times,log2(Table3[1,1:8]),ylab= "YBR042C Log2 Expression")
plot(times,log2(Table3[2,1:8]),ylab= "YDL243C Log2 Expression")
plot(times,log2(Table3[3,1:8]),ylab= "YER019W Log2 Expression")
hist(Table3$adjusted_p,xlab = "adjusted p",main= "Distribution of adjusted p-value", breaks = 50)
