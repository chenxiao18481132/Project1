library(caret)
library("rlist")
library(aricode)
set.seed(3456)

expression<-read.table("gene_expression.txt")
label<-read.table("labels.txt")
trainIndex <- createDataPartition(label$V1, p = .8, 
                                  list = FALSE, 
                                  times = 1)
head(trainIndex)
expressiontrain<-expression[trainIndex,]
expressiontest<-expression[-trainIndex,]
Y<-as.matrix(expression)
m<-fitmodel(Y,10,iter=100)
kmfull<-kmeans(Y,7,nstart=40)   
km<-kmeans(m$EZ,7,nstart = 40)   
list.save(m, 'list.json')
m=list.load('list.json')
clustComp(km$cluster,label$V1)
clustComp(kmfull$cluster,label$V1)
tablecluster=list()
for(i in 1:7){
tablecluster[[i]]=table(km$cluster[label$V1==i-1])
#sum_right=sum_right+max(tablecluster[[i]]
}
