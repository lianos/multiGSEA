## This was supplemental file 1 from the Q-Gen publication:
## http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0707-9
## http://static-content.springer.com/esm/art%3A10.1186%2Fs12859-015-0707-9/MediaObjects/12859_2015_707_MOESM1_ESM.txt
library(qusage)

#Q-Gen function

qusage_gen<-function(resids,labels,estimates,dof,std.errors,geneSets,var.equal=TRUE){

if(var.equal){labels<-rep("Resids",ncol(resids))}
if(nrow(resids)!=length(estimates)){return("Error: Number of rows in residual matrix do not equal length of estimate vectors")}
if(nrow(resids)!=length(dof)){return("Error: Number of rows in residual matrix do not equal length of dof vectors")}
if(nrow(resids)!=length(std.errors)){return("Error: Number of rows in residual matrix do not equal length of std,errors vectors")}

names(estimates)<-rownames(resids)
names(dof)<-rownames(resids)
names(std.errors)<-rownames(resids)

qlist<-list(mean=estimates,SD=std.errors,dof=dof,labels=labels)
results<-newQSarray(qlist)

cat("Aggregating gene data for gene sets.")
results<-aggregateGeneSet(results, geneSets,n.points=2^14)
cat("Done. \nCalculating VIF's on residual matrix.")
results<-calcVIF(resids,results,useCAMERA=FALSE)
cat("\nQ-Gen analysis complete.")
results
}





######################
#Modification of Original Qusage Example to illustrate Q-Gen
######################


##create example data - a set of 500 genes normally distributed across 20 patients
  eset = matrix(rnorm(500*20),500,20, dimnames=list(1:500,1:20))
  labels = c(rep("A",10),rep("B",10))

##create a number of gene sets with varying levels of differential expression.
  geneSets = list()
  for(i in 0:10){
    genes = ((30*i)+1):(30*(i+1))
    eset[genes,labels=="B"] = eset[genes,labels=="B"] + rnorm(1)

    geneSets[[paste("Set",i)]] = genes
  }


#Fitting a simple linear model for B vs A to get residual matrix
#Note these can come from conducting the models manually in R or they can simply come
#from result files from analyses conducted in other statistical software or packages

resids<-  t(apply(eset,1,function(x){lm(x~labels)$residuals }))
rownames(resids)<-1:500
colnames(resids)<-labels

#Calculating estimates,standarderrors,and dof
difference<-apply(eset,1,function(x){result=t.test(x[11:20],x[1:10],var.equal=T)$estimate
                                     return(result[1]-result[2])})
deg.of.freedom<-apply(eset,1,function(x){t.test(x[11:20],x[1:10],var.equal=T)$parameter})
standarderror<-difference/apply(eset,1,function(x){t.test(x[11:20],x[1:10],var.equal=T)$statistic})

##Generalized Qusage Analysis derived from any linear mixed model
##Results object is the standard qusage object and all qusage functions can be applied directly


results = qusage_gen(resids,labels,difference,deg.of.freedom,standarderror,geneSets,T)
plot(results)










