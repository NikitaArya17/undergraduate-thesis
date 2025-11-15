#Note: The code to create the output file "NB_data.csv" is missing.
#Therefore, this file has been copied from the authors' folder
#in order to continue running the code.

#In addition, the file has been misnamed here and was corrected:
df=read.csv("2_analysis_ML_visualisation/output_data/NB_data_REPRODUCTION.csv",fileEncoding = "latin1", header=TRUE) 

cul_f_last = tail(grep("Pert_",colnames(df)),n=1)
p = dim(df)[2]

mut = df[,(cul_f_last+2):p]

p1=dim(mut)[2]
n= dim(mut)[1]

PAB= matrix(nrow=p1,ncol=p1)

i=1

PA = as.numeric()

library(doParallel)

detectCores()
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
getDoParWorkers()

library(foreach)

ptm = proc.time()

for(i in 1:p1)
{

rs <- foreach(j = 1:p1, .combine = rbind) %dopar% {
  a = which(mut[,i]=="Presence")
  b= which(mut[,j]=="Presence")
  c=length(intersect(a,b))
  return(c(length(a)/n, c/n))
}
PA[i]= rs[1,1]
PAB[i,]=rs[,2]
}

time = proc.time() - ptm 

out=data.frame(PA,PAB)

stopCluster(cl) #An important coding practice was omitted.

n= length(PA)
dep = matrix(nrow=n,ncol=n)

i=1
for (i in 1:n)
{
j=1

for (j in 1:n)
{
dep[i,j] = PAB[i,j]/(PA[i]*PA[j])
}

}

diag(dep)=0

write.csv(dep,"2_analysis_ML_visualisation/output_data/dependence_repro.csv",row.names=FALSE)


