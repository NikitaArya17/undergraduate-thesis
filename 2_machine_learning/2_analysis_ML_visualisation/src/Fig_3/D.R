source("2_analysis_ML_visualisation/src/setpath.R")
library(dplyr)

d1= read.csv(file.path(data_path,"Database.csv"),header=TRUE,stringsAsFactors = FALSE) 
df = d1[d1[,1]=="",-2] #rewritten more simply - original code threw an error
gene_start = grep("rrrD",colnames(df))
stress_col <- grep("Stress", colnames(df)) #line added to ensure this column is not excluded from subsetting
p = dim(df)[2]
mut = df[, c(stress_col, gene_start:p)]
mut[, -1][mut[, -1] != mut[1, 2] & mut[, -1] != ""] = "Hit"
id_mut =data.frame(df[,1],mut)

n <- colnames(mut)
cols = c(3,4,5,6,7,8,9,10,11,13,18)
df$label_col =  apply( df[ , cols ] , 1 , paste , collapse = "-" )
mut_group = split(df,df[,"label_col"])
total = length(mut_group)
nums_stress = as.numeric()
nums = as.numeric()

for ( i in 1:total)
{
  tmp = mut_group[[i]]
  cols = c(13,15,16,17)
  
  nums_stress=c(nums_stress, length(unique(apply( tmp[ , cols ] , 1 , paste , collapse = "-" ))))
  
  nums=c(nums,nrow(tmp))
  
}

ncol(df)
anti = mut_group[[18]][,c(13,21:(ncol(df)-1))] #corrected logical flaw in subsetting
write.csv(mut_group[[18]],"2_analysis_ML_visualisation/output_data/anti_group.csv",row.names = FALSE)

check = function(x)
{
  num_hit = length(which(x!=""))
  if (num_hit/length(x)!=0)
    return(1)
  else
    return(0)
}

n <- colnames(anti)[-1]
f <- as.formula(paste("cbind(",paste(n, collapse = " ,"),")","~Stress",sep=""))

categories=read.csv("2_analysis_ML_visualisation/ref/anti_category.csv",header=FALSE,stringsAsFactors = FALSE)
mut <- mut %>%
  filter(Stress %in% categories[,1])

#THIS STEP WAS CORRECTED:
################
anti_each <- aggregate(.~Stress, data=mut, FUN =check)
labels <- anti_each[, 1]
num_anti_each <- anti_each[, -1]
res <- scale(num_anti_each[, -13])
row.names(res) <- labels
d <- dist(res)
#################

methods = c("average","ward.D2","single","complete","median","centroid")

colors = c("#6ED9B0","#FA6c6C","#42135F","#FACB18","#7D4FFE")

types = c("DNA,protein and cell wall synthesis","Protein synthesis",
          "Folic acid synthesis","Cell wall synthesis","DNA synthesis")


assignColor <- function(label)
{
  index = which(label==categories[,1])
  return(colors[categories[index,2]])
}

## function to set label color
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label") 
    ## set label color to red for A and B, to blue otherwise
    attr(x, "nodePar") <-  list(lab.col=assignColor(label),pch = 19,col="white")
  }
  return(x)
}

for (method in methods)
{
  hc = hclust(d,method=method)
  
  
  
  ## apply labelCol on all nodes of the dendrogram
  de <- dendrapply(as.dendrogram(hc), labelCol)
  
  
  pdf(paste("2_analysis_ML_visualisation/output_figures/binary/",method,"_binary.pdf",sep=""), 7,6);
  par(mar=c(4,2, 2,7)+1)
  
  
  plot(de,horiz=TRUE,xlab="distance")#,main=method)
  
  
  legend("topleft", legend = types, fill = colors, title = "", box.col = "transparent")
  
  dev.off()
}