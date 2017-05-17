### change into current files
working_dir = "F:\\Projects\\Wang\\Trees\\02DEseq\\FPKM"

file_in =  "FPKM_tree2.txt"


##### 
# The format of input file:
#  1st column gene names
#  other columns expression data of genes in each sample


#### set working dir
setwd(working_dir)

### read input file
read.table(file_in,header=TRUE,row.names=1,sep="\t")->affydata;  

affydat_num = affydata
affydat_num[affydat_num<2] = 0
affydat_num[affydat_num>=2] = 1

sum_num = apply(affydat_num,1,sum)

affydata_sel = affydata[sum_num>=3,]

### remove the first columns
#data = as.matrix(affydata[,-c(1)])
data = as.matrix(affydata_sel)
data_log = log2(data+1)
### box plot of data
boxplot(data_log)

### tree plot
plot(hclust(dist(t(data_log))))


#######################
## Clustering using Pearson
hclust2 <- function(x, method="average", ...)
	hclust(x, method=method, ...)
dist2 <- function(x, ...)
	as.dist(1-cor(t(x), method="pearson"))

plot(hclust2(dist2(t(data_log ))))

