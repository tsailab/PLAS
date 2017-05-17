
workingDir<-"E:\\Project_2016\\20160712RNAseq\\01DEseq\\01FPKM"
file<-"FPKM2_YY_T2016.txt"
file_out<-"FPKM_YY_tension2016_out.csv"


### In the input file, the first column is the gene ID,
### all other columns are FPKM for each gene in samples

######
## set working dir
setwd(workingDir);
## read table
read.table(file,header=TRUE,sep="\t")->affydata;
## remove the first column
#data<-affydata[,-c(1)]


affydat_num = affydata
affydat_num[affydat_num<3] = 0
affydat_num[affydat_num>=3] = 1

sum_num = apply(affydat_num,1,sum)

affydata_sel = affydata[sum_num>=6,]



data = affydata_sel 

names<- colnames(data)

### run PCC
x <- data
y <- data
bl <- lapply(x, function(u){
			lapply(y, function(v){
						cor(u,v,method="pearson") # Function with column from x and column from y as inputs
					})
		})
out = matrix(unlist(bl), ncol=ncol(y), byrow=T)

### change names
rownames(out)= names
colnames(out)= names

write.csv(out,file_out)








