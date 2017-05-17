
#update GO.db
source("http://www.bioconductor.org/biocLite.R")
#biocLite("GO.db")
library(topGO)

GO_space = c("BP","MF","CC")

go_group<-GO_space[3]
#BP MF  CC

#read map
#geneID2GO <- readMappings(file = "D:\\Database\\topGO\\Ptri\\Ptr_GO_GOSLIM.map")
geneID2GO <- readMappings(file = "F:\\Projects\\Wang\\Trees\\10Gene_Function\\topGO\\PopEuph_topGo.map")

#setwd("D:\\Dropbox\\RNAseq_TUA_YY\\201507_tension\\ver3\\GO\\P0.01FC1.5")

setwd("F:\\Projects\\Wang\\Trees\\02DEseq\\GO\\list")


str(head(geneID2GO))

geneNames <- names(geneID2GO)
head(geneNames)




mygo<-function(gene_list){
	
#file<-"Glyma01g25690_0.85.txt"
#file<-"up_NT_WT_FD10.txt"
	file<-gene_list
			
#read select genes
read.table(file,header=FALSE,sep="\t")->data
gene<-as.character(data[[1]])




myInterestingGenes <- gene
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)


BPterms <- ls(GOBPTerm)
#head(BPterms)

	                           # MF BP CC
GOdata <- new("topGOdata", ontology = go_group, allGenes = geneList,nodeSize=3, annot = annFUN.gene2GO,
		 gene2GO = geneID2GO)

#graph(GOdata)


resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")

#resultFisher <- getSigGroups(GOdata, test.stat)
#test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
	num<-length(usedGO(GOdata))
allRes <- GenTable(GOdata, classic = resultFis,
		 ranksOf = "classic", topNodes = num)
allRes

file_out= paste("out_",go_group,"_",file,sep="")
write.table(allRes,file=file_out,sep="\t",quote = FALSE,row.names = FALSE,col.names =TRUE)
#dev.new()
#showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 20, useInfo = "all")
#dev.off()
}


my_list<-list.files(path =getwd(), pattern =".txt$", all.files = TRUE,
		full.names = FALSE, recursive = FALSE)

my_list2 = my_list

for(fff in my_list2){    
	#mygo(fff,"func")
	#mygo(fff,"comp")
	print (fff)
	mygo(fff)
}






