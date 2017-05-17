setwd("F:\\Projects\\Wang\\Trees\\02DEseq")
file_all = c("Tree_design.txt")

data_dir = "F:\\Projects\\Wang\\Trees\\02DEseq\\data"
gtf_file = "GCF_000495115.1_PopEup_1.0_genomic.gtf"


library(DESeq2)

design_file = file_all 
sampleTable0 = read.table(design_file, head=TRUE)

file_pair = "Tree_compare.txt"
compare = read.table(file_pair , head=TRUE)




gtf2GRangesList <- function(myfile="my.gff") {
	gtf <- read.delim(myfile, header=FALSE)
	colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame",      
			"attributes")
	#chronly <- c(1:22, "X", "Y", "MT")
	#gtf <- gtf[as.character(gtf$seqname) %in% chronly, ] # Cleanup to remove non-chromosome rows
	
	gene_ids <-  gsub(".*gene_id (.*?);.*", "\\1", gtf$attributes) #get gene_id from attributes column
	transcript_ids<- gsub(".*transcript_id (.*?);.*", "\\1", gtf$attributes) #get transcript_id from attributes column
	
	index<-gene_ids!="" #skip those have no value
	index2<-transcript_ids!=""
	
	gene.gr<-GRanges(seqnames=gtf$seqname[index],
			ranges=IRanges(gtf$start[index],gtf$end[index]),
			strand=gtf$strand[index],
			tx_id=transcript_ids[index],
			gene_id=gene_ids[index])
	gene.gr.list<-split(gene.gr,gene_ids[index])
	
	transcript.gr<-GRanges(seqnames=gtf$seqname[index2],
			ranges=IRanges(gtf$start[index2],gtf$end[index2]),
			strand=gtf$strand[index2],
			tx_id=transcript_ids[index2],
			gene_id=gene_ids[index2])
	transcript.gr.list<-split(transcript.gr,transcript_ids[index2])
	
	r<-list()
	gene<-gene.gr.list
	#r$transcript<-transcript.gr.list
	return(gene)
}



gene_range = gtf2GRangesList(gtf_file)


total_num = dim(compare)[1]
for (index_num in c(1:total_num)){
	test = as.character(compare[index_num,1])
	print(test)
	file_out = paste(test,"_out.txt",sep="")
	g1 = as.character(compare[index_num,2])
	g2 = as.character(compare[index_num,3])
	
	sampleTable1 = sampleTable0[sampleTable0$Sample == g1,]
	sampleTable2 = sampleTable0[sampleTable0$Sample == g2,]
	
	sampleTable = rbind(sampleTable1,sampleTable2)
	
	################################################################
	## Single factors
	
	sampleTable_sel = sampleTable
	
	
	ddsSF = DESeqDataSetFromHTSeqCount(sampleTable_sel, directory=data_dir,design= ~ Sample)
	design(ddsSF) <- formula(~ Sample)
	ddsSF <- DESeq(ddsSF)
	resultsNames(ddsSF)
	
	resSFtreatment <- results(ddsSF, cooksCutoff=FALSE, contrast=c("Sample",g2,g1))
	
	
	## filter basing on FPKM
	design_in = "Sample"
	## FPKM
	
	rowRanges(ddsSF) = 	gene_range 
	fpkm_out = fpkm(ddsSF)
	
	## filter basing on FPKM
	
	group_summary  = as.character(sampleTable_sel[[design_in]])
	sample_summary = as.character(sampleTable_sel$sample_id)
	group_summary_uni = unique(group_summary)
	
	group_num = length(group_summary_uni)
	
	if(group_num<2){
		cat ("Check the desing region")
	}
	
	fpkm_filter_index  = 0	
	for(group_test in group_summary_uni ){
		fpkm_filter_index = fpkm_filter_index+1
		sample_test = sample_summary[group_summary==group_test]
		fpkm_test = fpkm_out[, sample_test]
		head(fpkm_test)
		fpkm_test_0 = fpkm_test
		fpkm_test_0 [fpkm_test_0 >=1]   = 1
		fpkm_test_0 [fpkm_test_0 <1 ]   = 0
		sample_test_num = dim(fpkm_test_0)[2]
		fpkm_test_0_sum = apply(fpkm_test_0,1,sum)
		fpkm_test_0_sum[fpkm_test_0_sum < sample_test_num] = 0
		fpkm_test_0_sum[fpkm_test_0_sum >0 ] = 1
		if(fpkm_filter_index ==1){
			fpmk_checking_out = fpkm_test_0_sum
		}else{
			fpmk_checking_out = fpmk_checking_out + fpkm_test_0_sum
		}
	}
	
	out = as.data.frame(resSFtreatment)
	out_with_fpkm = cbind(fpkm_out,out)
	
	not_na = ! ( is.na(out$pvalue))
	out_sel=out_with_fpkm[not_na & fpmk_checking_out > 0,]	

	write.table(out_sel,file=file_out,sep="\t")
	

}


