setwd("C:\\Users\\Shawn\\Documents\\GitHub\\Scripts\\RNAseq\\02DEseq")
file_all = c("Tree_design.txt") ##Figure this out

data_dir = "C:\\Users\\Shawn\\Documents\\GitHub\\Scripts\\RNAseq\\02DEseq"
gtf_file = "PdeltoidesWV94_445.gtf" ##Figure this out

library(Biostrings)
library(IRanges)
library(GenomicRanges)
library(DESeq2)

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





dire = data_dir
design_file = file_all[1]
sampleTable = read.table(design_file, head=TRUE)

DESeqDataSetFromHTSeqCount_lxue = function (sampleTable, directory = "", design, ...)
{
	if (missing(design)) {
		stop("design is missing")
	}
	l <- lapply(as.character(sampleTable[, 2]), function(fn) read.table(file.path(directory,
								fn)))
	if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1))))
		stop("Gene IDs (first column) differ between files.")
	tbl <- sapply(l, function(a) a$V2)
	rownames(tbl) <- l[[1]]$V1
	rownames(sampleTable) <- sampleTable[, 1]
	specialRows <- rownames(tbl) %in% c("no_feature", "ambiguous","too_low_aQual", "not_aligned", "alignment_not_unique",
			"__no_feature", "__ambiguous","__too_low_aQual", "__not_aligned", "__alignment_not_unique")
	tbl <- tbl[!specialRows, ]
	dds <- DESeqDataSetFromMatrix(countData = tbl, colData = sampleTable[,
					-(1:2), drop = FALSE], design = design, ...)
	return(dds)
}


xxx=strsplit(design_file,"\\_")[[1]]
tissue = xxx[length(xxx)]

tissue = "pollen"
#########################################################
## Multiple factors

dds_raw = DESeqDataSetFromHTSeqCount_lxue(sampleTable,directory=dire, design= ~ 1)
dds = dds_raw

rowRanges(dds) = gene_range

fpkm_out = fpkm(dds)

file_out = paste("FPKM_",tissue,".txt",sep="")
write.table(as.data.frame(fpkm_out),	file=file_out,sep="\t")


