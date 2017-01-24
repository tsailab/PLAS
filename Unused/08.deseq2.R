###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################

rm(list=ls())
gc()

###########################################################################################
###  Load necessary R libraries for later usage.										###
###########################################################################################

cat("\n\n########## Loading Libraries ##########\n", sep="")
Biostrings = library("Biostrings", logical.return = TRUE)
IRanges = library("IRanges", logical.return = TRUE)
GenomicRanges = library("GenomicRanges", logical.return = TRUE)
DESeq2 = library("DESeq2", logical.return = TRUE)

if(!Biostrings){
	source("http://bioconductor.org/biocLite.R")
	biocLite("Biostrings")
	library("Biostrings")
}
if(!IRanges){
	source("http://bioconductor.org/biocLite.R")
	biocLite("IRanges")
	library("IRanges")
}
if(!GenomicRanges){
	source("http://bioconductor.org/biocLite.R")
	biocLite("GenomicRanges")
	library("GenomicRanges")
}
if(!DESeq2){
	source("http://bioconductor.org/biocLite.R")
	biocLite("DESeq2")
	library("DESeq2")
}

args = commandArgs(TRUE)
gtf_file = args[1] #"../00.genome/Pta717s_v1.1_merge.gtf"
file_all = args[2] #"../02.fastq/01.data/tissue_record.txt"
dire = args[3] #"../11.htseq/01.data/02.uniq.hit"
tgtfolder = args[4]

gtf2GRangesList <- function(myfile="my.gff", skip=0, format="gtf") {
	gtf <- read.delim(myfile, header=FALSE, skip=skip, as.is=TRUE)
	colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
	#chronly <- c(1:22, "X", "Y", "MT")
	#gtf <- gtf[as.character(gtf$seqname) %in% chronly, ] # Cleanup to remove non-chromosome rows
	
	if(format == "gtf"){
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
	}else if(format == "gff3"){
		transcript_ids<- gsub(".*Parent=(.*?);.*", "\\1", gtf$attributes[gtf$feature!="gene" & gtf$feature!="mRNA"]) #get transcript_id from attributes column
		transcript2gene = cbind(gsub(".*ID=(.*?);.*", "\\1", gtf$attributes[gtf$feature=="mRNA"]),
								gsub(".*Parent=(.*?)", "\\1", gtf$attributes[gtf$feature=="mRNA"]))
		loc = match(transcript_ids, transcript2gene[,1])
		gene_ids = transcript2gene[loc,2]
		index<-gene_ids!="" #skip those have no value
		index2<-transcript_ids!=""
		
		gene.gr<-GRanges(seqnames=gtf$seqname[gtf$feature!="gene" & gtf$feature!="mRNA"][index],
				ranges=IRanges(gtf$start[gtf$feature!="gene" & gtf$feature!="mRNA"][index],gtf$end[gtf$feature!="gene" & gtf$feature!="mRNA"][index]),
				strand=gtf$strand[gtf$feature!="gene" & gtf$feature!="mRNA"][index],
				tx_id=transcript_ids[index],
				gene_id=gene_ids[index])
		gene.gr.list<-split(gene.gr,gene_ids[index])
		
		transcript.gr<-GRanges(seqnames=gtf$seqname[gtf$feature!="gene" & gtf$feature!="mRNA"][index2],
				ranges=IRanges(gtf$start[gtf$feature!="gene" & gtf$feature!="mRNA"][index2],gtf$end[gtf$feature!="gene" & gtf$feature!="mRNA"][index2]),
				strand=gtf$strand[gtf$feature!="gene" & gtf$feature!="mRNA"][index2],
				tx_id=transcript_ids[index2],
				gene_id=gene_ids[index2])
		transcript.gr.list<-split(transcript.gr,transcript_ids[index2])
	}
	
	
	
	
	r<-list()
	gene<-gene.gr.list
	#r$transcript<-transcript.gr.list
	return(gene)
}

gene_range = gtf2GRangesList(gtf_file, skip=2, format="gff3")

design_file = file_all[1]
sampleTable = read.table(design_file, header=FALSE, as.is=TRUE)
Condition = matrix(unlist(strsplit(sampleTable[,3], split="_")), ncol=3, byrow=TRUE)
Condition = Condition[, -dim(Condition)[2]]
Condition = paste(Condition[,1], Condition[,2], sep="_")
sampleTable = cbind(sampleTable, Condition)
colnames(sampleTable) = c("Experiment", "Sample_ID", "Sample_Name", "QualityScheme", "Condition")
design = ~Condition

DESeqDataSetFromHTSeqCount_lxue = function (sampleTable, directory = "", design, col=1, ...)
{
	if (missing(design)) {
		stop("design is missing")
	}
	l <- lapply(paste(sampleTable[,col], "/", sampleTable[,col], rep(".counts.txt", length(sampleTable[,1])), sep=""),
			function(fn) read.table(file.path(directory, fn)))
	if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1))))
		stop("Gene IDs (first column) differ between files.")
	tbl <- sapply(l, function(a) a$V2)
	rownames(tbl) <- l[[1]]$V1
	rownames(sampleTable) <- sampleTable[, col]
	specialRows <- rownames(tbl) %in% c("no_feature", "ambiguous","too_low_aQual", "not_aligned", "alignment_not_unique",
			"__no_feature", "__ambiguous","__too_low_aQual", "__not_aligned", "__alignment_not_unique")
	tbl <- tbl[!specialRows, ]
	dds <- DESeqDataSetFromMatrix(countData = tbl, colData = sampleTable[,
					-(1:col), drop = FALSE], design = design)
	return(dds)
}


###########################################################################################
### Multiple factors																	###
###########################################################################################
## DESeqDataSetFromHTSeqCount cannot be used directly, because the special rows need to be removed first.
dds_raw = DESeqDataSetFromHTSeqCount_lxue(sampleTable,directory=dire, design= design, col=2)

dds = dds_raw

rowRanges(dds) = gene_range

fpkm_out = fpkm(dds) # must use R version 3.1.1 or above

system(paste("mkdir -p ", tgtfolder, sep=""))
file_out = paste(tgtfolder, "genes.fpkm_tracking",sep="/")
write.table(as.data.frame(fpkm_out), file=file_out, sep="\t", quote=FALSE)

###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################
rm(list=ls())
gc()


