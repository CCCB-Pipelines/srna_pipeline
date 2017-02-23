# Input: Rscript --vanilla srna_processing.R [PROJECT_DIR] [OUTPUT_DIR](optional)
# 
# This script assumes the PROJECT_DIR contains at least two subdirectories:
#   1) pipeline_results/  # containing sRNABench result output
#   2) Annotations/       # contain sample and subject annotation files
#
# This program takes directories of sRNAWorkBench Output and generate count matrices
# for miRBase_main.txt, 
#     mature_sense_nonRed.grouped
#     Homo_sapiens_trna_sense_nonRed.grouped
#     pirna_sense_nonRed.grouped
#     snorna_sense_nonRed.grouped
#
#

library(reshape2)
library(xlsx)


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Input: Rscript --vanilla srna_processing.R [PROJECT_DIR] [OUTPUT_DIR](optional)\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "matrices"
}

if(args[1] == '-h'){
  help='
Input: Rscript --vanilla srna_processing.R [PROJECT_DIR] [OUTPUT_DIR](optional)

This script assumes the PROJECT_DIR contains at least two subdirectories:
  1) pipeline_results/  # containing sRNABench result output AS IS
  2) Annotations/       # contain sample (*_Samples.tsv) and subject (*_Subjects.tsv) annotation files
     - Sample Annotation File (*_Samples.tsv): Require 4 columns,  Project, Subject, Sample_ID, SeqFile
                Batch:     String [:alnum:] and _
                Subject:   Subject that the sample was derived from
                Sample_ID: Unique Label for the sample sequenced
                SeqFile:   Sequence file name used to generate sRNABench pipeline data without file extension (.fastq)
     - Subject Annotation File (*_Subjects.tsv):  Require 2 columns: Project, Subject
                Project:   Project name using String [::alnum:] and _
                Subject:   Subject enrolled into this project.
     NOTE: If a *_SampleSheet.xls exists withint the directory,  *_SampleSheet.xls will be used and assume there are
           two sheets: Samples and Subjects

This program takes directories of sRNAWorkBench Output and generate count matrices
for miRBase_main.txt,
    mature_sense_nonRed.grouped
    Homo_sapiens_trna_sense_nonRed.grouped
    pirna_sense_nonRed.grouped
    snorna_sense_nonRed.grouped
'
  stop(help, call.=F)
}

HOME=args[1]
OUTDIR=paste0(HOME, args[2])

#HOME="~/scribbles/projects/IonitaGhiran/srna_hemoglobin/"
#OUTDIR=paste0(HOME, "tmp/")  

print(paste("start processing:", HOME))
print(paste("OUTDIR:", OUTDIR))

AnnotateFiles=list.files(paste0(HOME, "Annotations/"), full.names=T)
SampleAnnot=tryCatch({
  if(any(grep("SampleSheet.xlsx", AnnotateFiles))){
        SampleAnnotFile=AnnotateFiles[grep("SampleSheet.xlsx", AnnotateFiles)[1]]
        read.xlsx(SampleAnnotFile,  sheetName="Samples")
    }else if(any(grep("Samples.tsv", AnnotateFiles))){
        SampleAnnotFile=AnnotateFiles[grep("Samples.tsv", AnnotateFiles)[1]]
        read.table(SampleAnnotFile, sep="\t", header=T)
    }}, error = function(err) {
       stop(err)
    })


DATADIR=paste0(HOME,"/pipeline_results/")  # srna pipeline result directory

# DATADIR:  all batches processed by sRNAWorkBench should reside within the directory "pipeline_results"
#         the script will automatically search for the data files
#
#

# Get all miRBase_main files
FILES=list(
           #miRBase=list.files(path=DATADIR, pattern="miRBase_main.txt", 
           #      full.names=T, all.files=T, recursive=T),
           mature_sense_nonRed=list.files(path=DATADIR, pattern="mature_sense_nonRed.grouped", 
                                   full.names=T, all.files=T, recursive=T),
           tRNA_nonRed=list.files(path=DATADIR, pattern="Homo_sapiens_trna_sense_nonRed.grouped", 
                           full.names=T, all.files=T, recursive=T),
           pirna_nonRed=list.files(path=DATADIR, pattern="pirna_sense_nonRed.grouped", 
                           full.names=T, all.files=T, recursive=T),
           snorna_nonRed=list.files(path=DATADIR, pattern="snorna_sense_nonRed.grouped", 
                           full.names=T, all.files=T, recursive=T)
)
print(FILES)

srna_RC=lapply(names(FILES), function(X, FILES, SampleAnnot){
  MasterMtx=c()
  for(infile in FILES[[X]]){
    if(grepl("PlantAndVirus",infile)){next}  # skip ones in PlantAndVirus files
    SeqID=gsub(".*pipeline_results\\/\\/([[:alnum:]_]+)\\/.*\\.grouped", "\\1", infile, perl=T)    # Resolve file pattern to extract Sample names
  
    print(SeqID)  
  
    Sample_ID=subset(SampleAnnot, SeqFile==SeqID, Sample_ID)
    inputMtx=read.table(infile, sep="\t", comment.char = "", header=T)
    inputMtx=cbind(Sample_ID=Sample_ID, inputMtx)
    MasterMtx=rbind(MasterMtx, inputMtx)
  }
  
  if("RC" %in% colnames(MasterMtx)){
    count.mtx=dcast(MasterMtx, name ~ Sample_ID, value.var="RC", mean)
  }else{
    count.mtx=dcast(MasterMtx, name ~ Sample_ID, value.var="read.count.non.redundant.", mean)
  }
  features=count.mtx[,1]
  count.mtx=data.matrix(count.mtx[-1]) # remove name for the matrix
  count.mtx[is.nan(count.mtx)]<-0
  rownames(count.mtx)=features
    
  return(count.mtx)
}, FILES=FILES, SampleAnnot=SampleAnnot)
names(srna_RC)=names(FILES)



ifelse(!dir.exists(OUTDIR), dir.create(OUTDIR), FALSE)
for(srna_type in names(srna_RC)){
  write.table(srna_RC[[srna_type]], file=paste0(OUTDIR,"/",srna_type, "_raw_count_matrix.tsv"),
              sep="\t", quote=F, row.names=T, col.names=T)
}
