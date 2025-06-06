library(dplyr)
source("R/Count.reads.in.fastq.R")
source("R/Find.primer.seqs.in.fastq.R")

project <- "Reference.genome"
#primer.file <- "Illumina.old.smallRNA.primer"
primer.file <- "primer.sequences"
fastq.dir <- paste0("/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/Data/Fastq/", project)
primers <- read.csv(paste0("data-raw/", primer.file, ".csv"))
read.type <- "single" # single or paired

fastq.files <- list.files(path = fastq.dir, pattern = "*.fasta.gz")

# use these lines to limit the number of files and primers for testing purposes
#fastq.files <- fastq.files[25]
primers <- primers[1:3,]

if(read.type == "paired") {
  fastq.files <- sapply(strsplit(fastq.files, split = "_R"), function(f){f[1]}) %>%
    unique()
  fastq.files <- paste0(fastq.files, "_R1.fastq.gz")
}

date()
all.matches <- primer.seqs.in.fastq(fastq.dir, fastq.files, primers, read.type)
date()

tot.reads <- count.reads.in.fastq(fastq.dir, fastq.files)

rev <- all.matches %>% group_by(Sample) %>% summarise(rev = sum(rev))
sum.by.ind <- all.matches %>% group_by(Sample) %>% summarise(fwd = sum(fwd)) %>%
  left_join(rev) %>% left_join(tot.reads)

rev <- all.matches %>% group_by(locus) %>% summarise(rev = sum(rev))
sum.by.loc <- all.matches %>% group_by(locus) %>% summarise(fwd = sum(fwd)) %>%
  left_join(rev)
sum.by.loc$ratio <- sapply(1:nrow(sum.by.loc), function(l){
  max(sum.by.loc[l,2:3]) / min(sum.by.loc[l,2:3])
})

save(all.matches, sum.by.ind, sum.by.loc, file = paste0("data/Fastq.summary.", project, ".", primer.file,".rda"))


