# This function determines the number of reads in fastq.gz files. It counts the 
# number of lines in each file and divides by four (since there are four lines 
# per read in a fastq). It can be used with either single or paired reads. For 
# paired reads it counts the number of reads in the forward file.

count.reads.in.fastq <- function(fastq.dir, fastq.files){ #, read.type){
  
#  if(read.type == "single"){.  ### FOR SINGLE END FILES
    tot.reads <- do.call(rbind, lapply(1:length(fastq.files), function(f){
      cnt <- system2(command = "zcat",
                     args = c("<", paste0(fastq.dir, "/", fastq.files[f]), "|", "wc", "-l"),
                     stdout = TRUE, stderr = TRUE)
      return(as.numeric(cnt)/4)
    }))
    
#  } else if(read.type == "paired"){
    
#    tot.reads <- do.call(rbind, lapply(1:length(fastq.files), function(f){
#      fname <- fastq.files[f]
#      fwd.fastq <- paste0(fname, "_R1.fastq.gz")
#      fwd.reads <- system2(command = "zcat",
#                           args = c("<", paste0(fastq.dir, "/", fwd.fastq), "|", "wc", "-l"),
#                           stdout = TRUE, stderr = TRUE)
#      return(tot.reads = as.numeric(fwd.reads))
#    }))

#  } else(print("read.type must be single or paired"))
  
  return(data.frame(Sample = fastq.files, tot.reads = tot.reads))
}