#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#USAGE - Rscript VariantAnnotation.R [input.vcf] [reference.gb] [output.csv]


#check for number of arguments
if (length(args)<3) {
  stop("Three arguments must be supplied", call.=FALSE)
}

#load libraries
pacman::p_load(
      VariantAnnotation,
      genbankr,
      tidyverse, 
      GenomicFeatures,
      randomcoloR,
      gtools,
      GenomeInfoDb,
      dplyr
)


#arguments
VCF=args[1]
GB=args[2]
outputfile = args[3]



#import vcf generated from CovR MSA using jvarkit msa2vcf script
vcf <- readVcf(VCF)

#import reference genbank file - formated with NCBI
gb<-readGenBank(GB)

#import reference gff file that has been bgzipped and index
txdb <- makeTxDbFromGenBank(gb)
#txdb = makeTxDbFromGFF(file="NZ_CP049085.2.sort.gff3.gz")
#fa = open(FaFile("NZ_CP049085.2.fa"))

# Need seqlevels to be equivalent between query (vcf) and subject (txdb)
seqlevels(vcf) <- seqlevels(txdb) 

#predict amino acid change
coding <- predictCoding(vcf, txdb, gb)
#coding = predictCoding(vcf, txdb, fa)
coding2 <- coding %>% as_tibble()
#coding2 <- subset(coding2, select = -c(6, 8, 10))

# create matrix of alt alleles
GT <- as.matrix(geno(vcf)$GT)

# Remove /* from GT
GT <- sub("\\/.*", "", GT)

# transpose and create a dataframe
df1 <- as.data.frame(t(GT))

# TODO
#Add consequence to the column names
# First possibly concatenate first and second column of coding effects 
# paste(coding2$seqnames,coding2$start, sep=":")
# And then have partial match condition to append REFAA PROTEINLOC VARAA based on that match
# paste AA change to separate dataframe and remove whitespace
#df2<- paste0(coding2$REFAA,coding2$PROTEINLOC,coding2$VARAA, sep = "")%>%as_data_frame()
#df3 <- paste(coding2$CONSEQUENCE,df2$value, sep = " ")%>%as_data_frame()
#colnames(df1)<-paste(colnames(df1),df3$value)

# convert all values into numerics and create a row sums which will sum all mutations found per isolate
df1 <- df1 %>% 
      mutate_all(funs(as.numeric(.)))
df1 <- df1 %>% 
      mutate(mutations = rowSums(.!=0))

# Make rownames into first column
df1 <- rownames_to_column(df1, "sample_ID")

# Create vector that sums all non-zero values to get totals of how many times a mutation occurred at this locus
colSums <- colSums(df1!=0)
df1[nrow(df1) + 1, ] <- colSums

## Need to filter mutations for >= 1
df1<-filter(df1, mutations !=0)

#colnames(df1)=gsub("NZ_CP043530.","",colnames(df1))
colnames(df1)=gsub("^.*?\\:","",colnames(df1))

# Also going to substitute '.' with '>'
colnames(df1)=gsub("\\/",">",colnames(df1))

# Prior to any annotation write out dataframes
write_csv(coding2,"coding_changes.csv" )
write_csv(df1, outputfile)
