#title: "FoulkesGen"
#author: "Foulkes Lab"
#date: "11/19/2019"
#output: html_document
# Customize as needed for file locations

# Modify data.dir to indicate the location of the GWAStutorial files
# Intermediate data files will also be stored in this same location unless you set out.dir
data.dir <- '/Users/johnsonam/Documents/Independent Study/Evoke'

out.dir <- data.dir                     # may want to write to a separate dir to avoid clutter

# Download files
#urlSupport <- "https://www.mtholyoke.edu/courses/afoulkes/Data/GWAStutorial/GWASTutorial_Files.zip"
#zipSupport.fn <- sprintf("%s/GWAStutorial_Files.zip", data.dir) 

# Input files
gwas.fn <- lapply(c(bed='bed',bim='bim',fam='fam',gds='gds'), function(n) sprintf("%s/GWAStutorial.%s", data.dir, n))
clinical.fn <- sprintf("%s/GWAStutorial_clinical.csv", data.dir) 
onethou.fn <- lapply(c(info='info',ped='ped'), function(n) sprintf("%s/chr16_1000g_CEU.%s", data.dir, n))
protein.coding.coords.fname <- sprintf("%s/ProCodgene_coords.csv", data.dir)

# Output files
gwaa.fname <- sprintf("%s/GWAStutorialout.txt", out.dir)
gwaa.unadj.fname <- sprintf("%s/GWAStutorialoutUnadj.txt", out.dir)
impute.out.fname <- sprintf("%s/GWAStutorial_imputationOut.csv", out.dir)
CETP.fname <- sprintf("%s/CETP_GWASout.csv", out.dir)

# Working data saved between each code snippet so each can run independently.
# Use save(data, file=working.data.fname(num))
working.data.fname <- function(num) { sprintf("%s/working.%s.Rdata", out.dir, num) }

# Download and unzip data needed for this tutorial
library(downloader)

download(urlSupport, zipSupport.fn)
unzip(zipSupport.fn, exdir = data.dir)

library(snpStats)

# Read in PLINK files
geno <- read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))

# Obtain the SnpMatrix object (genotypes) table from geno list
# Note: Phenotypes and covariates will be read from the clinical data file, below
genotype <- geno$genotype
print(genotype)                  


#Obtain the SNP information from geno list
genoBim <- geno$map
colnames(genoBim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(genoBim))


# Remove raw file to free up memory
rm(geno)

# Read in clinical file
clinical <- read.csv(clinical.fn,
                     colClasses=c("character", "factor", "factor", rep("numeric", 4)))
rownames(clinical) <- clinical$FamID
print(head(clinical))


# Subset genotype for subject data
genotype <- genotype[clinical$FamID, ]
print(genotype)  # Tutorial: All 1401 subjects contain both clinical and genotype data

# Write genotype, genoBim, clinical for future use
save(genotype, genoBim, clinical, file = working.data.fname(1))

# Create SNP summary statistics (MAF, call rate, etc.)
snpsum.col <- col.summary(genotype)
print(head(snpsum.col))

# Setting thresholds
call <- 0.95
minor <- 0.01

# Filter on MAF and call rate
use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(ncol(genotype)-sum(use),"SNPs will be removed due to low MAF or call rate.\n") #203287 SNPs will be removed

# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]

print(genotype)                           
# 658186 SNPs remain

# Write subsetted genotype data and derived results for future use
save(genotype, snpsum.col, genoBim, clinical, file=working.data.fname(2))

# Sample level filtering

#source("globals.R")
#error here

# load data created in previous snippets
load(working.data.fname(2))

library(snpStats)

source("https://bioconductor.org/biocLite.R")
biocLite("SNPRelate")
library(SNPRelate)                     # LD pruning, relatedness, PCA
library(plyr)

# Create sample statistics (Call rate, Heterozygosity)
snpsum.row <- row.summary(genotype)

# Add the F stat (inbreeding coefficient) to snpsum.row
MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)

head(snpsum.row)

# Setting thresholds
sampcall <- 0.95    # Sample call rate cut-off
hetcutoff <- 0.1    # Inbreeding coefficient cut-off

sampleuse <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampcall & abs(hetF) <= hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE    # remove NA's as well
cat(nrow(genotype)-sum(sampleuse), "subjects will be removed due to low sample call rate or inbreeding coefficient.\n") #0 subjects removed

# Subset genotype and clinical data for subjects who pass call rate and heterozygosity crtieria
genotype <- genotype[sampleuse,]
clinical<- clinical[ rownames(genotype), ]

# Checking for Relatedness

ld.thresh <- 0.2    # LD cut-off
kin.thresh <- 0.1   # Kinship cut-off

# Create gds file, required for SNPRelate functions
snpgdsBED2GDS(gwas.fn$bed, gwas.fn$fam, gwas.fn$bim, gwas.fn$gds)

genofile <- openfn.gds(gwas.fn$gds, readonly = FALSE)

# Automatically added "-1" sample suffixes are removed
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)

#Prune SNPs for IBD analysis
set.seed(1000)
geno.sample.ids <- rownames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          sample.id = geno.sample.ids, # Only analyze the filtered samples
                          snp.id = colnames(genotype)) # Only analyze the filtered SNPs

snpset.ibd <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.ibd),"will be used in IBD analysis\n")  # Tutorial: expect 72812 SNPs

# Find IBD coefficients using Method of Moments procedure.  Include pairwise kinship.
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE,
                    sample.id = geno.sample.ids,
                    snp.id = snpset.ibd,
                    num.thread = 1)
ibdcoeff <- snpgdsIBDSelection(ibd)     # Pairwise sample comparison
head(ibdcoeff)

# Check if there are any candidates for relatedness
ibdcoeff <- ibdcoeff[ ibdcoeff$kinship >= kin.thresh, ]

# iteratively remove samples with high kinship starting with the sample with the most pairings
related.samples <- NULL
while ( nrow(ibdcoeff) > 0 ) {
  
  # count the number of occurrences of each and take the top one
  sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
  rm.sample <- sample.counts[1, 'x']
  cat("Removing sample", as.character(rm.sample), 'too closely related to', sample.counts[1, 'freq'],'other samples.\n')
  
  # remove from ibdcoeff and add to list
  ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample,]
  related.samples <- c(as.character(rm.sample), related.samples)
}

# filter genotype and clinical to include only unrelated samples
genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$FamID %in% related.samples), ]

geno.sample.ids <- rownames(genotype)

cat(length(related.samples), "similar samples removed due to correlation coefficient >=", kin.thresh,"\n") 

print(genotype)

# Checking for ancestry

# Find PCA matrix
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.ibd, num.thread=1)

# Create data frame of first two principal comonents
pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],    # the first eigenvector
                    PC2 = pca$eigenvect[,2],    # the second eigenvector
                    stringsAsFactors = FALSE)

# Plot the first two principal comonents
plot(pctab$PC2, pctab$PC1, xlab="Principal Component 2", ylab="Principal Component 1", main = "Ancestry Plot")

# Close GDS file
closefn.gds(genofile)

# Overwrite old genotype with new filtered version
save(genotype, genoBim, clinical, file=working.data.fname(3))

# Hardy-Weinberg SNP filtering on CAD controls

hardy <- 10^-6      # HWE cut-off

CADcontrols <- clinical[ clinical$CAD==0, 'FamID' ]
snpsum.colCont <- col.summary( genotype[CADcontrols,] )
HWEuse <- with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)

HWEuse[is.na(HWEuse)] <- FALSE          # Remove NA's as well
cat(ncol(genotype)-sum(HWEuse),"SNPs will be removed due to high HWE.\n")  # 1296 SNPs removed

# Subset genotype and SNP summary data for SNPs that pass HWE criteria
genotype <- genotype[,HWEuse]

print(genotype)                           # 656890 SNPs remain

# Overwrite old genotype with new filtered version
save(genotype, genoBim, clinical, file=working.data.fname(4))

#data generation
# Set LD threshold to 0.2
ld.thresh <- 0.2

set.seed(1000)
geno.sample.ids <- rownames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          sample.id = geno.sample.ids, # Only analyze the filtered samples
                          snp.id = colnames(genotype)) # Only analyze the filtered SNPs

snpset.pca <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.pca),"\n")  #72578 SNPs will be used in PCA analysis

pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.pca, num.thread=1)

# Find and record first 10 principal components
# pcs will be a N:10 matrix.  Each column is a principal component.
pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : 10],
                  stringsAsFactors = FALSE)
print(head(pcs))

# Close GDS file
closefn.gds(genofile)

# Store pcs for future reference with the rest of the derived data
save(genotype, genoBim, clinical, pcs, file=working.data.fname(5))

# Read in 1000g data for given chromosome 16
thougeno <- read.pedfile(onethou.fn$ped, snps = onethou.fn$info, which=1)

# Obtain genotype data for given chromosome
genoMatrix <- thougeno$genotypes

# Obtain the chromosome position for each SNP
support <- thougeno$map
colnames(support)<-c("SNP", "position", "A1", "A2")
head(support)

# Imputation of non-typed 1000g SNPs
presSnps <- colnames(genotype)

# Subset for SNPs on given chromosome
presSnps <- colnames(genotype)
presDatChr <- genoBim[genoBim$SNP %in% presSnps & genoBim$chr==16, ]
targetSnps <- presDatChr$SNP

# Subset 1000g data for our SNPs
# "missing" and "present" are snpMatrix objects needed for imputation rules
is.present <- colnames(genoMatrix) %in% targetSnps

missing <- genoMatrix[,!is.present]
print(missing)                          # Almost 400,000 SNPs

present <- genoMatrix[,is.present]
print(present)                          # Our typed SNPs

# Obtain positions of SNPs to be used for imputation rules
pos.pres <- support$position[is.present]
pos.miss <- support$position[!is.present]

# Calculate and store imputation rules using snp.imputation()
rules <- snp.imputation(present, missing, pos.pres, pos.miss)

# Remove failed imputations
rules <- rules[can.impute(rules)]
cat("Imputation rules for", length(rules), "SNPs were estimated\n")  # Imputation rules for 197888 SNPs were estimated

# Quality control for imputation certainty and MAF
# Set thresholds
r2threshold <- 0.7
minor <- 0.01

# Filter on imputation certainty and MAF
rules <- rules[imputation.r2(rules) >= r2threshold]

cat(length(rules),"imputation rules remain after imputations with low certainty were removed\n")  # 162565 imputation rules remain after imputations with low certainty were removed

rules <- rules[imputation.maf(rules) >= minor]
cat(length(rules),"imputation rules remain after MAF filtering\n")  # 162565 imputation rules remain after MAF filtering

# Obtain posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules, target, as.numeric=FALSE)
print(imputed)  # 162565 SNPs were imputed

rm(genoMatrix)
rm(missing)
rm(present)

# Add new imputed, target and rules data to saved results
save(genotype, genoBim, clinical, pcs, imputed, target, 
     rules, support, file=working.data.fname(6))

##Gene wide association analysis
# Genome-wide Association Analysis
# Parallel implementation of linear model fitting on each SNP

GWAA <- function(genodata=genotypes,  phenodata=phenotypes, family = gaussian, filename=NULL,
                 append=FALSE, workers=getOption("mc.cores",2L), flip=TRUE,
                 select.snps=NULL, hosts=NULL, nSplits=10)
{
  if (!require(doParallel)) { stop("Missing doParallel package") }
  
  #Check that a filename was specified
  if(is.null(filename)) stop("Must specify a filename for output.")
  
  #Check that the genotype data is of class 'SnpMatrix'
  if( class(genodata)!="SnpMatrix") stop("Genotype data must of class 'SnpMatrix'.")
  
  #Check that there is a variable named 'phenotype' in phenodata table
  if( !"phenotype" %in% colnames(phenodata))  stop("Phenotype data must have column named 'phenotype'")
  
  #Check that there is a variable named 'id' in phenodata table
  if( !"id" %in% colnames(phenodata)) stop("Phenotype data must have column named 'id'.")
  
  #If a vector of SNPs is given, subset genotype data for these SNPs
  if(!is.null(select.snps)) genodata<-genodata[,which(colnames(genodata)%in%select.snps)]
  
  #Check that there are still SNPs in 'SnpMatrix' object
  if(ncol(genodata)==0) stop("There are no SNPs in the 'SnpMatrix' object.")
  
  #Print the number of SNPs to be checked
  cat(paste(ncol(genodata), " SNPs included in analysis.\n"))
  
  #If append=FALSE than we will overwrite file with column names
  if(!isTRUE(append)) {
    columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")
    write.table(t(columns), filename, row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
  
  # Check sample counts
  if (nrow(phenodata) != nrow(genodata)) {
    warning("Number of samples mismatch.  Using subset found in phenodata.")
  }
  
  # Order genodata rows to be the same as phenodata
  genodata <- genodata[phenodata$id,]
  
  cat(nrow(genodata), "samples included in analysis.\n")
  
  # Change which allele is counted (major or minor)
  flip.matrix<-function(x) {
    zero2 <- which(x==0)
    two0 <- which(x==2)
    x[zero2] <- 2
    x[two0] <- 0
    return(x)
  }
  
  nSNPs <- ncol(genodata)
  genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset
  
  snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
  snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group
  
  if (is.null(hosts)) {
    # On Unix this will use fork and mclapply.  On Windows it
    # will create multiple processes on localhost.
    cl <- makeCluster(workers)
  } else {
    # The listed hosts must be accessible by the current user using
    # password-less ssh with R installed on all hosts, all 
    # packages installed, and "rscript" is in the default PATH.
    # See docs for makeCluster() for more information.
    cl <- makeCluster(hosts, "PSOCK")
  }
  show(cl)                            # report number of workers and type of parallel implementation
  registerDoParallel(cl)
  
  foreach (part=1:nSplits) %do% {
    # Returns a standar matrix of the alleles encoded as 0, 1 or 2
    genoNum <- as(genodata[,snp.start[part]:snp.stop[part]], "numeric")
    
    # Flip the numeric values of genotypes to count minor allele
    if (isTRUE(flip)) genoNum <- flip.matrix(genoNum)
    
    # For each SNP, concatenate the genotype column to the
    # phenodata and fit a generalized linear model
    rsVec <- colnames(genoNum)
    res <- foreach(snp.name=rsVec, .combine='rbind') %dopar% {
      a <- summary(glm(phenotype~ . - id, family=family, data=cbind(phenodata, snp=genoNum[,snp.name])))
      a$coefficients['snp',]
    }
    
    # write results so far to a file
    write.table(cbind(rsVec,res), filename, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)
    
    cat(sprintf("GWAS SNPs %s-%s (%s%% finished)\n", snp.start[part], snp.stop[part], 100*part/nSplits))
  }
  
  stopCluster(cl)
  
  return(print("Done."))
}

library(GenABEL)
source("GWAA.R")

# Merge clincal data and principal components to create phenotype table
phenoSub <- merge(clinical,pcs)      # data.frame => [ FamID CAD sex age hdl pc1 pc2 ... pc10 ]

# We will do a rank-based inverse normal transformation of hdl
phenoSub$phenotype <- rntransform(phenoSub$hdl, family="gaussian")

# Show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(phenoSub$hdl, main="Histogram of HDL", xlab="HDL")
hist(phenoSub$phenotype, main="Histogram of Tranformed HDL", xlab="Transformed HDL")

# Remove unnecessary columns from table
phenoSub$hdl <- NULL
phenoSub$ldl <- NULL
phenoSub$tg <- NULL
phenoSub$CAD <- NULL

# Rename columns to match names necessary for GWAS() function
phenoSub <- rename(phenoSub, replace=c(FamID="id"))

# Include only subjects with hdl data
phenoSub<-phenoSub[!is.na(phenoSub$phenotype),]
# 1309 subjects included with phenotype data

print(head(phenoSub))

# Run GWAS analysis
# Note: This function writes a file, but does not produce an R object
start <- Sys.time()
GWAA(genodata=genotype, phenodata=phenoSub, filename=gwaa.fname)

end <- Sys.time()
print(end-start)

# Add phenosub to saved results
save(genotype, genoBim, clinical, pcs, imputed, target, rules, phenoSub, support, file=working.data.fname(7))

# Carry out association testing for imputed SNPs using snp.rhs.tests()
rownames(phenoSub) <- phenoSub$id

imp <- snp.rhs.tests(phenotype ~ sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                     family = "Gaussian", data = phenoSub, snp.data = target, rules = rules)

# Obtain p values for imputed SNPs by calling methods on the returned GlmTests object.
results <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results <- results[!is.na(results$p.value),]

#Write a file containing the results
write.csv(results, impute.out.fname, row.names=FALSE)

# Merge imputation testing results with support to obtain coordinates
imputeOut<-merge(results, support[, c("SNP", "position")])
imputeOut$chr <- 16

imputeOut$type <- "imputed"

# Find the -log_10 of the p-values
imputeOut$Neg_logP <- -log10(imputeOut$p.value)

# Order by p-value
imputeOut <- arrange(imputeOut, p.value)
print(head(imputeOut))

# Returns the subset of SNPs that are within extend.boundary of gene
# using the coords table of gene locations.
map2gene <- function(gene, coords, SNPs, extend.boundary = 5000) {
  coordsSub <- coords[coords$gene == gene,] #Subset coordinate file for spcified gene
  
  coordsSub$start <- coordsSub$start - extend.boundary # Extend gene boundaries
  coordsSub$stop <- coordsSub$stop + extend.boundary
  
  SNPsub <- SNPs[SNPs$position >= coordsSub$start & SNPs$position <= coordsSub$stop &
                   SNPs$chr == coordsSub$chr,] #Subset for SNPs in gene
  
  return(data.frame(SNPsub, gene = gene, stringsAsFactors = FALSE))
}

# Read in file containing protein coding genes coords
genes <- read.csv(protein.coding.coords.fname, stringsAsFactors = FALSE)

# Subset for CETP SNPs
impCETP <- map2gene("CETP", coords = genes, SNPs = imputeOut)

# Filter only the imputed CETP SNP genotypes 
impCETPgeno <- imputed[, impCETP$SNP ]

save(genotype, genoBim, clinical, pcs, imputed, target, rules,
     phenoSub, support, genes, impCETP, impCETPgeno, imputeOut, file = working.data.fname(8))

## Post-analytic visualization and genomic interrogation
# Read in GWAS output that was produced by GWAA function
GWASout <- read.table(gwaa.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Find the -log_10 of the p-values
GWASout$Neg_logP <- -log10(GWASout$p.value)

# Merge output with genoBim by SNP name to add position and chromosome number
GWASout <- merge(GWASout, genoBim[,c("SNP", "chr", "position")])
rm(genoBim)

# Order SNPs by significance
GWASout <- arrange(GWASout, -Neg_logP)
print(head(GWASout))

# Combine typed and imputed
GWASout$type <- "typed"

GWAScomb<-rbind.fill(GWASout, imputeOut)
head(GWAScomb)
tail(GWAScomb)

# Subset for CETP SNPs
typCETP <- map2gene("CETP", coords = genes, SNPs = GWASout)

# Combine CETP SNPs from imputed and typed analysis
CETP <- rbind.fill(typCETP, impCETP)[,c("SNP","p.value","Neg_logP","chr","position","type","gene")]
print(CETP)

write.csv(CETP, CETP.fname, row.names=FALSE) # save for future use

save(genotype, clinical, pcs, imputed, target, rules, phenoSub, support, genes,
     impCETP, impCETPgeno, imputeOut, GWASout, GWAScomb, CETP, file=working.data.fname(9))

# Receives a data.frame of SNPs with Neg_logP, chr, position, and type.
# Plots Manhattan plot with significant SNPs highlighted.
GWAS_Manhattan <- function(GWAS, col.snps=c("black","gray"),
                           col.detected=c("black"), col.imputed=c("blue"), col.text="black",
                           title="GWAS Tutorial Manhattan Plot", display.text=TRUE,
                           bonferroni.alpha=0.05, bonferroni.adjustment=1000000,
                           Lstringent.adjustment=10000) {
  
  bonferroni.thresh <- -log10(bonferroni.alpha / bonferroni.adjustment)
  Lstringent.thresh <- -log10(bonferroni.alpha / Lstringent.adjustment)
  xscale <- 1000000
  
  manhat <- GWAS[!grepl("[A-z]",GWAS$chr),]
  
  #sort the data by chromosome and then location
  manhat.ord <- manhat[order(as.numeric(manhat$chr),manhat$position),]
  manhat.ord <- manhat.ord[!is.na(manhat.ord$position),]
  
  ##Finding the maximum position for each chromosome
  max.pos <- sapply(1:21, function(i) { max(manhat.ord$position[manhat.ord$chr==i],0) })
  max.pos2 <- c(0, cumsum(max.pos))                  
  
  #Add spacing between chromosomes
  max.pos2 <- max.pos2 + c(0:21) * xscale * 10
  
  #defining the positions of each snp in the plot
  manhat.ord$pos <- manhat.ord$position + max.pos2[as.numeric(manhat.ord$chr)]
  
  # alternate coloring of chromosomes
  manhat.ord$col <- col.snps[1 + as.numeric(manhat.ord$chr) %% 2]
  
  # draw the chromosome label roughly in the middle of each chromosome band
  text.pos <- sapply(c(1:22), function(i) { mean(manhat.ord$pos[manhat.ord$chr==i]) })
  
  # Plot the data
  plot(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
       pch=20, cex=.3, col= manhat.ord$col[manhat.ord$type=="typed"], xlab=NA,
       ylab="Negative Log P-value", axes=F, ylim=c(0,max(manhat$Neg_logP)+1))
  #Add x-label so that it is close to axis
  mtext(side = 1, "Chromosome", line = 1.25)
  
  points(manhat.ord$pos[manhat.ord$type=="imputed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="imputed"],
         pch=20, cex=.4, col = col.imputed)
  
  points(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
         pch=20, cex=.3, col = manhat.ord$col[manhat.ord$type=="typed"])
  
  axis(2)
  abline(h=0)
  
  SigNifSNPs <- as.character(GWAS[GWAS$Neg_logP > Lstringent.thresh & GWAS$type=="typed", "SNP"])
  
  #Add legend
  legend("topright",c("Bonferroni corrected threshold (p = 5E-8)", "Candidate threshold (p = 5E-6)"),
         border="black", col=c("gray60", "gray60"), pch=c(0, 0), lwd=c(1,1),
         lty=c(1,2), pt.cex=c(0,0), bty="o", cex=0.6)
  
  #Add chromosome number
  text(text.pos/xscale, -.3, seq(1,22,by=1), xpd=TRUE, cex=.8)
  
  #Add bonferroni line
  abline(h=bonferroni.thresh, untf = FALSE, col = "gray60")
  
  #Add "less stringent" line
  abline(h=Lstringent.thresh, untf = FALSE, col = "gray60", lty = 2 )
  
  #Plotting detected genes
  #Were any genes detected?
  if (length(SigNifSNPs)>0){
    
    sig.snps <- manhat.ord[,'SNP'] %in% SigNifSNPs
    
    points(manhat.ord$pos[sig.snps]/xscale,
           manhat.ord$Neg_logP[sig.snps],
           pch=20,col=col.detected, bg=col.detected,cex=0.5)
    
    text(manhat.ord$pos[sig.snps]/xscale,
         manhat.ord$Neg_logP[sig.snps],
         as.character(manhat.ord[sig.snps,1]), col=col.text, offset=1, adj=-.1, cex=.5)
  }
}

# Create Manhattan Plot
GWAS_Manhattan(GWAScomb)

# Rerun the GWAS using unadjusted model
phenoSub2 <- phenoSub[,c("id","phenotype")] # remove all extra factors, leave only phenotype

GWAA(genodata=genotype, phenodata=phenoSub2, filename=gwaa.unadj.fname)

GWASoutUnadj <- read.table(gwaa.unadj.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Create QQ plots for adjusted and unadjusted model outputs
par(mfrow=c(1,2))
lambdaAdj <- estlambda(GWASout$t.value^2,plot=TRUE,method="median")
lambdaUnadj <- estlambda(GWASoutUnadj$t.value^2,plot=TRUE,method="median")

cat(sprintf("Unadjusted lambda: %s\nAdjusted lambda: %s\n", lambdaUnadj$estimate, lambdaAdj$estimate))

# Calculate standardized lambda
lambdaAdj_1000<-1+(lambdaAdj$estimate-1)/nrow(phenoSub)*1000
lambdaUnadj_1000<-1+(lambdaUnadj$estimate-1)/nrow(phenoSub)*1000
cat(sprintf("Standardized unadjusted lambda: %s\nStandardized adjusted lambda: %s\n", lambdaUnadj_1000, lambdaAdj_1000))

library(LDheatmap)
library(rtracklayer)

# Add "rs247617" to CETP
CETP <- rbind.fill(GWASout[GWASout$SNP == "rs247617",], CETP)

# Combine genotypes and imputed genotypes for CETP region
subgen <- cbind(genotype[,colnames(genotype) %in% CETP$SNP], impCETPgeno)     # CETP subsets from typed and imputed SNPs

# Subset SNPs for only certain genotypes
certain <- apply(as(subgen, 'numeric'), 2, function(x) { all(x %in% c(0,1,2,NA)) })
subgen <- subgen[,certain]

# Subset and order CETP SNPs by position
CETP <- CETP[CETP$SNP %in% colnames(subgen),]
CETP <- arrange(CETP, position)
subgen <- subgen[, order(match(colnames(subgen),CETP$SNP)) ]

# Create LDheatmap
ld <- ld(subgen, subgen, stats="R.squared") # Find LD map of CETP SNPs

ll <- LDheatmap(ld, CETP$position, flip=TRUE, name="myLDgrob", title=NULL)

# Add genes, recombination
llplusgenes <- LDheatmap.addGenes(ll, chr = "chr16", genome = "hg19", genesLocation = 0.01)

# Add plot of -log(p)
library(ggplot2)

plot.new()
llQplot2<-LDheatmap.addGrob(llplusgenes, rectGrob(gp = gpar(col = "white")),height = .34)
pushViewport(viewport(x = 0.483, y= 0.76, width = .91 ,height = .4))

grid.draw(ggplotGrob({
  qplot(position, Neg_logP, data = CETP, xlab="", ylab = "Negative Log P-value", xlim = range(CETP$position),
        asp = 1/10, color = factor(type), colour=c("#000000", "#D55E00")) + 
    theme(axis.text.x = element_blank(),
          axis.title.y = element_text(size = rel(0.75)), legend.position = "none", 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_color_manual(values = c("red", "black"))
}))

# Create regional association plot
# Create data.frame of most significant SNP only
library(postgwas)

snps<-data.frame(SNP=c("rs1532625"))

# Change column names necessary to run regionalplot function
GWAScomb <- rename(GWAScomb, c(p.value="P", chr="CHR", position="BP"))


# Edit biomartConfigs so regionalplot function
# pulls from human genome build 37/hg19

myconfig <- biomartConfigs$hsapiens
myconfig$hsapiens$gene$host <- "grch37.ensembl.org"
myconfig$hsapiens$gene$mart <- "ENSEMBL_MART_ENSEMBL"
myconfig$hsapiens$snp$host <- "grch37.ensembl.org"
myconfig$hsapiens$snp$mart <- "ENSEMBL_MART_SNP"


# Run regionalplot using HAPMAP data (pop = CEU)
regionalplot(snps, GWAScomb, biomart.config = myconfig, window.size = 400000, draw.snpname = data.frame(
  snps = c("rs1532625", "rs247617"), 
  text = c("rs1532625", "rs247617"),
  angle = c(20, 160),
  length = c(1, 1), 
  cex = c(0.8)
),
ld.options = list(
  gts.source = 2, 
  max.snps.per.window = 2000, 
  rsquare.min = 0.8, 
  show.rsquare.text = FALSE
),
out.format = list(file = NULL, panels.per.page = 2))


