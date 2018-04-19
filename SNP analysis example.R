#Example for handling SNPs
#Reference book: Bioinformatics using R (Acorn)
#Reference site: http://davinci.crg.es/estivill_lab/tools/SNPassoc/SupplementaryMaterial.pdf

#file directory: 'C:/Users/BISYN/Documents/R/win-library/3.4/'
#pacakge directory: C:\Users\BISYN\AppData\Local\Temp\RtmpQvfZ5g\downloaded_packages

#Notice: If library is not implement, please restart the R-Studio.

#Download SNPassoc and SNPs data
#source("http://bioconductor.org/biocLite.R")
biocLite("SNPassoc")
#install.packages("SNPassoc")
#remove.packages("SNPassoc")
library(SNPassoc)
#install.packages("haplo.stats")
library(haplo.stats)
#install.packages("zoo")
library(zoo)

#1 Phenotype SNP Association Analysis of each SNP
data("SNPs")
head(SNPs) #id, casco(case or control), sex(Femal, Male), blood.pre, protein, snp10001 ~ snpXXXXX(key)
head(SNPs.info.pos) #snp(key), chr, pos

SNPs[1:10,1:9]

#setup SNP data to column-6:40 in all data
mySNP <- setupSNP(SNPs, 6:40, sep="")

#association analysis
#dominant - 우성, codominant - 공우성, recessive - 열성, overdominant - 초우성
#casco: case or control / snp10001 - particular SNP
myres <- association(casco~sex+snp10001+blood.pre, data=mySNP,model.interaction = c("dominant","codominant"))
myres

#=================================================

#2 SNP Association Analysis of All SNPs
myres_all <- WGassociation(protein~1, data=mySNP, model="all")
myres_all <- WGassociation(protein, data=mySNP, model="all")
myres_all

#extract the dominant, revessive column in table
dominant(myres_all)
recessive(myres_all)

#Whole Gene Statistic: 
#row:Codominant, Dominant, Recessive, overdominant,log-Additive
#column: n, me, se, dif, lower, upper, p-value, AIC
WGstats(myres_all)

summary(myres_all)
plot(myres_all)
#extract the table : row-SNPs / column-comments, log-Addictive
resHapMap <- WGassociation(protein, data=mySNP, model="log")

#================================================

#3 SNP Association Analysis of Large HapMap data
data("HapMap")
str(HapMap)
str(HapMap.SNPs.pos)
MyHapMap <- setupSNP(HapMap, colSNPs=3:9307, sort=TRUE,info=HapMap.SNPs.pos,sep="")
#head(MyHapMap)
myHapMapres <- WGassociation(group, data=MyHapMap,model="dominant")
head(myHapMapres)
print(myHapMapres)
plot(myHapMapres, whole=TRUE)

#================================================

#4 Load the data for GWAS of Plink (TBD)

#================================================

#5 Data handling using GWASTools

#remove.packages("RSQLite")
source("http://bioconductor.org/biocLite.R")
biocLite("GWASTools")
biocLite("GWASdata")
biocLite("gdsfmt")
#install.packages("RSQLite")
#install.packages("gdsfmt")
#install.packages("ncdf4")
library(ncdf4)
library(RSQLite)
library(GWASTools)
library(GWASdata)
library(gdsfmt)

data(affy_scan_annot)
data(affy_snp_annot)

#set the affy_geno.nc file route in windows directory
file <- system.file("extdata", "affy_geno.nc", package="GWASdata")
file

nc <- NcdfGenotypeReader (file)

#check the number of scan and snp
nscan(nc) #47
nsnp(nc) #3300

geno <- getGenotype(nc, snp=c(1,20), scan=c(1,5))

data("affyScanADF") #ScanAnnotationDataFrame
data("affySnpADF") #SnpAnnotationDataFrame

#set up for extracting phenotype in snp
varMetadata(affySnpADF)

MyGenoData <- GenotypeData(nc, snpAnnot = affySnpADF, scanAnnot = affyScanADF)

str (MyGenoData)

mcr_chr <- missingGenotypeByScanChrom(MyGenoData)
head(mcr_chr$missing.counts)

mcr_sex <- missingGenotypeBySnpSex(MyGenoData)

par(mfrow = c(1,2))
hist(mcr_sex$missing.fraction, xlab="Fraction", main="Histogram for the missing fraction for every sex")
hist(mcr_chr$missing.fraction, xlab="Fraction", main="Histogram for the missing fraction along chromosomes")

sum(mcr_sex$missing.fraction>0.05)

#================================================

#6 Data handling of different type data for GWAS (TBD)

dir <- "<path/to/desired/directory>"
#install.packages("GenABEL")
library(GenABEL)
convert.snp.affymetrix(dir, map, outfile, skipaffym)

#================================================

#7 Annotation of SNP (TBD)
#site: https://cran.r-project.org/src/contrib/Archive/NCBI2R/

biocLite("NCBI2R")
#install.packages("NCBI2R")
library(NCBI2R)

#================================================

#8 Hardy-Weinberg Equilibrium (TBD)

remove.packages("GWASExactHW")
source("http://bioconductor.org/biocLite.R")
biocLite("GWASExactHW")
library(GWASExactHW)

pA <- runif(1)
pAA <- pA^2
pAa <- 2*pA*(1-pA)
paa <- (1-pA)^2
myCounts <- rmultinom(100,500,c(pAA,pAa,paa))

genotypes <- data.frame(t(myCounts))
colnames(genotypes) = c("nAA","nAa","naa")

hwPvalues <- HWExact(genotypes)

load("path/to/code/file/dir/Alzheimers.rda")

#================================================

#9 GWAS using CNV data (TBD)
#site: https://cran.r-project.org/src/contrib/Archive/CNVassoc/

#NOT available to CNVassoc packages
#install.packages("CNVassoc")
#biocLite("CNVassoc")

library(CNVassoc)

#================================================

#10 Visualization of GWAS

library(GWASTools)
library(SNPassoc)
myres <- WGassociation(protein,data=mySNP,model="all")
pvals <- dominant(myres)
qqPlot(pvals)

n <- 1000
pvals <- sample(-log10((1:n)/n),n,replace=TRUE)
chromosome <- c(rep(1,100), rep(2,150), rep(3,80), rep(4,90),
                rep(5,100), rep(6,60), rep(7,70), rep(8,70),
                rep(9,70), rep(10,50), rep("X",110), rep("Y",50))

manhattanPlot(pvals, chromosome)

pvals = dominant (myHapMapres)
pvals = -log10(pvals)

chromosome <- HapMap.SNPs.pos$chromosome

manhattanPlot(pvals, chromosome, signif=1e-5)

#install.packages("gap")
#install.packages("gap.datasets")
library(gap)
library(gap.datasets)
data(CDKN)

head(CDKNlocus)

asplot (CDKNlocus, CDKNmap, CDKNgenes, best.pval = 5.4e-8,sf=c(3,6))
