## module add R/4.2.2-foss-2022b

#install.packages("devtools")
#library(devtools)
#install_github("jrs95/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = TRUE)
library(hyprcoloc)

args = commandArgs(trailingOnly=TRUE)
loci<-paste0("LOC",args[1])
directory<-args[2]

message(loci)
message(directory)

betas<-read.table(file = paste0(directory,loci,".beta.txt"),
                   header=T,
                  row.names = 1)

#head(betas)

SEs<-read.table(file = paste0(directory,loci,".se.txt"),
                        header=T,
                        row.names = 1)

#head(SEs)

traits<-colnames(betas)
snps<-rownames(betas)

# Colocalisation analyses
results <- hyprcoloc(as.matrix(betas), as.matrix(SEs), 
                     trait.names=traits, snp.id=snps)

#print(results)

write.table(x=results$results,file=paste0(directory,loci,".out.txt"), quote=F, sep="\t", row.names=F)
