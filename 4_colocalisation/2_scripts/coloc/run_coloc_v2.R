library(coloc)

args = commandArgs(trailingOnly=TRUE)
loc<-args[1]
trait1<-args[2]
trait2<-args[3]

print(paste(loc,trait1,trait2))

t1<-read.table(file=paste0(trait1,".LOC",loc,".txt"), header=T, sep="\t")
t1<-t1[!duplicated(t1$snp),]
head(t1)

t2<-read.table(file=paste0(trait2,".LOC",loc,".txt"), header=T, sep="\t")
t2<-t2[!duplicated(t2$snp),]
head(t2)

t1.list<-list(beta=t1$beta,
               varbeta=t1$varbeta,
               snp=t1$snp,
               position=t1$position,
               type=unique(t1$type),
               N=max(t1$N),
               MAF=t1$MAF)

t2.list<-list(beta=t2$beta,
               varbeta=t2$varbeta,
               snp=t2$snp,
               position=t2$position,
               type=unique(t2$type),
               N=max(t2$N),
               MAF=t2$MAF)

my.res <- coloc.abf(dataset1=t1.list,
                    dataset2=t2.list)

saveRDS(object=my.res, file=paste0(trait1,".",trait2,".LOC",loc,".rds"))

## #names(my.res)
