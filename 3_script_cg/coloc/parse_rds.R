args = commandArgs(trailingOnly=TRUE)
loc<-args[1]
trait1<-args[2]
trait2<-args[3]

rds<-readRDS(paste0(trait1,".",trait2,".LOC",loc,".rds"))

#print(rds$priors)

write.table(x=rds$results,
            file=paste0(paste0(trait1,".",trait2,".LOC",loc),".results.txt"),
            quote=F, sep="\t",
            row.names=F, col.names=T)

write.table(x=rds$summary,
            file=paste0(paste0(trait1,".",trait2,".LOC",loc),".summary.txt"),
            quote=F, sep="\t",
            row.names=F, col.names=T)
