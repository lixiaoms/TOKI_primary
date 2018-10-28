args<-commandArgs(T)
length<-ncol(as.matrix(read.table(args[1])))
fimo<-read.csv(args[2],sep='\t')
bin<-as.numeric(args[3])*1000
centro<-as.numeric(strsplit(args[4],',')[[1]])
outdir<-args[5]

options(digits=3)

ctcfplus=rep(0,length)
ctcfminus=rep(0,length)
for (i in 1:(nrow(fimo)-3)){
    pos=ceiling((as.numeric(strsplit(as.character(fimo[i,3]),'[:,-]')[[1]][3])+1)/bin)
    strand=as.character(fimo[i,6])
    if (strand=='+'){ctcfplus[pos]=ctcfplus[pos]+1}
    if (strand=='-'){ctcfminus[pos]=ctcfminus[pos]+1}
}

ctcfplus=ctcfplus/mean(c(ctcfplus,ctcfminus))
ctcfminus=ctcfminus/mean(c(ctcfplus,ctcfminus))

pair=matrix(0,length,length)
M=matrix(0,length,length)

for (k in 1:(2*ceiling(600000/bin))){
for (i in 1:(length-k)){
    pair[i,i+k]=max(sum(ctcfplus[i:(i+k-1)]),ctcfminus[i]+pair[i+1,i+k])
    M[i,i+k]=1/(pair[i,i+k]+1)
}}

M[lower.tri(M)]=t(M)[lower.tri(t(M))]
ctcf=matrix(0,length,length)

if (centro[1]>1){
    ctcf[1:(centro[1]-1),1:(centro[1]-1)]=M[1:(centro[1]-1),1:(centro[1]-1)]
    }

if (centro[2]<length){
    ctcf[(centro[2]+1):length,(centro[2]+1):length]=M[(centro[2]+1):length,(centro[2]+1):length]
    }

write.table(ctcf,paste0(outdir,'/motif_matrix'),row.names=F,col.names=F)