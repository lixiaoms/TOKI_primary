args<-commandArgs(T)
length<-ncol(as.matrix(read.table(args[1])))
tab<-read.csv(args[2],sep='\t',header=FALSE)
bin<-as.numeric(args[3])
centro<-as.numeric(strsplit(args[4],',')[[1]])
outdir<-args[5]

options(digits=3)


inter=bin/5
exp=rep(0,length)

for (j in 1:length){
    mean=mean(log(tab[,4]+1,10))+2.5
    exp[j]=mean(tab[(inter*(j-1)+1):min(inter*j,nrow(tab)),4])+1
    exp[j]=min(max(0,log(exp,10)-mean),1)}

E=matrix(0,length,length)
    for(i in 1:(length-1)){
    for(j in (i+1):min(length,i+2*ceiling(600/bin))){
    E[i,j]=(1-(exp[i]-exp[j])**2)*(min(exp[i],exp[j])+1)*(j-i)**(-1)
}}

E[lower.tri(E)]=t(E)[lower.tri(t(E))]
exp=matrix(0,length,length)

if (centro[1]>1){
    exp[1:(centro[1]-1),1:(centro[1]-1)]=E[1:(centro[1]-1),1:(centro[1]-1)]
    }

if (centro[2]<length){
    exp[(centro[2]+1):length,(centro[2]+1):length]=E[(centro[2]+1):length,(centro[2]+1):length]
    }

write.table(exp,paste0(outdir,'/expression_matrix'),row.names=F,col.names=F)