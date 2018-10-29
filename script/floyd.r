args<-commandArgs(T)
hic<-as.matrix(read.table(args[1]))
max_inter<-2*ceiling(600/as.numeric(args[2]))
centro<-as.numeric(strsplit(args[3],',')[[1]])
ncore<-as.numeric(args[4])
outdir<-args[5]

suppressMessages(suppressWarnings(library(foreach)))
suppressMessages(suppressWarnings(library(doParallel)))
options(digits=3)

Floyd_single=function(W,max_inter){
D=1/W
for (i in 1:ncol(W)){D[i,i]=0}
a=length(which(W[col(W)==row(W)+1]>0))/(ncol(W)-1)
for (i in 1:(ncol(W)-1)){
                        D[i+1,i]=min(1/a,D[i+1,i])
                       D[i,i+1]=min(1/a,D[i,i+1])}
for (k in 1:ncol(W)){for (i in max(k-max_inter,1):min(k+max_inter,ncol(W))){for (j in max(k-max_inter,1):min(k+max_inter,ncol(W))){D[i,j]=min(D[i,j],D[i,k]+D[j,k])}}}
A=1/D
for (i in 1:ncol(W)){A[i,i]=W[i,i]}
return(A)}

Floyd_multi=function(W,ncore,max_inter){
D=1/W
for (i in 1:ncol(W)){D[i,i]=0}
a=length(which(W[col(W)==row(W)+1]>0))/(ncol(W)-1)
for (i in 1:(ncol(W)-1)){
                        D[i+1,i]=min(1/a,D[i+1,i])
                       D[i,i+1]=min(1/a,D[i,i+1])}
cl=makeCluster(ncore)
registerDoParallel(cl, cores=ncore)
size=ceiling((ncol(W)-max_inter)/ncore)+max_inter
system.time(
  D1 <- foreach(i=1:ncore,.combine='cbind') %dopar%
  { 
     res=D[((i-1)*(size-max_inter)+1):min((i-1)*(size-max_inter)+size,ncol(W)),((i-1)*(size-max_inter)+1):min((i-1)*(size-max_inter)+size,ncol(W))]
     for (k in (1:ncol(res))){for (i in max(k-max_inter,1):min(k+max_inter,ncol(res))){for (j in max(k-max_inter,1):min(k+max_inter,ncol(res))){res[i,j]=min(res[i,j],res[i,k]+res[j,k])}}}
     res=rbind(res,matrix(0,nrow=size-nrow(res),ncol=ncol(res)))
  res
  }
)
stopImplicitCluster()
stopCluster(cl)
for (i in 1:ncore){
    for (j in 1:min(size,ncol(D1)-(i-1)*size)){
        for (k in j:min(j+max_inter,size,ncol(D1)-(i-1)*size)){
            D[(i-1)*(size-max_inter)+j,(i-1)*(size-max_inter)+k]=min(D[(i-1)*(size-max_inter)+j,(i-1)*(size-max_inter)+k],D1[k,(i-1)*size+j])
            D[(i-1)*(size-max_inter)+k,(i-1)*(size-max_inter)+j]=min(D[(i-1)*(size-max_inter)+j,(i-1)*(size-max_inter)+k],D1[k,(i-1)*size+j])
        }
    }
}
A=1/D
for (i in 1:ncol(W)){A[i,i]=W[i,i]}
return(A)}

Floyd=function(W,ncore,max_inter){
    if (ncol(W)>2000 & ncore>2){
        return(Floyd_multi(W,ncore,max_inter))
    }else{
        return(Floyd_single(W,max_inter))}}

if (centro[1]>1){
    W1=hic[1:(centro[1]-1),1:(centro[1]-1)]
    hic[1:(centro[1]-1),1:(centro[1]-1)]=Floyd(W1,ncore,max_inter)
    }

if (centro[2]<ncol(hic)){
    W2=hic[(centro[2]+1):ncol(hic),(centro[2]+1):ncol(hic)]
    hic[(centro[2]+1):ncol(hic),(centro[2]+1):ncol(hic)]=Floyd(W2,ncore,max_inter)
    }

write.table(hic,paste0(outdir,'/distance_matrix'),row.names=F,col.names=F)
    