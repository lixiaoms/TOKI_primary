args<-commandArgs(T)
hic<-as.matrix(read.table(args[1]))
max_inter<-ceiling(600/as.numeric(args[2]))
centro<-as.numeric(strsplit(args[3],',')[[1]])
outdir<-args[4]
floyd<-as.matrix(read.table(paste0(outdir,'/distance_matrix')))
exp<-matrix(0,ncol(hic),ncol(hic))
ctcf<-matrix(0,ncol(hic),ncol(hic))

if (file.exists(paste0(outdir,'/expression_matrix'))){
    exp<-as.matrix(read.table(paste0(outdir,'/expression_matrix')))}
if (file.exists(paste0(outdir,'/motif_matrix'))){
    ctcf<-as.matrix(read.table(paste0(outdir,'/motif_matrix')))}

suppressMessages(suppressWarnings(library(pROC)))
options(digits=4)

bound_dist=function(window,bound){
    dist=c()
    for (i in window){
        dist=c(dist,min(abs(bound-i)))
    }
return(dist)}

AUC=function(window,boundary,bound){
    index=rep(0,length(window))
    index[which(window %in% boundary)]=1
    data=data.frame(label=index,dist=bound_dist(window,bound))
    return(auc(roc(as.factor(data$label), data$dist)))
}

IS=function(A,slide=max_inter,delta=3){
    signal_list=c()
    for (i in (slide+1):(ncol(A)-slide)){
      signal_list=c(signal_list,mean(A[(i-slide):(i-1),(i+1):(i+slide)])/(mean(A[(i-slide):(i-1),(i-slide):(i-1)])+mean(A[(i+1):(i+slide),(i+1):(i+slide)])))
    }
    if (length(which(is.na(signal_list)))+length(which(is.infinite(signal_list)))>0){signal_list=log(signal_list/mean(signal_list[-c(which(is.na(signal_list)),which(is.infinite(signal_list)))]),2)}
    if (length(which(is.na(signal_list)))+length(which(is.infinite(signal_list)))==0){signal_list=log(signal_list/mean(signal_list),2)}
    return(signal_list)}

IS_delta=function(A,slide=max_inter,delta=3){
    delta_list=c()
    signal_list=IS(A,slide,delta)
    for (i in (delta+1):(length(signal_list)-delta)){
        delta_list=c(delta_list,mean(signal_list[(i-delta):(i-1)])-mean(signal_list[(i+1):(i+delta)]))
    }
    return(delta_list)
}

valley=function(x,y,delta=3){
    delete=c(which(is.na(x)),which(is.na(y))+delta,which(is.na(y))+delta-1,length(y)+delta+3)
    pos=c()
    for (i in 1:(length(y)-1)){
    if (min(abs(i+delta-delete))>0 && y[i]>=0 && y[i+1]<0 && x[i+delta]<0){pos=c(pos,i)}}
    return(pos)}

TOKI=function(S,A,M,E){
    ref=3+max_inter+valley(IS(S),IS_delta(S))
    if (length(ref)==0){return(c())}
    t=sum(sum(S))
    s=1
    l=1
    m=1
    e=1
    for (i in 0:(2*max_inter)){
        s=s+sum(S[col(S)-row(S)==i])
        l=l+sum(A[col(A)-row(A)==i])
        m=m+sum(M[col(M)-row(M)==i])
        e=e+sum(E[col(E)-row(E)==i])}
    A=A*s/(l*t)
    M=M*s/(m*t)
    E=E*s/(e*t)
    auc=0
    for (j in c(1/3,1/2,1,2,3)){for (k in c(1/3,1/2,1,2,3)){
        auc1=AUC(1:ncol(S),ref,3+max_inter+valley(IS(A*l/(m+e)+j*M+k*E),IS_delta(A*l/(m+e)+j*M+k*E)))
    if (auc1 > auc){
        b=j
        c=k
        auc=auc1
    }}}
    possi=(A*l/(m+e)+b*M+c*E)/(l/(m+e)+b+c)
    bayes=matrix(0,ncol(possi),ncol(possi))
    for (i in 1:ncol(bayes)){
    for (j in 0:min(2*max_inter,ncol(bayes)-i)){
        possi[i,i+j]=max(0,(possi[i,i+j]-t*(S[i,i+j]/t-possi[i,i+j])**2+((possi[i,i+j]-t*(S[i,i+j]/t-possi[i,i+j])**2)**2+4*S[i,i+j]*(S[i,i+j]/t-possi[i,i+j])**2)**(1/2))/2)
    }}
    possi[lower.tri(possi)]=t(possi)[lower.tri(t(possi))]
    return(3+max_inter+valley(IS(possi),IS_delta(possi)))
}

TAD=c()
length=ncol(hic)

if (centro[1]>1){
    TAD=c(TOKI(hic[1:(centro[1]-1),1:(centro[1]-1)],floyd[1:(centro[1]-1),1:(centro[1]-1)],ctcf[1:(centro[1]-1),1:(centro[1]-1)],exp[1:(centro[1]-1),1:(centro[1]-1)]),centro[1])
    }

if (centro[2]<length){
    TAD=c(TAD,centro[2]+1,centro[2]+TOKI(hic[(centro[2]+1):length,(centro[2]+1):length],floyd[(centro[2]+1):length,(centro[2]+1):length],ctcf[(centro[2]+1):length,(centro[2]+1):length],exp[(centro[2]+1):length,(centro[2]+1):length]))
    }

write.table(TAD,paste0(outdir,'/TAD'),row.names=F,col.names=F)

1