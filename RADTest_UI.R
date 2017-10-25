## test of association by rarity of pathogenicity (tarp)

#install.packages("doParallel")
#install.packages("foreach")

library(foreach)
library(doParallel)

qq<-function(dat,txt){
  plot(-log10((length(dat):1)/length(dat)),-log10(sort(dat,decreasing = T)),xlab = "-log10 expected",ylab="-log10 observed ",main=txt)
  abline(0,1,lty=2,col="gray")
  abline(h=-log10(0.05),lty=2,col="gray")
  abline(h=-log10(0.05/length(dat)),col="gray",lty=2)
}


minpvalue<-function(pvec,N){
  
  s = sort(pvec) 
  k = length(s)
  
  ## get min p value
  p = 1
  len=0;
  fix=0;
  for (i in 1:k) {
    sim_pvalue<-poisson.test(i, N *s[i], alternative="greater")$p.value
    if(sim_pvalue<=p){
      p =  sim_pvalue
      len=i
      fix=s[i]
    }
  }
  #cat("min p-value =", p, "\n")
  vec<-c(p,len,fix)
  names(vec)<-c("minp","#casesIncutoff","cutoff")
  return(vec)
}

permutation<-function(permu,N,cores){
  totalCores = detectCores()
  ncores = min(cores, totalCores - 1)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  bigN = min(permu * N, 1e8)
  
  y = runif(bigN)  # under the null, p-value follows a uniform dist. 
  permut_pvalue<-c()
  permut_pvalue = foreach(j=1:permu, .combine=cbind) %dopar% {
    sp = sort(sample(y, N, replace=T))  # bootstrap sampling
    sp = sp[sp<0.05]
    k = length(sp)
    arr<-rep(1,N);
    if (k > 0) {
      for (i in 1:k) {
        ## keep all permutation result ##
        arr[i]<-poisson.test(i, N * sp[i], alternative="greater")$p.value 
      }
    }	
   
    min(arr)
  }
  stopCluster(cl)
  return(permut_pvalue)
}

groupDat<-function(dat,diseaseModel,N){
  gene_pvalue<-c()
  dat<-dat[which(dat$Dz.Model==diseaseModel),]
  if(dim(dat)[1]==0){return(NULL)}
  for (gene in unique(dat$GeneName)){
    pvec<-dat$popScore[which(dat$GeneName==gene)]
    pvec<-as.numeric(pvec)
    # print(length(pvec))
    minp<-minpvalue(pvec,N)
    gene_pvalue<-rbind(gene_pvalue,c(gene,minp))
    #set_minp<-c(set_minp,minp)
  }
  gene_pvalue<-data.frame(gene_pvalue,stringsAsFactors = F)
  names(gene_pvalue)<-c("GeneName","min_pvalue","cohortNumInThreshold","popscore_threshold")
  return(gene_pvalue)
}

fixedT_test<-function (dat,N,diseaseModel,T=0.001){
 # tests <- c()
  dat<-dat[which(dat$Dz.Model==diseaseModel),]
  rs<-c()
  for (g in unique(dat$GeneName)){
    o<-length(dat$ProbandName[which(dat$GeneName==g)])
    pt<-poisson.test(o,N*T);
    rs<-rbind(rs,c(g,o,pt$p.value))
  }
  rs<-data.frame(rs,stringsAsFactors = F)
  names(rs)<-c("Gene","Ncase","pvalue");
  return(rs)
}

##"PSAP/Data_Method/",
runRADT<-function(dat,prefix,N,permu,cores,freq=0.001){
  if(missing(freq)){freq=0.001;print("default population frequency is 0.1%")}
  Ftest<-fixedT_test(dat,N,diseaseModel = "DOM-het",T=freq)
  write.csv(Ftest,file = paste(prefix,"fixedThreshold_PAH.csv",sep=""),row.names = F,quote = F)
  png(paste(prefix,"Fixed_threshold.png",sep="_"),width=7,height = 7, units = 'in', res = 500)
  qq(as.numeric(Ftest$pvalue),txt = "Fixed threshold")
  dev.off()
  het_minp<-groupDat(dat,"DOM-het",N)
  homo_minp<-groupDat(dat,"REC-hom",N)
  chet_minp<-groupDat(dat,"REC-chet",N)
  
  ### get the min-pvalue list for every gene
  mat_permut<-permutation(permu,N,cores)
  ## assess the significance of min-p based on empirical distribution
  ## traverse the mat_permut and compared with min pvalue
  mat_permut<-sort(mat_permut)
  if(dim(het_minp)[1]>0){
    het_minp$empirical<-unlist(lapply(1:dim(het_minp)[1],function(x) length(which(mat_permut<=as.numeric(as.character(het_minp$min_pvalue[x]))))))+1
    het_minp$empirical<-het_minp$empirical/(permu+1)
  }
   if((!is.na(chet_minp) && !is.null(chet_minp)) &&dim(chet_minp)[1]>0){
     chet_minp$empirical<-unlist(lapply(as.numeric(chet_minp$min_pvalue),function(x) length(which(mat_permut<=x))))+1
     chet_minp$empirical<-chet_minp$empirical/(permu+1)
   }
   if((!is.na(homo_minp)&&!is.null(homo_minp) )&&dim(homo_minp)[1]>0){
     homo_minp$empirical<-unlist(lapply(as.numeric(homo_minp$min_pvalue),function(x) length(which(mat_permut<=x))))+1
     homo_minp$empirical<-homo_minp$empirical/(permu+1)
   }
  pdf(paste(prefix,"_QQ_empirical.pdf",sep=""))
  if(dim(het_minp)[1]>0){
    qq(het_minp$empirical,"DOM-het")
    
    write.csv(het_minp,file = paste(prefix,"_Dom_het.association.csv",sep=""),row.names = F)
  }
   if(exists(chet_minp) && dim(chet_minp)[1]>0){
     qq(chet_minp$empirical,"REC-Chet")
     write.csv(chet_minp,file = paste(prefix,"_REC_chet.association.csv",sep=""),row.names = F)
   }
   if(exists(homo_minp) && dim(homo_minp)[1]>0){
     qq(homo_minp$empirical,"REC-Hom")
     write.csv(homo_minp,file = paste(prefix,"_REC_hom.association.csv",sep=""),row.names = F)
   }
  dev.off()
  
  
}

tarpTest <- function(file,N, permu = 1e6, cores = 10, prefix="cohorts",folder,freq ) {
  popscores<-read.table(file,header=1,comment.char = "",stringsAsFactors = F,check.names = F,sep="\t",nrows = 2);
  if(dim(popscores)[2]<2){
    popscores<-read.csv(file,header=1,comment.char = "",stringsAsFactors = F,check.names = F);
  }else{
    popscores<-read.table(file,header=1,comment.char = "",stringsAsFactors = F,check.names = F,sep="\t");
  }
 
  runRADT(popscores,prefix,N,permu,cores,freq )
  
}

print("required format is Rscript tarp.R file N pop_freq permutation  cores   output_prefix")
## file: the psap result, it has to include columnName: GeneName, Dz.Model, ProbandName, popScore
## before runing the script, pleas make sure you have installed R library: doParallel and foreach 
args<-commandArgs(trailingOnly=T) ## args[6] = family - must be annovar annotated (avinput.hg19_multianno.txt) and have a separate header file for the vcf colums(header); args[7] = family pedigree file
file=args[1] ## the PSAP result, merged result of all samples together, required columnName: GeneName, Dz.Model, ProbandName, popScore
             ## supposed the qc have been applied, such as VQSR, inheritance violation, mappablity and so on.
N=as.numeric(args[2]) ## The total number of samples in the cohort, including the ones do not have PSAP output
freq=args[3] ## the population frequency (<=freq) to be tested
permu=as.numeric(args[4]) ## the simulation times.
cores=as.numeric(args[5]) ## it could work parrally with multiple cores, set then number of cores used in the test
prefix=args[6] ### the prefix name for output file
if(length(args)<6){print("required format : file N pop_freq permutation  cores   output_prefix") ;stop();}
tarpTest(file,N,permu,cores,prefix,folder,as.numeric(freq))
