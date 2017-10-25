
source("Pipeline/NA_script/R/qq_unif.R")

call_pvalue<-function(case,case_total,control,control_total,side){
  if(missing(side)){side="two.sided"}
  return(binom.test(case,control+case,p=case_total/(case_total+control_total),alternative = side,conf.level = 0.95))
}
call_enrichment<-function(case,control,len_case,len_control){
  return (length(case)*len_control/(len_case*length(control)))
  
}

COV_filter <- function(data){
  if(length(grep("GnomAD_Genome_cov",names(data)))>0){
    data<-data[which(as.numeric(data$GnomAD_Genome_cov10)>0.85),]
  }
  # if(length(grep("AC_PCGC",names(data)))>0 && length(grep("AC_xGen",names(data)))>0 && length(grep("AC_VCR",names(data)))>0){
  #   data <- data[which(2*data$AC_PCGC/data$AN_PCGC >0.8 & 2*data$AC_xGen/data$AN_xGen >0.8 & 2*(data$AC_VCR/data$AN_VCR) >0.8  &data$GnomAD_Genome_cov10>0.8 ),]
  # }
  # if(length(grep("AC_Control",names(data)))>0 && length(grep("AC_VU",names(data)))>0){
  #   data <- data[which(2*data$AC_Control/data$AN_Control >0.8 & 2*data$AC_VU/data$AN_VU >0.8 ),]
  # }
  return(data)
}

VQSR_filter <- function(data){
  
  #  if(length(grep("VQSLOD",names(data)))>0){
  #    data <- data[which(data$VQSLOD > -0.5),]
  #  }
  
  if(length(grep("FILTER",names(data)))>0){
    filter_ignore <- c("VQSRTrancheSNP99.90to100.00",
                       "VQSRTrancheSNP99.80to99.90","VQSRTrancheSNP99.70to99.80",
                       "VQSRTrancheSNP99.60to99.70",#"VQSRTrancheSNP99.50to99.60",
                       #"VQSRTrancheINDEL99.50to99.60",
                       "VQSRTrancheINDEL99.60to99.70",                 
                       "VQSRTrancheINDEL99.70to99.80",  
                       "VQSRTrancheINDEL99.80to99.90", "VQSRTrancheINDEL99.90to100.00")
    
    data <- data[which(!data$FILTER %in%  filter_ignore |data$FILTER=="."),]
  }
  
  return(data)
}


map_filter<-function(data){
  index<-grep("mappability",names(data),ignore.case = T)
  if(length(index)>0){ 
    data <- data[which(as.numeric(data[,index])==1),]
  }
  
  if(length(grep("genomicSuperDups",names(data)))>0){
    index<-grep("Score",data$genomicSuperDups)
    
    dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-|=|;"))[2])))>0.95)]
    if(length(dup_indexs)>0){
      data <- data[-dup_indexs,]
    }
  }
  return (data)
}

basic_RADtest<-function(pvec,N){
  
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
    #  print ( poisson.test(i, N * s[i], alternative="greater")$p.value)
  }
  #cat("min p-value =", p, "\n")
  vec<-c(p,len,fix)
  names(vec)<-c("minp","#casesIncutoff","cutoff")
  return(vec)
}



format_AAchange<-function(dat){
  h<-which(names(dat)=="AAChange")
  if(length(h)>0){
    dat$Transcript<-unlist(lapply(dat$AAChange,FUN = function(x)unlist(strsplit(x,split = ":"))[2] ))
    dat$Exon<-unlist(lapply(dat$AAChange,FUN = function(x)unlist(strsplit(x,split = ":"))[3] ))
    
    dat$Nucleiotide_Change<-unlist(lapply(dat$AAChange,FUN = function(x)unlist(strsplit(x,split = ":"))[4] ))
  dat$Protein_Change<-unlist(lapply(dat$AAChange,FUN = function(x)unlist(strsplit(x,split = ":"))[5] ))

  
 
  dat$Mutation_Type<-""
  dat$Mutation_Type[which(dat$VarClass=="stopgain")]<-"Nonsense_Mutation"
  dat$Mutation_Type[which(dat$VarClass=="frameshiftinsertion")]<-"Frame_Shift_ins"
  dat$Mutation_Type[which(dat$VarClass=="nonframeshiftdeletion")]<-"In_Frame_del"
  dat$Mutation_Type[which(dat$VarClass=="nonframeshiftinsertion")]<-"In_Frame_ins"
  dat$Mutation_Type[which(dat$VarClass=="frameshiftdeletion")]<-"Frame_Shift_del"
  dat$Mutation_Type[which(dat$VarClass=="nonsynonymousSNV")]<-"Missense_Mutation"
  dat$Mutation_Type[which(dat$VarFunc=="splicing")]<-"Splice_Site"
  }else{
    h<-which(names(dat)=="AAChange.refGene")
    if(length(h)>0){
      dat$Protein_Change<-unlist(lapply(dat$AAChange.refGene,FUN = function(x)unlist(strsplit(x,split = ":"))[5] ))
      dat$Mutation_Type<-""
      dat$Mutation_Type[which(dat$ExonicFunc.refGene=="stopgain")]<-"Nonsense_Mutation"
      dat$Mutation_Type[which(dat$ExonicFunc.refGene=="frameshift_deletion")]<-"Frame_Shift_del"
      dat$Mutation_Type[which(dat$ExonicFunc.refGene=="frameshift_insertion")]<-"Frame_Shift_ins"
      dat$Mutation_Type[which(dat$ExonicFunc.refGene=="nonframeshift_deletion")]<-"In_Frame_del"
      dat$Mutation_Type[which(dat$ExonicFunc.refGene=="nonframeshift_insertion")]<-"In_Frame_ins"
      dat$Mutation_Type[which(dat$ExonicFunc.refGene=="nonsynonymous_SNV")]<-"Missense_Mutation"
      dat$Mutation_Type[which(dat$Func.refGene=="splicing")]<-"Splice_Site"
    } 
    
  }
  return(dat)
  
}
get_lgd<-function(dat){
  index1<-grep("stop",dat$VarClass)
  index2<-grep("^frame",dat$VarClass)
  index3<-grep("splicing",dat$VarFunc)
  return(unique(c(index1,index2,index3)))
}
get_indel<-function(dat){
  len_ref<-unlist(lapply(dat$REF,FUN = function(x) nchar(x)))
  len_alt<-unlist(lapply(dat$ALT,FUN=function(x)nchar(x)))
  if(length(len_ref)!=length(len_alt)){stop("Error: the number of ref and alt are not equal!"); }
  return(which(len_ref!=len_alt))
}
get_snv<-function(dat){
  len_ref<-unlist(lapply(dat$REF,FUN = function(x) nchar(x)))
  len_alt<-unlist(lapply(dat$ALT,FUN=function(x)nchar(x)))
  if(length(len_ref)!=length(len_alt)){stop("Error: the number of ref and alt are not equal!"); }
  return(which(len_ref==len_alt))
}



formatFreq<-function(data){
  if(length(grep("1KGfreq",names(data)))>0){
  data$`1KG.afr.freq`[which(is.na(as.numeric(data$`1KG.afr.freq`)))]=0
  data$`1KG.eas.freq`[which(is.na(as.numeric(data$`1KG.eas.freq`)))]=0
  data$`1KG.amr.freq`[which(is.na(as.numeric(data$`1KG.amr.freq`)))]=0
  data$`1KG.eur.freq`[which(is.na(as.numeric(data$`1KG.eur.freq`)))]=0
  data$`1KG.sas.freq`[which(is.na(as.numeric(data$`1KG.sas.freq`)))]=0
  data$`1KGfreq`[which(is.na(as.numeric(data$`1KGfreq`)))]=0
  }
  if(length(grep("ESPfreq",names(data)))>0){
    data$ESPfreq[which(is.na(as.numeric(data$ESPfreq)))]=0
    data$ESP.aa.freq[which(is.na(as.numeric(data$ESP.aa.freq)))]=0
    data$ESP.ea.freq[which(is.na(as.numeric(data$ESP.ea.freq)))]=0
  }
  if(length(grep("ExACfreq",names(data)))>0){
    data$ExACfreq[which(is.na(as.numeric(data$ExACfreq)))]=0
    data$ExAC.eas.freq[which(is.na(as.numeric(data$ExAC.eas.freq)))]=0
    data$ExAC.afr.freq[which(is.na(as.numeric(data$ExAC.afr.freq)))]=0
    data$ExAC.amr.freq[which(is.na(as.numeric(data$ExAC.amr.freq)))]=0
    data$ExAC.nfe.freq[which(is.na(as.numeric(data$ExAC.nfe.freq)))]=0
    data$ExAC.sas.freq[which(is.na(as.numeric(data$ExAC.sas.freq)))]=0
    data$ExAC.fin.freq[which(is.na(as.numeric(data$ExAC.fin.freq)))]=0
    data$ExAC.oth.freq[which(is.na(as.numeric(data$ExAC.oth.freq)))]=0
  }
  if(length(grep("gnomAD_Exome_AF",names(data)))>0){
    data$gnomAD_Exome_AF[which(is.na(as.numeric(data$gnomAD_Exome_AF)))]=0
    data$gnomAD_Exome_AF_ASJ[which(is.na(as.numeric(data$gnomAD_Exome_AF_ASJ)))]=0
    data$gnomAD_Exome_AF_AMR[which(is.na(as.numeric(data$gnomAD_Exome_AF_AMR)))]=0
    data$gnomAD_Exome_AF_AFR[which(is.na(as.numeric(data$gnomAD_Exome_AF_AFR)))]=0
    data$gnomAD_Exome_AF_EAS[which(is.na(as.numeric(data$gnomAD_Exome_AF_EAS)))]=0
    data$gnomAD_Exome_AF_SAS[which(is.na(as.numeric(data$gnomAD_Exome_AF_SAS)))]=0
    data$gnomAD_Exome_AF_OTH[which(is.na(as.numeric(data$gnomAD_Exome_AF_OTH)))]=0
    data$gnomAD_Exome_AF_NFE[which(is.na(as.numeric(data$gnomAD_Exome_AF_NFE)))]=0
    data$gnomAD_Exome_AF_FIN[which(is.na(as.numeric(data$gnomAD_Exome_AF_FIN)))]=0
  }
  if(length(grep("gnomAD_Genome_AF",names(data)))>0){
    data$gnomAD_Genome_AF[which(is.na(as.numeric(data$gnomAD_Genome_AF)))]=0
    data$gnomAD_Genome_AF_AFR[which(is.na(as.numeric(data$gnomAD_Genome_AF_AFR)))]=0
    data$gnomAD_Genome_AF_ASJ[which(is.na(as.numeric(data$gnomAD_Genome_AF_ASJ)))]=0
    data$gnomAD_Genome_AF_AMR[which(is.na(as.numeric(data$gnomAD_Genome_AF_AMR)))]=0
    data$gnomAD_Genome_AF_AFR[which(is.na(as.numeric(data$gnomAD_Genome_AF_AFR)))]=0
    data$gnomAD_Genome_AF_EAS[which(is.na(as.numeric(data$gnomAD_Genome_AF_EAS)))]=0
    data$gnomAD_Genome_AF_OTH[which(is.na(as.numeric(data$gnomAD_Genome_AF_OTH)))]=0
    data$gnomAD_Genome_AF_NFE[which(is.na(as.numeric(data$gnomAD_Genome_AF_NFE)))]=0
    data$gnomAD_Genome_AF_FIN[which(is.na(as.numeric(data$gnomAD_Genome_AF_FIN)))]=0
    data$gnomAD_Genome_AF_POPMAX[which(is.na(as.numeric(data$gnomAD_Genome_AF_POPMAX)))]=0
  }
  return(data)
}


formatFreq_new<-function(data){
  if(length(grep("1000g2015aug_all",names(data)))>0){
    data$`1000g2015aug_afr`[which(is.na(as.numeric(data$`1000g2015aug_afr`)))]=0
    data$`1000g2015aug_eas`[which(is.na(as.numeric(data$`1000g2015aug_eas`)))]=0
    data$`1000g2015aug_amr`[which(is.na(as.numeric(data$`1000g2015aug_amr`)))]=0
    data$`1000g2015aug_eur`[which(is.na(as.numeric(data$`1000g2015aug_eur`)))]=0
    data$`1000g2015aug_sas`[which(is.na(as.numeric(data$`1000g2015aug_sas`)))]=0
    data$`1000g2015aug_all`[which(is.na(as.numeric(data$`1000g2015aug_all`)))]=0
  }
  if(length(grep("ExAC_ALL",names(data)))>0){
    data$ExAC_ALL[which(is.na(as.numeric(data$ExAC_ALL)))]=0
    data$ExAC_EAS[which(is.na(as.numeric(data$ExAC_EAS)))]=0
    data$ExAC_AFR[which(is.na(as.numeric(data$ExAC_AFR)))]=0
    data$ExAC_AMR[which(is.na(as.numeric(data$ExAC_AMR)))]=0
    data$ExAC_NFE[which(is.na(as.numeric(data$ExAC_NFE)))]=0
    data$ExAC_SAS[which(is.na(as.numeric(data$ExAC_SAS)))]=0
    data$ExAC_FIN[which(is.na(as.numeric(data$ExAC_FIN)))]=0
    data$ExAC_OTH[which(is.na(as.numeric(data$ExAC_OTH)))]=0
  }
  if(length(grep("esp6500siv2_all",names(data)))>0){
    data$esp6500siv2_all[which(is.na(as.numeric(data$esp6500siv2_all)))]=0
    data$esp6500siv2_aa[which(is.na(as.numeric(data$esp6500siv2_aa)))]=0
    data$esp6500siv2_ea[which(is.na(as.numeric(data$esp6500siv2_ea)))]=0
  }
  if(length(grep("gnomAD_exome_ALL",names(data)))>0){
    data$gnomAD_exome_ALL[which(is.na(as.numeric(data$gnomAD_exome_ALL)))]=0
    data$gnomAD_exome_ASJ[which(is.na(as.numeric(data$gnomAD_exome_ASJ)))]=0
    data$gnomAD_exome_AMR[which(is.na(as.numeric(data$gnomAD_exome_AMR)))]=0
    data$gnomAD_exome_AFR[which(is.na(as.numeric(data$gnomAD_exome_AFR)))]=0
    data$gnomAD_exome_EAS[which(is.na(as.numeric(data$gnomAD_exome_EAS)))]=0
    data$gnomAD_exome_FIN[which(is.na(as.numeric(data$gnomAD_exome_FIN)))]=0
    data$gnomAD_exome_SAS[which(is.na(as.numeric(data$gnomAD_exome_SAS)))]=0
    data$gnomAD_exome_OTH[which(is.na(as.numeric(data$gnomAD_exome_OTH)))]=0
    data$gnomAD_exome_NFE[which(is.na(as.numeric(data$gnomAD_exome_NFE)))]=0
  }
  if(length(grep("gnomAD_genome_ALL",names(data)))>0){
    data$gnomAD_genome_ALL[which(is.na(as.numeric(data$gnomAD_genome_ALL)))]=0
    data$gnomAD_genome_AFR[which(is.na(as.numeric(data$gnomAD_genome_AFR)))]=0
    data$gnomAD_genome_ASJ[which(is.na(as.numeric(data$gnomAD_genome_ASJ)))]=0
    data$gnomAD_genome_AMR[which(is.na(as.numeric(data$gnomAD_genome_AMR)))]=0
    data$gnomAD_genome_FIN[which(is.na(as.numeric(data$gnomAD_genome_FIN)))]=0
    data$gnomAD_genome_AFR[which(is.na(as.numeric(data$gnomAD_genome_AFR)))]=0
    data$gnomAD_genome_EAS[which(is.na(as.numeric(data$gnomAD_genome_EAS)))]=0
    data$gnomAD_genome_OTH[which(is.na(as.numeric(data$gnomAD_genome_OTH)))]=0
    data$gnomAD_genome_NFE[which(is.na(as.numeric(data$gnomAD_genome_NFE)))]=0
   # data$gnomAD_genome_POPMAX[which(is.na(as.numeric(data$gnomAD_genome_POPMAX)))]=0
  }
  return(data)
}



filter_allfreq <- function(data,freq_avg,freq_max){
#  freq <- 0.001
#  freq2 <- 0.001
  l<-length(grep("ExACfreq",names(data)))
  
  if(l<1) return( data);
  
  data <- data[which(na.pass(as.numeric(data$ExACfreq)< freq_avg)
                     &na.pass(as.numeric(data$ExAC.amr.freq)< freq_max)
                     &as.numeric(data$ExAC.afr.freq)< freq_max
                     &as.numeric(data$ExAC.nfe.freq)< freq_max
                     &as.numeric(data$ExAC.sas.freq)< freq_max
                     &as.numeric(data$ExAC.eas.freq)< freq_max
                     &as.numeric(data$ExAC.oth.freq)< freq_max
                     #  & as.numeric(data$`1KGfreq`) < freq2
                     #  &as.numeric(data$ESPfreq)< freq2
                     #  &as.numeric(data$gnomAD_Genome_AF)< freq2
  ),]
#  data <- data[which( as.numeric(data$AC)< 25
 #                     &as.numeric(data$AB)>0.2
 # ),]
  if(length(grep("Mappability",names(data)))>0){ 
    data <- data[which(data$Mappability==1),]
  }
  return(data)
}


load_lof_CADDscore<-function(gene,dataset){
  ### /home/local/ARCS/nz2274/Application/psap/lookups/full.lof.pCADD.gencodeV19.allsites.txt.gz
  record<-0
  len<-length(which(dataset$V1==gene))
  if(len>0)  record<-dataset[which(dataset$V1==gene),3]
  return(record)  
}

get_popscore<- function(gene,score,dataset){
  ##/home/local/ARCS/nz2274/Application/psap/lookups/full.het.pCADD.gencodeV19.allsites.txt.gz
  popscore<-1;
  scale<-seq(0,70,0.05)
  len<-length(which(dataset$V1==gene))
  if(len>0){
   vec<-as.numeric(dataset[which(dataset$V1==gene),2:dim(dataset)[2]])
    popscore<-vec[findInterval(score,scale)+1]
  }
  return(popscore)
}

filter_allfreq_new <- function(data,freq_avg,freq_max){
  
  data <- data[which(na.pass(as.numeric(data$ExAC_ALL)< freq_avg)
                     &na.pass(as.numeric(data$ExAC_AMR)< freq_max)
                     &as.numeric(data$ExAC_AFR)< freq_max
                     &as.numeric(data$ExAC_NFE)< freq_max
                     &as.numeric(data$ExAC_FIN)< freq_max
                     &as.numeric(data$ExAC_SAS)< freq_max
                     &as.numeric(data$ExAC_EAS)< freq_max
                     &as.numeric(data$ExAC_OTH)< freq_max
                     &as.numeric(data$gnomAD_genome_EAS)<freq_max
                     &as.numeric(data$gnomAD_genome_NFE)<freq_max
                     &as.numeric(data$gnomAD_genome_FIN)<freq_max
                     &as.numeric(data$gnomAD_genome_OTH)<freq_max
                     &as.numeric(data$gnomAD_genome_ASJ)<freq_max
                     &as.numeric(data$gnomAD_genome_AMR)<freq_max
                     &as.numeric(data$gnomAD_genome_ALL)<freq_avg
                     &as.numeric(data$gnomAD_exome_ALL)<freq_avg
                     &as.numeric(data$gnomAD_exome_EAS)<freq_max
                     &as.numeric(data$gnomAD_exome_NFE)<freq_max
                     &as.numeric(data$gnomAD_exome_FIN)<freq_max
                     &as.numeric(data$gnomAD_exome_OTH)<freq_max
                     &as.numeric(data$gnomAD_exome_ASJ)<freq_max
                     &as.numeric(data$gnomAD_exome_AMR)<freq_max
                     #  & as.numeric(data$`1KGfreq`) < freq2
                     #  &as.numeric(data$ESPfreq)< freq2
                     #  &as.numeric(data$gnomAD_Genome_AF)< freq2
  ),]
 # data <- data[which( as.numeric(data$AC)< 25
#                      &as.numeric(data$AB)>0.2
#  ),]
  if(length(grep("Mappability",names(data)))>0){ 
    data <- data[which(data$Mappability==1),]
  }
  if(length(grep("genomicSuperDups",names(data)))>0){
    index<-grep("Score",data$genomicSuperDups)
    as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-"))[2])))
    dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-"))[2])))>0.9)]
    if(length(dup_indexs)>0){
      data <- data[-dup_indexs,]
    }
  }
  return(data)
}



load_dataset<-function(){
  resources="Resources/"
 # fexac<-paste(resources,"ExAC_r03_0316_z_pli_rec_null_data.bed",sep="")
  fheart<-paste(resources,"GeneExpression/mousebrain_Heart.csv",sep="")
  flung<-paste(resources,"GeneExpression/Lung_rank_RNAseq Asselin-Labat-GPL13112_human.csv",sep="")
  fdiaphram<-paste(resources,"GeneExpression/diaphragm_rank.csv",sep="")
 # exac<-read.table(fexac,stringsAsFactors = F,header = 1,check.names = F,comment.char = "")
  heart_exp<-read.csv(fheart,sep = ",",header = 1,check.names = F)
  lung_exp<-read.csv(flung,sep=",",header = 1,check.names = F)
  diaphragm_exp<-read.csv(fdiaphram,sep=",",header = 1,check.names = F)
  
 
  return(list(heart_exp=heart_exp,lung_exp=lung_exp,diaph_exp=diaphragm_exp))
}
load_exac<-function(){
  resources="Resources/"
  fexac<-paste(resources,"ExAC_r03_0316_z_pli_rec_null_data.bed.gz",sep="")
  exac<-read.table(fexac,stringsAsFactors = F,header = 1,check.names = F,comment.char = "")
  return (exac)
}
load_roadmap_RKPM<-function(){
  resources="~/server/Resources/"
 
  file<-paste(resources,"GeneExpression/57epigenomes.RPKM.pc.csv.gene.txt",sep="")
  exp<-read.table(file,header = 1,sep = "\t",check.names = F,stringsAsFactors = F,row.names = NULL)
  #flung<-paste(resources,"GeneExpression/Lung_rank_RNAseq Asselin-Labat-GPL13112_human.csv",sep="")
  fdiaphram<-paste(resources,"GeneExpression/diaphragm_rank.csv",sep="")
  aorta<-"E065"
  lung<-"E096"
  heart_exp<-exp[,c("GeneName",aorta)]
  lung_exp<-exp[,c("GeneName",lung)]
  diaphragm_exp<-read.csv(fdiaphram,sep=",",header = 1,check.names = F,strip.white = T)
  
  return(list(heart_exp=heart_exp,lung_exp=lung_exp,diaph_exp=diaphragm_exp))
}



format_consensus<-function(data){
  nms<-names(data)
  nms[which(nms=="Gene.refGene")]="GeneName"
  nms[which(nms=="ExonicFunc.refGene")]="VarClass"
  nms[which(nms=="Func.refGene")]="VarFunc"
  nms[which(nms=="ExAC_ALL")]="ExACfreq"
  nms[which(nms=="ExAC_AMR")]="ExAC.amr.freq"
  nms[which(nms=="ExAC_EAS")]="ExAC.eas.freq"
  nms[which(nms=="ExAC_NFE")]="ExAC.nfe.freq"
  nms[which(nms=="ExAC_AFR")]="ExAC.afr.freq"
  nms[which(nms=="ExAC_FIN")]="ExAC.fin.freq"
  nms[which(nms=="1000g2015aug_afr")]="1KG.afr.freq"
  nms[which(nms=="1000g2015aug_all")]="1KG.all.freq"
  nms[which(nms=="1000g2015aug_amr")]="1KG.amr.freq"
  nms[which(nms=="1000g2015aug_eas")]="1KG.eas.freq"
  nms[which(nms=="1000g2015aug_eur")]="1KG.eur.freq"
  nms[which(nms=="1000g2015aug_sas")]="1KG.sas.freq"
  nms[which(nms=="AAChange.refGene")]="AAChange"
#  nms[which(nms=="CADD13_phred")]="CADDphred"
  nms[which(nms=="CADD_phred")]="CADDphred"
  nms[which(nms=="CADDraw")]="CADD13_raw"
 # nms[which(nms=="VarFunc")]="Func.refGene"
  nms[which(nms=="MetaLR_pred")]="MetaLRprd"
  nms[which(nms=="MetaSVM_pred")]="MetaSVMprd"
  nms[which(nms=="MetaSVMscr")]="MetaSVM_score"
  nms[which(nms=="MutAprd")]="MutationAssessor_pred"
  nms[which(nms=="MutAscr")]="MutationAssessor_score"
  nms[which(nms=="MutTscr")]="MutationTaster_score" 
  nms[which(nms=="MutTprd")]="MutationTaster_pred" 
  nms[which(nms=="Polyphen2_HVAR_pred")]="PP2.hvar.prd" 
  nms[which(nms=="Polyphen2_HDIV_pred")]="PP2.hdiv.prd" 
  #nms[which(nms=="proband")]="proband_GT"
  if(length(which(nms=="proband_GT"))==0){
    data$proband_GT<-data$proband
  }
  data$proband<-data$ProbandName
  #nms[which(nms=="ProbandName")]="proband"
  
  names(data)<-nms
  return(data)
  
}

gene_exp_diaphragm<-function(data,lib){
  for(id in 1:dim(data)[1]){
    gid<-which(names(data)=="Gene.refGene"|names(data)=="GeneName")
    index<-grep(data[id,gid],lib$gene)[1]
    if(length(index)>0){
      data$diaphragm_rank[id]=lib$rank[index]
      data$diaphragm_E115[id]=lib$E11.5[index]
    }else{
      data$diaphragm_rank[id]="."
      data$diaphragm_E115[id]="."
      
    }
  }
  return(data)
}

gene_exp_Heart<-function(data,lib){
  for(id in 1:dim(data)[1]){
    gid<-which(names(data)=="Gene.refGene"|names(data)=="GeneName")
    index<-grep(data[id,gid],lib$human.External.Gene.Name)[1]
    if(length(index)>0){
      data$HEART_EXP[id]=lib$e14.5_rank[index]
      #  data$diaphragm_E115[id]=lib$E11.5[index]
    }else{
      data$HEART_EXP[id]="."
      #data$diaphragm_E115[id]="."
      
    }
  }
  print("Heart Expression, the higher the more related")
  return(data)
}



gene_exp_Lung<-function(data,lib){ ## reversed
  for(id in 1:dim(data)[1]){
    gid<-which(names(data)=="Gene.refGene"|names(data)=="GeneName")[1]
    
    index<-grep(data[id,gid],lib$human_gene)[1]
    if(length(index)>0){
      data$LUNG_EXP[id]=lib$`Control-Stroma rank`[index]
      #  data$diaphragm_E115[id]=lib$E11.5[index]
    }else{
      data$LUNG_EXP[id]="."
      #data$diaphragm_E115[id]="."
      
    }
  }
  print("Lung Expression, the higher the more related")
  return(data)
}


gene_zscore<-function(data,lib){
  for(id in 1:dim(data)[1]){
    gid<-which(names(data)=="Gene.refGene"|names(data)=="GeneName")[1]
    index<-grep(data[id,gid],lib$gene)[1]
    if(length(index)>0){
      data$pLI[id]=lib$pLI[index]
      data$mis_z[id]=lib$mis_z[index]
      data$lof_z[id]=lib$lof_z[index]
    }else{
      data$mis_z[id]="."
      data$pLI[id]="."
      data$lof_z[id]="."
    }
  }
  return(data)
}




seperate_proband_info<-function(data){
  a<-data$proband
  if(grepl(pattern = "(",x = a,fixed = T)==F||grepl(pattern = ":",x = a,fixed = T)==F){ return(data)}
  probs<-mapply(a,FUN = function(x) unlist(strsplit(unlist(strsplit (x,split=":"))[1],split = "(",fixed = T))[1] )
  gts<-mapply(a,FUN = function(x) unlist(strsplit(unlist(strsplit (x,split=":"))[1],split = "(",fixed = T))[2] )
  ads<-mapply(a,FUN = function(x) unlist(strsplit (x,split=":"))[2])
  dps<-mapply(a,FUN = function(x) as.numeric(unlist(strsplit (x,split=":"))[3]))
  gqs<-mapply(a,FUN = function(x) unlist(strsplit (x,split=":"))[4])
  
  ab<-c();
  
  for(i in 1:length(gts)){
    gt=gts[i];
    
    alt=max(as.numeric(unlist(strsplit(gt,split = "/",fixed = T))),na.rm = T)
    
    ad=as.numeric(unlist(strsplit(ads[i],split = ","))[alt+1])
    ab<-c(ab,ad/dps[i])
    
  }
  max_ad<-mapply(ads,FUN = function(x) max(unlist(strsplit(x,split = ",")) ))
  data$ProbandName<-probs
  data$GT<-gts
  data$AD<-ads
  data$Ind_DP<-dps
  data$AB<-ab
  data$Max_AD<-max_ad
  data$GQ<-gqs
  return (data)
}


addPhenotype2<-function(data,pheno){
  for(id in 1:length(data$ProbandName)){
    key=data$ProbandName[id];
    item<-which(pheno$ID==key)
    nms<-names(pheno)
    nms<-nms[which(nms!="")]
    if(length(item)>0){
      for(m in nms){
        data[id,as.character(m)]<-pheno[item,as.character(m)]
      }
    }
  }
  return (list(data=data,items=nms))
}

addPhenotype<-function(data,pheno){
  for(id in 1:length(data$ProbandName)){
    key=unlist(strsplit(data$ProbandName[id],split="_"))[2];
    item<-which(pheno$ID==key)
    nms<-names(pheno)
    nms<-nms[which(nms!="")]
    if(length(item)>0){
      for(m in nms){
        data[id,as.character(m)]<-pheno[item,as.character(m)]
      }
    }
  }
  return (list(data=data,items=nms))
}

output_candidates_v0<-function(denovo6,file2,items){
  if(missing(items)){items<-c()}
  denovo6$VarType="."
  denovo6$VarType[which(denovo6$MetaSVMprd=="D")]<-"DMIS"
  
  denovo6$VarType[which(denovo6$CADDphred>15 &denovo6$PP2.hdiv.prd=="D")]<-"DMIS"
  denovo6$VarType[which(denovo6$CADDphred>15 &denovo6$MetaSVMprd=="T")]<-"PDMIS"
  
  denovo6$VarType[which(denovo6$VarClass %in%c("stopgain","stoploss","frameshiftinsertion","frameshiftdeletion"))]<-"LGD"
  
  denovo6$VarType[which(denovo6$CADDphred<15 &denovo6$VarClass=="nonsynonymousSNV")]<-"MIS"
  
  aa<-mapply(denovo6$AAChange,FUN = function(x) unlist(strsplit(x,":"))[5])
  denovo6$AAC<-aa
  Nucleo<-mapply(denovo6$AAChange,FUN = function(x) unlist(strsplit(x,":"))[4])
  denovo6$Nucleotide<-Nucleo
  
  trans<-mapply(denovo6$AAChange,FUN = function(x) unlist(strsplit(x,":"))[2])
  denovo6$Transcript<-trans
  outsets<-c("proband","parents","CHROM","ID","POS","REF","ALT","FILTER","VarType","GeneName","VarClass","Transcript","Nucleotide","AAC","CADDphred","MetaSVMprd",
             "PP2.hdiv.prd","ExACfreq","1KGfreq","ESPfreq","75bp_mappability","AD","Ind_DP","GQ","IGV",
             "pLI","mis_z","LUNG_EXP","HEART_EXP","state","Gender","AGE","DIAG_AGE","pheno","group","population","AF","AC","AN",items)
  outset<-intersect(outsets,names(denovo6))
  write.csv(denovo6[,outset],
              file =file2,quote = T,row.names = F)
  
  return(denovo6)
}




output_candidates_v1<-function(denovo6,file2,items){
  if(missing(items)){items<-c()}
  denovo6$VarType="."
  denovo6$VarType[which(denovo6$MetaSVM_pred=="D")]<-"DMIS"
  
  denovo6$VarType[which(denovo6$CADD13_PHRED>15 &denovo6$Polyphen2_HVAR_pred=="D")]<-"DMIS"
  denovo6$VarType[which(denovo6$CADD13_PHRED>15 &denovo6$MetaSVM_pred=="T")]<-"PDMIS"
  
  denovo6$VarType[which(denovo6$ExonicFunc.refGene %in%c("stopgain","stoploss","frameshift_insertion","frameshift_deletion"))]<-"LGD"
  
  denovo6$VarType[which(denovo6$CADD13_PHRED<15 &denovo6$ExonicFunc.refGene=="nonsynonymous_SNV")]<-"MIS"
  
  aa<-mapply(denovo6$AAChange.refGene,FUN = function(x) unlist(strsplit(x,":"))[5])
  denovo6$AAC<-aa
  Nucleo<-mapply(denovo6$AAChange.refGene,FUN = function(x) unlist(strsplit(x,":"))[4])
  denovo6$Nucleotide<-Nucleo
  
  trans<-mapply(denovo6$AAChange.refGene,FUN = function(x) unlist(strsplit(x,":"))[2])
  denovo6$Transcript<-trans
  outsets<-c("ProbandName","VarType","ExonicFunc.refGene","Gene.refGene","Transcript","Nucleotide","AAC","pLI","lof_z","mis_z","HEART_EXP","LUNG_EXP",
             "MetaSVM_pred", "Polyphen2_HDIV_pred","CADD13_PHRED",
             "ExAC_ALL","DIAG_AGE","pheno","SIFT_pred","MCAP","ExAC_AFR","ExAC_NFE","ExAC_AMR","ExAC_EAS","ExAC_OTH",
             "1000g2015aug_all","1000g2015aug_eas","1000g2015aug_eur",
             "AF","AC","AB","Ind_DP","QUAL","proband","state","Gender","AGE","CHROM","ID","POS","REF","ALT","FILTER","parents",items)
  outsets2<-c("proband","GeneName","VarClass","CADDphred","MetaSVMprd",
             "PP2.hdiv.prd","ExACfreq","1KGfreq","ESPfreq","Mappability","pheno","group","population")
  outsets<-c(outsets,outsets2)
  outset<-intersect(outsets,names(denovo6))
  write.csv(denovo6[,outset],
            file =file2,quote = T,row.names = F)
  
  return(denovo6)
}

panel.qqconf<-function(data,maxp, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
  n=length(data)
  require(grid)
  conf.points =n
#  conf.points = min(conf.points, n-1);
  mpts<-matrix(nrow=conf.points*2, ncol=2)
  for(i in seq(from=1, to=conf.points)) {
    mpts[i,1]<- -log10((i-.5)/n) -log10(maxp) #data[i]  #
    mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i)) -log10(maxp)
    mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n) -log10(maxp)
    mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i)) -log10(maxp)
  }
  grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
}


qq_uniform<-function(data,min_dist,max_dist){
#  data<-c(data,runif(n=2490000,min=0,max=1))
  
  #x<--log10(sort(runif(n = length(unique(data)),min=min_dist,max=max_dist),decreasing = F))
 # y<--log10(sort(unique(as.numeric(data)),decreasing = F))
  min_dist<-min(data)
  max_dist<-max(data)
  x<-as.numeric(format(-log10(sort(runif(n = length((data)),min=min_dist,max=max_dist),decreasing = F)),digits = 4))
  y<-as.numeric(format(-log10(sort((as.numeric(data)),decreasing = F)),digits = 4))
  plot(x,y,ylab="-log10 (observed)",xlab="-log10(expected)")
  abline(0,1,col="gray",lty=2)
  # xyplot(y~x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
  #        prepanel=prepanel, scales=list(axs="i"), pch=pch,ylim=c(min(y)-0.5,max(y)+0.5),xlim=c(min(x)-0.5,max(x)+1),
  #        panel = function(x, y) {
  #          if (draw.conf) {
  #            panel.qqconf(x,max_dist, conf.points=conf.points, 
  #                         conf.col=conf.col, conf.alpha=conf.alpha)
  #          };
  #          panel.xyplot(x,y);
  #          panel.abline(0,1);
  #        }
  # )
}

