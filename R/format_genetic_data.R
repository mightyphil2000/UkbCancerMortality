create_gs<-function(dat=NULL,gs_dat=NULL,BETA=NULL){
	Dat2<-dat
	# Dat3<-Dat
	# gs_dat<-gs_dat[!duplicated(gs_dat$rsid),]
	snps<-gs_dat$rsid
	# j<-1
	if(any(duplicated(snps))) stop("duplicate SNPs present")
	for(j in 1:length(snps)){
		print(j)
		# Beta<-gs_dat[,BETA][gs_dat$rsid == snps[j]]
		# Dat2[1,snps[4]]*Beta
		# (2-Dat3[1,snps[4]])*Beta*-1
		# sum(Dat2[,snps[j]][1:10]*Beta)
		# sum((2-Dat3[,snps[j]][1:10])*Beta*-1)
		# sum(Dat2[,snps[j]][1:10])
		Beta<-gs_dat[,BETA][gs_dat$rsid == snps[j]]
		Test<-sign(Beta)
		if(Test<0){
			Beta<-Beta*-1
			# Dat2[,snps[j]]<-Dat2[,snps[j]]*gs_dat[,BETA][gs_dat$rsid == snps[j]]
			Dat2[,snps[j]]<-(2-Dat2[,snps[j]])*Beta		
			# Beta<-gs_dat[,BETA][gs_dat$rsid == snps[j]]*-1
			# Dat3[,snps[j]]<-(2-Dat3[,snps[j]])*Beta	
		}
		if(Test>0){
			Dat2[,snps[j]]<-Dat2[,snps[j]]*Beta		
		}

		# check allele frequencies same in snptest and data files
		# sum(Dat2[,snps[4]])/(nrow(Dat2)*2)
		# snptest[snptest$rsid == snps[4],"alleleB_frequency"]
		# print(sum(Dat2[,snps[j]]))
		# head(Dat2)
	}
	if(length(snps)>1){
		GS<-rowSums(Dat2[,names(Dat2) %in% snps])			
		# GS2<-rowSums(Dat3[,names(Dat3) %in% snps])			
		# plot(GS1,GS)
	}
	if(length(snps)==1){
		GS<-Dat2[,names(Dat2) %in% snps]
	}
	return(GS)	
}

create_gs_unweighted<-function(dat=NULL,gs_dat=NULL,BETA=NULL){
	Dat2<-dat
	gs_dat<-gs_dat[!duplicated(gs_dat$rsid),]
	snps<-gs_dat$rsid
	# j<-4
	if(any(duplicated(snps))) stop("duplicate SNPs present")
	for(j in 1:length(snps)){
		print(j)
		Test<-sign(gs_dat[,BETA][gs_dat$rsid == snps[j]])
		# make effect alleles reflect COX increasing allele
		if(Test<0){
			Dat2[,snps[j]]<-2-Dat2[,snps[j]]	
		}
		# Dat2[,snps[j]]<-Dat2[,snps[j]]*sign(gs_dat[,BETA][gs_dat$rsid == snps[j]])
		# print(sum(Dat2[,snps[j]]))
		# head(Dat2)
	}
	if(length(snps)>1){
		GS<-rowSums(Dat2[,names(Dat2) %in% snps]) #number of COX increasing alleles			
	}
	if(length(snps)==1){
		GS<-Dat2[,names(Dat2) %in% snps]
	}
	return(GS)	
}


flip_strand<-function(Dat=NULL){
	Pos<-Dat$Effect.Allele.x!=Dat$Effect.Allele.y	
	strand1<-c("A","T","G","C")
	strand2<-c("T","A","C","G")
	# lnor.y<-Dat$lnor.y[Pos]*-1
	# Dat$lnor.y[Pos]<-lnor.y
	ea<-Dat$Effect.Allele.y[Pos]
	oa<-Dat$Other.Allele[Pos]
	Dat$Effect.Allele.y[Pos]<-strand2[match(ea,strand1)]
	Dat$Other.Allele[Pos]<-strand2[match(oa,strand1)]				
	return(Dat)
}
		
harmonise_effect_allele<-function(Dat_harmonise=NULL,ea.x=NULL,ea.y=NULL,oa.x=NULL,oa.y=NULL,eaf.x=NULL,eaf.y=NULL,b.y=NULL){
	Pos<-Dat_harmonise[,ea.x]!=Dat_harmonise[,ea.y]
	if(sum(Pos)==0){
		all(Dat_harmonise[,ea.x]==Dat_harmonise[,ea.y] & Dat_harmonise[,oa.x]==Dat_harmonise[,oa.y])
		print("effect and other alleles in datasets x and y are already the same")
		return(Dat_harmonise)
	}
	if(sum(Pos)!=0){
	# Dat_harmonise[,c(ea.x,ea.y,oa.x,oa.y)]
		# all(Dat_harmonise[,ea.x]!=Dat_harmonise[,ea.y] & Dat_harmonise[,oa.x]!=Dat_harmonise[,oa.y])
		# print("effect and other alleles in datasets x and y are already the same")
		beta.y<-Dat_harmonise[,b.y][Pos]*-1
		Dat_harmonise[,b.y][Pos]<-beta.y
		oa<-Dat_harmonise[,ea.y][Pos]
		ea<-Dat_harmonise[,oa.y][Pos]
		Dat_harmonise[,ea.y][Pos]<-ea
		Dat_harmonise[,oa.y][Pos]<-oa
		eaf<-1-Dat_harmonise[,eaf.y][Pos]
		Dat_harmonise[,eaf.y][Pos]<-eaf	
		Test<-all(Dat_harmonise[,ea.x]==Dat_harmonise[,ea.y] & Dat_harmonise[,oa.x]==Dat_harmonise[,oa.y])
		if(Test) print("effect and other alleles in datasets x and y are now the same")
		if(!Test) print("attempted harmonisation failed probably because the studies are on different strands")
		return(Dat_harmonise)
	}
}
# head(cpd[order(cpd$PVALUE),])
cpd_genetic_score<-function(dat=NULL){
# https://conservancy.umn.edu/bitstream/handle/11299/201564/README.txt?sequence=29&isAllowed=y
	cpd<-utils::read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cpd_no23andme_sig_clump_relaxed.txt",sep="\t",head=TRUE,stringsAsFactors=F)
	snptest<-utils::read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/ukb_cpd_snp_stats_all.txt",sep=" ",head=T,stringsAsFactors=F)
	cpd_snptest<-merge(snptest,cpd,by.x="rsid",by.y="RSID")	
	cpd_snptest<-harmonise_effect_allele(Dat_harmonise=cpd_snptest,ea.x="alleleB",ea.y="ALT",oa.x="alleleA",oa.y="REF",eaf.x="alleleB_frequency",eaf.y=NULL,b.y="BETA")            
	dat$GS_cpd<-create_gs(dat=dat,gs_dat=cpd_snptest,BETA="BETA")
	return(dat)
}

csi_genetic_score<-function(dat=NULL){
	csi<-utils::read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_clumped.txt",sep=" ",head=TRUE,stringsAsFactors=F)
	csi2<-csi[csi$clump_r2 == 0.3,]
	snptest<-utils::read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/ukb_snp_stats_all.txt",sep=" ",head=T,stringsAsFactors=F)
	# snptest2<-snptest[,c("rsid","alleleA","alleleB","alleleB_frequency")]
	# names(snptest)[names(snptest) %in% c("alleleA","alleleB","alleleB_frequency")]<-c("Other.allele")
	csi_snptest<-merge(snptest,csi2,by.y="SNP",by.x="rsid")	
	csi_snptest<-harmonise_effect_allele(Dat_harmonise=csi_snptest,ea.x="alleleB",ea.y="EFFECT_ALLELE",oa.x="alleleA",oa.y="OTHER_ALLELE",eaf.x="alleleB_frequency",eaf.y="EAF",b.y="BETA")
	# now effect allele for csi beta is same as coded allele in UKB 
	# Pal<-paste0(Test$EFFECT_ALLELE,Test$OTHER_ALLELE)	
	# Test[which(Pal %in% c("GC","CG","TA","AT")),c("EAF","alleleB_frequency")]
	dat$GS_csi<-create_gs(gs_dat=csi_snptest,BETA="BETA")		
	return(dat)
}


