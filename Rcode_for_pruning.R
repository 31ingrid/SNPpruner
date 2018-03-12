#read in the LD and missing data files.

dataMISS=read.table("/path/to/data/miss_stat.lmiss",header=TRUE)
dataLD=read.table("/path/to/data/plink.ld",header=TRUE)
mapfile=read.table("/path/to/data/data.map")
#Make a new data frame (called aligns) that includes only pairwise comparisons 
#with r2 greater than some threshold (defined as R2_threshold)
R2_threshold=0.8;

aligns= dataLD[which(dataLD$R2>R2_threshold),]
misses=seq(1,nrow(aligns),1)#sequence of numbers equal to the number of pairwise comparisons 

 
aligns$SNP_A2 <- gsub('_', '.', aligns$SNP_A)#converting the SNP names to a number because it makes it easier to search for patterns.
aligns$SNP_A2 <- as.numeric(aligns$SNP_A2)#all numbers are unique - replaced the underscore with a period.
aligns$SNP_B2 <- gsub('_', '.', aligns$SNP_B)#same for second column of SNPs
aligns$SNP_B2 <= as.numeric(aligns$SNP_B2)
mapfile$V5 <= gsub('_', '.', mapfile$V2)

#Convert SNPs in the missing data file to numbers as well (underscore to decimal)
 dataMISS$SNP1 <- gsub('_', '.',  dataMISS$SNP)#converting the SNP names to a number. 
 dataMISS$SNP1=as.numeric( dataMISS$SNP1)#replaced the underscore with a period.

lst3=vector();#list of all the rows checked. 
lst4=list();#a list of each cluster found. 
clustersize=vector()#vector of the size of each cluster

k=1 #initiate this as a counter
for (i in 1:length(misses)){  #first make sure your i is not in lst3   length(misses)
 print(i)
 if (i%in%lst3==FALSE){ #lst3 is all the rows you have already examined. This allows each cluster to be distinct because a row will not be examined twice.
  c=vector();d=vector(); #temp vectors to hold rows of aligns with SNPs in LD for a single cluster
  lst1=vector()#rows of aligns in the current cluster - start fresh for each cluster
  lst2=vector()#names of aligns in current cluster - start fresh each cluster
  lst1=i#lst1 is for storing row numbers of SNPs in LD
  
  lst2=c(aligns$SNP_A2[i],aligns$SNP_B2[i])#put your aligns$SNP_A2 and your aligns$SNP_B2 in lst2
  
  lst2_len=vector()#length of the vector containing linked snps. If it is still increasing then continue this loop until all linked loci have been found
  foundall=1; #a check to see whether more linked SNPs are found upon additional repetitions of this loop
  while(foundall>0){
  for(m in 1:2) {#start by doing this loop 2 times. do not stop until you have identified the entire cluster (until length lst2 stops increasing)
   for (j in 1:length(lst2)){
   if((lst2[j]%in%aligns$SNP_B2)){c=c(c,which(aligns$SNP_B2==lst2[j]))}#lst2 is your first set of linked markers
   if((lst2[j]%in%aligns$SNP_A2)){d=c(d,which(aligns$SNP_A2==lst2[j]))}#see if there are other SNPS that are linked to these
  }
   if (length(c)>0){lst1=c(lst1,c[!is.na(c)])}#add these indices to list1 if they are new not already in lst1
   if (length(d)>0){lst1=c(lst1,d[!is.na(d)])}
   
   #now lst1 should have all the c and ds.
   lst1=unique(lst1[lst1!=0])
   lst2=unique(c(lst2,aligns$SNP_A2[lst1],aligns$SNP_B2[lst1]))
   lst2_len[m]=length(lst2)
   print(length(lst2))
  }
  foundall=lst2_len[length(lst2_len)]-lst2_len[length(lst2_len)-1]
  }
  
   lst4[[k]]=lst2 #keep track of the list of linked loci
   lst3=unique(c(lst3,lst1)); 
   k=k+1;
   all=c(all,lst2)
   clustersize[k]=length(lst2)
 }
}

#now find the SNP in each list with the lowest F_Miss (least missing data)
remove_list=list(); #list of SNP in each cluster with less missing data
keeps=vector()#vector of SNPs you are keeping
removes=vector()#vector of SNPs to remove
nclust=length(clustersize[!is.na(clustersize)]) 

for(i in 1:nclust){
 lst2=lst4[[i]]
 fmiss=vector()
 for(j in 1:length(lst2)){
   fmiss[j]=  dataMISS$F_MISS[which( dataMISS$SNP1==lst2[j])]
   }
   removes=c(removes,(lst2[-which(fmiss==min(fmiss))[1]]))
   if((clustersize[!is.na(clustersize)][i])!=length(lst2[-(which(fmiss==min(fmiss))[1])])+length(lst2[which(fmiss==min(fmiss))[1]])){print(i);print("discrepancy")}
   keeps=c(keeps,lst2[which(fmiss==min(fmiss))[1]])
   remove_list[[i]]=lst2[-which(fmiss==min(fmiss))[1]]
}

#These are some checks to make sure everything ran smoothly  
length(unique(keeps))+length(unique(removes)) #these should sum up to
length(unique(c(aligns$SNP_A2,aligns$SNP_B2))) #the same as the number of loci you started with
#also the list of loci should be the same number of SNPs
length(unique(all)) 

max(clustersize[!is.na(clustersize)])#what is the largest group of linked loci
hist(clustersize[!is.na(clustersize)])#here is a plot of the cluster sizes, probably mostly all pairs of 2

#data outputs
#lst4 - this is a list of all clusters of linked SNPs at the specified threshold
#remove_list - this is a list of all loci for each cluster with more missing data that has been selected for removal.
#note -  you will see that remove_list contains all but 1 SNP in the lst4, the one with the most data was not included.

#removes is the vector of SNPs to prune, and it can be used as a "blacklist" in stacks software (note that remove_list was a list).
#set 4th column of the MAP file to negative for the removed SNPs.

for (i in 1:length(removes)){
 mapfile$V4[which(mapfile$V5==removes[i])]=-1
}   
write.table(mapfile[,1:4],"/path/to/data/datapruned.map",quote=FALSE,col.names=F,row.names=F)

#then test to make sure the new map file works. If the new map file has a different name you can read it in separately using this code.
#plink --ped data.ped --map datapruned.map