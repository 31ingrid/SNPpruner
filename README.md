# SNPpruner
This R code complements PLINK output to identify groups of linked SNPs and retain a single SNP from each group with the least missing data.

Although it is common to find software that will provide estimate of LD among pairs of SNPs, it may also be important to identify groups of linked SNPs, which is not commonly provided. This could occur if SNP X is linked to SNP Y, and SNP Y is also linked to SNP Z. If you simply remove one SNP from each pair, you might be removing SNPs that you do not intend, or leaving ones you should be pruning. If SNPs X, Y, and Z are all in linkage disequilibrium, it would make sense to discard all but the SNP with the least missing data. Knowing that SNPs X, Y, and Z are all linked may be important if you find that certain SNPs appear to be under selection, or are particularly interesting for some reason. You can go back to the list of linked SNPs and try aligning SNPs that were pruned to the reference genome. This functino may be most relevant for non model organisms.

SNPpruner.R is a function written in R that takes output from PLINK and provides groups of linked SNP loci over a specified r2 value. It provides a list of SNPs with the least missing data from each group, and a list of SNPs to prune.  These can be used to create a blacklist in STACKS. It will also provide an updated .map file with a -1 in the 4th column so that PLINK will ignore those SNPs that are pruned.

For this analysis, use PLINK to generate a list of SNPs in linkage disequilibrium, and to generate the amount of missing data for each SNP. PLINK is recommended for its computing speed and the SNPpruner function uses PLINK formatted output files.

Once your data is loaded into PLINK and you create a BED file (http://zzz.bwh.harvard.edu/plink/tutorial.shtml), PLINK will analyze all pairs of linked SNPs in your data using this command:

plink --bfile datafilename --r2 #(output is plink.ld)

It will also provide the amount of missing data for each SNP

plink --bfile datafilename --missing --out miss_stat  #(output is miss_stat.miss)

#To run the R code, first download it and open it in R. Then set your working directory to the location where your plink.ld, miss_stat.miss, and data.map files are. The data.map file is the initial file you read into PLINK with information for each SNP locus.

Read in the SNPpruner function, then run by simply type:

SNPpruner(R2_threshold)

Note: R2_threshold is the R2 value {0<R2_threshold<1} above which you want to filter out SNPs.

Output Files: The function will save several files to your working directory: 

LinkedSNPS.txt - this is a list of all clusters of linked SNPs at the specified threshold. You can also read the list in R using this command: 

readRDS("LinkedSNPs.RData").

 SNPS_to_retain.txt - this is a vector of all loci to keep (the SNP in each cluster with the least missing data). This should probably not be used as a whitelist in STACKS because PLINK does not include pairwise comparisons with very low R2. 
 SNPS_to_prune.txt – this is a vector of all SNPs to prune (only the SNP with the most data was retained for each cluster of linked SNPs). It can be used as a "blacklist" in stacks software.
 
If the new map file has a different name than the original you can read it in separately using this code.
#plink --ped data.ped --map datapruned.map 

Here are a few references that discuss SNP pruning and use 0.8 as the R2 threshold.

Jasonowicz, A.J., Goetz, F.W., Goetz, G.W. and Nichols, K.M., 2016. Love the one you’re with: genomic evidence of panmixia in the sablefish (Anoplopoma fimbria). Canadian Journal of Fisheries and Aquatic Sciences, 74(3), pp.377-387.

Larson, W.A., Seeb, L.W., Everett, M.V., Waples, R.K., Templin, W.D. and Seeb, J.E., 2014. Genotyping by sequencing resolves shallow population structure to inform conservation of Chinook salmon (Oncorhynchus tshawytscha). Evolutionary Applications, 7(3), pp.355-369.

Three example data files are included here (data.map, plink.ld, miss_stat.miss).
