
#Differential abundance testing via "Deseq" using unrarified data from the Compiled2018.rmd file. 
#Filter out 18S, 16S, Mitochondria and Chloroplasts.
des_18S = subset_taxa(unrare_18S, Kingdom == "Eukaryota")
des_16S = subset_taxa(unrare_16S, (Kingdom != "Eukaryota")&(Class != "Chloroplast")&(Family != "Mitochondria"))
#Bin 16S samples by date collections, One= August 2017, Two= October 2017, Three=May 2018, Four=July 2018, Five=October 2018.
des_CB16LS_4= subset_samples(des_16S, Date=="d_Four")
des_CB16LS_3= subset_samples(des_16S, Date=="c_Three")
des_CB16LS_5= subset_samples(des_16S, Date=="e_Five")
des_CB16LS_2= subset_samples(des_16S, Date=="b_Two")
des_CB16LS_1= subset_samples(des_16S, Date=="a_One")
#Bin 18S samples by date collections. 
des_CB18LS_4= subset_samples(des_18S, Date=="d_Four")
des_CB18LS_3= subset_samples(des_18S, Date=="c_Three")
des_CB18LS_5= subset_samples(des_18S, Date=="e_Five")
des_CB18LS_2= subset_samples(des_18S, Date=="b_Two")
des_CB18LS_1= subset_samples(des_18S, Date=="a_One")

####16S, compare Control VS Lodgepole
LS1_cont_lodge = subset_samples(des_CB16LS_4, (Needle != "Red")&(Needle != "Green"))
LS1_cont_lodge_Genus <- tax_glom(LS1_cont_lodge, taxrank="Genus")
LS1_cont_lodge_Genus_rel = transform_sample_counts(LS1_cont_lodge_Genus, function(x) 100 * x/sum(x))
dat <- psmelt(LS1_cont_lodge_Genus_rel)
write.csv(dat, file='des1_LS4_Genus_rel_lodge.csv')
#Make DESeqDatasets using Genus level. 
ddsL.fam =phyloseq_to_deseq2(LS1_cont_lodge_Genus , ~TreeSpecies)
dds = DESeq(ddsL.fam, test="Wald", fitType="parametric") #holds normalized counts for each sample
res = results(dds, cooksCutoff = FALSE) # holds log2fold change and base means
ddsm = assay(dds) #convert ddsL.pt into a datatable
ddsdt = as.data.frame(ddsm)
restab = cbind(as(res, "data.frame"), as(tax_table(LS1_cont_lodge_Genus )[rownames(res), ], "matrix"))
restabCounts = cbind(as(restab, "data.frame"), as(ddsdt [rownames(restab),],"data.frame"))
#write out deseq2 results to csv file (P value <0.05)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(LS1_cont_lodge_Genus)[rownames(sigtab), ], "matrix"))
write.csv(sigtab, file = "/FilePath/LS4_cont_lodge_Genus.csv")
#________________________________________________________________________________________
###16S, compare Control VS Red Spruce
LS1_cont_red = subset_samples(des_CB16LS_4, (Needle != "Lodge")&(Needle != "Green"))
LS1_cont_red_Genus <- tax_glom(LS1_cont_red, taxrank="Genus")
LS1_cont_red_Genus_rel = transform_sample_counts(LS1_cont_red_Genus, function(x) 100 * x/sum(x))
dat <- psmelt(LS1_cont_red_Genus_rel)
write.csv(dat, file='des1_LS4_Genus_rel_red.csv')
#Make DESeqDatasets
ddsL.fam =phyloseq_to_deseq2(LS1_cont_red_Genus , ~TreeSpecies)
dds = DESeq(ddsL.fam, test="Wald", fitType="parametric") #holds normalized counts for each sample
res = results(dds, cooksCutoff = FALSE) # holds log2fold change and base means
ddsm = assay(dds) #convert ddsL.pt into a datatable
ddsdt = as.data.frame(ddsm)
restab = cbind(as(res, "data.frame"), as(tax_table(LS1_cont_red_Genus )[rownames(res), ], "matrix"))
restabCounts = cbind(as(restab, "data.frame"), as(ddsdt [rownames(restab),],"data.frame"))
#write out deseq2 results to csv file (P value <0.05)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(LS1_cont_red_Genus)[rownames(sigtab), ], "matrix"))
write.csv(sigtab, file = "/FilePath/LS4_cont_red_Genus.csv")
#________________________________________________________________________________________
###16S, compare Control VS Green Spruce
LS16_cont_green = subset_samples(des_CB16LS_4, (Needle != "Lodge")&(Needle != "Red"))
LS16_cont_green_Genus <- tax_glom(LS16_cont_green, taxrank="Genus")
LS16_cont_green_Genus_rel = transform_sample_counts(LS16_cont_green_Genus, function(x) 100 * x/sum(x))
dat <- psmelt(LS16_cont_green_Genus_rel)
write.csv(dat, file='des16_LS4_Genus_rel_green.csv')
#Make DESeqDatasets
ddsL.fam =phyloseq_to_deseq2(LS16_cont_green_Genus , ~TreeSpecies)
dds = DESeq(ddsL.fam, test="Wald", fitType="parametric") #holds normalized counts for each sample
res = results(dds, cooksCutoff = FALSE) # holds log2fold change and base means
ddsm = assay(dds) #convert ddsL.pt into a datatable
ddsdt = as.data.frame(ddsm)
restab = cbind(as(res, "data.frame"), as(tax_table(LS16_cont_green_Genus )[rownames(res), ], "matrix"))
restabCounts = cbind(as(restab, "data.frame"), as(ddsdt [rownames(restab),],"data.frame"))
#write out deseq2 results to csv file (P value <0.05)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(LS16_cont_green_Genus)[rownames(sigtab), ], "matrix"))
write.csv(sigtab, file = "/FilePath/LS4_cont_green_Genus.csv")
#________________________________________________________________________________________
##16S, compare Lodgepole VS Spruce- Tree species only
LS16_Lodge_Spruce = subset_samples(des_CB16LS_1, (TreeSpecies != "Control"))
LS16_Lodge_Spruce_Genus <- tax_glom(LS16_Lodge_Spruce, taxrank="Genus")
LS16_Lodge_Spruce_Genus_rel = transform_sample_counts(LS16_Lodge_Spruce_Genus, function(x) 100 * x/sum(x))
dat <- psmelt(LS16_Lodge_Spruce_Genus_rel)
write.csv(dat, file='des16_LS1_Genus_rel_Lodge_Spruce.csv')
#Make DESeqDatasets
ddsL.fam =phyloseq_to_deseq2(LS16_Lodge_Spruce_Genus , ~TreeSpecies)
dds = DESeq(ddsL.fam, test="Wald", fitType="parametric") #holds normalized counts for each sample
res = results(dds, cooksCutoff = FALSE) # holds log2fold change and base means
ddsm = assay(dds) #convert ddsL.pt into a datatable
ddsdt = as.data.frame(ddsm)
restab = cbind(as(res, "data.frame"), as(tax_table(LS16_Lodge_Spruce_Genus )[rownames(res), ], "matrix"))
restabCounts = cbind(as(restab, "data.frame"), as(ddsdt [rownames(restab),],"data.frame"))
#write out deseq2 results to csv file (P value <0.05)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(LS16_Lodge_Spruce_Genus)[rownames(sigtab), ], "matrix"))
write.csv(sigtab, file = "/FilePath/LS1_Lodge_Spruce_Genus.csv")
#________________________________________________________________________________________
###16S Red Spruce Versus Green Spruce
LS16_red_green = subset_samples(des_CB16LS_4, (Needle != "Lodge")&(Needle != "Control")&(Needle != "Shade"))
LS16_red_green_Genus <- tax_glom(LS16_red_green, taxrank="Genus")
LS16_red_green_Genus_rel = transform_sample_counts(LS16_red_green_Genus, function(x) 100 * x/sum(x))
dat <- psmelt(LS16_red_green_Genus_rel)
write.csv(dat, file='des16_LS1_Genus_rel_green.csv')
#Make DESeqDatasets
ddsL.fam =phyloseq_to_deseq2(LS16_red_green_Genus , ~Needle)
dds = DESeq(ddsL.fam, test="Wald", fitType="parametric") #holds normalized counts for each sample
res = results(dds, cooksCutoff = FALSE) # holds log2fold change and base means
ddsm = assay(dds) #convert ddsL.pt into a datatable
ddsdt = as.data.frame(ddsm)
restab = cbind(as(res, "data.frame"), as(tax_table(LS16_red_green_Genus )[rownames(res), ], "matrix"))
restabCounts = cbind(as(restab, "data.frame"), as(ddsdt [rownames(restab),],"data.frame"))
#write out deseq2 results to csv file (P value <0.05)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(LS16_red_green_Genus)[rownames(sigtab), ], "matrix"))
write.csv(sigtab, file = "/FilePath/LS1_red_green_Genus.csv")
