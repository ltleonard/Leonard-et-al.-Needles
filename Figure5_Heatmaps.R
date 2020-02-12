#This script converts phyloseq objects to ampvis2 data, readable by  commands we will be making to create the heatmaps for Figure 5. 
#Based on my PCOA analysis I determined that dates "one" and "four" were significant as a function of needle samples. I then proceeded with only these dates for further analysis.

#Combine OTU abundance table and taxonomy table from the 16S phyloseq object of date "four."
LS4 <- CB16LS_4
t_otu <- t(data.frame(otu_table(LS4)))
otutable_4ampvis2 <- data.frame(OTU = rownames(t_otu@.Data),
                                t_otu@.Data,
                                phyloseq::tax_table(LS4)@.Data,
                                check.names = FALSE
)

dim(otutable_4ampvis2)

otutable_4ampvis2[22] <- "NA"
colnames(otutable_4ampvis2)[22] <- "Species"
#Extract metadata from the phyloseq object:
metadata_4ampvis2 <- data.frame(phyloseq::sample_data(LS4), 
                                check.names = FALSE
)
#Load the data with amp_load:Warning message:
#Only 63 of 67 unique sample names match between metadata and otutable. The following unmatched samples have been removed:
LS4_ampvis <- amp_load(otutable_4ampvis2, metadata_4ampvis2)

#_______________________________________________________________________________________________
#Combine OTU abundance table and taxonomy table from the 16S phyloseq object of date "one."
LS1 <- CB16LS_1
t_otu <- t(data.frame(otu_table(LS1)))
otutable_4ampvis2 <- data.frame(OTU = rownames(t_otu@.Data),
                                t_otu@.Data,
                                phyloseq::tax_table(LS1)@.Data,
                                check.names = FALSE
)

dim(otutable_4ampvis2)

otutable_4ampvis2[20] <- "NA"
colnames(otutable_4ampvis2)[20] <- "Species"
#Extract metadata from the phyloseq object:
metadata_4ampvis2 <- data.frame(phyloseq::sample_data(LS1), 
                                check.names = FALSE
)
#Load the data with amp_load:Warning message:
#Only 63 of 67 unique sample names match between metadata and otutable. The following unmatched samples have been removed:
LS1_ampvis <- amp_load(otutable_4ampvis2, metadata_4ampvis2)
#_______________________________________________________________________________________________
#Combine OTU abundance table and taxonomy table from the 18S phyloseq object of date "four."
LS4_18 <- CB18LS_4
t_otu <- t(data.frame(otu_table(LS4_18)))
otutable_4ampvis2 <- data.frame(OTU = rownames(t_otu@.Data),
                                t_otu@.Data,
                                phyloseq::tax_table(LS4_18)@.Data,
                                check.names = FALSE
)

dim(otutable_4ampvis2)

otutable_4ampvis2[19] <- "NA"
colnames(otutable_4ampvis2)[19] <- "Species"
#Extract metadata from the phyloseq object:
metadata_4ampvis2 <- data.frame(phyloseq::sample_data(LS4_18), 
                                check.names = FALSE
)
#Load the data with amp_load:
LS4_18_ampvis <- amp_load(otutable_4ampvis2, metadata_4ampvis2)
#_______________________________________________________________________________________________
#Combine OTU abundance table and taxonomy table from the 18S phyloseq object of date "one."
LS1_18 <- CB18LS_1
t_otu <- t(data.frame(otu_table(LS5_18)))
otutable_4ampvis2 <- data.frame(OTU = rownames(t_otu@.Data),
                                t_otu@.Data,
                                phyloseq::tax_table(LS5_18)@.Data,
                                check.names = FALSE
)

dim(otutable_4ampvis2)

otutable_4ampvis2[21] <- "NA"
colnames(otutable_4ampvis2)[21] <- "Species"
#Extract metadata from the phyloseq object:
metadata_4ampvis2 <- data.frame(phyloseq::sample_data(LS5_18), 
                                check.names = FALSE
)
#Load the data with amp_load:
LS1_18_ampvis <- amp_load(otutable_4ampvis2, metadata_4ampvis2)

#_______________________________________________________________________________________________
#Plot the heatmap for date "one"
Heatmap1 <- amp_heatmap(data = LS1_ampvis ,
                      tax_aggregate = "Genus", 
                      tax_add = "Phylum", 
                      group_by = "TreeSpecies",
                      tax_show = 15, 
                      tax_empty = "remove", 
                      plot_values = TRUE, 
                      plot_values_size = 7,
                      plot_legendbreaks = c(0.1,1.0,10.0,50.0), 
                      max_abundance = 12, 
                      min_abundance = .1) +
  theme(axis.text.x = element_text(size = 20, color = "black", hjust = 0.5, angle = 0)) + 
  theme(axis.text.y = element_text(size = 20, color = "black", angle = 0))
  theme(strip.text.x = element_text(size = 20, color = 'Black')) + 
  theme(strip.background.x = element_rect(fill = "White")) + 
mapped$labels$x <- "Date of Sampling" #change x axis label
mapped$labels$fill <- "% Relative\nAbundance" #change legend label
print(Heatmap1)

#Plot the 16S heatmap for date "one"
Heatmap1 <- amp_heatmap(data = LS1_ampvis ,
                        tax_aggregate = "Genus", 
                        tax_add = "Phylum", 
                        group_by = "TreeSpecies",
                        tax_show = 15, 
                        tax_empty = "remove", 
                        plot_values = TRUE, 
                        plot_values_size = 7,
                        plot_legendbreaks = c(0.1,1.0,10.0,50.0), 
                        max_abundance = 12, 
                        min_abundance = .1) +
  theme(axis.text.x = element_text(size = 20, color = "black", hjust = 0.5, angle = 0)) + 
  theme(axis.text.y = element_text(size = 20, color = "black", angle = 0))
theme(strip.text.x = element_text(size = 20, color = 'Black')) + 
  theme(strip.background.x = element_rect(fill = "White")) + 
  mapped$labels$x <- "Date of Sampling" #change x axis label
mapped$labels$fill <- "% Relative\nAbundance" #change legend label
print(Heatmap1)
#_______________________________________________________________________________________________
#Plot the 16S heatmap for date "four"
Heatmap4 <- amp_heatmap(data = LS4_ampvis ,
                        tax_aggregate = "Genus", 
                        tax_add = "Phylum", 
                        group_by = "TreeSpecies",
                        tax_show = 15, 
                        tax_empty = "remove", 
                        plot_values = TRUE, 
                        plot_values_size = 7,
                        plot_legendbreaks = c(0.1,1.0,10.0,50.0), 
                        max_abundance = 12, 
                        min_abundance = .1) +
  theme(axis.text.x = element_text(size = 20, color = "black", hjust = 0.5, angle = 0)) + 
  theme(axis.text.y = element_text(size = 20, color = "black", angle = 0))
theme(strip.text.x = element_text(size = 20, color = 'Black')) + 
  theme(strip.background.x = element_rect(fill = "White")) + 
  mapped$labels$x <- "Date of Sampling" #change x axis label
mapped$labels$fill <- "% Relative\nAbundance" #change legend label
print(Heatmap4)
#_______________________________________________________________________________________________
#Plot the 18S heatmap for date "four"
Heatmap18S1 <- amp_heatmap(data = LS1_18_ampvis ,
                           tax_aggregate = "Phylum", 
                           #tax_add = "Family", 
                           group_by = "Needle",
                           tax_show = 15, 
                           tax_empty = "remove", 
                           plot_values = TRUE, 
                           plot_legendbreaks = c(0.1,1.0,10.0,50.0), 
                           #max_abundance = 10.5, 
                           min_abundance = .1) +
  theme(axis.text.x = element_text(size = 20, color = "black", hjust = 0.5, angle = 0)) + 
  theme(axis.text.y = element_text(size = 20, color = "black", angle = 0)) 
#theme(legend.text = element_text(size = 13, colour = 'black')) + 
#theme(legend.title = element_text(size = 15, color = 'black')) + 
#theme(strip.text.x = element_text(size = 15, color = 'Black')) + 
#theme(strip.background.x = element_rect(fill = "White")) + 
#theme(legend.key.size = unit(1.2, "cm"))+
#scale_x_discrete("Dates", labels = c("Control", "Lodgepole", "Spruce Combined")) 
mapped$labels$x <- "Date of Sampling" #change x axis label
mapped$labels$fill <- "% Relative\nAbundance" #change legend label
print(Heatmap18S1)
#_______________________________________________________________________________________________
#Plot the 18S heatmap for date "four"
Heatmap18S4 <- amp_heatmap(data = LS4_18_ampvis ,
                      tax_aggregate = "Phylum", 
                      #tax_add = "Family", 
                      group_by = "Needle",
                      tax_show = 15, 
                      tax_empty = "remove", 
                      plot_values = TRUE, 
                      plot_legendbreaks = c(0.1,1.0,10.0,50.0), 
                      #max_abundance = 10.5, 
                      min_abundance = .1) +
  theme(axis.text.x = element_text(size = 20, color = "black", hjust = 0.5, angle = 0)) + 
  theme(axis.text.y = element_text(size = 20, color = "black", angle = 0)) 
#theme(legend.text = element_text(size = 13, colour = 'black')) + 
#theme(legend.title = element_text(size = 15, color = 'black')) + 
#theme(strip.text.x = element_text(size = 15, color = 'Black')) + 
#theme(strip.background.x = element_rect(fill = "White")) + 
#theme(legend.key.size = unit(1.2, "cm"))+
#scale_x_discrete("Dates", labels = c("Control", "Lodgepole", "Spruce Combined")) 
mapped$labels$x <- "Date of Sampling" #change x axis label
mapped$labels$fill <- "% Relative\nAbundance" #change legend label
print(Heatmap18S4)
