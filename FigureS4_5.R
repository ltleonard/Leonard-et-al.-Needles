#Using the "Figure5_Heatmaps" code to convert the phyloeq object to ampvis objects, boxplots can now be created
#Again, focusing only on two dates that were shown to have significant clustering in the PCOA plots

#16S boxplot of date "one"- August 2017
boxplot1 <- amp_boxplot(data=LS1_ampvis,
                        group_by = "Needle",
                        tax_empty = "remove",
                        tax_aggregate = "Family",
                        tax_show = 15,
                        tax_add = "Phylum"
)+
  scale_color_manual(values=c("darkgreen","navajowhite4","red4"), name="Level of Impact")+  
  scale_y_continuous(limit=c(0,14), breaks=c(0,2,4,6,8,10,12,14))+
  #theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 18, color = "black", angle = 0))+ 
  theme(axis.text.x = element_text(size = 18, color = "black", hjust = 0.5, angle = 0))+ 
  theme(
    panel.grid.major = element_line(size = 0.5, linetype = "solid"))
print(boxplot1)
####################################################################
#16S boxplot of date "four"- July 2018
boxplot4 <- amp_boxplot(data=LS4_ampvis,
                        group_by = "Needle",
                        tax_empty = "remove",
                        tax_aggregate = "Family",
                        tax_show = 15,
                        tax_add = "Phylum"
)+
  scale_color_manual(values=c("darkgreen","navajowhite4","red4"), name="Level of Impact")+  
  scale_y_continuous(limit=c(0,14), breaks=c(0,2,4,6,8,10,12,14))+
  #theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 18, color = "black", angle = 0))+ 
  theme(axis.text.x = element_text(size = 18, color = "black", hjust = 0.5, angle = 0))+ 
  theme(
    panel.grid.major = element_line(size = 0.5, linetype = "solid"))
print(boxplot4)
####################################################################
#18S boxplot of date "one"- August 2017
boxplot18_1 <- amp_boxplot(data= LS1_18_ampvis,
                           group_by = "Needle",
                           tax_empty = "remove",
                           tax_aggregate = "Phylum",
                           tax_show = 15
)+
  scale_color_manual(values=c("black","navajowhite4","darkorange","blue","grey"), name="Level of Impact")+  
  scale_y_continuous(limit=c(0,40), breaks=c(0,10,20,30,40))+
  #theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 18, color = "black", angle = 0))+ 
  theme(axis.text.x = element_text(size = 18, color = "black", hjust = 0.5, angle = 0))+ 
  theme(
    panel.grid.major = element_line(size = 0.5, linetype = "solid"))
print(boxplot_18_1)
####################################################################
#18S boxplot of date "four"- July 2018
boxplot18_4 <- amp_boxplot(data= LS4_18_ampvis,
                           group_by = "Needle",
                           tax_empty = "remove",
                           tax_aggregate = "Phylum",
                           tax_show = 5
)+
  scale_color_manual(values=c("black","darkgreen","navajowhite4","red4","blue"), name="Level of Impact")+  
  #scale_y_continuous(limit=c(0,7), breaks=c(0,1,2,3,4,5,6,7))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 18, color = "black", angle = 0))+ 
  theme(axis.text.x = element_text(size = 18, color = "black", hjust = 0.5, angle = 0))+ 
  theme(
    panel.grid.major = element_line(size = 0.5, linetype = "solid"))
print(boxplot_18_4)