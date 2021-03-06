library(PMCMR)
#Let's plot the observed alpha diversity based on the phyloseq objects loaded into the environment from the preprocessing steps. 

#16S alpha diversity for 16S date "one"
LS16_1=plot_richness(CB16LS_1, x = "Needle", measures = "Observed") +
  geom_point(colour = "white") + 
  #Add tabs to the whiskers on the plot
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(fill=factor(Needle)))  + 
  facet_grid(Location~.) +
  theme(legend.position=c("none"),plot.title=element_text(size=16, hjust=0.5),
        strip.text.y = element_text(size=16), axis.text.x = element_text(size=16, angle = 0, vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), 
        axis.title.y = element_text(size=16), strip.text.x = element_text(size=16)) + 
  #Change the fill colors of the boxplots
  scale_fill_manual(values = c("grey","darkgreen","navajowhite4","red4")) +
  scale_y_continuous(limit=c(0,650), breaks=c(0,100,200,300,400,500,600))+
  ylab("Observed OTUs")
#Change the Labels of the x-axis ticks
#scale_x_discrete(breaks=c('Control','Green', 'Lodge','Red'),
# labels=c("Control","Lodgepole","Green Spruce","Red Spruce"))

#Alpha diversity for 16S date "four"
LS16_4=plot_richness(CB16LS_4, x = "Needle", measures = "Observed") +
  geom_point(colour = "white") + 
  #Add tabs to the whiskers on the plot
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(fill=factor(Needle)))  + 
  facet_grid(Location~.) +
  theme(legend.position=c("none"),plot.title=element_text(size=16, hjust=0.5),
        strip.text.y = element_text(size=16), axis.text.x = element_text(size=16, angle = 0, vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), 
        axis.title.y = element_text(size=16), strip.text.x = element_text(size=16)) + 
  #Change the fill colors of the boxplots
  scale_fill_manual(values = c("grey","darkgreen","navajowhite4","red4","blue")) +
  scale_y_continuous(limit=c(0,650), breaks=c(0,100,200,300,400,500,600))+
  ylab("Observed OTUs")
#Change the Labels of the x-axis ticks
#scale_x_discrete(breaks=c('Control','Green', 'Lodge','Red'),
# labels=c("Control","Lodgepole","Green Spruce","Red Spruce"))

#______________________________________________________________________________________________
#Alpha diversity for 18S date "one"
LS18_1=plot_richness(CB18LS_1, x = "Needle", measures = "Observed") +
  geom_point(colour = "white") + 
  #Add tabs to the whiskers on the plot
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(fill=factor(Needle)))  + 
  facet_grid(Location~.) +
  theme(legend.position=c("none"),plot.title=element_text(size=16, hjust=0.5),
        strip.text.y = element_text(size=16), axis.text.x = element_text(size=16, angle = 0, vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), 
        axis.title.y = element_text(size=16), strip.text.x = element_text(size=16)) + 
  #Change the fill colors of the boxplots
  scale_fill_manual(values = c("grey","darkgreen","navajowhite4","red4")) +
  scale_y_continuous(limit=c(0,75), breaks=c(0,25,50,75))+
  ylab("Observed OTUs")
#Change the Labels of the x-axis ticks
#scale_x_discrete(breaks=c('Control','Green', 'Lodge','Red'),
# labels=c("Control","Lodgepole","Green Spruce","Red Spruce"))

#Alpha diversity for 18S date "four"
LS18_4=plot_richness(CB18LS_4, x = "Needle", measures = "Observed") +
  geom_point(colour = "white") + 
  #Add tabs to the whiskers on the plot
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(aes(fill=factor(Needle)))  + 
  facet_grid(Location~.) +
  theme(legend.position=c("none"),plot.title=element_text(size=16, hjust=0.5),
        strip.text.y = element_text(size=16), axis.text.x = element_text(size=16, angle = 0, vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), 
        axis.title.y = element_text(size=16), strip.text.x = element_text(size=16)) + 
  #Change the fill colors of the boxplots
  scale_y_continuous(limit=c(0,75), breaks=c(0,25,50,75))+
  scale_fill_manual(values = c("grey","darkgreen","navajowhite4","red4","blue")) +
  ylab("Observed OTUs")
#Change the Labels of the x-axis ticks
#scale_x_discrete(breaks=c('Control','Green', 'Lodge','Red'),
# labels=c("Control","Lodgepole","Green Spruce","Red Spruce"))