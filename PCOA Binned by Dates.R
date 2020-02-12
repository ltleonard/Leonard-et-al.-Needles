ordu_we_LS16_1 = ordinate(CB16LS_1, "PCoA", "unifrac", weighted = TRUE)
ordu_we_LS16_2 = ordinate(CB16LS_2, "PCoA", "unifrac", weighted = TRUE)
ordu_we_LS16_3 = ordinate(CB16LS_3, "PCoA", "unifrac", weighted = TRUE)
ordu_we_LS16_4 = ordinate(CB16LS_4, "PCoA", "unifrac", weighted = TRUE)
ordu_we_LS16_5 = ordinate(CB16LS_5, "PCoA", "unifrac", weighted = TRUE)

sample_sums(CB18LS_4)

#Remember you didn't rename the ordu_we's so you need to rerun everything each time
one=plot_ordination(CB16LS_1,  ordu_we_LS16_1, color = "TreeSpecies") + 
  geom_point(size = 10) + 
  #facet_grid(~Snowmelt) +
  ggtitle("08/02/2017 P=0.026*, R=0.379") +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4"), name="Level of Impact") + 
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_1,  ordu_we_LS16_1, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_1  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_1)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

two=plot_ordination(CB16LS_2,  ordu_we_LS16_2, color = "TreeSpecies") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  theme(legend.position = "none")+
  ggtitle("10/14/2017") +
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4"), name="Level of Impact") + 
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_2,  ordu_we_LS16_2, type="scree")

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_2  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_2)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

three=plot_ordination(CB16LS_3,  ordu_we_LS16_3, color = "TreeSpecies") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  theme(legend.position = "none")+
  ggtitle("05/27/2018") +
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4"), name="Level of Impact") + 
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_3,  ordu_we_LS16_3, type="scree")

#ADONIS test
groupsig = transform_sample_counts(CB16LS_3  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_3)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

four=plot_ordination(CB16LS_4,  ordu_we_LS16_4, color = "TreeSpecies") + 
  geom_point(size = 10) + 
  #facet_grid(~Snowmelt) +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("7/19/2018 P=0.004**, R=0.530") +
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") + 
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_4,  ordu_we_LS16_4, type="scree")

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_4  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_4)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

five=plot_ordination(CB16LS_5,  ordu_we_LS16_5, color = "TreeSpecies") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  theme(legend.position = "none")+
  ggtitle("10/25/2018") +
  scale_color_manual(values=c("red4","grey","darkgreen","navajowhite4","blue"), name="Level of Impact") + 
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_5,  ordu_we_LS16_5, type="scree")

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_5  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_5)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

plot_grid(one, ncol=1)

#____________________________18S___________________________________________
ordu_we_LS18_1 = ordinate(CB18LS_1, "PCoA", "unifrac", weighted = TRUE)
ordu_we_LS18_2 = ordinate(CB18LS_2, "PCoA", "unifrac", weighted = TRUE)
ordu_we_LS18_3 = ordinate(CB18LS_3, "PCoA", "unifrac", weighted = TRUE)
ordu_we_LS18_4 = ordinate(CB18LS_4, "PCoA", "unifrac", weighted = TRUE)
ordu_we_LS18_5 = ordinate(CB18LS_5, "PCoA", "unifrac", weighted = TRUE)

#Remember you didn't rename the ordu_we's so you need to rerun everything each time
one=plot_ordination(CB18LS_1,  ordu_we_LS18_1, color = "TreeSpecies") + 
  geom_point(size = 10) + 
  #facet_grid(~Snowmelt) +
  ggtitle("08/02/2017 P=0.018*, R=0.382") +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4"), name="Level of Impact") + 
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))

#ADONIS test 
groupsig = transform_sample_counts(CB18LS_1  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB18LS_1)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

two=plot_ordination(CB18LS_2,  ordu_we_LS18_2, color = "TreeSpecies") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  theme(legend.position = "none")+
  ggtitle("10/14/2017") +
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4"), name="Level of Impact") + 
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))

#ADONIS test 
groupsig = transform_sample_counts(CB18LS_2  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB18LS_2)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

three=plot_ordination(CB18LS_3,  ordu_we_LS18_3, color = "TreeSpecies") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  theme(legend.position = "none")+
  ggtitle("05/27/2018") +
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4"), name="Level of Impact") + 
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))

#ADONIS test 
groupsig = transform_sample_counts(CB18LS_3  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB18LS_3)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

four=plot_ordination(CB18LS_4,  ordu_we_LS18_4, color = "TreeSpecies") + 
  geom_point(size = 10) + 
  #facet_grid(~Snowmelt) +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("7/19/2018") +
  scale_color_manual(values=c("grey","navajowhite4","darkgreen","red4","blue"), name="Level of Impact") + 
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
##########Look at q~!!!! I changed it to species
#ADONIS test 
groupsig = transform_sample_counts(CB18LS_4  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB18LS_4)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis


five=plot_ordination(CB18LS_5,  ordu_we_LS18_5, color = "TreeSpecies") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  theme(legend.position = "none")+
  ggtitle("10/25/2018") +
  scale_color_manual(values=c("red4","grey","darkgreen","navajowhite4","blue"), name="Level of Impact") + 
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))

#ADONIS test 
groupsig = transform_sample_counts(CB18LS_5  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB18LS_5)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

plot_grid(one, two, three, four, five, ncol=3)

#__________________________________________________________________________
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(CB18LS_5, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Needle", title="Bray NMDS")

ps.Ord <- ordinate(ps.prop, method = "PCoA", distance = "unifrac")
ord=plot_ordination(ps.prop, ps.Ord, color="Needle")
ord+geom_point(size=3.5)