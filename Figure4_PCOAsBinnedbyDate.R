#This code is written with the libraries and mapping file Compiled2018MappingFile.txt are already loaded from the Compiled2018_Prerocessing.Rmd and are recognized within the R environment. 
#____________________________ Ordination for 16S___________________________________________
#Object to ordinate August 2017 DNA samples
ordu_we_LS16_1 = ordinate(CB16LS_1, "PCoA", "unifrac", weighted = TRUE)
#Object to ordinate October 2017 DNA samples
ordu_we_LS16_2 = ordinate(CB16LS_2, "PCoA", "unifrac", weighted = TRUE)
#Object to ordinate May 2018 DNA samples
ordu_we_LS16_3 = ordinate(CB16LS_3, "PCoA", "unifrac", weighted = TRUE)
#Object to ordinate July 2018 DNA samples
ordu_we_LS16_4 = ordinate(CB16LS_4, "PCoA", "unifrac", weighted = TRUE)
#Object to ordinate October 2018 DNA samples
ordu_we_LS16_5 = ordinate(CB16LS_5, "PCoA", "unifrac", weighted = TRUE)
#Let's set a publishable theme:
theme_set(theme_bw())

#Let's create a PCOA plot for DNA samples binned into Date "one"
one=plot_ordination(CB16LS_1,  ordu_we_LS16_1, color = "Needle", shape = "Needle") + 
  geom_point(size = 10) + 
  ggtitle("08/02/2017 P=0.026*, R=0.379") +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#Let's assign colors to the samples. #009E73 and #D55EE00 are color-blind friendly colors for "red" and "green."
  scale_color_manual(values=c("grey","#009E73","navajowhite4","#D55E00"), name="Level of Impact") + 
#Let's also assign shapes to each sample.  
  scale_shape_manual(values=c(19,17,15,17))+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_1,  ordu_we_LS16_1, type="scree") 
#ADONIS test for ordination "one"
groupsig = transform_sample_counts(CB16LS_1  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_1)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis

#Let's create a PCOA plot for DNA samples binned into Date "two"
two=plot_ordination(CB16LS_2,  ordu_we_LS16_2, color = "Needle") + 
  geom_point(size = 4) + 
  theme(legend.position = "none")+
  ggtitle("10/14/2017") +
  scale_color_manual(values=c("grey","#009E73","navajowhite4","#D55E00"), name="Level of Impact") + 
#Let's also assign shapes to each sample.  
  scale_shape_manual(values=c(19,17,15,17))+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_2,  ordu_we_LS16_2, type="scree")
#ADONIS test for ordination "two"
groupsig = transform_sample_counts(CB16LS_2  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_2)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis

#Let's create a PCOA plot for DNA samples binned into Date "three"
three=plot_ordination(CB16LS_3,  ordu_we_LS16_3, color = "Needle") + 
  geom_point(size = 4) + 
  theme(legend.position = "none")+
  ggtitle("05/27/2018") +
  scale_color_manual(values=c("grey","#009E73","navajowhite4","#D55E00"), name="Level of Impact") + 
#Let's also assign shapes to each sample.  
  scale_shape_manual(values=c(19,17,15,17))+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_3,  ordu_we_LS16_3, type="scree")
#ADONIS test for ordination "three"
groupsig = transform_sample_counts(CB16LS_3  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_3)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis

#Let's create a PCOA plot for DNA samples binned into Date "four"
four=plot_ordination(CB16LS_4,  ordu_we_LS16_4, color = "Needle") + 
  geom_point(size = 10) + 
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("7/19/2018 P=0.004**, R=0.530") +
  scale_color_manual(values=c("grey","#009E73","navajowhite4","#D55E00","blue"), name="Level of Impact") + 
 #Let's also assign shapes to each sample.  
  scale_shape_manual(values=c(19,17,15,17,19))+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_4,  ordu_we_LS16_4, type="scree")
#ADONIS test for ordination "four"
groupsig = transform_sample_counts(CB16LS_4  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_4)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis

#Let's create a PCOA plot for DNA samples binned into Date "five"
five=plot_ordination(CB16LS_5,  ordu_we_LS16_5, color = "Needle") + 
  geom_point(size = 4) + 
  theme(legend.position = "none")+
  ggtitle("10/25/2018") +
  scale_color_manual(values=c("#D55E00","grey","#009E73","navajowhite4","blue"), name="Level of Impact") + 
 #Let's also assign shapes to each sample.  
  scale_shape_manual(values=c(19,17,15,17,19))+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_5,  ordu_we_LS16_5, type="scree")
#ADONIS test for ordination "five" 
groupsig = transform_sample_counts(CB16LS_5  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_5)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis

#____________________________Ordination for 18S___________________________________________
#Object to ordinate August 2017 DNA samples
ordu_we_LS18_1 = ordinate(CB18LS_1, "PCoA", "unifrac", weighted = TRUE)
#Object to ordinate October 2017 DNA samples                                   
ordu_we_LS18_2 = ordinate(CB18LS_2, "PCoA", "unifrac", weighted = TRUE)
#Object to ordinate May 2018 DNA samples                                   
ordu_we_LS18_3 = ordinate(CB18LS_3, "PCoA", "unifrac", weighted = TRUE)
#Object to ordinate July 2018 DNA samples                                   
ordu_we_LS18_4 = ordinate(CB18LS_4, "PCoA", "unifrac", weighted = TRUE)
#Object to ordinate October 2018 DNA samples                                   
ordu_we_LS18_5 = ordinate(CB18LS_5, "PCoA", "unifrac", weighted = TRUE)

#Let's create a PCOA plot for DNA samples binned into Date "one"
one=plot_ordination(CB18LS_1,  ordu_we_LS18_1, color = "TreeSpecies") + 
  geom_point(size = 10) + 
  ggtitle("08/02/2017 P=0.018*, R=0.382") +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values=c("grey","#009E73","navajowhite4","#D55E00"), name="Level of Impact") + 
#Let's also assign shapes to each sample.  
  scale_shape_manual(values=c(19,17,15,17))+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
#ADONIS test 
groupsig = transform_sample_counts(CB18LS_1  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB18LS_1)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

#Let's create a PCOA plot for DNA samples binned into Date "two"
two=plot_ordination(CB18LS_2,  ordu_we_LS18_2, color = "Needle") + 
  geom_point(size = 4) + 
  theme(legend.position = "none")+
  ggtitle("10/14/2017") +
  scale_color_manual(values=c("grey","#009E73","navajowhite4","#D55E00"), name="Level of Impact") + 
#Let's also assign shapes to each sample.  
  scale_shape_manual(values=c(19,17,15,17))+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
#ADONIS test 
groupsig = transform_sample_counts(CB18LS_2  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB18LS_2)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis

#Let's create a PCOA plot for DNA samples binned into Date "three"
three=plot_ordination(CB18LS_3,  ordu_we_LS18_3, color = "Needle") + 
  geom_point(size = 4) + 
  theme(legend.position = "none")+
  ggtitle("05/27/2018") +
  scale_color_manual(values=c("grey","#009E73","navajowhite4","#D55E00"), name="Level of Impact") + 
#Let's also assign shapes to each sample.  
  scale_shape_manual(values=c(19,17,15,17))+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
#ADONIS test 
groupsig = transform_sample_counts(CB18LS_3  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB18LS_3)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis

#Let's create a PCOA plot for DNA samples binned into Date "four"
four=plot_ordination(CB18LS_4,  ordu_we_LS18_4, color = "Needle") + 
  geom_point(size = 10) + 
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("7/19/2018 P=0.004**, R=0.530") +
  scale_color_manual(values=c("grey","#009E73","navajowhite4","#D55E00","blue"), name="Level of Impact") + 
#Let's also assign shapes to each sample.  
  scale_shape_manual(values=c(19,17,15,17,19))+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
#ADONIS test **Note that the q~object here is "TreeSpecies". This is due to n<3 for the spruce samples after rarifying. 
groupsig = transform_sample_counts(CB18LS_4  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB18LS_4)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis

#Let's create a PCOA plot for DNA samples binned into Date "five"
five=plot_ordination(CB18LS_5,  ordu_we_LS18_5, color = "Needle") + 
  geom_point(size = 4) + 
  theme(legend.position = "none")+
  ggtitle("10/25/2018") +
  scale_color_manual(values=c("#D55E00","grey","#009E73","navajowhite4","blue"), name="Level of Impact") + 
#Let's also assign shapes to each sample.  
  scale_shape_manual(values=c(19,17,15,17,19))+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
#ADONIS test 
groupsig = transform_sample_counts(CB18LS_5  , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB18LS_5)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis
