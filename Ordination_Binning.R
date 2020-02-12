
CB16LS_Red= subset_samples(CB16_LS_no2red_bot_dubs, Needle=="Red")
CB16LS_Green= subset_samples(CB16_LS_no2red_bot_dubs, Needle=="Green")
CB16LS_Lodge= subset_samples(CB16_LS_no2red_bot_dubs, Needle=="Lodge")
CB16LS_Control= subset_samples(CB16_LS_no2red_bot_dubs, Needle=="Control")

CB16LS_2018= subset_samples(CB16_LS_no2red_bot_dubs, Year=="b_2018")
CB16LS_2017= subset_samples(CB16_LS_no2red_bot_dubs, Year=="a_2017")

CB16LS_Late= subset_samples(CB16_LS_no2red_bot_dubs, Season=="Late")
CB16LS_Mid= subset_samples(CB16_LS_no2red_bot_dubs, Season=="Middle")
CB16LS_Early= subset_samples(CB16_LS_no2red_bot_dubs, Season=="Early")
##Only Lodgepole V Spruce
CB16LS_Needles= subset_samples(CB16_LS_no2red_bot_dubs)
Needles_2018= subset_samples(CB16LS_Needles, Year=="b_2018")
Needles_2017= subset_samples(CB16LS_Needles, Year=="a_2017")
Needles_Late= subset_samples(CB16LS_Needles, Season=="Late")
Needles_Mid= subset_samples(CB16LS_Needles, Season=="Middle")
Needles_Early= subset_samples(CB16LS_Needles, Season=="Early")


ordu_we_LS16_regular = ordinate(CB16_LS_no2red_bot_dubs, "PCoA", "unifrac", weighted = TRUE)

#Run with all 5 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16_LS_no2red_bot_dubs,  ordu_we_LS16_regular, color = "Date", shape="Needle") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("2017 & 2018, P=0.001**, R=0.300") +
  #theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16_LS_no2red_bot_dubs,  ordu_we_LS16_regular, type="scree") 

#ADONIS test 
##########LOOK AT THE q~ VALUE!!!!!!!!
groupsig = transform_sample_counts(CB16_LS_no2red_bot_dubs , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16_LS_no2red_bot_dubs)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis
##########################################################################
ordu_we_LS16_2018 = ordinate(CB16LS_2018, "PCoA", "unifrac", weighted = TRUE)

#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_2018,  ordu_we_LS16_2018, color = "Needle") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("2018 P=0.001**, R=0.319") +
  #theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_2018,  ordu_we_LS16_2018, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_2018 , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_2018)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis
##########################################################################
ordu_we_LS16_2017 = ordinate(CB16LS_2017, "PCoA", "unifrac", weighted = TRUE)

#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_2017,  ordu_we_LS16_2017, color = "Needle") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("2017") +
  #theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_2017,  ordu_we_LS16_2017, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_2017 , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_2017)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis
##########################################################################
ordu_we_LS16_Late = ordinate(CB16LS_Late, "PCoA", "unifrac", weighted = TRUE)

#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_Late,  ordu_we_LS16_Late, color = "Needle") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("October, P=0.001***, R=0.409") +
  #theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_Late,  ordu_we_LS16_Late, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_Late , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_Late)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis
##########################################################################
ordu_we_LS16_Mid = ordinate(CB16LS_Mid, "PCoA", "unifrac", weighted = TRUE)

#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_Mid,  ordu_we_LS16_Mid, color = "Needle", shape="Date") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("October, P=0.005**, R=0.124") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_Mid,  ordu_we_LS16_Mid, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_Mid , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_Mid)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis
##########################################################################
ordu_we_LS16_Early = ordinate(CB16LS_Early, "PCoA", "unifrac", weighted = TRUE)

#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_Early,  ordu_we_LS16_Early, color = "Needle") + 
geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("October, P=0.005**, R=0.124") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_Early,  ordu_we_LS16_Early, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_Early , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_Early)) 
adonis<-adonis(q~Needle, data=d,permutations=999)
adonis
##########################################################################
##########################################################################
##########################################################################
#Do it all as just function of type now. 
#Run with all 5 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16_LS_no2red_bot_dubs,  ordu_we_LS16_regular, color = "Type") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("2017 & 2018, P=0.001***, R=0.300") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16_LS_no2red_bot_dubs,  ordu_we_LS16_regular, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16_LS_no2red_bot_dubs , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16_LS_no2red_bot_dubs)) 
adonis<-adonis(q~Type, data=d,permutations=999)
adonis
##########################################################################
#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_2018,  ordu_we_LS16_2018, color = "Type") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("2018 P=0.001**, R=0.319") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_2018,  ordu_we_LS16_2018, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_2018 , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_2018)) 
adonis<-adonis(q~Type, data=d,permutations=999)
adonis
##########################################################################
#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_2017,  ordu_we_LS16_2017, color = "Type") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("2017") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_2017,  ordu_we_LS16_2017, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_2017 , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_2017)) 
adonis<-adonis(q~Type, data=d,permutations=999)
adonis
##########################################################################
#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_Late,  ordu_we_LS16_Late, color = "Type", shape="Species") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("October, P=0.001***, R=0.409") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_Late,  ordu_we_LS16_Late, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_Late , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_Late)) 
adonis<-adonis(q~Type, data=d,permutations=999)
adonis
##########################################################################
#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_Mid,  ordu_we_LS16_Mid, color = "Type", shape="Species") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("October, P=0.005**, R=0.124") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_Mid,  ordu_we_LS16_Mid, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_Mid , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_Mid)) 
adonis<-adonis(q~Type, data=d,permutations=999)
adonis
##########################################################################
#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_Early,  ordu_we_LS16_Early, color = "Type", shape="Date") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) + 
  ggtitle("October, P=0.005**, R=0.124") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_Early,  ordu_we_LS16_Early, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_Early , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_Early)) 
adonis<-adonis(q~Type, data=d,permutations=999)
adonis
##########################################################################
##########################################################################
##########################################################################
ordu_we_LS16_Needles = ordinate(CB16LS_Needles, "PCoA", "unifrac", weighted = TRUE)

#Do it all as just function of type now. 
#Run with all 5 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(CB16LS_Needles,  ordu_we_LS16_Needles, color = "TreeSpecies", shape="Needle") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("2017 & 2018, P=0.001***, R=0.325") +
  #theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(CB16LS_Needles,  ordu_we_LS16_Needles, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(CB16LS_Needles , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(CB16LS_Needles)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis
##########################################################################
ordu_we_Needles_2018 = ordinate(Needles_2018, "PCoA", "unifrac", weighted = TRUE)

#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(Needles_2018,  ordu_we_Needles_2018, color = "TreeSpecies", shape="Needle") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("2018 P=0.001***, R=0.318") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(Needles_2018,  ordu_we_Needles_2018, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(Needles_2018 , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(Needles_2018)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis
##########################################################################
ordu_we_Needles_2017 = ordinate(Needles_2017, "PCoA", "unifrac", weighted = TRUE)

#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(Needles_2017,  ordu_we_Needles_2017, color = "TreeSpecies", shape="Needle") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("2017, P=0.034*, R=0.115") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(Needles_2017,  ordu_we_Needles_2017, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(Needles_2017 , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(Needles_2017)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis
##########################################################################
ordu_we_Needles_Late = ordinate(Needles_Late, "PCoA", "unifrac", weighted = TRUE)

#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(Needles_Late,  ordu_we_Needles_Late, color = "TreeSpecies", shape="Needle") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("October, P=0.001***, R=0.418") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(Needles_Late,  ordu_we_Needles_Late, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(Needles_Late , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(Needles_Late)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis
##########################################################################
ordu_we_Needles_Mid = ordinate(Needles_Mid, "PCoA", "unifrac", weighted = TRUE)

#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(Needles_Mid,  ordu_we_Needles_Mid, color = "TreeSpecies", shape="Needle") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) +
  ggtitle("October, P=0.001***, R=0.222") +
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(Needles_Mid,  ordu_we_Needles_Mid, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(Needles_Mid , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(Needles_Mid)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis
##########################################################################
ordu_we_Needles_Early = ordinate(Needles_Early, "PCoA", "unifrac", weighted = TRUE)

#Run with all 3 2018 dates together for the Lower Subalpine Regular as a function of needle type
Regular=plot_ordination(Needles_Early,  ordu_we_Needles_Early, color = "TreeSpecies", shape="Needle") + 
  geom_point(size = 4) + 
  #facet_grid(~Snowmelt) + 
  ggtitle("October, P=0.001***, R=0.222") +
  #theme(legend.position = "none")+
  scale_color_manual(values=c("grey","darkgreen","navajowhite4","red4","blue"), name="Level of Impact") +
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), strip.text.x = element_text(size=16))
plot_ordination(Needles_Early,  ordu_we_Needles_Early, type="scree") 

#ADONIS test 
groupsig = transform_sample_counts(Needles_Early , function(x) 100* x/ sum(x))
ufrac<-UniFrac(groupsig, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)
d<-data.frame(sample_data(Needles_Early)) 
adonis<-adonis(q~TreeSpecies, data=d,permutations=999)
adonis