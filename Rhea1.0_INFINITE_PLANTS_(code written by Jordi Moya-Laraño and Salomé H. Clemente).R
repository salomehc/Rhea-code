#!/usr/bin/Rscript

setwd("")
getwd()

##setwd("")

version="NULL" ###"SCATTERED" OR "EMPIRICAL"  
###"NULL" means a null model without male damage from interference
###"SCATTERED" means a version in which the damage is scattered among all possible matings (order nor timing is not relevant)
###"EMPIRICAL" based on published paper

scenario="FORCED" ###"EXCLUSION" no culling, "FORCED" forced coexistence from species-bias in culling

n_tevansi=150   
n_turticae=150
phi_m_dmg=0.01   ##plug 0.999 for simulations without evolution of RI
steps<-300       ##days of simulation
##tuner_for_dispersal=10

source("all_functions_INFINITE_PLANTS_24_02_2016bis.r")     ##

sim_num=2
folder=paste("",sim_num,sep="")
dir.create(folder,recursive=T)
setwd(folder)
getwd()

A

   #############################################
   #### GENETIC FUNCTIONS ######################
   #############################################
   n_loci=20 ### number of loci per trait - pseudovalues for interpolation will need
   ##to be changed if a different number of loci is used
   n_traits=5 ##includes the pseudo-trait drift
   n_modules=2 ###the module  n_traits%%n_moudles corresponds to the drift trait
                ###the round(n_traits/n_modules) no the number of traits per module

  ###for recombination build realistic chromosome in which all the loci from all traits
  ###are randomly scattered around a single matrix
  ###instead of renaming all the allel IDs and things, we just assing a position to each
  ###locus for each trait and then randomly reassing the order in a 20x5 matrix
library(reshape)

##source("chrom.r")

tu_chromosome<-chrom()  ###assings random positions to the alleles as in true quantitive genetics
te_chromosome<-chrom()



  library(psych) ##to test if correlations are right

      n_alleles=10 ####alleles available in the population
      rho1=0.0  ###corr dispersal vs assim efficiency - mod1 - play with rho=0 and rho=-0.99
      rho2=0.0  ###corr male damage vs sex ratio - mod2
      mod_type_tu="neg"##whether there are two negative correlations among traits or are all positive "pos"
      mod_type_te="neg"

      ini_sex_ratio=0.75##0.75

    ###FOR 3-TRAIT MODULES
    #### this creates a chromosome type from which each individual samples its genetics for each of its 2 chromosomes


    tevansi<-as.data.frame(matrix(nrow=n_tevansi),ncol=1)
    turticae<-as.data.frame(matrix(nrow=n_turticae),ncol=1)

    names(tevansi)<-"sex"
    names(turticae)<-"sex"

    tevansi<-sex_det(tevansi,ini_sex_ratio)
    turticae<-sex_det(turticae,ini_sex_ratio)

    table(tevansi$sex)
    table(turticae$sex)
  ################################

    #######################################################
    #### tevansi and turticae TRAITS ######################
    #######################################################

    ##for simplification the loci values are the same for all traits
    loci<-crea_loci(n_loci,n_alleles)



    ###PREY MOD1 - calls functions to originate the phenotypes for module 1
    mod_tu1<-crea_module(turticae,mod_type=mod_type_tu,loci,n_loci,n_alleles,rho1)
    mod_te1<-crea_module(tevansi,mod_type=mod_type_te,loci,n_loci,n_alleles,rho1)

    ###PREY MOD2 - calls functions to originate the phenotypes for module 2
    mod_tu2<-crea_module(turticae,mod_type=mod_type_tu,loci,n_loci,n_alleles,rho2)
    mod_te2<-crea_module(tevansi,mod_type=mod_type_te,loci,n_loci,n_alleles,rho2)


   # names(positions)

    drift_te<-crea_trait(tevansi,loci,n_loci,n_alleles)
    drift_tu<-crea_trait(turticae,loci,n_loci,n_alleles)


  #################################################################################################################################################
  #################################################################################################################################################
  #################################################################################################################################################
  #################################################################################################################################################
  #################################################################################################################################################
  #################################################################################################################################################
  #################################################################################################################################################
  ###PHENOTYPES AND MAIN LOOP START HERE


  ##################
        library(reshape)
        
        #ini_sex_ratio=2
        tevansi<-as.data.frame(matrix(nrow=n_tevansi),ncol=1)
        turticae<-as.data.frame(matrix(nrow=n_turticae),ncol=1)
        names(tevansi)<-"sex"
        names(turticae)<-"sex"
        density="uncrowded"

        phi_ass=1     ##the phi value in equations A3 and A4 in Moya-Laraño et al. 2012
        phi_R_P=1     ##higher values mean high additive phenotypic variation due to QG
        phi_sexr=1

        ###assign same sex as in the genetics
        tevansi$sex<-mod_te1$sex
        turticae$sex<-mod_tu1$sex




  mod_te1$sp=rep("te",n_tevansi)
  mod_tu1$sp=rep("tu",n_turticae)

  mod_te2$sp=rep("te",n_tevansi)
  mod_tu2$sp=rep("tu",n_turticae)

  mod_te1$ID=paste("te",mod_te1$ID,sep="_")
  mod_tu1$ID=paste("tu",mod_tu1$ID,sep="_")

  mod_te2$ID=paste("te",mod_te2$ID,sep="_")
  mod_tu2$ID=paste("tu",mod_tu2$ID,sep="_")


  #both1 creates a list with all individuals, of both species, and with the genetics of each one,
  #it can be thus used for trait values( and for reproduction further on )
  ##for module 1
  both1 <- list(sp=c(mod_te1$sp,mod_tu1$sp) ,ID=c(mod_te1$ID,mod_tu1$ID),sex=c(mod_te1$sex,mod_tu1$sex),
  trt1_cr1_ID=rbind(mod_te1$trt1_cr1_ID,mod_tu1$trt1_cr1_ID), trt1_cr1_val=rbind(mod_te1$trt1_cr1_val,
  mod_tu1$trt1_cr1_val), trt1_cr2_ID=rbind(mod_te1$trt1_cr2_ID,mod_tu1$trt1_cr2_ID),
  trt1_cr2_val=rbind(mod_te1$trt1_cr2_val,mod_tu1$trt1_cr2_val),trt2_cr1_ID=rbind(mod_te1$trt2_cr1_ID,
  mod_tu1$trt2_cr1_ID), trt2_cr1_val=rbind(mod_te1$trt2_cr1_val,mod_tu1$trt2_cr1_val),
  trt2_cr2_ID= rbind(mod_te1$trt2_cr2_ID,mod_tu1$trt2_cr2_ID),trt2_cr2_val= rbind(mod_te1$trt2_cr2_val,
  mod_tu1$trt2_cr2_val),phenotypes=rbind(mod_te1$phenotypes,mod_tu1$phenotypes),
  expressing1=rbind(mod_te1$expressing1,mod_tu1$expressing1),
  expressing2=rbind(mod_te1$expressing2,mod_tu1$expressing2))

  ##for module 2
  both2 <- list(sp=c(mod_te2$sp,mod_tu2$sp) ,ID=c(mod_te2$ID,mod_tu2$ID),sex=c(mod_te2$sex,mod_tu2$sex),
  trt1_cr1_ID=rbind(mod_te2$trt1_cr1_ID,mod_tu2$trt1_cr1_ID), trt1_cr1_val=rbind(mod_te2$trt1_cr1_val,
  mod_tu2$trt1_cr1_val), trt1_cr2_ID=rbind(mod_te2$trt1_cr2_ID,mod_tu2$trt1_cr2_ID),
  trt1_cr2_val=rbind(mod_te2$trt1_cr2_val,mod_tu2$trt1_cr2_val),trt2_cr1_ID=rbind(mod_te2$trt2_cr1_ID,
  mod_tu2$trt2_cr1_ID), trt2_cr1_val=rbind(mod_te2$trt2_cr1_val,mod_tu2$trt2_cr1_val),
  trt2_cr2_ID= rbind(mod_te2$trt2_cr2_ID,mod_tu2$trt2_cr2_ID),trt2_cr2_val= rbind(mod_te2$trt2_cr2_val,
  mod_tu2$trt2_cr2_val),phenotypes=rbind(mod_te2$phenotypes,mod_tu2$phenotypes),
  expressing1=rbind(mod_te2$expressing1,mod_tu2$expressing1),
  expressing2=rbind(mod_te2$expressing2,mod_tu2$expressing2))

  drift_te$sp=rep("te",n_tevansi)
  drift_tu$sp=rep("tu",n_turticae)
  drift_te$ID<-paste("te",drift_te$ID,sep="_")
  drift_tu$ID<-paste("tu",drift_tu$ID,sep="_")


  ## for the single trait-drift thing
  both3 <- list(sp=c(drift_te$sp,drift_tu$sp) ,ID=c(drift_te$ID,drift_tu$ID),
  sex=c(drift_te$sex,drift_tu$sex),trt1_cr1_ID=rbind(drift_te$trt1_cr1_ID,
  drift_tu$trt1_cr1_ID),
  trt1_cr1_val=rbind(drift_te$trt1_cr1_val,drift_tu$trt1_cr1_val),
  trt1_cr2_ID=rbind(drift_te$trt1_cr2_ID,drift_tu$trt1_cr2_ID),
  trt1_cr2_val=rbind(drift_te$trt1_cr2_val,drift_tu$trt1_cr2_val),
  phenotypes=rbind(drift_te$phenotypes,drift_tu$phenotypes),
  expressing1=rbind(drift_te$expressing1,drift_tu$expressing1))

  genotypes<-list(both1=both1,both2=both2,both3=both3)

######################################################
############PHENOTYPES START HERE ####################
######################################################

    bichinhos<-  as.data.frame(matrix(nrow=(n_tevansi+n_turticae),ncol=11,))
    names(bichinhos)[c(1:11)]<-c("instar","sex","sp","In_ID","order","eggs","drift","ass","R_P","m_damage","g_sex_ratio")


    bichinhos$sex[1:n_tevansi]<- (tevansi$sex)                                     #puts sex
    bichinhos$sex[(n_tevansi+1):(n_turticae+n_tevansi)]<- (turticae$sex)

    bichinhos$sp[1:n_tevansi]<- "te"                                               #puts species
    bichinhos$sp[(n_tevansi+1):(n_turticae+n_tevansi)]<- "tu"

    bichinhos$In_ID[1:n_tevansi]<-1:n_tevansi                                      #puts individual identification
    bichinhos$In_ID[(n_tevansi+1):(n_turticae+n_tevansi)]<-1:n_turticae ##(n_tevansi+1):(n_turticae+n_tevansi)

    bichinhos[,1]<- as.integer(round(runif((nrow(bichinhos)), min=1,max=5))) #attributes instar(between 1 and 5) to each individual)

    bichinhos$In_ID<-paste(bichinhos$sp,bichinhos$In_ID,sep="_")

    # we can do eventually an initialization with only adults

   ## corresponds to the number to the name of the instar, is this really necessary? - is better so we do not have to remember
    bichinhos[bichinhos$instar==1,1]<-"egg"
    bichinhos[bichinhos$instar==2,1]<-"larva"
    bichinhos[bichinhos$instar==3,1]<-"proto"
    bichinhos[bichinhos$instar==4,1]<-"deuto"
    bichinhos[bichinhos$instar==5,1]<-"adult"


    ## egg counter- each female is assigned 200 eggs in the beggining but males have 0
    ## but fecundity is ruled by energy, putting 200 allows for evolution of fecundity
    bichinhos$eggs<-ifelse(bichinhos$sex=="female",200,0)

  ###and Salomé's signature from when she loved science remains!
  ###adding up the numbers for each instar and species at inizialization
  ###this will have to be done during the simulation also - convert to generalized function
  bichinhos_all<-  as.data.frame(matrix(nrow=4,ncol=5,))
  names(bichinhos_all)[c(1:5)]<-c("egg","larva","proto","deuto","adult")
  row.names(bichinhos_all)<-c("male_tu", "female_tu","male_te", "female_te")
  bichinhos_all[1,1]<-sum( bichinhos$instar=="egg" & bichinhos$sex=="male" & bichinhos$sp=="tu" )
  bichinhos_all[2,1]<-sum( bichinhos$instar=="egg" & bichinhos$sex=="female" & bichinhos$sp=="tu" )
  bichinhos_all[3,1]<-sum( bichinhos$instar=="egg" & bichinhos$sex=="male" & bichinhos$sp=="te" )
  bichinhos_all[4,1]<-sum( bichinhos$instar=="egg" & bichinhos$sex=="female" & bichinhos$sp=="te" )

  bichinhos_all[1,2]<-sum( bichinhos$instar=="larva" & bichinhos$sex=="male" & bichinhos$sp=="tu" )
  bichinhos_all[2,2]<-sum( bichinhos$instar=="larva" & bichinhos$sex=="female" & bichinhos$sp=="tu" )
  bichinhos_all[3,2]<-sum( bichinhos$instar=="larva" & bichinhos$sex=="male" & bichinhos$sp=="te" )
  bichinhos_all[4,2]<-sum( bichinhos$instar=="larva" & bichinhos$sex=="female" & bichinhos$sp=="te" )

  bichinhos_all[1,3]<-sum( bichinhos$instar=="proto" & bichinhos$sex=="male" & bichinhos$sp=="tu" )
  bichinhos_all[2,3]<-sum( bichinhos$instar=="proto" & bichinhos$sex=="female" & bichinhos$sp=="tu" )
  bichinhos_all[3,3]<-sum( bichinhos$instar=="proto" & bichinhos$sex=="male" & bichinhos$sp=="te" )
  bichinhos_all[4,3]<-sum( bichinhos$instar=="proto" & bichinhos$sex=="female"& bichinhos$sp=="te" )

  bichinhos_all[1,4]<-sum( bichinhos$instar=="deuto" & bichinhos$sex=="male" & bichinhos$sp=="tu"  )
  bichinhos_all[2,4]<-sum( bichinhos$instar=="deuto" & bichinhos$sex=="female" & bichinhos$sp=="tu" )
  bichinhos_all[3,4]<-sum( bichinhos$instar=="deuto" & bichinhos$sex=="male"& bichinhos$sp=="te"   )
  bichinhos_all[4,4]<-sum( bichinhos$instar=="deuto" & bichinhos$sex=="female" & bichinhos$sp=="te" )

  bichinhos_all[1,5]<-sum( bichinhos$instar=="adult" & bichinhos$sex=="male" & bichinhos$sp=="tu"  )
  bichinhos_all[2,5]<-sum( bichinhos$instar=="adult" & bichinhos$sex=="female" & bichinhos$sp=="tu"   )
  bichinhos_all[3,5]<-sum( bichinhos$instar=="adult" & bichinhos$sex=="male" & bichinhos$sp=="te"  )
  bichinhos_all[4,5]<-sum( bichinhos$instar=="adult" & bichinhos$sex=="female"& bichinhos$sp=="te"  )

  ###from Mitchell 1973
  #we are assuming diff in deutonynf females (between crowded and uncrowded) are already related with energy acumulated for the eggs #so we are using only
  growth_mass<-as.data.frame(matrix(nrow=4,ncol=5,))
  names(growth_mass)[c(1:5)]<-c("egg","larva","proto","deuto","adult")
  row.names(growth_mass)<-c("male_u", "female_u","male_c", "female_c")
  growth_mass[1,] <- c( 0.49,  2.11, -0.65, 0.1039175,  0.5260825)
  growth_mass[2,] <- c(0.49,  2.11,  7.12, 2.4636364, 11.0863636)  ####adult column refers to gain during 1st 6 days as adults(after emergence)
  growth_mass[3,] <- c(0.49,  2.11,  0.22, 0.1105155,  0.5594845)
  growth_mass[4,] <- c(0.49,  2.11,  1.49, 2.8418182, 12.7881818)

  active_days<-as.data.frame(matrix(nrow=4,ncol=5,))
  names(active_days)[c(1:5)]<-c( "egg","larva","proto","deuto","adult")
  row.names(active_days)<-c("male_u", "female_u","male_c", "female_c")
  active_days[1,] <- c(1, 2.1, 1.5, 1.6, 6)
  active_days[2,] <- c(1, 2.1, 1.5, 1.8, 6)
  active_days[3,] <- c(1, 2.1, 1.5, 1.6, 6)
  active_days[4,] <- c(1, 2.1, 1.5, 1.8, 6)

  quiescent_days<-as.data.frame(matrix(nrow=4,ncol=5,))
  names(quiescent_days)[c(1:5)]<-c("egg","larva","proto","deuto","adult")
  row.names(quiescent_days)<-c("male_u", "female_u","male_c", "female_c")
  quiescent_days[1,] <- c(7, 1.6, 1.5, 1.5,0)
  quiescent_days[2,] <- c(7, 1.6, 1.5, 1.7,0)
  quiescent_days[3,] <- c(7, 1.6, 1.5, 1.5,0)
  quiescent_days[4,] <- c(7, 1.6, 1.5, 1.7,0)

  mass_day<-as.data.frame(matrix(nrow=4,ncol=7,))
  names(mass_day)[c(1:7)]<-c("sex", "density","egg","larva","proto","deuto","adult")
  #row.names(mass_day)<-c("male_u", "female_u","male_c", "female_c")
  mass_day[1,]<- c("male", "uncrowded",growth_mass[1,(1:5)]/active_days[1,(1:5)])
  mass_day[2,]<- c("female", "uncrowded",growth_mass[2,(1:5)]/active_days[2,(1:5)])
  mass_day[3,]<- c("male", "crowded",growth_mass[3,(1:5)]/active_days[3,(1:5)])
  mass_day[4,]<- c("female", "crowded",growth_mass[4,(1:5)]/active_days[4,(1:5)])
  mass_day$sex<-as.character(mass_day$sex)
  mass_day$density<-as.character(mass_day$density)
  mass_day$egg<-as.numeric(mass_day$egg)
  mass_day$larva<-as.numeric(mass_day$larva)
  mass_day$proto<-as.numeric(mass_day$proto)
  mass_day$deuto<-as.numeric(mass_day$deuto)
  mass_day$adult<-as.numeric(mass_day$adult)

  ### Oliveira et al. Oecologia
  ass_eff<- as.data.frame(matrix(nrow=2,ncol=4,))
  names(ass_eff)[c(1:4)]<-c("inf evansi",	"clean",	"coinfect",	"inf urticae")
  row.names(ass_eff)<-c("evansi", "urticae")
  ass_eff[1,] <- c(0.723903106,	0.615286425,	0.564057885,	0.480798595)
  ass_eff[2,] <- c(0.707579888,	0.627503111,	0.576787937,	0.45)



  frac_LxUx_ass=0.15 ## this is a number added or substracted to the mean to define Lx and Ux for equations A3 and A4 in Moya-Lara^no et al 2012



  ###applying equations A3 and A4 in Moya-Laraño et al. 2012
  ranges2<-as.data.frame(matrix(nrow=8,ncol=2,))       #
   names (ranges2)<-c("Lx","Ux")
  row.names (ranges2)<-c ("evansi_coinf", "urticae_coinf", "evansi_urticae","urticae_urticae","evansi_evansi","urticae_evansi", "evansi_clean", "urticae_clean")
  ranges2[1,] <- ranges(ass_eff[1,3]-frac_LxUx_ass,ass_eff[1,3]+frac_LxUx_ass,phi_ass)
  ranges2[2,] <- ranges(ass_eff[2,3]-frac_LxUx_ass,ass_eff[2,3]+frac_LxUx_ass,phi_ass)
  ranges2[3,] <- ranges(ass_eff[1,4]-frac_LxUx_ass,ass_eff[1,4]+frac_LxUx_ass,phi_ass)
  ranges2[4,] <- ranges(ass_eff[2,4]-frac_LxUx_ass,ass_eff[2,4]+frac_LxUx_ass,phi_ass)
  ranges2[5,] <- ranges(ass_eff[1,1]-frac_LxUx_ass,ass_eff[1,1]+frac_LxUx_ass,phi_ass)
  ranges2[6,] <- ranges(ass_eff[2,1]-frac_LxUx_ass,ass_eff[2,1]+frac_LxUx_ass,phi_ass)
  ranges2[7,] <- ranges(ass_eff[1,2]-frac_LxUx_ass,ass_eff[1,2]+frac_LxUx_ass,phi_ass)
  ranges2[8,] <- ranges(ass_eff[2,2]-frac_LxUx_ass,ass_eff[2,2]+frac_LxUx_ass,phi_ass)


  ############ WORLD ################
  surface=18*20*50  ##number of leaflets * width * length  (in mm2)
  cell_value=22.58 ## value in ug of a 1mm2 which we call cell for historical reasons of Rhea
  R<-surface*cell_value/10 ###in ug- Info in Rant 2004 - file mite_growth_curves.xls 1st sheet in Salomé's computer
  n_plants=40  ##plants are merely in a row, so dispersal goes necessarily from plant 1 to 2 and so on.
  plant<-as.data.frame(matrix(nrow=n_plants,ncol=3))
  names(plant)<-c("plant_ID","R","P")
  plant$plant_ID<-1:n_plants
  plant$R<-rep(R,n_plants) ##now biomass is R and R gets updated as food is taken from the environment
  plant$P<-rep(0,n_plants)
  plant$age<-rep(0,n_plants)##age is considered from the time of infestation
  plant$infect_date=0
  plant$K_date=0
  plant$switching="on"

  ###################################


  ####building the dispersal trait
  ####minimal food necessary for an adult female divided by assimilation efficiency for
  ####each type of environment - this gives the energy taken by one female
  ###this number 11.0863636 comes from growth_mass above
  min_P<- c(11.0863636/ 0.7239031, 11.0863636/0.6275031, 11.0863636/0.5767879,11.0863636/0.4807986)
  min_P<-as.data.frame(t(min_P))
  names(min_P)<-names(ass_eff)

  max_R_P<-c(R- min_P [1,1],R- min_P[1,2] ,R- min_P[1,3],R- min_P[1,4])
  max_R_P<-  as.data.frame(t(max_R_P))
  names(max_R_P)<-names(ass_eff)


  ####to establish how close to R-P=0 is the mean of the distribution for the disperal R-P trait
  pro_range<- 1/4


  mean_R_P<-c(max_R_P[1,1]*pro_range, max_R_P[1,2]*pro_range,max_R_P[1,3]*pro_range,
  max_R_P[1,4]*pro_range)
  mean_R_P<- as.data.frame(t(mean_R_P))
  names(mean_R_P)<-names(ass_eff)
  upperlevel<- c(mean_R_P[1,1]*2, mean_R_P[1,2]*2,mean_R_P[1,3]*2,mean_R_P[1,4]*2)
  upperlevel<- as.data.frame(t(upperlevel))
  names(upperlevel)<-names(ass_eff)


  ranges3<-as.data.frame(matrix(nrow=4,ncol=2,))       #
   names (ranges3)<-c("Lx","Ux")
  row.names (ranges3)<-c ("evansi","clean", "coinf", "urticae")

  ###applying eq A3 and A4 - minimum R-P is zero... assume animals won't stay in a patch w/o resources
  ranges3[1,] <- ranges(0,upperlevel[,1],phi_R_P)##*tuner_for_dispersal
  ranges3[2,] <- ranges(0,upperlevel[,2],phi_R_P)##*tuner_for_dispersal
  ranges3[3,] <- ranges(0,upperlevel[,3],phi_R_P)##*tuner_for_dispersal
  ranges3[4,] <- ranges(0,upperlevel[,4],phi_R_P)##*tuner_for_dispersal


  ###pseudo-values (genetics) obtained for a single trait in 50000 individuals
  ###- Arale did it!!! this is for 20 loci, changing the number of loci involves
  ### recalculating these values - so, this is a parameter with which we cannot play
  Min_te<-7.12
  Max_te<-13.68
  Min_tu<-7.12
  Max_tu<-13.68

  ###traits per module
  #both1 - R_P vs assimm eff - play with rho=0 vs rho=0.99 neg
  #both2 - male damage vs sex ratio - rho=0


  ###drift - a single non-functional trait to monitor drift
  bichinhos$drift<-0
  bichinhos$ass<-0
  bichinhos$R_P<-0
  bichinhos$m_damage<-0
  bichinhos$g_sex_ratio<-0



  ####interpolating pseudo-values (genetics) to ecological traits (phenotypes)
  for(i in 1:(length(bichinhos$sp))){
  bichinhos[i,"drift"]<- both3$phenotypes$phenotype1[i] ##assigns drift from single trait
  if (bichinhos$sp[i]=="te" && sum( bichinhos_all[1:2,])>0 && sum(bichinhos_all[3:4,])>0){     # coinfection
    bichinhos[i,"ass"]<- interpol2(genotypes$both1$phenotypes$phenotype1[i],Min_te, Max_te,ranges2[1,1],ranges2[1,2] )
    bichinhos[i,"R_P"]<- interpol2(genotypes$both1$phenotypes$phenotype2[i],Min_te, Max_te,ranges3[3,1],ranges3[3,2])
     }
  if (bichinhos$sp[i]=="tu" && sum(bichinhos_all[1:2,])>0 && sum(bichinhos_all[3:4,])>0){
    bichinhos[i,"ass"]<- interpol2(genotypes$both1$phenotypes$phenotype1[i],Min_tu, Max_tu,ranges2[2,1],ranges2[2,2] )
    bichinhos[i,"R_P"]<- interpol2(genotypes$both1$phenotypes$phenotype2[i],Min_tu, Max_tu,ranges3[3,1],ranges3[3,2] )
    }
  if (bichinhos$sp[i]=="te" && sum(bichinhos_all[1:2,])>0 && sum(bichinhos_all[3:4,])==0){     #infection urticae
    bichinhos[i,"ass"]<- interpol2(genotypes$both1$phenotype$phenotype1[i],Min_te, Max_te,ranges2[3,1],ranges2[3,2] )
    bichinhos[i,"R_P"]<- interpol2(genotypes$both1$phenotypes$phenotype2[i],Min_te, Max_te,ranges3[4,1],ranges3[4,2] )
    }
   if (bichinhos$sp[i]=="tu" && sum(bichinhos_all[1:2,])>0 && sum(bichinhos_all[3:4,])==0){
    bichinhos[i,"ass"]<-interpol2(genotypes$both1$phenotypes$phenotype1[i],Min_tu, Max_tu,ranges2[4,1],ranges2[4,2] )
    bichinhos[i,"R_P"]<- interpol2(genotypes$both1$phenotypes$phenotype2[i],Min_tu, Max_tu,ranges3[4,1],ranges3[4,2] )
    }
  if (bichinhos$sp[i]=="te" && sum(bichinhos_all[1:2,])==0 && sum(bichinhos_all[3:4,])>0){        #infection evansi
    bichinhos[i,"ass"]<- interpol2(genotypes$both1$phenotypes$phenotype1[i],Min_te, Max_te,ranges2[5,1],ranges2[5,2] )
    bichinhos[i,"R_P"]<- interpol2(genotypes$both1$phenotypes$phenotype2[i],Min_te, Max_te,ranges3[1,1],ranges3[1,2] )
    }
   if( bichinhos$sp[i]=="tu" && sum(bichinhos_all[1:2,])==0 && sum(bichinhos_all[3:4,])>0){
    bichinhos[i,"ass"]<- interpol2(genotypes$both1$phenotypes$phenotype1[i],Min_tu, Max_tu,ranges2[6,1],ranges2[6,2] )
    bichinhos[i,"R_P"]<- interpol2(genotypes$both1$phenotypes$phenotype2[i],Min_tu, Max_tu,ranges3[1,1],ranges3[1,2])
    }
   if( bichinhos$sp[i]=="te" && sum(bichinhos_all[1:2,])==0 && sum(bichinhos_all[3:4,])==0){        #clean
    bichinhos[i,"ass"]<- interpol2(genotypes$both1$phenotypes$phenotype1[i],Min_te, Max_te,ranges2["ass",1],ranges2["ass",2])
    bichinhos[i,"R_P"]<- interpol2(genotypes$both1$phenotypes$phenotype2[i],Min_te, Max_te,ranges3[2,1],ranges3[2,2]  )
    }
   if (bichinhos$sp[i]=="tu" && sum(bichinhos_all[1:2,])==0 && sum(bichinhos_all[3:4,])==0){
    bichinhos[i,"ass"]<- interpol2(genotypes$both1$phenotypes$phenotype1[i],Min_tu, Max_tu,ranges2[8,1],ranges2[8,2] )
    bichinhos[i,"R_P"]<- interpol2(genotypes$both1$phenotypes$phenotype2[i],Min_tu, Max_tu,ranges3[2,1],ranges3[2,2] )
    }
  }


  ###male damage from "Figuras e Tabela v2_Jordi.xls"
  ###if for Te the last male in the FIRST DAY is a Tu then 36% decrease in the remaining fecundity
  ###if for Tu the last male in the SECOND DAY is a Te then 53% decrease in the sex ratio
  ###of the remaining offspring
  
  if(version=="EMPIRICAL"|version=="NULL"){
  ranges4<-ranges(0.31,0.41,phi_m_dmg) ##damage male Tu on female Te - fecundity -> we model effect on available body mass
  ranges5<-ranges(0.43,0.63,phi_m_dmg) ##damage male Te on female Tu - sex ratio
  }
  
  if(version=="SCATTERED"){
  ranges4<-ranges(0.08,0.18,phi_m_dmg) ##damage male Tu on female Te - fecundity -> we model effect on available body mass
  ranges5<-ranges(0.20,0.30,phi_m_dmg) ##damage male Te on female Tu - sex ratio
  }
  
  
  bichinhos$m_damage<-ifelse(bichinhos$sp=="tu",interpol2(genotypes$both2$phenotypes$phenotype1,Min_te, Max_te,
  ranges4[1],ranges4[2]),interpol2(genotypes$both2$phenotypes$phenotype1,Min_te, Max_te,
  ranges5[1],ranges5[2]))


  ###sex ratio with genetic basis
  ranges6<-ranges(0.6,0.9,phi_sexr) ##probability of being female - for no variation we can do (0.75,0.75)
  
  bichinhos$g_sex_ratio<-interpol2(genotypes$both2$phenotypes$phenotype2,Min_te, Max_te,ranges6[1],ranges6[2])
  mean(bichinhos$g_sex_ratio)


  #initialitations of counters and states
  #############################################################
  #############################################################
  #############################################################
  #############################################################
  bichinhos$position<-1
  bichinhos$day_food<-0
  bichinhos$mass<-0
  bichinhos$quiesc[bichinhos$instar!="egg"]<-0
  n_eggs<-nrow(bichinhos[bichinhos$instar=="egg",])
  bichinhos$quiesc[bichinhos$instar=="egg"]<-round(runif(n_eggs,min=0.51,max=7.49)) ##so eggs do not have to wait systematically 7 days.
  bichinhos$adult_age=1 ###ifelse(bichinhos$instar=="adult"&bichinhos$sex=="female",1,0)
  bichinhos$alive=1
  bichinhos$sex_ratio<-bichinhos$g_sex_ratio
  bichinhos$gen_dam=0
  bichinhos$generation=1
  bichinhos$death_date=0
  bichinhos$migrate_date=0

  ##length_for_matings<-length(bichinhos) ##to assess the position of the first mating counter (1 to 10)

  #mating counters per day
  vectorM<-c(15,10,9,8,7,6,5,4,3,3)
  vectorF<-c(3,1,0,0,0,0,0,0,0,0)

  ###line for adding times of each mating for each female at adult age 1!!!
  bichinhos$t_mate1=0
  bichinhos$t_mate2=0
  bichinhos$t_mate3=0

  for(i in 1:nrow(bichinhos)){
    if(bichinhos$sex[i]=="female"){
      bichinhos[i,(length(bichinhos)-2):length(bichinhos)]<-sort(sample(1:12, 3, replace = FALSE))
    }
  }

  ##the following is to control the species identity of the mating sequence to calculate interference effects on females
  length_for_update<-length(bichinhos) ##to assess the position of the last variable included in bichinhos


  bichinhos$mate1_SP_ID<-as.character(NA)
  bichinhos$mate2_SP_ID<-as.character(NA)
  bichinhos$mate3_SP_ID<-as.character(NA)
  bichinhos$mate4_SP_ID<-as.character(NA)

  bichinhos$mate_ID<-as.character(NA)
  ##############################################################
  ##############################################################
  ##############################################################
  ##############################################################


  ###initialize mass using matrix of masses...
  mass_ini<- as.data.frame(matrix(nrow=10,ncol=4,))
  names(mass_ini)[c(1:4)]<-c("sex",	"instar",	"crowded",	"uncrowded")
  mass_ini[1,] <- c("female", "egg", 1.23,1.23)
  mass_ini[2,] <- c("male", "egg", 1.23,1.23)
  mass_ini[3,] <- c("female", "larva", 1.72,1.72)
  mass_ini[4,] <- c("male", "larva", 1.72,1.72)
  mass_ini[5,] <- c("female", "proto", 3.83,3.83)
  mass_ini[6,] <- c("male", "proto", 3.83,3.83)
  mass_ini[7,] <- c("female", "deuto", 4.32,10.95)##we take this as the mass of an adult female that cannot lay eggs yet
  mass_ini[8,] <- c("male", "deuto", 4.05,3.18)
  mass_ini[9,] <- c("female", "adult", 20.95,	24.5)
  mass_ini[10,] <- c("male", "adult", 4.72,	3.81)
  mass_ini$sex<-as.character(mass_ini$sex)
  mass_ini$instar<-as.character(mass_ini$instar)
  mass_ini$crowded<-as.numeric(mass_ini$crowded)
  mass_ini$uncrowded<-as.numeric(mass_ini$uncrowded)

  ###assigning initial masses
  for(i in 1:nrow(bichinhos)){
          bichinhos$mass[i]<- mass_ini[mass_ini$sex==bichinhos$sex[i]&
          mass_ini$instar==bichinhos$instar[i],
          colnames(mass_ini)==density]
  }

  bichinhos$mass<-ifelse(bichinhos$instar=="egg",bichinhos$mass+0.49,
    bichinhos$mass)


  #### its the mass necessary to become adult/ass_eff  - data in Mitchell 1973 - table is on damage sheet in mite_growth_curves.xls
  new_damage<- as.data.frame(matrix(nrow=8,ncol=4,))       # only for coinfection situation
  names (new_damage)<-c("egg","larva","proto","deuto")
  row.names (new_damage)<-c ("male_u_evansi", "male_u_urticae",
  "female_u_evansi","female_u_urticae","male_c_evansi","male_c_urticae",
  "female_c_evansi" , "female_c_urticae")
  new_damage[1,] <- c(4.573999, 3.705293, 0.000000, 1.116907)
  new_damage[2,] <- c(4.473048, 3.623515, 0.000000, 1.092256)
  new_damage[3,] <- c(41.25463, 40.38593, 36.64518, 24.02236)
  new_damage[4,] <- c(40.34412, 39.49458, 35.83639, 23.49217)
  new_damage[5,] <- c(6.187308, 5.318603, 1.577852, 1.187821)
  new_damage[6,] <- c(6.050751, 5.201218, 1.543028, 1.161605)
  new_damage[7,] <- c(34.96095, 34.09225, 30.35149, 27.70992)
  new_damage[8,] <- c(34.18934, 33.33981, 29.68162, 27.09835)

  ##the following probabilities come from paper1 graphs.xls
  ##timings have been changed to probabilities
  p_matings<-as.data.frame(matrix(nrow=8,ncol=2))
  names(p_matings)<-c("t_0","t_24")
  row.names(p_matings)<-c("tu_tu_tu","tu_tu_te","tu_te_tu","tu_te_te","te_te_te","te_te_tu",
  "te_tu_te","te_tu_tu")
  p_matings[1,]<-c(0.48,0.23)
  p_matings[2,]<-c(0.39,0.33)
  p_matings[3,]<-c(0.70,0.59)
  p_matings[4,]<-c(0,0)
  p_matings[5,]<-c(0.38,0.26)
  p_matings[6,]<-c(0.41,0.27)
  p_matings[7,]<-c(1,0.94)
  p_matings[8,]<-c(0,0)





  ################ MAIN FUNCTION ###########################
  instars<-c("egg","larva","proto","deuto","adult") ##this is in the order of growth precisely for growth

  table_bich<-table(bichinhos$instar, bichinhos$sex, bichinhos$sp)


  ###this is for debugging purposes only
  genotypes2<-genotypes
  bichinhos2<-bichinhos
#  bichinhos<-bichinhos2
#  genotypes<-genotypes2
#

break_all="off"

culling_freq=30 ###each 30 days is approx. the as-infected-lifespan of an infected plant
total_cullings=round(steps/culling_freq)
culling_counter=1
next_culling=culling_freq

######main function
for(i in 1:steps){

    if(scenario=="EXCLUSION"){
      ## atributes the "order" to the individuals- to randomize the movement-
      ###because they are in bichinhos arranged by species
      bichinhos$order<-runif(nrow(bichinhos))
      #sorts the individuals by random order number instead of by In_ID
      bichinhos<-sort_df(bichinhos,"order")
    }else{   ##culling by species to study trait evolution scenario=="FORCED"
      if(next_culling==i){
        #
        ## atributes the "order" to the individuals- to randomize the movement-
        ###because they are in bichinhos arranged by species
        bichinhos$order<-runif(nrow(bichinhos))
        #sorts the individuals by random order number instead of by In_ID
        bichinhos<-sort_df(bichinhos,"order")
        
        kept_tu<-subset(bichinhos,sp=="tu"&sex=="female"&alive==1)
        kept_te<-subset(bichinhos,sp=="te"&sex=="female"&alive==1)
        
        
        ##females to keep  tu & te
        if((nrow(kept_tu)>=n_turticae)&(nrow(kept_te))>=(n_tevansi)){
         kept_tu<-kept_tu[1:n_turticae,]
         kept_te<-kept_te[1:n_tevansi,]
        }else{#if not enough keep what is left taking the minimum btw species so no bias
         el_min<-min(nrow(kept_tu),nrow(kept_te))
         kept_tu<-kept_tu[1:el_min,]
         kept_te<-kept_te[1:el_min,]
        }
        
 
        bichinhos$alive<-ifelse(bichinhos$In_ID %in% c(kept_tu$In_ID,kept_te$In_ID),1,0)
        
        culling_counter=culling_counter+1
        next_culling=culling_freq*culling_counter
      }
    }

#    
 ##active_plants... those with mites assuming they are colonized in a row - last plant colonized last in row
 active_plants<-as.integer(sort(unique(bichinhos$position[bichinhos$alive==1]))) ##careful, it did not come ordered by default

 ##plant_loops<-length(active_plants)
 
 plant$infect_date<-ifelse(plant$infect_date==0 & plant$plant_ID %in% active_plants,i,plant$infect_date)
 
 for(p in active_plants[1]:active_plants[length(active_plants)]){

    plant$age[p]=plant$age[p]+1 ##plant age refers to days of plant live since infection, to calculate approx. crowded conditions at the 17 according to Mitchell 1973

    if(plant$age[p]>16){ ##plant changes to "crowded" situation after 17 days of infection (Mitchell suggests 16-18).
        density="crowded"
    }else{
        density="uncrowded"
    }


    sub_bich<-subset(bichinhos,quiesc==0&alive==1&position==p)

    if(nrow(sub_bich)>0){

    #####################
    ##FORAGING ALGORITHM (includes death from starvation)
    #####################

    ##removing non-feeders to increase speed
    only_feeders<-subset(sub_bich,!instar %in% "egg")
    non_feeders<-subset(sub_bich,instar %in% "egg")

    if(nrow(only_feeders)>0&plant$R[p]>0){

        for(j in 1:nrow(only_feeders)){

           ###the dammage on the plant is the mass day (requirement) divided by assim efficiency
           if(plant$R[p]>0){
           only_feeders$day_food[j]<-mass_day[mass_day$sex %in% only_feeders$sex[j] &
           mass_day$density==density,colnames(mass_day)==only_feeders$instar[j]]/only_feeders$ass[j]

           ##update biomass (R) for the plant
           plant$R[p]=plant$R[p]-only_feeders$day_food[j]

           ##net gain after assimilation efficiency
           only_feeders$day_food[j]=only_feeders$day_food[j]*only_feeders$ass[j]
           only_feeders$mass[j]=only_feeders$mass[j]+only_feeders$day_food[j]
           }else{
              ##if just one does not feed because there is no food all die right the way "break"
              ##however, adult females at age one can still migrate so we recover them and will die the next day

              only_feeders$alive<-ifelse(only_feeders$sex %in% "female"&only_feeders$instar %in% "adult"&
              only_feeders$adult_age==1,1,0)
              non_feeders$alive=0
              only_feeders$death_date<-ifelse(only_feeders$alive==0,i,only_feeders$death_date)
              non_feeders$death_date=i
              plant$K_date[p]=i
              break

           }

        }

      }else{

      ##rekill to be sure that no left-overs are there
      ###if just one does not feed because there is no food all die right the way "break"
         ##however, adult females at age one can still migrate so we recover them and will die the next day
         if(plant$R[p]<=0){
         only_feeders$alive<-ifelse(only_feeders$sex %in% "female"&only_feeders$instar %in% "adult"&
         only_feeders$adult_age==1,1,0)
         non_feeders$alive=0
         only_feeders$death_date<-ifelse(only_feeders$alive==0,i,only_feeders$death_date)
         non_feeders$death_date=i
         plant$K_date[p]=i
         }

      }

    #recovering sub_bich with all animals for the day
    if(nrow(only_feeders)>0){
    sub_bich<-rbind(only_feeders,non_feeders)
    }else{
    sub_bich<-non_feeders
    }
    ################ END OF FORAGING ALGORITHM

    rm(only_feeders)
    rm(non_feeders)

  #################
  #MATING ALGORITHM
  #################
  ##remove momentarily non_matees from sub_bich as to speed up the code

  ##this is necessary for the reproductive algorithm and to allow removing the non-adults for speed up
  sub_bich$mass_eggs_day<-0

  sub_bich$mass_eggs_day<-ifelse(sub_bich$sex %in% "female"&sub_bich$instar %in% "adult",
  sub_bich$mass-10.95,sub_bich$mass_eggs_day)  ##10.95 ug = mass of female deutonymph

  sub_bich$male_eggs<-0
  sub_bich$female_eggs<-0

  sub_bich$eggs_day<-round(sub_bich$mass_eggs_day/1.23)

  ##split files to speed up
  non_matees<-subset(sub_bich,!instar %in% "adult")
  matees<-subset(sub_bich,instar %in% "adult")

  rm(sub_bich)

  if(nrow(matees)>0){
  male_matees<-matees[matees$sex %in% "male",]
  female_matees<-matees[matees$sex %in% "female"&matees$adult_age<3,]

  #the following only occurs if there are males and females available to mate
  #################################################################
  #################################################################
  if(nrow(male_matees)>0&nrow(female_matees)>0){
      male_matees$matings_today=0

      for(j in 1:nrow(male_matees)){
        male_matees$matings_today[j]<-vectorM[male_matees$adult_age[j]]
      }


      female_matees$matings_today=0

      for(j in 1:nrow(female_matees)){
        female_matees$matings_today[j]<-vectorF[female_matees$adult_age[j]]
      }

      male_vector<-as.vector(matrix(nrow=10000,ncol=1))
      vector_position=1

      for(j in 1:nrow(male_matees)){
          number<-male_matees$matings_today[j]
          for(k in 1:number){
             male_vector[vector_position]<-male_matees$In_ID[j]
             vector_position=vector_position+1
          }
      }

      male_vector<-male_vector[!is.na(male_vector)]
      male_vector

      female_vector<-as.vector(matrix(nrow=10000,ncol=1))
      vector_position=1

      for(j in 1:nrow(female_matees)){
          number<-female_matees$matings_today[j]
          for(k in 1:number){
             female_vector[vector_position]<-female_matees$In_ID[j]
             vector_position=vector_position+1
          }
      }

      female_vector<-female_vector[!is.na(female_vector)]
      female_vector

      male_vector<-as.data.frame(male_vector)
      male_vector$rand<-runif(nrow(male_vector))
      male_vector<-sort_df(male_vector,"rand")

      female_vector<-as.data.frame(female_vector)
      female_vector$rand<-runif(nrow(female_vector))
      female_vector<-sort_df(female_vector,"rand")

      if(nrow(male_vector)>nrow(female_vector)){
        pairs_matrix<-data.frame(female_vector[,1],male_vector[1:nrow(female_vector),1])
      }else{
        pairs_matrix<-data.frame(female_vector[1:nrow(male_vector),1],male_vector[,1])
      }

      names(pairs_matrix)<-c("female_ID","male_ID")

      ##print(pairs_matrix)
      ####end mating assingment #################
      ###########################################

      ##providing the male mating species and id to the female - updating sub_bich
      ############################################################################
      for(j in 1:nrow(matees)){

          if(matees$sex[j] %in% "female"&matees$instar[j] %in% "adult"&matees$adult_age[j]==1){

             sub_females<-pairs_matrix[pairs_matrix$female_ID %in% matees$In_ID[j],]

             if(nrow(sub_females)==3){
             matees[j,(length_for_update+1):(length_for_update+3)]<-
             as.character(sub_females$male_ID)
             }
             if(nrow(sub_females)==2){
             matees[j,(length_for_update+1):(length_for_update+2)]<-
             as.character(sub_females$male_ID)
             }
             if(nrow(sub_females)==1){
             matees[j,length_for_update+1]<-
             as.character(sub_females$male_ID)
             }


          }

          if(matees$sex[j] %in% "female"&matees$instar[j] %in% "adult"&matees$adult_age[j]==2){
             sub_females<-pairs_matrix[pairs_matrix$female_ID %in% matees$In_ID[j],]

             if(nrow(sub_females)==1){
             matees[j,length_for_update+4]<-
             as.character(sub_females$male_ID)
             }
          }

        }


   } ##end if(nrow(female_matees)>0)

  ######## re-doing sub_bich according to the probabilities of actually mating with
  ######## each pair, taking into account first and second day matings.

  matees<-reassign(matees)


  #######################################
  #######################################
  ####################################### END MATING ALGORITHM

  ###REPRODUCTIVE ALGORITHM
  ##According to Paper2... number of eggs per day depend also on available energy

  for(j in 1:nrow(matees)){
    if(matees$sex[j] %in% "female"&matees$instar[j] %in% "adult"&matees$eggs_day[j]>0){

      the_list<-mate_effect(matees,j)

        if(length(the_list)==1){
          matees[j,]<-the_list[[1]]  ###this is if mass was not enough to lay eggs
        }else{
          matees[j,]<-the_list[[1]]

          genotypes<-the_list[[2]]

          bichinhos<-rbind(bichinhos,the_list[[3]])

          matees$mass[j]<-matees$mass[j]-matees$mass_eggs_day[j]##only if the eggs are laid is the mass lost!!!

          dups<-nrow(bichinhos[duplicated(bichinhos$In_ID),])
          if(dups>0){
           print(c(i,j,dups))
          }

        }

    }
  }

  }##end if(nrow(matees)>0)

  ##recover non_matees in sub_bich
  if(nrow(matees)>0){
  sub_bich<-rbind(matees,non_matees)
  }else{
  sub_bich<-rbind(matees,non_matees)
  }
  ########


  rm(matees)
  rm(non_matees)

     ################################################
 ################################################
 ################################################
 ################################################

  ##updating the overall data table with all animals
  for(j in 1:nrow(bichinhos)){
          if(nrow(sub_bich[sub_bich$In_ID %in% bichinhos$In_ID[j],c("eggs","quiesc","instar","mass","adult_age","mate1_SP_ID","mate2_SP_ID","mate3_SP_ID","mate4_SP_ID","mate_ID","alive","position","death_date")])>0){            ##"eggs_day","male_eggs","female_eggs"

          bichinhos[j,c("eggs","quiesc","instar","mass","adult_age","mate1_SP_ID","mate2_SP_ID","mate3_SP_ID","mate4_SP_ID","mate_ID","alive","position","death_date")]<-
          sub_bich[sub_bich$In_ID %in% bichinhos$In_ID[j],c("eggs","quiesc","instar","mass","adult_age","mate1_SP_ID","mate2_SP_ID","mate3_SP_ID","mate4_SP_ID","mate_ID","alive","position","death_date")]                        ##"eggs_day","male_eggs","female_eggs"

          }
  }



    rm(sub_bich)
    }##if(nrow(sub_bich)>0)


 
    ##mass to calculate projected secondary production
    pseudo_bich<-bichinhos[bichinhos$position==p&bichinhos$alive==1,]

    

    table_ind<-table(pseudo_bich$instar, pseudo_bich$sex, pseudo_bich$sp)
    
   
    table_ind<-fix_table(table_ind,table_bich)


    ##calculating secondary production living out of that plant
    plant$P[p]<-P_function(table_ind)##projected secondary production in the system to assess dispersal

    
    rm(table_ind)

    

    if(plant$plant_ID[p]>1 & plant$K_date[p]==i & plant$switching[p]=="on"){
#      break_all="on"
#      break

      n_to_remove<-nrow(bichinhos[bichinhos$position==(plant$plant_ID[p]-1),])
      IDs_to_remove<-bichinhos$In_ID[bichinhos$position %in% (plant$plant_ID[p]-1)]
      print(c(plant$plant_ID[p],plant$plant_ID[p]-1,n_to_remove))
      result<-shorten(genotypes,bichinhos,p,n_to_remove,IDs_to_remove)
      
      genotypes<-result[[1]]
      bichinhos<-result[[2]]
      
      plant$switching[p]="off"
      
    }

#if(break_all=="on"){
#break
#}
##
    



  }##end of plant loop

#if(break_all=="on"){
#break
#}
##
  #################################
  ################################
  #MARE INDIVIDUALS GROW AND MOLT
  ###############################
  ###############################

  ##removing young eggs to increase speed
#  young_eggs<-subset(sub_bich,instar=="egg")  ##&quiesc>1
#  moltees<-subset(sub_bich,instar!="egg")     ##|(instar=="egg"&quiesc==1)

no_moltees<-subset(bichinhos,quiesc>0|alive==0)

juvs<-subset(bichinhos,!instar %in% "egg"&!instar %in% "adult"&quiesc==0&alive==1)

  eggs<-subset(bichinhos,instar %in% "egg"&quiesc==0&alive==1)

  if(nrow(eggs)>0){
  eggs$instar<-"larva"
  eggs$quiesc<-round(runif(1,min=1,max=2))
  }

adults<-subset(bichinhos,instar %in% "adult"&quiesc==0&alive==1)
adults$adult_age=adults$adult_age+1



  if(nrow(juvs)>0){
    for(j in 1:nrow(juvs)){

    the_instar<-juvs$instar[j]


        next_instar<-instars[which(instars==the_instar)+1]


        if (juvs$mass[j] > mass_ini[mass_ini$sex==juvs$sex[j] & mass_ini$instar==next_instar,
        colnames(mass_ini)==density]) {  ##this is assuming   uncrowdiness - test sensitivity or change to crowdiness dynamically
          juvs$instar[j] <- next_instar
          if(!next_instar %in% "adult"){
          juvs$quiesc[j]=round(runif(1,min=1,max=2))
          }else{
          juvs$quiesc[j]=0
          juvs$adult_age[j]=1
          }
        }



    }
  }



bichinhos<-rbind(juvs,adults,eggs,no_moltees)

rm(juvs)
rm(adults)
rm(eggs)
rm(no_moltees)



###kill males and females from aging -
  bichinhos$alive<-with(bichinhos,ifelse(sex %in% "male" & adult_age>=10,0,ifelse(sex %in% "female" & adult_age>=30,0,alive)))

##update quiescence for next day
  bichinhos$quiesc<-ifelse(bichinhos$quiesc>0,bichinhos$quiesc-1,bichinhos$quiesc)


####DISPERSAL ALGORITHM#####
###only mated females at age one migrate - Mitchell 1973
## age==2 is really age==1 because they just updated their age but had no chance to migrate
## at age==1 would not work because they do not have mates assigned, so it is just to fix this problem
fem_age_1<-subset(bichinhos,adult_age==2&sex %in% "female"&instar %in% "adult"&!is.na(mate_ID))
others<-subset(bichinhos,adult_age!=2|!sex %in% "female"|!instar %in% "adult"|is.na(mate_ID))

  if(nrow(fem_age_1)>0){
  for(j in 1:nrow(fem_age_1)){

       here=fem_age_1$position[j]
       right_spot="off"

#       print(c("age adult female to migrate",fem_age_1$In_ID[j],fem_age_1$adult_age[j]))
#       if(i>66&i<69){
#          write.table(c(i,"age adult female to migrate",fem_age_1$In_ID[j],fem_age_1$adult_age[j]),"file.txt",append=T)
#       }
       while(here <= active_plants[length(active_plants)] & right_spot %in% "off"){
#            la_k<-plant$R[plant$plant_ID %in% fem_age_1$position[j]]
#            la_p<-plant$P[plant$plant_ID %in% fem_age_1$position[j]]
#            print(c(la_k,la_p))
            eval<-fem_age_1$R_P[j]>(plant$R[plant$plant_ID %in% fem_age_1$position[j]]-
            plant$P[plant$plant_ID %in% fem_age_1$position[j]])
            ##a low trait R_P value means that females migrate very late and stay in competitive environments
            if(eval){
              fem_age_1$position[j]=fem_age_1$position[j]+1
              fem_age_1$migrate_date[j]<-i
            }else{
            right_spot="on"
            }


       ##print(c(here,right_spot,eval,fem_age_1$position[j]))
       here=here+1
       }

  }
  }

bichinhos<-rbind(fem_age_1,others)

rm(fem_age_1)
rm(others)


print(i)
} ##end of main function


saveRDS(genotypes,paste("rm_genotypes_",1000,".rds",sep=""))
write.table(bichinhos,paste("rm_bichinhos_",1000,".txt",sep=""))



