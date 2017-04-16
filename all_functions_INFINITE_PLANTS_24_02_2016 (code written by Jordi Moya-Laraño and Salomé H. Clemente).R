chrom<-function()
{

    correlosome<-matrix(1:(n_traits*n_loci),nrow=n_loci,ncol=n_traits)

    new_order<-runif(n_traits*n_loci)
    order<-1:(n_traits*n_loci)

    for_order<-data.frame(new_order,order)

    final_order<-sort_df(for_order,"new_order")

    chromosome<-correlosome    ###the positions coming from the correlosome will
                               ####  be maintained 1 to 100 loci in the entire genome
                               #### as genome is 4 traits + 1 drift trait

    for(i in 1:n_traits){
      if(i==1){
        chromosome[1:n_loci]<-final_order$order[1:n_loci]
      }else{
        chromosome[((i*n_loci)-n_loci+1):(i*n_loci)]<-final_order$order[((i*n_loci)-n_loci+1):(i*n_loci)]
      }
    }

return(chromosome)
}

crea_loci<-function(n_loci,n_alleles)
{

      crom<-list(loc_num=as.data.frame(matrix(nrow=n_loci,ncol=1,1:n_loci)),
      allel_ID=as.data.frame(matrix(nrow=n_loci,ncol=n_alleles)),
      allel_val=as.data.frame(matrix(nrow=n_loci,ncol=n_alleles)))
      names(crom[[1]])[1]<-"loc_num"   ## Para cambiar los nombres, primer elemento del data.frame y el 2o es el V1 (redundante pero necesario).


        for (i in 1:n_loci){   ##locus i
          for (j in 1:n_alleles){  ##allele j
          crom$allel_ID[i,j]<-paste("ID",i,j,sep="_")
          min1=runif(1,min=0,max=0.5)
          max1=runif(1,min=0.5,max=1)
          crom$allel_val[i,j]<-runif(1,min=min1,max=max1)
          }
        }

      return(crom)
}

  ##########################

sex_det<-function(data,sex_ratio)
{
  for (i in 1:nrow(data)){
    p_sex_det=runif(1,min=0,max=1)
      if(p_sex_det<sex_ratio){
       data$sex[i]<-"female"
      }else{
       data$sex[i]<-"male"
      }
  }
return(data)
}

male_func<-function(lista,i,rho,mod_type,n_loci)
{
      ####for males
      lista$phenotypes$phenotype1[i]<-sum(lista$trt1_cr1_val[i,])

      diffrho=1-rho ###to induce the desired effect and parallel this to genetic correlation

      if(rho==0){
        lista$phenotypes$phenotype2[i]<-sum(lista$trt2_cr1_val[i,])
      }

      if(rho>0){
        if (mod_type=="neg"){
        lista$phenotypes$phenotype2[i]<-sum(lista$trt2_cr1_val[i,1:round(diffrho*n_loci,0)])+
        sum(1-lista$trt1_cr1_val[i,(round(diffrho*n_loci,0)+1):n_loci])
        }

        if (mod_type=="pos"){
        lista$phenotypes$phenotype2[i]<-sum(lista$trt2_cr1_val[i,1:round(diffrho*n_loci,0)])+
        sum(lista$trt1_cr1_val[i,(round(diffrho*n_loci,0)+1):n_loci])
        }
      }

return(lista)
}###closing the male function

  ###Below the order for dominance follows the order of the IDs.
  ###Therefore, the expression comes always from the largest ID.
  chrom_expres<-function(lista,i)   #data,
  {
    ##  for(i in 1:nrow(data)){
        for(j in 1:n_loci){
        if(order(c(lista$trt1_cr1_ID[i,j],lista$trt1_cr2_ID[i,j]))[2]==2){
        lista$expressing1[i,j]=lista$trt1_cr2_val[i,j]
        }else{
        lista$expressing1[i,j]=lista$trt1_cr1_val[i,j]
        }
        if(order(c(lista$trt2_cr1_ID[i,j],lista$trt2_cr2_ID[i,j]))[2]==2){
        lista$expressing2[i,j]=lista$trt2_cr2_val[i,j]
        }else{
        lista$expressing2[i,j]=lista$trt2_cr1_val[i,j]
        }
        }
   ##   }
  return(lista)
  }

  female_func<-function(lista,i,rho,mod_type,n_loci)
  {
      ####for females
      lista$phenotypes$phenotype1[i]<-sum(lista$expressing1[i,])

      diffrho=1-rho ###to induce the desired effect and parallel this to genetic correlation

      if(rho==0){
        lista$phenotypes$phenotype2[i]<-sum(lista$expressing2[i,])
      }

      if(rho>0){
        if (mod_type=="neg"){
        lista$phenotypes$phenotype2[i]<-sum(lista$expressing2[i,1:round(diffrho*n_loci,0)])+
        sum(1-lista$expressing1[i,(round(diffrho*n_loci,0)+1):n_loci])
        }

        if (mod_type=="pos"){
        lista$phenotypes$phenotype2[i]<-sum(lista$expressing2[i,1:round(diffrho*n_loci,0)])+
        sum(lista$expressing1[i,(round(diffrho*n_loci,0)+1):n_loci])
        }
      }

  return(lista)
  }###closing the female function



   ##SAMPLING CHROMOSOMES FOR THE INDIVIDUALS FOR EACH 2-VARIABLE MODULE
  crea_module<-function(data,mod_type,loci,n_loci,n_alleles,rho)
   {


      lista<-list(ID=as.integer(1:nrow(data)),sex=data$sex,
      trt1_cr1_ID=matrix(nrow=nrow(data),ncol=n_loci),
      trt1_cr1_val=matrix(nrow=nrow(data),ncol=n_loci),
      trt1_cr2_ID=matrix(nrow=nrow(data),ncol=n_loci),
      trt1_cr2_val=matrix(nrow=nrow(data),ncol=n_loci),
      trt2_cr1_ID=matrix(nrow=nrow(data),ncol=n_loci),
      trt2_cr1_val=matrix(nrow=nrow(data),ncol=n_loci),
      trt2_cr2_ID=matrix(nrow=nrow(data),ncol=n_loci),
      trt2_cr2_val=matrix(nrow=nrow(data),ncol=n_loci),
      phenotypes=as.data.frame(matrix(nrow=nrow(data),ncol=2,0)),
      expressing1=matrix(nrow=nrow(data),ncol=n_loci),
      expressing2=matrix(nrow=nrow(data),ncol=n_loci))

      ###assigning alleles at random for each individual, locus and each of the three traits in the module
      for (i in 1:nrow(data)){
        for (j in 1:n_loci){
          lista$trt1_cr1_ID[i,j]<-loci$allel_ID[j,as.integer(round(runif(1,min=1,max=n_alleles)))]
          lista$trt1_cr2_ID[i,j]<-loci$allel_ID[j,as.integer(round(runif(1,min=1,max=n_alleles)))]
          lista$trt2_cr1_ID[i,j]<-loci$allel_ID[j,as.integer(round(runif(1,min=1,max=n_alleles)))]
          lista$trt2_cr2_ID[i,j]<-loci$allel_ID[j,as.integer(round(runif(1,min=1,max=n_alleles)))]

          ###assigning values according to the alleles selected
          lista$trt1_cr1_val[i,j]<-loci$allel_val[j,lista$trt1_cr1_ID[i,j]==loci$allel_ID[j,]]
          lista$trt1_cr2_val[i,j]<-loci$allel_val[j,lista$trt1_cr2_ID[i,j]==loci$allel_ID[j,]]
          lista$trt2_cr1_val[i,j]<-loci$allel_val[j,lista$trt2_cr1_ID[i,j]==loci$allel_ID[j,]]
          lista$trt2_cr2_val[i,j]<-loci$allel_val[j,lista$trt2_cr2_ID[i,j]==loci$allel_ID[j,]]

        }
      }




      names(lista$phenotypes)[1:2]<-c("phenotype1","phenotype2")


      for (i in 1:length(lista$ID)){
      if(!is.na(match(lista$sex[i],"male"))){
      lista<-male_func(lista,i,rho,mod_type,n_loci)
      }else{
      lista<-chrom_expres(lista,i)    ##data,
      lista<-female_func(lista,i,rho,mod_type,n_loci)
      }
      }


      print(corr.test(lista$phenotypes[lista$sex=="male",]))

      print(corr.test(lista$phenotypes[lista$sex=="female",]))

  return(lista)
  }


  #############################################################################
  #############################################################################
  #############################################################################
  #############################################################################
  #############################################################################
  #############################################################################
  #for single trait
  male_func_one_trait<-function(lista,i,mod_type,n_loci)
  {
      ####for males

      lista$phenotypes$phenotype1[i]<-sum(lista$trt1_cr1_val[i,])

  return(lista)
  }###closing the male function

  chrom_expres_trait<-function(lista,i)
  {
  #    for(i in 1:nrow(data)){
        for(j in 1:n_loci){
        if(order(c(lista$trt1_cr1_ID[i,j],lista$trt1_cr2_ID[i,j]))[2]==2){
        lista$expressing1[i,j]=lista$trt1_cr2_val[i,j]
        }else{
        lista$expressing1[i,j]=lista$trt1_cr1_val[i,j]
        }
       }
  #    }
  return(lista)
  }

  female_func_one_trait<-function(lista,i,mod_type,n_loci)
  {
      ####for females

      lista$phenotypes$phenotype1[i]<-sum(lista$expressing1[i,])


  return(lista)
  }###closing the female function

  crea_trait<-function(data,loci,n_loci,n_alleles)
      {

      lista<-list(ID=as.integer(1:nrow(data)),sex=data$sex,
      trt1_cr1_ID=matrix(nrow=nrow(data),ncol=n_loci),
      trt1_cr1_val=matrix(nrow=nrow(data),ncol=n_loci),
      trt1_cr2_ID=matrix(nrow=nrow(data),ncol=n_loci),
      trt1_cr2_val=matrix(nrow=nrow(data),ncol=n_loci),
      phenotypes=as.data.frame(matrix(nrow=nrow(data),ncol=1,0)),
      expressing1=matrix(nrow=nrow(data),ncol=n_loci))

      for (i in 1:nrow(data)){
        for (j in 1:n_loci){
          lista$trt1_cr1_ID[i,j]<-loci$allel_ID[j,as.integer(round(runif(1,min=1,max=n_alleles)))]
          lista$trt1_cr2_ID[i,j]<-loci$allel_ID[j,as.integer(round(runif(1,min=1,max=n_alleles)))]
          }
      }


      ###assigning values according to the alleles selected
      for (i in 1:nrow(data)){
        for (j in 1:n_loci){
          lista$trt1_cr1_val[i,j]<-loci$allel_val[j,lista$trt1_cr1_ID[i,j]==loci$allel_ID[j,]]
          lista$trt1_cr2_val[i,j]<-loci$allel_val[j,lista$trt1_cr2_ID[i,j]==loci$allel_ID[j,]]
         }
      }


     names(lista$phenotypes)<-"phenotype1"

     for (i in 1:length(lista$ID)){
      if(!is.na(match(lista$sex[i],"male"))){
      lista<-male_func_one_trait(lista,i,mod_type,n_loci)
      }else{
      lista<-chrom_expres_trait(lista,i)
      lista<-female_func_one_trait(lista,i,mod_type,n_loci)
      }
      }

      trait<-lista

      return(trait)
  }


ranges<-function(Lx,Ux,phi){
    #eq A3
    lx=Lx+phi*((Ux-Lx)/2)
    #eq A4
    ux=Ux-phi*((Ux-Lx)/2)
return(c(lx,ux))
}

###assining realistic phenotypes to genotypes
interpol2<-function(x,minx,maxx,miny,maxy)  ##interpolation for values
{
  y=miny+((x-minx)*((maxy-miny)/(maxx-minx)))
  return(y)
}

###to fill gaps in table_ind prior to calculating P
fix_table<-function(table_ind,table_bich)
{
      table_zeroes<-table_bich

      names_for_1<-rownames(table_bich)
      names_for_2<-colnames(table_bich)
      names_for_3<-c("te","tu")

      ###this creates a full table with zeroes
      for(for1 in 1:length(names_for_1)){
       for(for2 in 1:length(names_for_2)){
        for(for3 in 1:length(names_for_3)){

          table_zeroes[names_for_1[for1],names_for_2[for2],names_for_3[for3]]<-0

         }
        }
       }

      if(dim(table_ind)[1]!=dim(table_zeroes)[1]|dim(table_ind)[2]!=dim(table_zeroes)[2]|dim(table_ind)[3]!=dim(table_zeroes)[3]){

      ###this adds the values of the incomplete table to the table full of zeroes
      for(for1 in 1:length(names_for_1)){
       for(for2 in 1:length(names_for_2)){
        for(for3 in 1:length(names_for_3)){

       texto<-try(table_ind[names_for_1[for1],names_for_2[for2],names_for_3[for3]],silent=T)

       if(!strsplit(as.character(texto)," ")[[1]][1] %in% "Error"){

       table_zeroes[names_for_1[for1],names_for_2[for2],names_for_3[for3]]<-table_ind[names_for_1[for1] ,names_for_2[for2],names_for_3[for3]]

       }

        }
       }
      }
      table_ind<-table_zeroes
      }else{
      table_ind<-table_ind
      }

return(table_ind)
}

####how much energy the current individuals on the plant will take from the system
P_function<-function(table_ind)
{
    ##table_ind<-fix_table(table_ind)
    if(density=="uncrowded"){
      P<- table_ind[1,1,1]* growth_mass[2,5] / ass_eff [1,3] + table_ind[1,2,1]* growth_mass[1,5] /
      ass_eff [1,3] + table_ind[2,1,1]* new_damage[3,4] + table_ind[2,2,1]* new_damage [1,4] +
      table_ind[3,1,1]* new_damage[3,1] + table_ind[3,2,1]* new_damage [1,1] + table_ind[4,1,1]*
      new_damage[3,2] + table_ind[4,2,1]* new_damage [1,2] + table_ind[5,1,1]* new_damage[3,3] +
      table_ind[5,2,1]* new_damage [1,3] +
      table_ind[1,1,2]* growth_mass[2,5] / ass_eff [2,3] + table_ind[1,2,2]* growth_mass[1,5] /
      ass_eff [2,3]+ table_ind[2,1,2]* new_damage[4,4] + table_ind[2,2,2]* new_damage [2,4] +
      table_ind[3,1,2]* new_damage[4,1] + table_ind[3,2,2]* new_damage [2,1] + table_ind[4,1,2]*
      new_damage[4,2] + table_ind[4,2,2]* new_damage [2,2] + table_ind[5,1,2]* new_damage[4,3] +
      table_ind[5,2,2]* new_damage [2,3]
    }else{
      P<- table_ind[1,1,1]* growth_mass[4,5] / ass_eff [1,3] + table_ind[1,2,1]* growth_mass[3,5] /
      ass_eff [1,3] + table_ind[2,1,1]* new_damage[3+4,4] + table_ind[2,2,1]* new_damage [1+4,4] +
      table_ind[3,1,1]* new_damage[3+4,1] + table_ind[3,2,1]* new_damage [1+4,1] + table_ind[4,1,1]*
      new_damage[3+4,2] + table_ind[4,2,1]* new_damage [1+4,2] + table_ind[5,1,1]* new_damage[3+4,3] +
      table_ind[5,2,1]* new_damage [1+4,3] +
      table_ind[1,1,2]* growth_mass[4,5] / ass_eff [2,3] + table_ind[1,2,2]* growth_mass[3,5] /
      ass_eff [2,3]+ table_ind[2,1,2]* new_damage[4+4,4] + table_ind[2,2,2]* new_damage [2+4,4] +
      table_ind[3,1,2]* new_damage[4+4,1] + table_ind[3,2,2]* new_damage [2+4,1] + table_ind[4,1,2]*
      new_damage[4+4,2] + table_ind[4,2,2]* new_damage [2+4,2] + table_ind[5,1,2]* new_damage[4+4,3] +
      table_ind[5,2,2]* new_damage [2+4,3]
    }
return(P)
}



###matings occur according the the observed probability in Clemente et al. paper1
reassign<-function(sub_bich)
{
      for(j in 1:nrow(sub_bich)){
       if(sub_bich$sex[j]=="female"){
          if(sub_bich$adult_age[j]==1){
               first_sp<-strsplit(sub_bich[j,length_for_update+1],"_")[[1]][1]
               second_sp<-strsplit(sub_bich[j,length_for_update+2],"_")[[1]][1]
               third_sp<-strsplit(sub_bich[j,length_for_update+3],"_")[[1]][1]
               first_second_str<-paste(sub_bich$sp[j],first_sp,second_sp,sep="_")
               second_third_str<-paste(sub_bich$sp[j],second_sp,third_sp,sep="_")

               if(!is.na(first_sp)&!is.na(second_sp)){
               p_1_vs_2<-p_matings[rownames(p_matings)==first_second_str,1]
                 r1<-runif(1)
                 if(p_1_vs_2<r1){ ##if second mating does not occur

                    sub_bich[j,length_for_update+2]<-NA
                 }
               }

               if(!is.na(second_sp)&!is.na(third_sp)){
               p_2_vs_3<-p_matings[rownames(p_matings)==second_third_str,1]
                 r2<-runif(1)
                 if(p_2_vs_3<r2){ ##if third mating does not occur
                 sub_bich[j,length_for_update+3]<-NA
                 }

               }
          }

          if(sub_bich$adult_age[j]==2){ ##if mate for second day is already assigned

             ####determine which is the last male truly mating!
             a<-!is.na(sub_bich[j,(length_for_update+1):(length_for_update+3)])[1]
             b<-!is.na(sub_bich[j,(length_for_update+1):(length_for_update+3)])[2]
             c<-!is.na(sub_bich[j,(length_for_update+1):(length_for_update+3)])[3]

             if(!a==F&!b==F&!c==F){

             last<-max(which(!is.na(sub_bich[j,(length_for_update+1):(length_for_update+3)])))

               last_sp<-strsplit(sub_bich[j,length_for_update+last],"_")[[1]][1]
               fourth_sp<-strsplit(sub_bich[j,length_for_update+4],"_")[[1]][1]
               last_fourth_str<-paste(sub_bich$sp[j],last_sp,fourth_sp,sep="_")

               if(!is.na(last_sp)&!is.na(fourth_sp)){
                 p_last_vs_4<-p_matings[rownames(p_matings)==last_fourth_str,2]

                 r3<-runif(1)
                 if(p_last_vs_4<r3){ ##if last mating does not occur

                    sub_bich[j,length_for_update+3]<-NA
                 }
               }
             }

          }
       }
    }

return(sub_bich)
}

crlsm_to_chrom<-function(sub_bich,j)
{
    male_ID<-sub_bich$mate_ID[j]
    female_ID<-sub_bich$In_ID[j]
    male_sp<-bichinhos$sp[bichinhos$In_ID %in% male_ID]  ##use bichinhos because male can be dead at this point
    female_sp<-sub_bich$sp[sub_bich$In_ID %in% female_ID]


    if(female_sp=="tu"){
      to_order<-tu_chromosome
      correlosome1_ID<-to_order
      correlosome1_val<-to_order
      correlosome2_ID<-to_order
      correlosome2_val<-to_order
    }else{
      to_order<-te_chromosome
      correlosome1_ID<-to_order
      correlosome1_val<-to_order
      correlosome2_ID<-to_order
      correlosome2_val<-to_order
    }


    n_loops=round(n_traits/n_modules)+1

    ##females
    ##IDs
    for(k in 1:n_loops){
      if(k==1){ ##module 1
        correlosome1_ID[,k]<-genotypes[[k]]$trt1_cr1_ID[genotypes[[k]]$ID %in% female_ID]
        correlosome2_ID[,k]<-genotypes[[k]]$trt1_cr2_ID[genotypes[[k]]$ID %in% female_ID]
        correlosome1_ID[,k+1]<-genotypes[[k]]$trt2_cr1_ID[genotypes[[k]]$ID %in% female_ID]
        correlosome2_ID[,k+1]<-genotypes[[k]]$trt2_cr2_ID[genotypes[[k]]$ID %in% female_ID]
      }else{    ##module above one
        if(length(genotypes[[k]])>9){
          correlosome1_ID[,k*2-1]<-genotypes[[k]]$trt1_cr1_ID[genotypes[[k]]$ID %in% female_ID]
          correlosome2_ID[,k*2-1]<-genotypes[[k]]$trt1_cr2_ID[genotypes[[k]]$ID %in% female_ID]
          correlosome1_ID[,k*2]<-genotypes[[k]]$trt2_cr1_ID[genotypes[[k]]$ID %in% female_ID]
          correlosome2_ID[,k*2]<-genotypes[[k]]$trt2_cr2_ID[genotypes[[k]]$ID %in% female_ID]
        }else{
          correlosome1_ID[,k*2-1]<-genotypes[[k]]$trt1_cr1_ID[genotypes[[k]]$ID %in% female_ID]
          correlosome2_ID[,k*2-1]<-genotypes[[k]]$trt1_cr2_ID[genotypes[[k]]$ID %in% female_ID]
        }
      }
    }

    ##vals
    for(k in 1:n_loops){
      if(k==1){ ##module 1
        correlosome1_val[,k]<-genotypes[[k]]$trt1_cr1_val[genotypes[[k]]$ID %in% female_ID]
        correlosome2_val[,k]<-genotypes[[k]]$trt1_cr2_val[genotypes[[k]]$ID %in% female_ID]
        correlosome1_val[,k+1]<-genotypes[[k]]$trt2_cr1_val[genotypes[[k]]$ID %in% female_ID]
        correlosome2_val[,k+1]<-genotypes[[k]]$trt2_cr2_val[genotypes[[k]]$ID %in% female_ID]
      }else{    ##module above one
        if(length(genotypes[[k]])>9){
          correlosome1_val[,k*2-1]<-genotypes[[k]]$trt1_cr1_val[genotypes[[k]]$ID %in% female_ID]
          correlosome2_val[,k*2-1]<-genotypes[[k]]$trt1_cr2_val[genotypes[[k]]$ID %in% female_ID]
          correlosome1_val[,k*2]<-genotypes[[k]]$trt2_cr1_val[genotypes[[k]]$ID %in% female_ID]
          correlosome2_val[,k*2]<-genotypes[[k]]$trt2_cr2_val[genotypes[[k]]$ID %in% female_ID]
        }else{
          correlosome1_val[,k*2-1]<-genotypes[[k]]$trt1_cr1_val[genotypes[[k]]$ID %in% female_ID]
          correlosome2_val[,k*2-1]<-genotypes[[k]]$trt1_cr2_val[genotypes[[k]]$ID %in% female_ID]
        }
      }
    }


    ##transform correlosome in chromosome with randomly ordered genome - random but the same for all individuals of the species
    chromosome1_ID_ordered<-correlosome1_ID[order(to_order)]
    chromosome2_ID_ordered<-correlosome2_ID[order(to_order)]
    chromosome1_val_ordered<-correlosome1_val[order(to_order)]
    chromosome2_val_ordered<-correlosome2_val[order(to_order)]

    ##put in matrix form before cross-over
    dim(chromosome1_ID_ordered)<-c(n_loci,n_traits)
    dim(chromosome2_ID_ordered)<-c(n_loci,n_traits)
    dim(chromosome1_val_ordered)<-c(n_loci,n_traits)
    dim(chromosome2_val_ordered)<-c(n_loci,n_traits)

    catch_list<-crossover(chromosome1_ID_ordered,chromosome2_ID_ordered,chromosome1_val_ordered,
    chromosome2_val_ordered)

    ##transform it back to correlosome (ordered genome by trait)
    chromosome1_ID_reordered<-catch_list[[1]][order(order(to_order))]
    chromosome2_ID_reordered<-catch_list[[2]][order(order(to_order))]
    chromosome1_val_reordered<-catch_list[[3]][order(order(to_order))]
    chromosome2_val_reordered<-catch_list[[4]][order(order(to_order))]

    ##put in back in matrix form after cross-over
    dim(chromosome1_ID_reordered)<-c(n_loci,n_traits)
    dim(chromosome2_ID_reordered)<-c(n_loci,n_traits)
    dim(chromosome1_val_reordered)<-c(n_loci,n_traits)
    dim(chromosome2_val_reordered)<-c(n_loci,n_traits)

return(list(chromosome1_ID_reordered,chromosome2_ID_reordered,chromosome1_val_reordered,
chromosome2_val_reordered))
}

###this determines the gamete genetics from the males that mate with the females
###they are falsely 2n, as only cr1 is used throughout (it is just to make the code simpler)
###notice that they do not recombine
male_genetics<-function(sub_bich,j)
{


    male_ID<-sub_bich$mate_ID[j]
    male_sp<-bichinhos$sp[bichinhos$In_ID %in% male_ID]  ##use bichinhos because male can be dead at this point
    
    if(male_sp=="tu"){
      to_order<-tu_chromosome
      correlosome1_ID<-to_order
      correlosome1_val<-to_order
      correlosome2_ID<-to_order
      correlosome2_val<-to_order
    }else{
      to_order<-te_chromosome
      correlosome1_ID<-to_order
      correlosome1_val<-to_order
      correlosome2_ID<-to_order
      correlosome2_val<-to_order
    }


    n_loops=round(n_traits/n_modules)+1

    ##females
    ##IDs
    for(k in 1:n_loops){
      if(k==1){ ##module 1
        
#        esto<-genotypes[[k]]$trt1_cr1_ID[genotypes[[k]]$ID %in% male_ID]
#        print(esto)
        correlosome1_ID[,k]<-genotypes[[k]]$trt1_cr1_ID[genotypes[[k]]$ID %in% male_ID]
        correlosome2_ID[,k]<-genotypes[[k]]$trt1_cr2_ID[genotypes[[k]]$ID %in% male_ID]
        correlosome1_ID[,k+1]<-genotypes[[k]]$trt2_cr1_ID[genotypes[[k]]$ID %in% male_ID]
        correlosome2_ID[,k+1]<-genotypes[[k]]$trt2_cr2_ID[genotypes[[k]]$ID %in% male_ID]
      }else{    ##module above one
        if(length(genotypes[[k]])>9){
          correlosome1_ID[,k*2-1]<-genotypes[[k]]$trt1_cr1_ID[genotypes[[k]]$ID %in% male_ID]
          correlosome2_ID[,k*2-1]<-genotypes[[k]]$trt1_cr2_ID[genotypes[[k]]$ID %in% male_ID]
          correlosome1_ID[,k*2]<-genotypes[[k]]$trt2_cr1_ID[genotypes[[k]]$ID %in% male_ID]
          correlosome2_ID[,k*2]<-genotypes[[k]]$trt2_cr2_ID[genotypes[[k]]$ID %in% male_ID]
        }else{
          correlosome1_ID[,k*2-1]<-genotypes[[k]]$trt1_cr1_ID[genotypes[[k]]$ID %in% male_ID]
          correlosome2_ID[,k*2-1]<-genotypes[[k]]$trt1_cr2_ID[genotypes[[k]]$ID %in% male_ID]
        }
      }
    }

    ##vals
    for(k in 1:n_loops){
      if(k==1){ ##module 1
        correlosome1_val[,k]<-genotypes[[k]]$trt1_cr1_val[genotypes[[k]]$ID %in% male_ID]
        correlosome2_val[,k]<-genotypes[[k]]$trt1_cr2_val[genotypes[[k]]$ID %in% male_ID]
        correlosome1_val[,k+1]<-genotypes[[k]]$trt2_cr1_val[genotypes[[k]]$ID %in% male_ID]
        correlosome2_val[,k+1]<-genotypes[[k]]$trt2_cr2_val[genotypes[[k]]$ID %in% male_ID]
      }else{    ##module above one
        if(length(genotypes[[k]])>9){
          correlosome1_val[,k*2-1]<-genotypes[[k]]$trt1_cr1_val[genotypes[[k]]$ID %in% male_ID]
          correlosome2_val[,k*2-1]<-genotypes[[k]]$trt1_cr2_val[genotypes[[k]]$ID %in% male_ID]
          correlosome1_val[,k*2]<-genotypes[[k]]$trt2_cr1_val[genotypes[[k]]$ID %in% male_ID]
          correlosome2_val[,k*2]<-genotypes[[k]]$trt2_cr2_val[genotypes[[k]]$ID %in% male_ID]
        }else{
          correlosome1_val[,k*2-1]<-genotypes[[k]]$trt1_cr1_val[genotypes[[k]]$ID %in% male_ID]
          correlosome2_val[,k*2-1]<-genotypes[[k]]$trt1_cr2_val[genotypes[[k]]$ID %in% male_ID]
        }
      }
    }

return(list(correlosome1_ID,correlosome2_ID,correlosome1_val,correlosome2_val))
}



##for this function to be general we need a way to put the entire genome together
##from the lists and WARNING: remember is done for IDs only, need to have also the val thing
##done / need to use both1, both2, both3 instead and then new individuals are added to this
##list - open delete and save list as in the case of phenotypes... more important here!!!
crossover<-function(chromosome1_ID_ordered,chromosome2_ID_ordered,chromosome1_val_ordered,
    chromosome2_val_ordered)
{

    out_cr1_ID<-chromosome1_ID_ordered
    out_cr2_ID<-chromosome2_ID_ordered
    out_cr1_val<-chromosome1_val_ordered
    out_cr2_val<-chromosome2_val_ordered


    ##cross-over
    for(i in 1:n_traits){

        break_point=round(runif(1,min=0.51,max=n_loci+0.49)) ##to make tails really uniform after rounding


        if(break_point<n_loci){
          top_cr1_ID<-chromosome1_ID_ordered[1:break_point,i] ##top cr1
          top_cr2_ID<-chromosome2_ID_ordered[1:break_point,i] ##top cr2
          top_cr1_val<-chromosome1_val_ordered[1:break_point,i] ##top cr1
          top_cr2_val<-chromosome2_val_ordered[1:break_point,i] ##top cr2

          btm_cr1_ID<-chromosome1_ID_ordered[(break_point+1):n_loci,i] ##bottom cr1
          btm_cr2_ID<-chromosome2_ID_ordered[(break_point+1):n_loci,i] ##bottom cr2
          btm_cr1_val<-chromosome1_val_ordered[(break_point+1):n_loci,i] ##bottom cr1
          btm_cr2_val<-chromosome2_val_ordered[(break_point+1):n_loci,i] ##bottom cr2

          out_cr1_ID[,i]<-c(top_cr1_ID,btm_cr2_ID)
          out_cr2_ID[,i]<-c(top_cr2_ID,btm_cr1_ID)
          out_cr1_val[,i]<-c(top_cr1_val,btm_cr2_val)
          out_cr2_val[,i]<-c(top_cr2_val,btm_cr1_val)


        }else{
          out_cr1_ID[,i]<-chromosome1_ID_ordered[1:n_loci,i] ##bottom cr1
          out_cr2_ID[,i]<-chromosome2_ID_ordered[1:n_loci,i] ##bottom cr2
          out_cr1_val[,i]<-chromosome1_val_ordered[1:n_loci,i] ##bottom cr1
          out_cr2_val[,i]<-chromosome2_val_ordered[1:n_loci,i] ##bottom cr2

        }

    }



    ##and feed this back including vals to original module and trait structure...
    ##this will be the genetics of one gamete

return(list(out_cr1_ID,out_cr2_ID,out_cr1_val,out_cr2_val))
}

##male damaging effects plus egg laying depending on male identities
if(version=="EMPIRICAL"){
mate_effect<-function(sub_bich,j)
{

        
        first_sp<-strsplit(sub_bich[j,length_for_update+1],"_")[[1]][1]
        second_sp<-strsplit(sub_bich[j,length_for_update+2],"_")[[1]][1]
        third_sp<-strsplit(sub_bich[j,length_for_update+3],"_")[[1]][1]

        ######### if no males are available, they are forced to lay male eggs
        if(sub_bich$sex[j]=="female"&sub_bich$instar[j]=="adult"){
        if(is.na(first_sp)){
             
              sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
              if(sub_bich$male_eggs[j]>0){
                new_list<-reproduce(sub_bich,j)
                sub_bich[j,]<-new_list[[1]]
              }


        }else{
        

          if(sub_bich$adult_age[j]==1){


               if(!is.na(first_sp)){
               ## if first male is wrong species, they are forced to lay male eggs during the first hours
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }

                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }

                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }

                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }

                  }

               ## if first male is right species, both males and female eggs are laid
                  if(sub_bich$sp[j]=="te"&first_sp=="te"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


               }

               if(!is.na(first_sp)&!is.na(second_sp)){
                  ##if second males is still wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                 ##wrong male arrives, start interference effects
                 ##the tetetu combination, tu interferes with fecundity with females
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     p_hours=hours/12
                     male_damaging=sub_bich$mate2_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutute combination does not have interference in the first 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="te"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }

               ##when there is no second male mating in the second attemp, but there is a third male
               if(!is.na(first_sp)&is.na(second_sp)&!is.na(third_sp)){
                  ##if second males is still wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                 ##wrong male arrives, start interference effects
                 ##the tetetu combination, tu interferes with fecundity with females
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutute combination does not have interference in the first 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }


               ##effects of third male acting as actual third
               if(!is.na(first_sp)&!is.na(second_sp)&!is.na(third_sp)){
                  ##if third males is still wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                 ##wrong male arrives as third mate, start interference effects
                 ##the tetetetu combination, tu interferes with fecundity with females
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tututute combination does not have interference in the first 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }

          }


          ##AGE 2 or older - a fourth male is offered at day 2!!
          fourth_sp<-strsplit(as.character(sub_bich[j,length_for_update+4]),"_")[[1]][1]

          if(sub_bich$adult_age[j]>1&!is.na(fourth_sp)){ ##if mate for second day is already assigned

               ##if female has mated with three males previously
               if(!is.na(first_sp)&!is.na(second_sp)&!is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"&fourth_sp=="tu"){
                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  ###from here on we implement sex ratio with a random component
                  ###otherwise is impossible to go away from 2:1, and rather will approach 2:0 from rounding
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="te"){##&fourth_sp=="tu" fourth male is irrelevant
                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="tu"){##&third_sp=="tu" third is irrelevant
                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"){##&fourth_sp=="te" fourth is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){##&second_sp=="te"&third_sp=="tu" second and third irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"&fourth_sp=="te"){
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="te"){##&third_sp=="tu" third male is irrelevant
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){##&second_sp=="tu"&third_sp=="tu" second and third male is irrelevant
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  ##wrong male arrives as second mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                     male_damaging=sub_bich$mate2_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##wrong male arrives as third mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##wrong male arrives as third mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="tu"){
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


               }



               ##if female has mated with only two first males
               if(!is.na(first_sp)&!is.na(second_sp)&is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&fourth_sp=="tu"){
                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="tu"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrive at first mating and fourth has no effect
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){ ##&second_sp=="te" second is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){ ##&second_sp=="tu" second male is irrelevant
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="te"){  ## third male is irrelevant
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  ##wrong male arrives as second mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                     male_damaging=sub_bich$mate2_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

               }

               ##if female has mated with two males previously, but not second
               if(!is.na(first_sp)&is.na(second_sp)&!is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"&fourth_sp=="tu"){
                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){## &third_sp=="te" third is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="te"&fourth_sp=="te"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="tu"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }                  
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="te"){ ##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){  ##&third_sp=="tu" third male is irrelevant
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"&fourth_sp=="te"){  ## third male is irrelevant
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  ##wrong male arrives as third mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="tu"){
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

               }

               ##if female has mated only with the first male
               if(!is.na(first_sp)&is.na(second_sp)&is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&fourth_sp=="tu"){
                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"){##&fourth_sp="tu" fourth male is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }


          }##end of age 2 and fourth male

          if(sub_bich$adult_age[j]>1&is.na(fourth_sp)){ ##if fourth male is not available
                   ##if first male is available but there is no second nor third males
                   if(!is.na(first_sp)&is.na(second_sp)&is.na(third_sp)){
                    ## if first male is wrong species, they are forced to lay male eggs during the first hours
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                          if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                     ##if first male is the right species, both males and female eggs are laid
                      if(sub_bich$sp[j]=="te"&first_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }


                   }

                   ##if first male and second male available but there is no third male
                   if(!is.na(first_sp)&!is.na(second_sp)&is.na(third_sp)){
                      ##if second males is still wrong species, continue to lay male eggs
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##right male arrives, start laying females as well as males
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }


                     ##wrong male arrives, CONTINUE interference effects
                     ##the tetetu combination, tu interferes with fecundity with females
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                         male_damaging=sub_bich$mate2_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##the tutute combination does not have interference if mating occurs within the first 24 hours
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                   }

                   ##when there is no second male mating in the second attemp, but there is a third male
                   if(!is.na(first_sp)&is.na(second_sp)&!is.na(third_sp)){
                      ##if second males is still wrong species, CONTINUE to lay male eggs
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##right male arrived, CONTINUE laying females as well as males
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                     ##wrong male arrives, start interference effects
                     ##the tetetu combination, tu interferes with fecundity with females
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="tu"){
                         male_damaging=sub_bich$mate3_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##the tutute combination does not have interference in the first 24 hours
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                   }


                   ##effects of third male acting as actual third
                   if(!is.na(first_sp)&!is.na(second_sp)&!is.na(third_sp)){
                      ##if third males is still wrong species, continue to lay male eggs
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##right male arrives, start laying females as well as males
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }


                     ##wrong male arrives as third mate, start interference effects
                     ##the tetetetu combination, tu interferes with fecundity with females
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                         male_damaging=sub_bich$mate3_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="tu"){
                         male_damaging=sub_bich$mate3_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){##third male irrelevant
                         male_damaging=sub_bich$mate2_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      ##the tututute combination does not have interference in the first 24 hours
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&third_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="te"&third_sp=="tu"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                   }

              }##end of age 2 and NO fourth male

       

        }##end of else
      }##end of female and adult   

if(sub_bich$eggs_day[j]==0){
  return(list(sub_bich[j,]))
}

if(sub_bich$eggs_day[j]>0){
return(list(sub_bich[j,],new_list[[2]],new_list[[3]]))
}

}
}##END EMPIRICAL


if(version=="SCATTERED"){
##male damaging effects plus egg laying depending on male identities
mate_effect<-function(sub_bich,j)
{

        
        first_sp<-strsplit(sub_bich[j,length_for_update+1],"_")[[1]][1]
        second_sp<-strsplit(sub_bich[j,length_for_update+2],"_")[[1]][1]
        third_sp<-strsplit(sub_bich[j,length_for_update+3],"_")[[1]][1]

        ######### if no males are available, they are forced to lay male eggs
        if(sub_bich$sex[j]=="female"&sub_bich$instar[j]=="adult"){
        if(is.na(first_sp)){
             
              sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
              if(sub_bich$male_eggs[j]>0){
                new_list<-reproduce(sub_bich,j)
                sub_bich[j,]<-new_list[[1]]
              }


        }else{
        

          if(sub_bich$adult_age[j]==1){


               if(!is.na(first_sp)){
               ## if first male is wrong species, they are forced to lay male eggs during the first hours
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }

                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }

                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }

                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }

                  }

               ## if first male is right species, both males and female eggs are laid
                  if(sub_bich$sp[j]=="te"&first_sp=="te"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


               }

               if(!is.na(first_sp)&!is.na(second_sp)){
                  ##if second males is still wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                 ##wrong male arrives, start interference effects
                 ##the tetetu combination, tu interferes with fecundity with females
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     p_hours=hours/12
                     male_damaging=sub_bich$mate2_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutute combination does not have interference in the first 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="te"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate2_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }

               ##when there is no second male mating in the second attemp, but there is a third male
               if(!is.na(first_sp)&is.na(second_sp)&!is.na(third_sp)){
                  ##if second males is still wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                 ##wrong male arrives, start interference effects
                 ##the tetetu combination, tu interferes with fecundity with females
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutute combination does not have interference in the first 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }


               ##effects of third male acting as actual third
               if(!is.na(first_sp)&!is.na(second_sp)&!is.na(third_sp)){
                  ##if third males is still wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                 ##wrong male arrives as third mate, start interference effects
                 ##the tetetetu combination, tu interferes with fecundity with females
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tututute combination does not have interference in the first 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }

          }


          ##AGE 2 or older - a fourth male is offered at day 2!!
          fourth_sp<-strsplit(as.character(sub_bich[j,length_for_update+4]),"_")[[1]][1]

          if(sub_bich$adult_age[j]>1&!is.na(fourth_sp)){ ##if mate for second day is already assigned

               ##if female has mated with three males previously
               if(!is.na(first_sp)&!is.na(second_sp)&!is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"&fourth_sp=="tu"){
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  ###from here on we implement sex ratio with a random component
                  ###otherwise is impossible to go away from 2:1, and rather will approach 2:0 from rounding
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="te"){##&fourth_sp=="tu" fourth male is irrelevant
                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="tu"){##&third_sp=="tu" third is irrelevant
                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"){##&fourth_sp=="te" fourth is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){##&second_sp=="te"&third_sp=="tu" second and third irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"&fourth_sp=="te"){
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="te"){##&third_sp=="tu" third male is irrelevant
                     male_damaging=sub_bich$mate2_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){##&second_sp=="tu"&third_sp=="tu" second and third male is irrelevant
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  ##wrong male arrives as second mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                     male_damaging=sub_bich$mate2_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##wrong male arrives as third mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##wrong male arrives as third mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="tu"){
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


               }



               ##if female has mated with only two first males
               if(!is.na(first_sp)&!is.na(second_sp)&is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&fourth_sp=="tu"){
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="tu"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"){##&fourth_sp=="te" fourth is irrelevant here
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrive at first mating and fourth has no effect
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){ ##&second_sp=="te" second is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&fourth_sp=="tu"){
                     
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){ ##&second_sp=="tu" second male is irrelevant
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="te"){  ## third male is irrelevant
                     male_damaging=sub_bich$mate2_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  ##wrong male arrives as second mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                     male_damaging=sub_bich$mate2_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

               }

               ##if female has mated with two males previously, but not second
               if(!is.na(first_sp)&is.na(second_sp)&!is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"&fourth_sp=="tu"){
                  
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){## &third_sp=="te" third is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="te"&fourth_sp=="te"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="tu"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }                  
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="te"){ ##&fourth_sp=="te" fourth is irrelevant here
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){  ##&third_sp=="tu" third male is irrelevant
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"&fourth_sp=="te"){  ## third male is irrelevant
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  ##wrong male arrives as third mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="tu"){
                     male_damaging=sub_bich$mate3_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

               }

               ##if female has mated only with the first male
               if(!is.na(first_sp)&is.na(second_sp)&is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&fourth_sp=="tu"){
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"){##&fourth_sp="tu" fourth male is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&fourth_sp=="tu"){
                     
                     male_damaging=sub_bich$mate1_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){
                     male_damaging=sub_bich$mate4_SP_ID[j]
                     male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                     eggs<-sub_bich$eggs_day[j]
                     sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }


          }##end of age 2 and fourth male

          if(sub_bich$adult_age[j]>1&is.na(fourth_sp)){ ##if fourth male is not available
                   ##if first male is available but there is no second nor third males
                   if(!is.na(first_sp)&is.na(second_sp)&is.na(third_sp)){
                    ## if first male is wrong species, they are forced to lay male eggs during the first hours
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"){
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                          if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                     ##if first male is the right species, both males and female eggs are laid
                      if(sub_bich$sp[j]=="te"&first_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }


                   }

                   ##if first male and second male available but there is no third male
                   if(!is.na(first_sp)&!is.na(second_sp)&is.na(third_sp)){
                      ##if second males is still wrong species, continue to lay male eggs
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"){
                     
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##right male arrives, start laying females as well as males
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }


                     ##wrong male arrives, CONTINUE interference effects
                     ##the tetetu combination, tu interferes with fecundity with females
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                         male_damaging=sub_bich$mate2_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##the tutute combination does not have interference if mating occurs within the first 24 hours
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         
                         male_damaging=sub_bich$mate2_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                   }

                   ##when there is no second male mating in the second attemp, but there is a third male
                   if(!is.na(first_sp)&is.na(second_sp)&!is.na(third_sp)){
                      ##if second males is still wrong species, CONTINUE to lay male eggs
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"){
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##right male arrived, CONTINUE laying females as well as males
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                     ##wrong male arrives, start interference effects
                     ##the tetetu combination, tu interferes with fecundity with females
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="tu"){
                         male_damaging=sub_bich$mate3_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##the tutute combination does not have interference in the first 24 hours
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         
                         male_damaging=sub_bich$mate3_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                   }


                   ##effects of third male acting as actual third
                   if(!is.na(first_sp)&!is.na(second_sp)&!is.na(third_sp)){
                      ##if third males is still wrong species, continue to lay male eggs
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##right male arrives, start laying females as well as males
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)

                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }


                     ##wrong male arrives as third mate, start interference effects
                     ##the tetetetu combination, tu interferes with fecundity with females
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                         male_damaging=sub_bich$mate3_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="tu"){
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){##third male irrelevant
                         male_damaging=sub_bich$mate2_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$mass_eggs_day[j]<-sub_bich$mass_eggs_day[j]-sub_bich$mass_eggs_day[j]*male_damage
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      ##the tututute combination does not have interference in the first 24 hours
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         
                         male_damaging=sub_bich$mate3_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&third_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         
                         male_damaging=sub_bich$mate1_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="te"&third_sp=="tu"){
                         eggs<-sub_bich$eggs_day[j]
                         
                         male_damaging=sub_bich$mate2_SP_ID[j]
                         male_damage=bichinhos$m_damage[bichinhos$In_ID %in% male_damaging]
                         sub_bich$sex_ratio[j]=sub_bich$g_sex_ratio[j]-sub_bich$g_sex_ratio[j]*male_damage

                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                   }

              }##end of age 2 and NO fourth male

       

        }##end of else
      }##end of female and adult   

if(sub_bich$eggs_day[j]==0){
  return(list(sub_bich[j,]))
}

if(sub_bich$eggs_day[j]>0){
return(list(sub_bich[j,],new_list[[2]],new_list[[3]]))
}

}


}###END OF SCATTERED

if(version=="NULL"){
mate_effect<-function(sub_bich,j)
{

        
        first_sp<-strsplit(sub_bich[j,length_for_update+1],"_")[[1]][1]
        second_sp<-strsplit(sub_bich[j,length_for_update+2],"_")[[1]][1]
        third_sp<-strsplit(sub_bich[j,length_for_update+3],"_")[[1]][1]

        ######### if no males are available, they are forced to lay male eggs
        if(sub_bich$sex[j]=="female"&sub_bich$instar[j]=="adult"){
        if(is.na(first_sp)){
             
              sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
              if(sub_bich$male_eggs[j]>0){
                new_list<-reproduce(sub_bich,j)
                sub_bich[j,]<-new_list[[1]]
              }


        }else{
        

          if(sub_bich$adult_age[j]==1){


               if(!is.na(first_sp)){
               ## if first male is wrong species, they are forced to lay male eggs during the first hours
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }

                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }

                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }

                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }

                  }

               ## if first male is right species, both males and female eggs are laid
                  if(sub_bich$sp[j]=="te"&first_sp=="te"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"){
                     if(!is.na(second_sp)){
                     hours=sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate2[j]
                     }
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


               }

               if(!is.na(first_sp)&!is.na(second_sp)){
                  ##if second males is still wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                 ##wrong male arrives, start interference effects
                 ##the tetetu combination, tu interferes with fecundity with females
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     p_hours=hours/12
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutute combination does not have interference in the first 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="te"){
                     if(!is.na(third_sp)){
                     hours=sub_bich$t_mate2[j]-sub_bich$t_mate1[j]
                     }else{
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate1[j]
                     }

                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }

               ##when there is no second male mating in the second attemp, but there is a third male
               if(!is.na(first_sp)&is.na(second_sp)&!is.na(third_sp)){
                  ##if second males is still wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                 ##wrong male arrives, start interference effects
                 ##the tetetu combination, tu interferes with fecundity with females
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutute combination does not have interference in the first 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }


               ##effects of third male acting as actual third
               if(!is.na(first_sp)&!is.na(second_sp)&!is.na(third_sp)){
                  ##if third males is still wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$male_eggs[j]=round(p_hours*sub_bich$eggs_day[j])
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     sub_bich$mate_ID[j]=sub_bich$mate1_SP_ID[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                 ##wrong male arrives as third mate, start interference effects
                 ##the tetetetu combination, tu interferes with fecundity with females
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tututute combination does not have interference in the first 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                     hours=sub_bich$t_mate3[j]-sub_bich$t_mate2[j]
                     p_hours=hours/12
                     eggs<-p_hours*sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }

          }


          ##AGE 2 or older - a fourth male is offered at day 2!!
          fourth_sp<-strsplit(as.character(sub_bich[j,length_for_update+4]),"_")[[1]][1]

          if(sub_bich$adult_age[j]>1&!is.na(fourth_sp)){ ##if mate for second day is already assigned

               ##if female has mated with three males previously
               if(!is.na(first_sp)&!is.na(second_sp)&!is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"&fourth_sp=="tu"){
                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  ###from here on we implement sex ratio with a random component
                  ###otherwise is impossible to go away from 2:1, and rather will approach 2:0 from rounding
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="te"){##&fourth_sp=="tu" fourth male is irrelevant
                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="tu"){##&third_sp=="tu" third is irrelevant
                     sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"){##&fourth_sp=="te" fourth is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){##&second_sp=="te"&third_sp=="tu" second and third irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"&fourth_sp=="te"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="te"){##&third_sp=="tu" third male is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){##&second_sp=="tu"&third_sp=="tu" second and third male is irrelevant
                      eggs<-sub_bich$eggs_day[j]
                      if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  ##wrong male arrives as second mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##wrong male arrives as third mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##wrong male arrives as third mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="tu"){
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


               }



               ##if female has mated with only two first males
               if(!is.na(first_sp)&!is.na(second_sp)&is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&fourth_sp=="tu"){
                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="tu"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrive at first mating and fourth has no effect
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){ ##&second_sp=="te" second is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){ ##&second_sp=="tu" second male is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&fourth_sp=="te"){  ## third male is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  ##wrong male arrives as second mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

               }

               ##if female has mated with two males previously, but not second
               if(!is.na(first_sp)&is.na(second_sp)&!is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"&fourth_sp=="tu"){
                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){## &third_sp=="te" third is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="te"&fourth_sp=="te"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="tu"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }                  
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="te"){ ##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="te"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){  ##&third_sp=="tu" third male is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"&fourth_sp=="te"){  ## third male is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

                  ##wrong male arrives as third mate, CONTINUE interference effects
                  ##the tetetu combination has effect within the first 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="tu"){
                     sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                     eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }

               }

               ##if female has mated only with the first male
               if(!is.na(first_sp)&is.na(second_sp)&is.na(third_sp)){
                  ##if fourth male is still the wrong species, continue to lay male eggs
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&fourth_sp=="tu"){
                     sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&fourth_sp=="te"){
                    sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##right male arrives, start laying females as well as males
                  if(sub_bich$sp[j]=="te"&first_sp=="tu"&fourth_sp=="te"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"){##&fourth_sp=="te" fourth is irrelevant here
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  if(sub_bich$sp[j]=="te"&first_sp=="te"){##&fourth_sp="tu" fourth male is irrelevant
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }


                  if(sub_bich$sp[j]=="tu"&first_sp=="te"&fourth_sp=="tu"){
                     sub_bich$mate_ID[j]=sub_bich$mate4_SP_ID[j]
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                 ##wrong male arrives as fourth mate, start interference effects
                 ##the tetetetetu combination has no effect after 24h
                  if(sub_bich$sp[j]=="te"&first_sp=="te"&fourth_sp=="tu"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
                  ##the tutututute combination DOES HAVE an interference effect after 24 hours
                  if(sub_bich$sp[j]=="tu"&first_sp=="tu"&fourth_sp=="te"){
                     eggs<-sub_bich$eggs_day[j]
                     if(sub_bich$sex_ratio[j]>runif(1)){
                      sub_bich$female_eggs[j]=round(eggs)
                     }else{
                      sub_bich$male_eggs[j]=round(eggs)
                     }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                  }
               }


          }##end of age 2 and fourth male

          if(sub_bich$adult_age[j]>1&is.na(fourth_sp)){ ##if fourth male is not available
                   ##if first male is available but there is no second nor third males
                   if(!is.na(first_sp)&is.na(second_sp)&is.na(third_sp)){
                    ## if first male is wrong species, they are forced to lay male eggs during the first hours
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                          if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                     ##if first male is the right species, both males and female eggs are laid
                      if(sub_bich$sp[j]=="te"&first_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }

                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }


                   }

                   ##if first male and second male available but there is no third male
                   if(!is.na(first_sp)&!is.na(second_sp)&is.na(third_sp)){
                      ##if second males is still wrong species, continue to lay male eggs
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##right male arrives, start laying females as well as males
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate2_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }


                     ##wrong male arrives, CONTINUE interference effects
                     ##the tetetu combination, tu interferes with fecundity with females
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##the tutute combination does not have interference if mating occurs within the first 24 hours
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                   }

                   ##when there is no second male mating in the second attemp, but there is a third male
                   if(!is.na(first_sp)&is.na(second_sp)&!is.na(third_sp)){
                      ##if second males is still wrong species, CONTINUE to lay male eggs
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="tu"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##right male arrived, CONTINUE laying females as well as males
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                     ##wrong male arrives, start interference effects
                     ##the tetetu combination, tu interferes with fecundity with females
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&third_sp=="tu"){
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##the tutute combination does not have interference in the first 24 hours
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&third_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                   }


                   ##effects of third male acting as actual third
                   if(!is.na(first_sp)&!is.na(second_sp)&!is.na(third_sp)){
                      ##if third males is still wrong species, continue to lay male eggs
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                         sub_bich$male_eggs[j]=sub_bich$eggs_day[j]
                         if(sub_bich$male_eggs[j]>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      ##right male arrives, start laying females as well as males
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="te"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="tu"){
                         #sub_bich$mate_ID[j]=sub_bich$mate3_SP_ID[j]
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }


                     ##wrong male arrives as third mate, start interference effects
                     ##the tetetetu combination, tu interferes with fecundity with females
                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="te"&third_sp=="tu"){
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="te"&first_sp=="tu"&second_sp=="te"&third_sp=="tu"){
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="te"&first_sp=="te"&second_sp=="tu"){##third male irrelevant
                         sub_bich$eggs_day[j]<-round(sub_bich$mass_eggs_day[j]/1.23)
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      ##the tututute combination does not have interference in the first 24 hours
                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="tu"&third_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }
                      if(sub_bich$sp[j]=="tu"&first_sp=="te"&second_sp=="tu"&third_sp=="te"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                      if(sub_bich$sp[j]=="tu"&first_sp=="tu"&second_sp=="te"&third_sp=="tu"){
                         eggs<-sub_bich$eggs_day[j]
                         if(sub_bich$sex_ratio[j]>runif(1)){
                          sub_bich$female_eggs[j]=round(eggs)
                         }else{
                          sub_bich$male_eggs[j]=round(eggs)
                         }
                         if(round(eggs)>0){
                         new_list<-reproduce(sub_bich,j)
                         sub_bich[j,]<-new_list[[1]]
                         }
                      }

                   }

              }##end of age 2 and NO fourth male

       

        }##end of else
      }##end of female and adult   

if(sub_bich$eggs_day[j]==0){
  return(list(sub_bich[j,]))
}

if(sub_bich$eggs_day[j]>0){
return(list(sub_bich[j,],new_list[[2]],new_list[[3]]))
}

}

}###END OF NULL



make_both<-function(sp,last_ID,n_offs_m,n_offs_f,male_offs_male_gamete_genotype,
      female_offs_male_gamete_genotype,female_offs_female_gamete_genotype)
{
  n_loops=round(n_traits/n_modules)+1

  k=1
  while(k<=n_loops*2-1){
  if(k<n_loops*2-1){
   if(n_offs_f>0){
     list_f <- list(sp=rep(sp,n_offs_f),ID=paste(sp,(last_ID+n_offs_m+1):(last_ID+n_offs_m+n_offs_f),sep="_"),sex=rep("female",n_offs_f),
     trt1_cr1_ID=matrix(nrow=n_offs_f,ncol=n_loci),trt1_cr2_ID=matrix(nrow=n_offs_f,ncol=n_loci),
     trt1_cr1_val=matrix(nrow=n_offs_f,ncol=n_loci),trt1_cr2_val=matrix(nrow=n_offs_f,ncol=n_loci),
     trt2_cr1_ID=matrix(nrow=n_offs_f,ncol=n_loci),trt2_cr2_ID=matrix(nrow=n_offs_f,ncol=n_loci),
     trt2_cr1_val=matrix(nrow=n_offs_f,ncol=n_loci),trt2_cr2_val=matrix(nrow=n_offs_f,ncol=n_loci),
     phenotypes=as.data.frame(matrix(nrow=n_offs_f,ncol=2)),expressing1=matrix(nrow=n_offs_f,ncol=n_loci),
     expressing2=matrix(nrow=n_offs_f,ncol=n_loci))
     names(list_f$phenotypes)<-c("phenotype1","phenotype2")

     for(xx in 1:n_offs_f){
        list_f$trt1_cr1_ID[xx,]=female_offs_male_gamete_genotype[[xx]][[1]][,k]  ##male genotype always comes from 1st chromosome, second chromosome is aborted for real haplodiploidy
        list_f$trt1_cr2_ID[xx,]=female_offs_female_gamete_genotype[[xx]][[2]][,k]
        list_f$trt1_cr1_val[xx,]=female_offs_male_gamete_genotype[[xx]][[3]][,k]
        list_f$trt1_cr2_val[xx,]=female_offs_female_gamete_genotype[[xx]][[4]][,k]
        list_f$trt2_cr1_ID[xx,]=female_offs_male_gamete_genotype[[xx]][[1]][,k+1]  ##male genotype always comes from 1st chromosome, second chromosome is aborted for real haplodiploidy
        list_f$trt2_cr2_ID[xx,]=female_offs_female_gamete_genotype[[xx]][[2]][,k+1]
        list_f$trt2_cr1_val[xx,]=female_offs_male_gamete_genotype[[xx]][[3]][,k+1]
        list_f$trt2_cr2_val[xx,]=female_offs_female_gamete_genotype[[xx]][[4]][,k+1]

     }
   }

   if(n_offs_m>0){
     list_m <- list(sp=rep(sp,n_offs_m),ID=paste(sp,(last_ID+1):(last_ID+n_offs_m),sep="_"),sex=rep("male",n_offs_m),
     trt1_cr1_ID=matrix(nrow=n_offs_m,ncol=n_loci),trt1_cr2_ID=matrix(nrow=n_offs_m,ncol=n_loci),
     trt1_cr1_val=matrix(nrow=n_offs_m,ncol=n_loci),trt1_cr2_val=matrix(nrow=n_offs_m,ncol=n_loci),
     trt2_cr1_ID=matrix(nrow=n_offs_m,ncol=n_loci),trt2_cr2_ID=matrix(nrow=n_offs_m,ncol=n_loci),
     trt2_cr1_val=matrix(nrow=n_offs_m,ncol=n_loci),trt2_cr2_val=matrix(nrow=n_offs_m,ncol=n_loci),
     phenotypes=as.data.frame(matrix(nrow=n_offs_m,ncol=2)),expressing1=matrix(nrow=n_offs_m,ncol=n_loci),
     expressing2=matrix(nrow=n_offs_m,ncol=n_loci))
     names(list_m$phenotypes)<-c("phenotype1","phenotype2")
   

      for(xx in 1:n_offs_m){
          list_m$trt1_cr1_ID[xx,]=male_offs_male_gamete_genotype[[xx]][[1]][,k]  ##male genotype always comes from 1st chromosome, second chromosome is aborted for real haplodiploidy
          list_m$trt1_cr2_ID[xx,]=male_offs_male_gamete_genotype[[xx]][[2]][,k]
          list_m$trt1_cr1_val[xx,]=male_offs_male_gamete_genotype[[xx]][[3]][,k]
          list_m$trt1_cr2_val[xx,]=male_offs_male_gamete_genotype[[xx]][[4]][,k]
          list_m$trt2_cr1_ID[xx,]=male_offs_male_gamete_genotype[[xx]][[1]][,k+1]  ##male genotype always comes from 1st chromosome, second chromosome is aborted for real haplodiploidy
          list_m$trt2_cr2_ID[xx,]=male_offs_male_gamete_genotype[[xx]][[2]][,k+1]
          list_m$trt2_cr1_val[xx,]=male_offs_male_gamete_genotype[[xx]][[3]][,k+1]
          list_m$trt2_cr2_val[xx,]=male_offs_male_gamete_genotype[[xx]][[4]][,k+1]
      }
    
    }
  }else{
   if(n_offs_m>0){
     list_m <- list(sp=rep(sp,n_offs_m) ,ID=paste(sp,(last_ID+1):(last_ID+n_offs_m),sep="_"),sex=rep("male",n_offs_m),
     trt1_cr1_ID=matrix(nrow=n_offs_m,ncol=n_loci),trt1_cr2_ID=matrix(nrow=n_offs_m,ncol=n_loci),
     trt1_cr1_val=matrix(nrow=n_offs_m,ncol=n_loci),trt1_cr2_val=matrix(nrow=n_offs_m,ncol=n_loci),
     phenotypes=as.data.frame(matrix(nrow=n_offs_m,ncol=1)),expressing1=matrix(nrow=n_offs_m,ncol=n_loci))
     names(list_m$phenotypes)<-"phenotype1"
     
     for(xx in 1:n_offs_m){
        list_m$trt1_cr1_ID[xx,]=male_offs_male_gamete_genotype[[xx]][[1]][,k]  ##male genotype always comes from 1st chromosome, second chromosome is aborted for real haplodiploidy
        list_m$trt1_cr2_ID[xx,]=male_offs_male_gamete_genotype[[xx]][[2]][,k]
        list_m$trt1_cr1_val[xx,]=male_offs_male_gamete_genotype[[xx]][[3]][,k]
        list_m$trt1_cr2_val[xx,]=male_offs_male_gamete_genotype[[xx]][[4]][,k]
     }

   
   }
   if(n_offs_f>0){
     list_f <- list(sp=rep(sp,n_offs_f) ,ID=paste(sp,(last_ID+n_offs_m+1):(last_ID+n_offs_m+n_offs_f),sep="_"),sex=rep("female",n_offs_f),
     trt1_cr1_ID=matrix(nrow=n_offs_f,ncol=n_loci),trt1_cr2_ID=matrix(nrow=n_offs_f,ncol=n_loci),
     trt1_cr1_val=matrix(nrow=n_offs_f,ncol=n_loci),trt1_cr2_val=matrix(nrow=n_offs_f,ncol=n_loci),
     phenotypes=as.data.frame(matrix(nrow=n_offs_f,ncol=1)),expressing1=matrix(nrow=n_offs_f,ncol=n_loci))
     names(list_f$phenotypes)<-"phenotype1"
  
     for(xx in 1:n_offs_f){
        list_f$trt1_cr1_ID[xx,]=female_offs_male_gamete_genotype[[xx]][[1]][,k]  ##male genotype always comes from 1st chromosome, second chromosome is aborted for real haplodiploidy
        list_f$trt1_cr2_ID[xx,]=female_offs_female_gamete_genotype[[xx]][[2]][,k]
        list_f$trt1_cr1_val[xx,]=female_offs_male_gamete_genotype[[xx]][[3]][,k]
        list_f$trt1_cr2_val[xx,]=female_offs_female_gamete_genotype[[xx]][[4]][,k]
    }
   }
  }##else

   if(k==1){
      if(n_offs_m>0&n_offs_f>0){
        pre_geno1<-list(sp=c(list_m$sp,list_f$sp) ,ID=c(list_m$ID,list_f$ID),sex=c(list_m$sex,list_f$sex),
        trt1_cr1_ID=rbind(list_m$trt1_cr1_ID,list_f$trt1_cr1_ID), trt1_cr1_val=rbind(list_m$trt1_cr1_val,
        list_f$trt1_cr1_val), trt1_cr2_ID=rbind(list_m$trt1_cr2_ID,list_f$trt1_cr2_ID),
        trt1_cr2_val=rbind(list_m$trt1_cr2_val,list_f$trt1_cr2_val),trt2_cr1_ID=rbind(list_m$trt2_cr1_ID,
        list_f$trt2_cr1_ID), trt2_cr1_val=rbind(list_m$trt2_cr1_val,list_f$trt2_cr1_val),
        trt2_cr2_ID= rbind(list_m$trt2_cr2_ID,list_f$trt2_cr2_ID),trt2_cr2_val= rbind(list_m$trt2_cr2_val,
        list_f$trt2_cr2_val),phenotypes=rbind(list_m$phenotypes,list_f$phenotypes),
        expressing1=rbind(list_m$expressing1,list_f$expressing1),
        expressing2=rbind(list_m$expressing2,list_f$expressing2))
        if(length(pre_geno1$ID)<2){break}
      }
      if(n_offs_m==0&n_offs_f>0){
        pre_geno1<-list(sp=list_f$sp ,ID=list_f$ID,sex=list_f$sex,
        trt1_cr1_ID=list_f$trt1_cr1_ID, trt1_cr1_val=list_f$trt1_cr1_val, trt1_cr2_ID=list_f$trt1_cr2_ID,
        trt1_cr2_val=list_f$trt1_cr2_val,trt2_cr1_ID=list_f$trt2_cr1_ID, trt2_cr1_val=list_f$trt2_cr1_val,
        trt2_cr2_ID=list_f$trt2_cr2_ID,trt2_cr2_val=list_f$trt2_cr2_val,phenotypes=list_f$phenotypes,
        expressing1=list_f$expressing1,
        expressing2=list_f$expressing2)
      }
      if(n_offs_m>0&n_offs_f==0){
        pre_geno1<-list(sp=list_m$sp ,ID=list_m$ID,sex=list_m$sex,
        trt1_cr1_ID=list_m$trt1_cr1_ID, trt1_cr1_val=list_m$trt1_cr1_val, trt1_cr2_ID=list_m$trt1_cr2_ID,
        trt1_cr2_val=list_m$trt1_cr2_val,trt2_cr1_ID=list_m$trt2_cr1_ID, trt2_cr1_val=list_m$trt2_cr1_val,
        trt2_cr2_ID=list_m$trt2_cr2_ID,trt2_cr2_val=list_m$trt2_cr2_val,phenotypes=list_m$phenotypes,
        expressing1=list_m$expressing1,
        expressing2=list_m$expressing2)
      }

   }else{
      if(k<n_loops*2-1){
        if(n_offs_m>0&n_offs_f>0){
          pre_geno2<-list(sp=c(list_m$sp,list_f$sp) ,ID=c(list_m$ID,list_f$ID),sex=c(list_m$sex,list_f$sex),
          trt1_cr1_ID=rbind(list_m$trt1_cr1_ID,list_f$trt1_cr1_ID), trt1_cr1_val=rbind(list_m$trt1_cr1_val,
          list_f$trt1_cr1_val), trt1_cr2_ID=rbind(list_m$trt1_cr2_ID,list_f$trt1_cr2_ID),
          trt1_cr2_val=rbind(list_m$trt1_cr2_val,list_f$trt1_cr2_val),trt2_cr1_ID=rbind(list_m$trt2_cr1_ID,
          list_f$trt2_cr1_ID), trt2_cr1_val=rbind(list_m$trt2_cr1_val,list_f$trt2_cr1_val),
          trt2_cr2_ID= rbind(list_m$trt2_cr2_ID,list_f$trt2_cr2_ID),trt2_cr2_val= rbind(list_m$trt2_cr2_val,
          list_f$trt2_cr2_val),phenotypes=rbind(list_m$phenotypes,list_f$phenotypes),
          expressing1=rbind(list_m$expressing1,list_f$expressing1),
          expressing2=rbind(list_m$expressing2,list_f$expressing2))
        }
        if(n_offs_m==0&n_offs_f>0){
          pre_geno2<-list(sp=list_f$sp ,ID=list_f$ID,sex=list_f$sex,
          trt1_cr1_ID=list_f$trt1_cr1_ID, trt1_cr1_val=list_f$trt1_cr1_val, trt1_cr2_ID=list_f$trt1_cr2_ID,
          trt1_cr2_val=list_f$trt1_cr2_val,trt2_cr1_ID=list_f$trt2_cr1_ID, trt2_cr1_val=list_f$trt2_cr1_val,
          trt2_cr2_ID=list_f$trt2_cr2_ID,trt2_cr2_val=list_f$trt2_cr2_val,phenotypes=list_f$phenotypes,
          expressing1=list_f$expressing1,
          expressing2=list_f$expressing2)
        }
        if(n_offs_m>0&n_offs_f==0){
          pre_geno2<-list(sp=list_m$sp ,ID=list_m$ID,sex=list_m$sex,
          trt1_cr1_ID=list_m$trt1_cr1_ID, trt1_cr1_val=list_m$trt1_cr1_val, trt1_cr2_ID=list_m$trt1_cr2_ID,
          trt1_cr2_val=list_m$trt1_cr2_val,trt2_cr1_ID=list_m$trt2_cr1_ID, trt2_cr1_val=list_m$trt2_cr1_val,
          trt2_cr2_ID=list_m$trt2_cr2_ID,trt2_cr2_val=list_m$trt2_cr2_val,phenotypes=list_m$phenotypes,
          expressing1=list_m$expressing1,
          expressing2=list_m$expressing2)
        }
      }else{
        if(n_offs_m>0&n_offs_f>0){
          pre_geno3<-list(sp=c(list_m$sp,list_f$sp) ,ID=c(list_m$ID,list_f$ID),sex=c(list_m$sex,list_f$sex),
          trt1_cr1_ID=rbind(list_m$trt1_cr1_ID,list_f$trt1_cr1_ID), trt1_cr1_val=rbind(list_m$trt1_cr1_val,
          list_f$trt1_cr1_val), trt1_cr2_ID=rbind(list_m$trt1_cr2_ID,list_f$trt1_cr2_ID),
          trt1_cr2_val=rbind(list_m$trt1_cr2_val,list_f$trt1_cr2_val),phenotypes=rbind(list_m$phenotypes,list_f$phenotypes),
          expressing1=rbind(list_m$expressing1,list_f$expressing1))
        }
        if(n_offs_m==0&n_offs_f>0){
          pre_geno3<-list(sp=list_f$sp ,ID=list_f$ID,sex=list_f$sex,
          trt1_cr1_ID=list_f$trt1_cr1_ID, trt1_cr1_val=list_f$trt1_cr1_val, trt1_cr2_ID=list_f$trt1_cr2_ID,
          trt1_cr2_val=list_f$trt1_cr2_val,phenotypes=list_f$phenotypes,expressing1=list_f$expressing1)
        }
        if(n_offs_m>0&n_offs_f==0){
          pre_geno3<-list(sp=list_m$sp ,ID=list_m$ID,sex=list_m$sex,
          trt1_cr1_ID=list_m$trt1_cr1_ID, trt1_cr1_val=list_m$trt1_cr1_val, trt1_cr2_ID=list_m$trt1_cr2_ID,
          trt1_cr2_val=list_m$trt1_cr2_val,phenotypes=list_m$phenotypes,expressing1=list_m$expressing1)
        }

      }
   }

     k=k+2
     if(k==n_loops*2){
        k=n_loops*2-1
     }
   }###end of while


    ###all these below is a fast solution, to generalize it we must add the fact that they can be any number of pre_geno's
    ##assign phenotypes
    for(ss in 1:(n_offs_m+n_offs_f)){
      if(pre_geno1$sex[ss]=="male"){
        pre_geno1<-male_func(pre_geno1,ss,rho1,"neg",n_loci)
      }else{
        pre_geno1<-chrom_expres(pre_geno1,ss)
        pre_geno1<-female_func(pre_geno1,ss,rho1,"neg",n_loci)
      }
      if(pre_geno1$sex[ss]=="male"){
        pre_geno1<-male_func(pre_geno1,ss,rho1,"neg",n_loci)
      }else{
        pre_geno1<-chrom_expres(pre_geno1,ss)
        pre_geno1<-female_func(pre_geno1,ss,rho1,"neg",n_loci)
      }
    }


    for(ss in 1:(n_offs_m+n_offs_f)){
      if(pre_geno2$sex[ss]=="male"){
        pre_geno2<-male_func(pre_geno2,ss,rho2,"neg",n_loci)
      }else{
        pre_geno2<-chrom_expres(pre_geno2,ss)
        pre_geno2<-female_func(pre_geno2,ss,rho2,"neg",n_loci)
      }
      if(pre_geno2$sex[ss]=="male"){
        pre_geno2<-male_func(pre_geno2,ss,rho2,"neg",n_loci)
      }else{
        pre_geno2<-chrom_expres(pre_geno2,ss)
        pre_geno2<-female_func(pre_geno2,ss,rho2,"neg",n_loci)
      }
    }

    for(ss in 1:(n_offs_m+n_offs_f)){
      if(pre_geno3$sex[ss]=="male"){
        pre_geno3<-male_func_one_trait(pre_geno3,ss,"neg",n_loci)
      }else{
        pre_geno3<-chrom_expres_trait(pre_geno3,ss)
        pre_geno3<-female_func_one_trait(pre_geno3,ss,"neg",n_loci)
      }
    }


    pre_geno<-list(pre_geno1,pre_geno2,pre_geno3)

    names(pre_geno)[[1]]<-"both1"
    names(pre_geno)[[2]]<-"both2"
    names(pre_geno)[[3]]<-"both3"

return(pre_geno)
}

phenotypes<-function(new_bich,new_offs)
{
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


  ####interpolating pseudo-values (genetics) to ecological traits (phenotypes)
  for(i in 1:(length(new_bich$sp))){
  new_bich[i,"drift"]<- both3$phenotypes$phenotype1[i] ##assigns drift from single trait
  if (new_bich$sp[i]=="te" && sum( bichinhos_all[1:2,])>0 && sum(bichinhos_all[3:4,])>0){     # coinfection
    new_bich[i,"ass"]<- interpol2(new_offs$both1$phenotypes$phenotype1[i],Min_te, Max_te,ranges2[1,1],ranges2[1,2] )
    new_bich[i,"R_P"]<- interpol2(new_offs$both1$phenotypes$phenotype2[i],Min_te, Max_te,ranges3[3,1],ranges3[3,2])
     }
  if (new_bich$sp[i]=="tu" && sum(bichinhos_all[1:2,])>0 && sum(bichinhos_all[3:4,])>0){
    new_bich[i,"ass"]<- interpol2(new_offs$both1$phenotypes$phenotype1[i],Min_tu, Max_tu,ranges2[2,1],ranges2[2,2] )
    new_bich[i,"R_P"]<- interpol2(new_offs$both1$phenotypes$phenotype2[i],Min_tu, Max_tu,ranges3[3,1],ranges3[3,2] )
    }
  if (new_bich$sp[i]=="te" && sum(bichinhos_all[1:2,])>0 && sum(bichinhos_all[3:4,])==0){     #infection urticae
    new_bich[i,"ass"]<- interpol2(new_offs$both1$phenotype$phenotype1[i],Min_te, Max_te,ranges2[3,1],ranges2[3,2] )
    new_bich[i,"R_P"]<- interpol2(new_offs$both1$phenotypes$phenotype2[i],Min_te, Max_te,ranges3[4,1],ranges3[4,2] )
    }
   if (new_bich$sp[i]=="tu" && sum(bichinhos_all[1:2,])>0 && sum(bichinhos_all[3:4,])==0){
    new_bich[i,"ass"]<-interpol2(new_offs$both1$phenotypes$phenotype1[i],Min_tu, Max_tu,ranges2[4,1],ranges2[4,2] )
    new_bich[i,"R_P"]<- interpol2(new_offs$both1$phenotypes$phenotype2[i],Min_tu, Max_tu,ranges3[4,1],ranges3[4,2] )
    }
  if (new_bich$sp[i]=="te" && sum(bichinhos_all[1:2,])==0 && sum(bichinhos_all[3:4,])>0){        #infection evansi
    new_bich[i,"ass"]<- interpol2(new_offs$both1$phenotypes$phenotype1[i],Min_te, Max_te,ranges2[5,1],ranges2[5,2] )
    new_bich[i,"R_P"]<- interpol2(new_offs$both1$phenotypes$phenotype2[i],Min_te, Max_te,ranges3[1,1],ranges3[1,2] )
    }
   if( new_bich$sp[i]=="tu" && sum(bichinhos_all[1:2,])==0 && sum(bichinhos_all[3:4,])>0){
    new_bich[i,"ass"]<- interpol2(new_offs$both1$phenotypes$phenotype1[i],Min_tu, Max_tu,ranges2[6,1],ranges2[6,2] )
    new_bich[i,"R_P"]<- interpol2(new_offs$both1$phenotypes$phenotype2[i],Min_tu, Max_tu,ranges3[1,1],ranges3[1,2])
    }
   if( new_bich$sp[i]=="te" && sum(bichinhos_all[1:2,])==0 && sum(bichinhos_all[3:4,])==0){        #clean
    new_bich[i,"ass"]<- interpol2(new_offs$both1$phenotypes$phenotype1[i],Min_te, Max_te,ranges2["ass",1],ranges2["ass",2])
    new_bich[i,"R_P"]<- interpol2(new_offs$both1$phenotypes$phenotype2[i],Min_te, Max_te,ranges3[2,1],ranges3[2,2]  )
    }
   if (new_bich$sp[i]=="tu" && sum(bichinhos_all[1:2,])==0 && sum(bichinhos_all[3:4,])==0){
    new_bich[i,"ass"]<- interpol2(new_offs$both1$phenotypes$phenotype1[i],Min_tu, Max_tu,ranges2[8,1],ranges2[8,2] )
    new_bich[i,"R_P"]<- interpol2(new_offs$both1$phenotypes$phenotype2[i],Min_tu, Max_tu,ranges3[2,1],ranges3[2,2] )
    }
  }


  new_bich$m_damage<-ifelse(new_bich$sp=="tu",interpol2(new_offs$both2$phenotypes$phenotype1,Min_te, Max_te,
  ranges4[1],ranges4[2]),interpol2(new_offs$both2$phenotypes$phenotype1,Min_te, Max_te,
  ranges5[1],ranges5[2]))


  ###sex ratio with genetic basis
  new_bich$g_sex_ratio<-interpol2(new_offs$both2$phenotypes$phenotype2,Min_te, Max_te,ranges6[1],ranges6[2])
  
  ###establish EVOLUTIONARY LIMITS
  new_bich$ass<-ifelse(new_bich$ass<0.25,0.25,ifelse(new_bich$ass>0.95,0.95,new_bich$ass))
  new_bich$R_P<-ifelse(new_bich$R_P<100,100,ifelse(new_bich$R_P>300000,300000,new_bich$R_P))
  new_bich$m_damage<-ifelse(new_bich$m_damage<0.25,0.25,ifelse(new_bich$m_damage>0.8,0.8,new_bich$m_damage))
  new_bich$g_sex_ratio<-ifelse(new_bich$g_sex_ratio<0.5,0.5,ifelse(new_bich$g_sex_ratio>0.98,0.98,new_bich$g_sex_ratio))



  #initialitations of counters and states
  #############################################################
  #############################################################
  #############################################################
  #############################################################
  new_bich$day_food<-0
  new_bich$mass<-1.23+0.49 ##egg plus "mass gain" from wherever it comes?
  new_bich$quiesc<-7
  new_bich$adult_age=0
  new_bich$alive=1
  new_bich$sex_ratio<-new_bich$g_sex_ratio
 
  

  ##length_for_matings<-length(new_bich) ##to assess the position of the first mating counter (1 to 10)

  ###line for adding times of each mating for each female at adult age 1!!!
  new_bich$t_mate1=0
  new_bich$t_mate2=0
  new_bich$t_mate3=0

  for(i in 1:nrow(new_bich)){
    if(new_bich$sex[i]=="female"){
      new_bich[i,c("t_mate1","t_mate2","t_mate3")]<-sort(sample(1:12, 3, replace = FALSE))
    }
  }

return(new_bich)
}


reproduce<-function(sub_bich,j)
{
    
    n_offs_f=sub_bich$female_eggs[j]
    n_offs_m=sub_bich$male_eggs[j]

    ###GENOTYPES#########################################################################
    
    ##females
    if(n_offs_f>0){
      ##the list below will be filled with as many items as female offspring will 
      ##be born in this round
      female_offs_female_gamete_genotype<-list(NULL)
      
      for(k in 1:n_offs_f){
        female_offs_female_gamete_genotype[[k]]<-crlsm_to_chrom(sub_bich,j) ##transforms correlosome to chromosome and call crossover
      }
      
      female_offs_male_gamete_genotype<-list(NULL)
  
      for(k in 1:n_offs_f){
        female_offs_male_gamete_genotype[[k]]<-male_genetics(sub_bich,j)
        ##print(female_offs_male_gamete_genotype[[k]][[1]][,1]) ##indeed males do not recombine
      }
    }
    ##males
    if(n_offs_m>0){
       male_offs_male_gamete_genotype<-list(NULL)
  
      for(k in 1:n_offs_m){
        male_offs_male_gamete_genotype[[k]]<-crlsm_to_chrom(sub_bich,j)
      }
    }

    
    last_tu_ID<-suppressWarnings(max(as.integer(unlist(strsplit(bichinhos$In_ID[bichinhos$sp=="tu"],"_"))),na.rm=T))
    last_te_ID<-suppressWarnings(max(as.integer(unlist(strsplit(bichinhos$In_ID[bichinhos$sp=="te"],"_"))),na.rm=T))

    if(sub_bich$sp[j]=="tu"){
        last_ID=last_tu_ID
        }else{
        last_ID=last_te_ID
    }

    
    if(n_offs_m>0 & n_offs_f>0){
      new_offs<-make_both(sub_bich$sp[j],last_ID,n_offs_m,n_offs_f,male_offs_male_gamete_genotype,
      female_offs_male_gamete_genotype,female_offs_female_gamete_genotype)
    }
    if(n_offs_m>0 & n_offs_f==0){
      new_offs<-make_both(sub_bich$sp[j],last_ID,n_offs_m,n_offs_f,male_offs_male_gamete_genotype,
      NA,NA)
    }
    if(n_offs_m==0 & n_offs_f>0){
      new_offs<-make_both(sub_bich$sp[j],last_ID,n_offs_m,n_offs_f,NA,
      female_offs_male_gamete_genotype,female_offs_female_gamete_genotype)
    }
    
    ###how to append lists, rbind or c will work differently for different types of variables, so...
    what1<-suppressWarnings(Map(function(x,y) Map(rbind,x,y) , genotypes, new_offs)) ###
    what2<-Map(function(x,y) Map(c,x,y) , genotypes, new_offs) ###
    
    m1<-c(what2[[1]][1:3],what1[[1]][4:14])
    m2<-c(what2[[2]][1:3],what1[[2]][4:14])
    m3<-c(what2[[3]][1:3],what1[[3]][4:9])
    new_geno<-list(m1,m2,m3)

    names(new_geno)[[1]]<-"both1"
    names(new_geno)[[2]]<-"both2"
    names(new_geno)[[3]]<-"both3"

    #str(new_geno)
   ###END GENOTYPES

  
    ###PHENOTYPES############################################################################
    new_bich<-bichinhos[1:(n_offs_m+n_offs_f),]
    new_bich[,]<-NA
    new_bich$instar="egg"
    new_bich$sex<-new_offs$both1$sex
    new_bich$sp<-new_offs$both1$sp
    new_bich$In_ID<-new_offs$both1$ID
    new_bich$eggs<-200
    new_bich$drift<-new_offs$both3$phenotypes$phenotype1
    new_bich$position=sub_bich$position[j]
    new_bich$gen_dam=sub_bich$generation[j]
    new_bich$generation=sub_bich$generation[j]+1   ##adding one generation to that of the dam
    new_bich$death_date=0
    new_bich$migrate_date=0

    new_bich<-phenotypes(new_bich,new_offs)
 
    ##update egg load
    sub_bich[j,"eggs"]=sub_bich[j,"eggs"]-sub_bich[j,"eggs_day"]

    
return(list(sub_bich[j,],new_geno,new_bich))##
}
  

shorten<-function(genotypes,bichinhos,p,n_to_remove,IDs_to_remove)
{
  how_far=5 ##number of plants to look ahead for male genotype or male damage to ensure not removing those males!!

  in_plant<-bichinhos[bichinhos$position %in% (plant$plant_ID[p]-1),]
  males_in_plant<-in_plant[in_plant$sex %in% "male"&in_plant$instar %in% "adult",]
  
  for(xxxx in (plant$plant_ID[p]:(plant$plant_ID[p]+how_far))){
      if(xxxx== plant$plant_ID[p]){
         females_in_plants2<-bichinhos[bichinhos$sex=="female" & bichinhos$instar=="adult"
         & bichinhos$position %in% xxxx,(length_for_update+1):(length_for_update+how_far)]
         
      }else{
         pre<-bichinhos[bichinhos$sex=="female" & bichinhos$instar=="adult"
         & bichinhos$position %in% xxxx,(length_for_update+1):(length_for_update+how_far)]
         
         females_in_plants2<-rbind(females_in_plants2,pre)
      }
  }
  
  fem_vector<-na.exclude(as.vector(as.matrix(females_in_plants2)))
  
#  for(xxxx in (plant$plant_ID[p]:(plant$plant_ID[p]+how_far))){
#      if(xxxx== plant$plant_ID[p]){
#         a<-females_in_plants2[females_in_plants2$position %in% xxxx,(length(females_in_plants2)-how_far):length(females_in_plants2)]     ## ]
#         a<-na.exclude(a)
#      }else{
#         pre<-females_in_plants2[females_in_plants2$position %in% xxxx,(length(females_in_plants2)-how_far):length(females_in_plants2)]
#         pre<-na.exclude(pre)
#         a<-c(a,pre)
#      }
#  } 
# 
  
#  females_in_plants<-bichinhos[bichinhos$sex=="female" & bichinhos$instar=="adult"
#  & bichinhos$position %in% (plant$plant_ID[p]:(plant$plant_ID[p]+how_far)),]
#  
#  a<-females_in_plants[females_in_plants$position %in% (plant$plant_ID[p]:(plant$plant_ID[p]+how_far)),             ##(plant$plant_ID[p]-1):(plant$plant_ID[p])
#  length_for_update+1]     ## :(length_for_update+5)
#  b<-females_in_plants[females_in_plants$position %in% (plant$plant_ID[p]:(plant$plant_ID[p]+how_far)),             ##(plant$plant_ID[p]-1):(plant$plant_ID[p])
#  length_for_update+2]     ## :(length_for_update+5)
#  c<-females_in_plants[females_in_plants$position %in% (plant$plant_ID[p]:(plant$plant_ID[p]+how_far)),             ##(plant$plant_ID[p]-1):(plant$plant_ID[p])
#  length_for_update+3]     ## :(length_for_update+5)
#  d<-females_in_plants[females_in_plants$position %in% (plant$plant_ID[p]:(plant$plant_ID[p]+how_far)),             ##(plant$plant_ID[p]-1):(plant$plant_ID[p])
#  length_for_update+4]     ## :(length_for_update+5)
#  e<-females_in_plants[females_in_plants$position %in% (plant$plant_ID[p]:(plant$plant_ID[p]+how_far)),             ##(plant$plant_ID[p]-1):(plant$plant_ID[p])
#  length_for_update+5]     ## :(length_for_update+5)
#  
#  a<-na.exclude(a)
#  b<-na.exclude(b)
#  c<-na.exclude(c)
#  d<-na.exclude(d)
#  e<-na.exclude(e)
#  
  a<-fem_vector
  a2<-a[a %in% males_in_plant$In_ID] ###a %in% males_in_plant$In_ID
#  b2<-b[b %in% males_in_plant$In_ID]
#  c2<-c[c %in% males_in_plant$In_ID]
#  d2<-d[d %in% males_in_plant$In_ID]
#  e2<-e[e %in% males_in_plant$In_ID]
# 
#  males_pre_keep<-unique(c(a2,b2,c2,d2,e2))
  
  males_to_keep<-unique(a2)
  
  ##those males that are on target plant
  ##males_to_keep<-males_in_plant[males_in_plant$In_ID %in% males_pre_keep,"In_ID"]

###genotypes to remove
  both1 <- list(sp=genotypes$both1$sp[genotypes$both1$ID %in% IDs_to_remove] ,
  ID=genotypes$both1$ID[genotypes$both1$ID %in% IDs_to_remove ],
  sex=genotypes$both1$sex[genotypes$both1$ID %in% IDs_to_remove ],
  trt1_cr1_ID=genotypes$both1$trt1_cr1_ID[genotypes$both1$ID %in% IDs_to_remove, ], 
  trt1_cr1_val=genotypes$both1$trt1_cr1_val[genotypes$both1$ID %in% IDs_to_remove, ],
  trt1_cr2_ID=genotypes$both1$trt1_cr2_ID[genotypes$both1$ID %in% IDs_to_remove, ],
  trt1_cr2_val=genotypes$both1$trt1_cr2_val[genotypes$both1$ID %in% IDs_to_remove, ],
  trt2_cr1_ID=genotypes$both1$trt2_cr1_ID[genotypes$both1$ID %in% IDs_to_remove, ],
  trt2_cr1_val=genotypes$both1$trt2_cr1_val[genotypes$both1$ID %in% IDs_to_remove, ],
  trt2_cr2_ID=genotypes$both1$trt2_cr2_ID[genotypes$both1$ID %in% IDs_to_remove, ],
  trt2_cr2_val= genotypes$both1$trt2_cr2_val[genotypes$both1$ID %in% IDs_to_remove, ],
  phenotypes=data.frame(phenotype1=genotypes$both1$phenotypes$phenotype1[genotypes$both1$ID %in% IDs_to_remove ],
  phenotype2=genotypes$both1$phenotypes$phenotype2[genotypes$both1$ID %in% IDs_to_remove ]),
  expressing1=genotypes$both1$expressing1[genotypes$both1$ID %in% IDs_to_remove, ],
  expressing2=genotypes$both1$expressing2[genotypes$both1$ID %in% IDs_to_remove, ])

  both1 <- list(sp=both1$sp[!both1$ID %in% males_to_keep ] ,
  ID=both1$ID[!both1$ID %in% males_to_keep ],
  sex=both1$sex[!both1$ID %in% males_to_keep ],
  trt1_cr1_ID=both1$trt1_cr1_ID[!both1$ID %in% males_to_keep, ], 
  trt1_cr1_val=both1$trt1_cr1_val[!both1$ID %in% males_to_keep, ],
  trt1_cr2_ID=both1$trt1_cr2_ID[!both1$ID %in% males_to_keep, ],
  trt1_cr2_val=both1$trt1_cr2_val[!both1$ID %in% males_to_keep, ],
  trt2_cr1_ID=both1$trt2_cr1_ID[!both1$ID %in% males_to_keep, ],
  trt2_cr1_val=both1$trt2_cr1_val[!both1$ID %in% males_to_keep, ],
  trt2_cr2_ID=both1$trt2_cr2_ID[!both1$ID %in% males_to_keep, ],
  trt2_cr2_val= both1$trt2_cr2_val[!both1$ID %in% males_to_keep, ],
  phenotypes=data.frame(phenotype1=both1$phenotypes$phenotype1[!both1$ID %in% males_to_keep ],
  phenotype2=both1$phenotypes$phenotype2[!both1$ID %in% males_to_keep ]),
  expressing1=both1$expressing1[!both1$ID %in% males_to_keep, ],
  expressing2=both1$expressing2[!both1$ID %in% males_to_keep, ])

  both2 <- list(sp=genotypes$both2$sp[genotypes$both2$ID %in% IDs_to_remove] ,
  ID=genotypes$both2$ID[genotypes$both2$ID %in% IDs_to_remove ],
  sex=genotypes$both2$sex[genotypes$both2$ID %in% IDs_to_remove ],
  trt1_cr1_ID=genotypes$both2$trt1_cr1_ID[genotypes$both2$ID %in% IDs_to_remove, ], 
  trt1_cr1_val=genotypes$both2$trt1_cr1_val[genotypes$both2$ID %in% IDs_to_remove, ],
  trt1_cr2_ID=genotypes$both2$trt1_cr2_ID[genotypes$both2$ID %in% IDs_to_remove, ],
  trt1_cr2_val=genotypes$both2$trt1_cr2_val[genotypes$both2$ID %in% IDs_to_remove, ],
  trt2_cr1_ID=genotypes$both2$trt2_cr1_ID[genotypes$both2$ID %in% IDs_to_remove, ],
  trt2_cr1_val=genotypes$both2$trt2_cr1_val[genotypes$both2$ID %in% IDs_to_remove, ],
  trt2_cr2_ID=genotypes$both2$trt2_cr2_ID[genotypes$both2$ID %in% IDs_to_remove, ],
  trt2_cr2_val= genotypes$both2$trt2_cr2_val[genotypes$both2$ID %in% IDs_to_remove, ],
  phenotypes=data.frame(phenotype1=genotypes$both2$phenotypes$phenotype1[genotypes$both2$ID %in% IDs_to_remove ],
  phenotype2=genotypes$both2$phenotypes$phenotype2[genotypes$both2$ID %in% IDs_to_remove ]),
  expressing1=genotypes$both2$expressing1[genotypes$both2$ID %in% IDs_to_remove, ],
  expressing2=genotypes$both2$expressing2[genotypes$both2$ID %in% IDs_to_remove, ])

  both2 <- list(sp=both2$sp[!both2$ID %in% males_to_keep ] ,
  ID=both2$ID[!both2$ID %in% males_to_keep ],
  sex=both2$sex[!both2$ID %in% males_to_keep ],
  trt1_cr1_ID=both2$trt1_cr1_ID[!both2$ID %in% males_to_keep, ], 
  trt1_cr1_val=both2$trt1_cr1_val[!both2$ID %in% males_to_keep, ],
  trt1_cr2_ID=both2$trt1_cr2_ID[!both2$ID %in% males_to_keep, ],
  trt1_cr2_val=both2$trt1_cr2_val[!both2$ID %in% males_to_keep, ],
  trt2_cr1_ID=both2$trt2_cr1_ID[!both2$ID %in% males_to_keep, ],
  trt2_cr1_val=both2$trt2_cr1_val[!both2$ID %in% males_to_keep, ],
  trt2_cr2_ID=both2$trt2_cr2_ID[!both2$ID %in% males_to_keep, ],
  trt2_cr2_val= both2$trt2_cr2_val[!both2$ID %in% males_to_keep, ],
  phenotypes=data.frame(phenotype1=both2$phenotypes$phenotype1[!both2$ID %in% males_to_keep ],
  phenotype2=both2$phenotypes$phenotype2[!both2$ID %in% males_to_keep ]),
  expressing1=both2$expressing1[!both2$ID %in% males_to_keep, ],
  expressing2=both2$expressing2[!both2$ID %in% males_to_keep, ])

  both3 <- list(sp=genotypes$both3$sp[genotypes$both3$ID %in% IDs_to_remove ] ,
  ID=genotypes$both3$ID[genotypes$both3$ID %in% IDs_to_remove ],
  sex=genotypes$both3$sex[genotypes$both3$ID %in% IDs_to_remove ],
  trt1_cr1_ID=genotypes$both3$trt1_cr1_ID[genotypes$both3$ID %in% IDs_to_remove, ], 
  trt1_cr1_val=genotypes$both3$trt1_cr1_val[genotypes$both3$ID %in% IDs_to_remove, ],
  trt1_cr2_ID=genotypes$both3$trt1_cr2_ID[genotypes$both3$ID %in% IDs_to_remove, ],
  trt1_cr2_val=genotypes$both3$trt1_cr2_val[genotypes$both3$ID %in% IDs_to_remove, ],
  phenotypes=data.frame(phenotype1=genotypes$both3$phenotypes$phenotype1[genotypes$both3$ID %in% IDs_to_remove ]),
  expressing1=genotypes$both3$expressing1[genotypes$both3$ID %in% IDs_to_remove, ])
  
  both3 <- list(sp=both3$sp[!both3$ID %in% males_to_keep] ,
  ID=both3$ID[!both3$ID %in% males_to_keep],sex=both3$sex[!both3$ID %in% males_to_keep],
  trt1_cr1_ID=both3$trt1_cr1_ID[!both3$ID %in% males_to_keep,], 
  trt1_cr1_val=both3$trt1_cr1_val[!both3$ID %in% males_to_keep,],
  trt1_cr2_ID=both3$trt1_cr2_ID[!both3$ID %in% males_to_keep,],
  trt1_cr2_val=both3$trt1_cr2_val[!both3$ID %in% males_to_keep,],
  phenotypes=data.frame(phenotype1=both3$phenotypes$phenotype1[!both3$ID %in% males_to_keep]),
  expressing1=both3$expressing1[!both3$ID %in% males_to_keep,])


  rm_genotypes<-list(both1=both1,both2=both2,both3=both3)

#  sort(rm_genotypes$both1$ID)
#  
#  sort(males_to_keep)
#
IDs_removed<-both1$ID

###genotypes to keep
  both1 <- list(sp=genotypes$both1$sp[!genotypes$both1$ID %in% IDs_removed] ,
  ID=genotypes$both1$ID[!genotypes$both1$ID %in% IDs_removed],
  sex=genotypes$both1$sex[!genotypes$both1$ID %in% IDs_removed],
  trt1_cr1_ID=genotypes$both1$trt1_cr1_ID[!genotypes$both1$ID %in% IDs_removed,], 
  trt1_cr1_val=genotypes$both1$trt1_cr1_val[!genotypes$both1$ID %in% IDs_removed,],
  trt1_cr2_ID=genotypes$both1$trt1_cr2_ID[!genotypes$both1$ID %in% IDs_removed,],
  trt1_cr2_val=genotypes$both1$trt1_cr2_val[!genotypes$both1$ID %in% IDs_removed,],
  trt2_cr1_ID=genotypes$both1$trt2_cr1_ID[!genotypes$both1$ID %in% IDs_removed,],
  trt2_cr1_val=genotypes$both1$trt2_cr1_val[!genotypes$both1$ID %in% IDs_removed,],
  trt2_cr2_ID=genotypes$both1$trt2_cr2_ID[!genotypes$both1$ID %in% IDs_removed,],
  trt2_cr2_val= genotypes$both1$trt2_cr2_val[!genotypes$both1$ID %in% IDs_removed,],
  phenotypes=data.frame(phenotype1=genotypes$both1$phenotypes$phenotype1[!genotypes$both1$ID %in% IDs_removed],
  phenotype2=genotypes$both1$phenotypes$phenotype2[!genotypes$both1$ID %in% IDs_removed]),
  expressing1=genotypes$both1$expressing1[!genotypes$both1$ID %in% IDs_removed,],
  expressing2=genotypes$both1$expressing2[!genotypes$both1$ID %in% IDs_removed,])

  both2 <- list(sp=genotypes$both2$sp[!genotypes$both2$ID %in% IDs_removed] ,
  ID=genotypes$both2$ID[!genotypes$both2$ID %in% IDs_removed],
  sex=genotypes$both2$sex[!genotypes$both2$ID %in% IDs_removed],
  trt1_cr1_ID=genotypes$both2$trt1_cr1_ID[!genotypes$both2$ID %in% IDs_removed,], 
  trt1_cr1_val=genotypes$both2$trt1_cr1_val[!genotypes$both2$ID %in% IDs_removed,],
  trt1_cr2_ID=genotypes$both2$trt1_cr2_ID[!genotypes$both2$ID %in% IDs_removed,],
  trt1_cr2_val=genotypes$both2$trt1_cr2_val[!genotypes$both2$ID %in% IDs_removed,],
  trt2_cr1_ID=genotypes$both2$trt2_cr1_ID[!genotypes$both2$ID %in% IDs_removed,],
  trt2_cr1_val=genotypes$both2$trt2_cr1_val[!genotypes$both2$ID %in% IDs_removed,],
  trt2_cr2_ID=genotypes$both2$trt2_cr2_ID[!genotypes$both2$ID %in% IDs_removed,],
  trt2_cr2_val= genotypes$both2$trt2_cr2_val[!genotypes$both2$ID %in% IDs_removed,],
  phenotypes=data.frame(phenotype1=genotypes$both2$phenotypes$phenotype1[!genotypes$both2$ID %in% IDs_removed],
  phenotype2=genotypes$both2$phenotypes$phenotype2[!genotypes$both2$ID %in% IDs_removed]),
  expressing1=genotypes$both2$expressing1[!genotypes$both2$ID %in% IDs_removed,],
  expressing2=genotypes$both2$expressing2[!genotypes$both2$ID %in% IDs_removed,])

  both3 <- list(sp=genotypes$both3$sp[!genotypes$both3$ID %in% IDs_removed] ,
  ID=genotypes$both3$ID[!genotypes$both3$ID %in% IDs_removed],
  sex=genotypes$both3$sex[!genotypes$both3$ID %in% IDs_removed],
  trt1_cr1_ID=genotypes$both3$trt1_cr1_ID[!genotypes$both3$ID %in% IDs_removed,], 
  trt1_cr1_val=genotypes$both3$trt1_cr1_val[!genotypes$both3$ID %in% IDs_removed,],
  trt1_cr2_ID=genotypes$both3$trt1_cr2_ID[!genotypes$both3$ID %in% IDs_removed,],
  trt1_cr2_val=genotypes$both3$trt1_cr2_val[!genotypes$both3$ID %in% IDs_removed,],
  phenotypes=data.frame(phenotype1=genotypes$both3$phenotypes$phenotype1[!genotypes$both3$ID %in% IDs_removed]),
  expressing1=genotypes$both3$expressing1[!genotypes$both3$ID %in% IDs_removed,])
  
  kept_genotypes<-list(both1=both1,both2=both2,both3=both3)
  
#  sort(kept_genotypes$both1$ID)
#  
#  sort(males_to_keep)
#



  ###remove phenotypes
  rm_bichinhos<-bichinhos[bichinhos$In_ID %in% rm_genotypes$both1$ID,]
  ###kept phenotypes
  kept_bichinhos<-bichinhos[!bichinhos$In_ID %in% rm_genotypes$both1$ID,]
  rownames(kept_bichinhos)<-as.character(1:nrow(kept_bichinhos))

  saveRDS(rm_genotypes,paste("rm_genotypes_",plant$plant_ID[p]-1,".rds",sep=""))
  write.table(rm_bichinhos,paste("rm_bichinhos_",plant$plant_ID[p]-1,".txt",sep=""))

return(list(kept_genotypes,kept_bichinhos))
}




