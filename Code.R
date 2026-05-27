
#### S0: Auxiliary functions ####

lt<-function(x){sign(x)*log1p(abs(x))} #log(abs(x)+1)
et<-function(x){sign(x)*expm1(abs(x))} #exp(abs(x))-1
nmr<-function(x){as.numeric(x)}
chr<-function(x){as.character(x)}
fct<-function(x){factor(x,levels=sort(unique(x)))}
mdl<-function(x){
  ux<-unique(x)
  tab<-tabulate(match(x, ux))
  ux[tab==max(tab)]
}
plt<-function(x,y,...){
  plot(x,y,pch=16,col=rgb(0,0,0,0.01),...)
}
GetPts<-function(x){vect(x,geom=c('x','y'),crs='epsg:4326',keepgeom=T)}
ws<-function(x,p){
  av<-quantile(x,p);x[x<av]<-av
  av<-quantile(x,1-p);x[x>av]<-av
  x
}

PrepareData<-function(aqrs='20km',datars='20km',type,scn=NULL,scope){
  
  if (scope=='China'){
    dir<-"F:/China_data/Code files/Datasets"
    dmm_pop<-readRDS("F:/China_data/Code files/Boundaries/PopPtsDist_10km.rds")
    aqrsbffr<-'20km'
  }
  if (scope=='Global'){
    dir<-"F:/Global_data/Code files/Datasets"
    dmm_pop<-NULL #Too big
    aqrsbffr<-'50km'
  }
  
  library(terra)
  library(dplyr)
  library(tidyr)
  
  ## Get data
  #Filter by resolution and name
  ad<-list.files(dir,datars,full.names=T);ad<-ad[!grepl('OffshrAq',ad,fixed=T)];ad1<-ad #Non-aquaculture data
  ad<-list.files(dir,aqrs,full.names=T);ad<-ad[grepl('OffshrAq_',ad,fixed=T)];ad2<-ad #Aquaculture data (except buffer) (using generic name, actual resolutino is pixel)
  ad<-list.files(dir,aqrsbffr,full.names=T);ad<-ad[grepl('OffshrAqBuffer',ad,fixed=T)];ad3<-ad #Aquaculture data (only buffer)
  ad<-sort(c(ad1,ad2,ad3))
  #Filter by category
  # ad1<-ad[!grepl('proj',ad,perl=T) & grepl('[0-9]{2}km',ad,perl=T)] #For historical
  ad2<-ad[grepl('_proj_bs_',ad,perl=T)] #For baseline
  if (!is.null(scn)){ad3<-ad[grepl(paste0('_proj_',scn,'_'),ad,perl=T)]} #For projections
  
  ## Prepare data from files
  {
    ## Generate raw data frame
    if (type=='base'){ad<-ad2}
    if (type=='proj'){ad<-ad3}
    vn<-vector('character',length(ad));for (i in 1:length(vn)){vn[i]<-strsplit(strsplit(ad[i],'\\/',perl=T)[[1]][5],'_..km',perl=T)[[1]][1]} #Variable names
    if (type!='proj'){
      av<-readRDS(ad[which(grepl('Bth',vn,fixed=T))])
      nobs<-nrow(readRDS(ad[which(grepl('Bth',vn,fixed=T))]))
    } ##Data frame template for PtID, x, y
    if (type=='proj'){
      av<-readRDS(ad[which(grepl('Population',vn,fixed=T))])
      av<-av[!duplicated(av$PtID),]
      nobs<-nrow(av)
    } ##Data frame template for PtID, x, y
    if (type=='base' & scope=='China'){tt<-seq(2005,2020,5)}
    if (type=='base' & scope=='Global'){tt<-seq(2015,2020,5)}
    if (type=='proj'){tt<-seq(2025,2100,5)}
    am0<-expand.grid(PtID=sort(unique(av$PtID)),x=-999999,y=-999999,Period=tt,Var=chr(vn),Val=nmr(NA),ValChr=chr(NA));am0$ValChr<-chr(am0$ValChr)
    am0<-data.frame(rows_update(am0,av[,c('PtID','x','y')],'PtID'))
    for (i in 1:length(vn)){
      
      am<-readRDS(ad[i])
      #print(dim(am))
      
      #print(unique(am$Period))
      am$Period[is.element(am$Period,2003:2005)]<-2005
      am$Period[is.element(am$Period,2015:2017)]<-2015
      am$Period[is.element(am$Period,2019:2023)]<-2020
      if (any(!is.element(am$Period,seq(2005,2100,5)))){stop(paste0('Year format error ',vn[i]))}
      #am$Var<-gsub('GDPP','GDPC',am$Var,fixed=T) #####################
      #am$Var<-gsub('VHM0|sws','vhm0',am$Var,perl=T) #####################
      #am$Var<-gsub('_projbs_','_proj_bs',am$Var,fixed=T) #####################
      if (is.element(vn[i],c('KD490','HDI','Price_AN','Price_PL'))){next()} #####################
      
      print(ad[i])
      # print(dim(am))
      
      #nobs based on single year
      if (nrow(am)!=nobs & class(am$Val)=='numeric'){am0<-data.frame(rows_update(am0,am[,c('Period','PtID','Var','Val')],c('Period','PtID','Var')));next()}
      if (class(am$Val)=='character'){names(am)[6]<-'ValChr';am0<-data.frame(rows_update(am0,am[,c('Period','PtID','Var','ValChr')],c('Period','PtID','Var')));next()}
      if (nrow(am)==nobs){am0<-data.frame(rows_update(am0,am[,c('PtID','Var','Val')],c('PtID','Var')));next()} #For variables that are constant across time
      #View(am0[am0$Var==vn[i],])
      
    }
    #Expand to columns
    am1<-data.frame(pivot_wider(am0[,-c(7)],id_cols=names(am0)[1:4],names_from='Var',values_from='Val'))
    am2<-data.frame(pivot_wider(am0[,-c(6)],id_cols=names(am0)[1:4],names_from='Var',values_from='ValChr'))
    am0<-am1
    #Add missing columns for projections
    if (type=='base'){am0$Subregion_proj_bs<-am2$Subregion_proj_bs}
    if (type=='proj'){
      
      for (vr in c('Subregion')){
        am<-readRDS(paste0(dir,"/",vr,"_proj_bs_20km.rds"));names(am)[6]<-vr
        am<-am[am$Period>=2019,]
        am0[,vr]<-chr(NA);am0[,vr]<-chr(am0[,vr])
        am0<-data.frame(rows_update(am0,am[,c('PtID',vr)],c('PtID')))
      }
      
      for (vr in c('Bth','ShrLength','ProtAreas')){
        am<-readRDS(paste0(dir,"/",vr,"_proj_bs_20km.rds"));names(am)[6]<-vr
        am<-am[am$Period>=2019,]
        am0[,vr]<-nmr(NA)
        am0<-data.frame(rows_update(am0,am[,c('PtID',vr)],c('PtID')))
      }
      
      for (vr in c('OffshrAq_AN','OffshrAq_PL','OffshrAqBuffer_AN','OffshrAqBuffer_PL')){
        am0[,vr]<-nmr(NA)
      }
      
    }
    #Other adjustments to format
    am0$Period<-fct(am0$Period)
    am0$Subregion<-fct(am0$Subregion)
    if (type=='base'){names(am0)<-gsub('_proj_bs','',names(am0),fixed=T)}
    if (type=='proj'){names(am0)<-gsub(paste0('_proj_',scn),'',names(am0),perl=T)}
    am0<-vect(am0,c("x","y"),crs='epsg:4326')
    am0$x<-geom(am0)[,3]
    am0$y<-geom(am0)[,4]
    if (scope=='China'){
      am0<-am0[,c('PtID','x','y',"Period",
                  "Population",
                  "GDPC",
                  "Fsh","PriSect","SecSect","TerSect",
                  "LUC_Crp","LUC_Frs","LUC_Urb","LUC_Wtr",
                  "thetao","so","chl",'ph','sws',
                  "Bth","ShrLength",
                  "ProtAreas","Subregion",
                  "OffshrAqBuffer_AN","OffshrAqBuffer_PL",
                  "OffshrAq_AN","OffshrAq_PL"
      )]
    }
    if (scope=='Global'){
      am0<-am0[,c('PtID','x','y',"Period",
                  "Population",
                  "GDPC",
                  "Fsh",
                  "TradeNonNorm_Consumption","TradeNonNorm_Exports","TradeNonNorm_Imports","TradeNonNorm_Production",
                  "TradeNorm_Consumption","TradeNorm_Exports","TradeNorm_Imports","TradeNorm_Production",
                  "LUC_Crp","LUC_Frs","LUC_Urb","LUC_Wtr",
                  "thetao","so","chl",'ph','sws',
                  "Bth","ShrLength",
                  "ProtAreas","Subregion",
                  "OffshrAqBuffer_AN","OffshrAqBuffer_PL",
                  "OffshrAq_AN","OffshrAq_PL"
      )]
    }
    str(data.frame(am0))
    summary(am0)
    #!is.unsorted(am0$Period)
    #for (k in unique(am0$Period)){print(!is.unsorted(am0$PtID[am0$Period==k]))}
    #all(spt$PtID==sort(spt$PtID))
  }
  
  ## Impute and amend data
  {
    #plot(am0[!is.na(am0$TradeNonNorm_Consumption),]);plot(am0[is.na(am0$TradeNonNorm_Consumption),],add=T,col='red')
    SptImputation<-function(am0,prp,type,fx=mean){
      
      if (type=='base'){
        tt<-c(2020)
        ddm<-ddm_pop
        ddm<-as.matrix(ddm);colnames(ddm)<-1:dim(ddm)[2];row.names(ddm)<-1:dim(ddm)[1]
      }
      if (type=='proj'){
        tt<-c(2025)
        ddm<-ddm_pop
        ddm<-as.matrix(ddm);colnames(ddm)<-1:dim(ddm)[2];row.names(ddm)<-1:dim(ddm)[1]
      }
      
      for (yr in tt){
        
        amm<-data.frame(am0[am0$Period==yr,])
        av<-!is.finite(amm[,prp])
        if (!any(av)){next()}
        am1<-amm[av,] #Missing data
        am2<-amm[!av,] #Available data
        
        ddd<-ddm[am1$PtID,am2$PtID]
        for (i in 1:nrow(am1)){
          u<-which(nmr(row.names(ddd))==am1$PtID[i])
          v<-which(ddd[u,]==min(ddd[u,]))
          av<-am2[is.element(am2$PtID,nmr(colnames(ddd)[v])),prp]
          # if (any(is.na(av))){stop()}
          # if (!is.finite(mean(av[!is.na(av)]))){stop()}
          amm[amm$PtID==am1$PtID[i],prp]<-fx(av[!is.na(av)])
        }
        am0[am0$Period==yr,prp]<-amm[,prp]
      }
      
      am0
    } #Use data from closest points
    SptImputationAppr<-function(am0,prp,tt,ctff=NULL){ #cutoff in degrees, ~200 km
      
      for (yr in tt){
        
        amm<-data.frame(am0[am0$Period==yr,])
        av<-!is.finite(amm[,prp])
        if (!any(av)){next()}
        am1<-amm[av,] #Missing data
        am2<-amm[!av,] #Available data
        
        if (is.null(ctff)){
          for (i in 1:nrow(am1)){
            u<-abs(am2[,'x']-am1[i,'x'])+abs(am2[,'y']-am1[i,'y'])
            amm[amm$PtID==am1$PtID[i],prp]<-am2[which.min(u),prp] #Unconditional
          }
        }
        if (!is.null(ctff)){
          for (i in 1:nrow(am1)){
            u<-abs(am2[,'x']-am1[i,'x'])+abs(am2[,'y']-am1[i,'y'])
            amm[amm$PtID==am1$PtID[i],prp]<-am2[which.min(u),prp] #Unconditional
            if (which.min(u)>ctff){amm[amm$PtID==am1$PtID[i],prp]<-0} #Add distance cutoff
          }
        }
        
        am0[am0$Period==yr,prp]<-amm[,prp]
        
      }
      
      am0
    }# Approximation for many points, use xy distance
    # Functions assume PtID is sorted
    InterpolateDataCont1<-function(am0,cl){
      
      yr<-chr(seq(2035,2095,10))
      am0<-data.frame(am0)
      
      for (i in yr){
        am0[am0$Period==i,cl]<-(am0[am0$Period==chr(nmr(i)+5),cl]+am0[am0$Period==chr(nmr(i)-5),cl])/2
      }
      
      am0<-vect(am0,geom=c('x','y'),crs='epsg:4326',keepgeom=T)
      am0
    }
    InterpolateDataCont2<-function(am0,cl){
      
      am0<-data.frame(am0)
      
      for (i in chr(seq(2035,2045,5))){
        a<-am0[am0$Period=='2030',cl]
        b<-am0[am0$Period=='2050',cl]
        am0[am0$Period==i,cl]<-a+(nmr(i)-nmr('2030'))*(b-a)/20
      }
      for (i in chr(seq(2055,2065,5))){
        a<-am0[am0$Period=='2050',cl]
        b<-am0[am0$Period=='2070',cl]
        am0[am0$Period==i,cl]<-a+(nmr(i)-nmr('2050'))*(b-a)/20
      }
      for (i in chr(seq(2075,2095,5))){
        a<-am0[am0$Period=='2070',cl]
        b<-am0[am0$Period=='2100',cl]
        am0[am0$Period==i,cl]<-a+(nmr(i)-nmr('2070'))*(b-a)/30
      }
      
      am0<-vect(am0,geom=c('x','y'),crs='epsg:4326',keepgeom=T)
      am0
    }
    ExtrapolateDataCont1<-function(am0,cl,scope){
      
      am0<-data.frame(am0)
      
      yr<-'2025'
      for (k in cl){
        # vct1<-readRDS(paste0("F:/China_data/Code files/Datasets/",k,"_proj_bs_20km.rds"))
        # ###Imputation
        # nn<-which(is.na(vct1$Val));if (length(nn)>0){for (l in nn){
        #   av<-abs(vct1$x-vct1$x[l])<=0.5 & abs(vct1$y-vct1$y[l])<=0.5 #Points within given distance
        #   vct1$Val[l]<-mean(vct1$Val[av],na.rm=T)
        # }}
        # ##
        # vct1<-vct1[,'Val']
        
        vct1<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_','base','_','0','.rds'));vct1<-vct1[vct1$Period=='2020',k]
        vct2<-am0[am0$Period=='2030',k]
        lc<-cbind(vct1,vct2)
        am0[am0$Period==yr,k]<-rowMeans(lc)
      }
      
      #Not needed for now
      # yr<-c('2095','2100')
      # for (i in yr){
      #   am0[am0$Period==i,cl]<-am0[am0$Period==chr(nmr(i)-5),cl]+
      #     (am0[am0$Period==chr(nmr(i)-5),cl]-am0[am0$Period==chr(nmr(i)-10),cl])
      # }
      
      am0<-vect(am0,geom=c('x','y'),crs='epsg:4326',keepgeom=T)
      am0
    }
    ExtrapolateDataCont2<-function(am0,cl,scope){
      
      am0<-data.frame(am0)
      
      yr<-'2025'
      for (k in cl){
        # vct1<-readRDS(paste0("F:/China_data/Code files/Datasets/",k,"_proj_bs_20km.rds"))
        # ###Imputation
        # nn<-which(is.na(vct1$Val));if (length(nn)>0){for (l in nn){
        #   av<-abs(vct1$x-vct1$x[l])<=0.25 & abs(vct1$y-vct1$y[l])<=0.25 #Points within given distance
        #   vct1$Val[l]<-mean(vct1$Val[av],na.rm=T)
        # }}
        # ##
        # vct1<-vct1[,'Val']
        
        vct1<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_','base','_','0','.rds'));vct1<-vct1[vct1$Period=='2020',k]
        vct2<-am0[am0$Period=='2030',k]
        lc<-cbind(vct1,vct2)
        am0[am0$Period==yr,k]<-rowMeans(lc)
      }
      
      am0<-vect(am0,geom=c('x','y'),crs='epsg:4326',keepgeom=T)
      am0
    }
    #WQ and LUC
    if (type!='proj'){
      
      if (scope=='China'){at<-seq(2005,2020,5)}
      if (scope=='Global'){at<-c(2015,2020)}
      
      #WQ and terrain
      am0<-SptImputationAppr(am0,'thetao',at)
      am0<-SptImputationAppr(am0,'so',at)
      am0<-SptImputationAppr(am0,'ph',at)
      am0<-SptImputationAppr(am0,'chl',at)
      am0<-SptImputationAppr(am0,'sws',at)
      
      if (scope=='Global'){
        am0<-SptImputationAppr(am0,'TradeNonNorm_Consumption',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNonNorm_Exports',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNonNorm_Imports',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNonNorm_Production',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNorm_Consumption',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNorm_Exports',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNorm_Imports',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNorm_Production',at,ctff=2)
      }
      
    }
    if (type=='proj'){
      # Unavailable data
      # WQ seq(2030,2090,10)
      # LUC c(2030,2050,2070,2100)
      
      #Fill data not available
      am0<-InterpolateDataCont1(am0,c('thetao','so','chl','ph','sws'))
      am0<-InterpolateDataCont2(am0,paste0('LUC_',c("Crp","Frs","Urb","Wtr")))
      am0<-ExtrapolateDataCont1(am0,c('thetao','so','chl','ph','sws'),scope)
      am0<-ExtrapolateDataCont2(am0,paste0('LUC_',c("Crp","Frs","Urb","Wtr")),scope)
      #length(which(is.na(xx[,'thetao'])))
      #length(which(is.na(xx[,'LUC_Crp'])))
      
      at<-seq(2025,2050,5)
      
      am0<-SptImputationAppr(am0,'chl',at)
      am0<-SptImputationAppr(am0,'thetao',at)
      am0<-SptImputationAppr(am0,'so',at)
      am0<-SptImputationAppr(am0,'ph',at)
      am0<-SptImputationAppr(am0,'sws',at)
      
      if (scope=='Global'){
        am0<-SptImputationAppr(am0,'TradeNonNorm_Consumption',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNonNorm_Exports',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNonNorm_Imports',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNonNorm_Production',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNorm_Consumption',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNorm_Exports',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNorm_Imports',at,ctff=2)
        am0<-SptImputationAppr(am0,'TradeNorm_Production',at,ctff=2)
      }
      
      if (scope=='China'){
        am0<-SptImputationAppr(am0,'PriSect',at,ctff=1)
        am0<-SptImputationAppr(am0,'SecSect',at,ctff=1)
        am0<-SptImputationAppr(am0,'TerSect',at,ctff=1)
      }
      
      am0<-data.frame(am0)
    }
    #Bathymetry
    #summary(am0$Bth)
    #am0$Bth<-am0$Bth+1000 #+3500 for global, not needed anymore
    #LUC
    {
      # 1 #Cropland
      # 2 #Forest
      # 3 #Shrubs and other sparse vegetation
      # 4 #Herbaceous wetland; Moss and lichen
      # 5 #Urban
      # 6, #Water
      # 7 #Ice and snow
      # av<-readRDS('F:/China_data/Code files/Datasets/LUC_20km.rds')
      # table(av$Val);length(which(is.na(av$Val)))
      # av<-sapply(chr(av$Val),function(x){u<-strsplit(x,'-')[[1]];u<-u[!is.na(u) & u!=6];if (length(u)>0){return(paste(sort(u),collapse='-'))}else{return(NA)}}) #chr(am0$LUC)
      # av[av=='NA']<-NA
      # av[!is.element(av,chr(1:7))]<-NA
      # table(av);length(which(is.na(av)))
      # am0$LUC<-factor(av,levels=chr(1:5))
      # am0<-SptImputation1(am0,'LUC',mdl)
    }
    #Subregion
    am0$Subregion<-chr(am0$Subregion)
    am0$Subregion[is.na(am0$Subregion)]<-'Other'
    am0$Subregion<-fct(am0$Subregion)
    #Protected areas
    # if (scope=='Global'){
    #   #hist(am0$ProtAreas)
    #   am0$ProtAreas[am0$ProtAreas<0.5]<-0
    #   am0$ProtAreas[am0$ProtAreas>=0.5]<-1
    #   am0$ProtAreas<-fct(am0$ProtAreas)
    # }
  }
  
  ## Add difference variables
  {  
    #Calculate difference variables
    av<-unique(am0$Period)
    dd0<-vector('list',length(av)-1)
    for (i in 1:(length(av)-1)){
      
      #Numeric variables
      an<-sapply(1:ncol(am0),function(x){class(data.frame(am0)[,x])});an<-an=='numeric'
      dd0[[i]]<-data.frame(am0[am0$Period==av[i+1],an])-data.frame(am0[am0$Period==av[i],an])
      
      if (type=='base'){tt<-'2020'}
      if (type=='proj'){tt<-'2025'}
      #Factor variables
      for (j in c('PtID','x','y','Subregion')){
        dd0[[i]][,j]<-data.frame(am0[am0$Period==tt,j])
      }
      
      #Add Period
      dd0[[i]]$Period<-paste0(av[i+1],'-',av[i])
      
    }
    dd0<-dplyr::bind_rows(dd0)
    dd0$Period<-fct(dd0$Period)
    dd0<-dd0[,names(am0)]
    
    #Merge data frames
    #av<-sapply(chr(dd0$Period),function(x){strsplit(x,'-',fixed=T)[[1]][1]}) #Later year from difference period (for notation)
    av<-sapply(chr(dd0$Period),function(x){strsplit(x,'-',fixed=T)[[1]][2]}) #Earlier year from difference period (for notation; e.g. 2015 denotes 2020-2015 difference)
    am<-dd0[,5:length(dd0)];names(am)<-sapply(names(am),function(x){paste0(x,'_D')})
    am$PtID<-dd0$PtID;am$Period<-av
    df<-left_join(data.frame(am0),am,join_by('PtID','Period'))
    df<-df[,!is.element(names(df),c("Bth_D","ShrLength_D","Subregion_D"))] #No change, no projections, label variable
    df$Period<-fct(df$Period)
    am0<-df
  }
  
  ## Add derived variables
  {
    if (type!='base'){
      am0$OffshrAq_AN_F<-am0$OffshrAq_AN+am0$OffshrAq_AN_D
      am0$OffshrAq_PL_F<-am0$OffshrAq_PL+am0$OffshrAq_PL_D
    }
    if (type=='base'){
      am0$OffshrAq_AN_F<-am0$OffshrAq_AN+am0$OffshrAq_AN_D#nmr(NA)
      am0$OffshrAq_PL_F<-am0$OffshrAq_PL+am0$OffshrAq_PL_D#nmr(NA)
    }
    #Add absence/presence variable
    for (spc in c('AN','PL','AN_F','PL_F')){
      astr<-paste0('OffshrAq___',spc)
      astr1<-paste0('OffshrAq_',spc)
      am0[,astr]<-NA
      am0[,astr][am0[,astr1]==0]<-'0'
      am0[,astr][am0[,astr1]>0]<-'1'
      am0[,astr]<-fct(am0[,astr])
    }
    #Add buffer variables
    for (ac in c('','_D')){
      if (type=='base' & ac=='_D'){next()}
      for (spc in c('AN','PL')){
        astr<-paste0('OffshrAqBuffer___',spc,ac)
        astr1<-paste0('OffshrAqBuffer_',spc,ac)
        am0[,astr]<-NA
        am0[,astr][am0[,astr1]<0]<-'0' #Not relevant, focus only on area in the buffer
        am0[,astr][am0[,astr1]==0]<-'0'
        am0[,astr][am0[,astr1]>0]<-'1'
        am0[,astr]<-fct(am0[,astr])
      }
    }
    #Add previous periods as explanatory variable. Up to 3 periods
    if (type=='hist'){
      for (spc in c('AN','PL')){
        am0[,paste0('NumPeriods_',spc)]<-NA
        for (i in unique(am0$PtID)){
          am0[am0$PtID==i & am0$Period=='2015',paste0('NumPeriods_',spc)]<-length(which(am0[am0$PtID==i & is.element(am0$Period,c('2005','2010','2015')),paste0('OffshrAq___',spc)]=='1'))
          #Should be 2020, left like that because for formatting issues in FarmArea Model
        }
      }
    }
    if (type!='hist'){
      am0$OffshrAqBuffer___AN_D<-nmr(NA)
      am0$OffshrAqBuffer___PL_D<-nmr(NA)
      am0$NumPeriods_AN<-nmr(NA)
      am0$NumPeriods_PL<-nmr(NA)
    }
    
    # #Add change patterns
    # for (spc in c('AN','PL')){
    #   astr<-paste0('ChgPat_',spc)
    #   astr1<-paste0('OffshrAq_',spc)
    #   astr2<-paste0('OffshrAq_',spc,'_F')
    #   am0[,astr]<-NA
    #   am0[,astr][am0[,astr1]==0 & am0[,astr2]==0]<-'00'
    #   am0[,astr][am0[,astr1]==0 & am0[,astr2]>0]<-'01'
    #   am0[,astr][am0[,astr1]>0 & am0[,astr2]==0]<-'10'
    #   am0[,astr][am0[,astr1]>0 & am0[,astr2]>0]<-'11'
    #   am0[,astr]<-fct(am0[,astr])
    # }
  }

  #Recalculate GDPC
  if (type=='base'){
    am0$GDPC<-am0$GDPC/am0$Population
    am0$GDPC[!is.finite(am0$GDPC)]<-0 #Omitting infinites length(am$GDPC[is.infinite(am$GDPC)])
    if (scope=='China'){TMP<-chr(seq(2005,2015,5))}
    if (scope=='Global'){TMP<-chr(2015)}
    for (yr in TMP){
      am0$GDPC_D[am0$Period==chr(yr)]<-am0$GDPC[am0$Period==chr(nmr(yr)+5)]-am0$GDPC[am0$Period==chr(nmr(yr))]
    }
  }
  
  if (is.null(scn)){ac<-'0'}else{ac<-scn}
  saveRDS(am0,paste0('F:/',scope,'_data/Code files/Modeling/Model_df_',type,'_',ac,'.rds'))
  #av<-vect(am0,geom=c('x','y'),crs='epsg:4326');plot(av);plot(av[is.element(av$PtID,gg)],add=T,col='red')
  
} #Includes recalculation of GDP
PrepareProjData<-function(scn,scope){
  
  df<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_','base','_','0','.rds'))
  pj<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_','proj','_',scn,'.rds'))
  
  ## Create and arrange columns
  #print(names(df)[which(!is.element(names(df),names(pj)))])
  ac<-names(pj)[which(!is.element(names(pj),names(df)))]
  for (i in ac){df[,ac]<-NA}
  df<-df[,names(pj)]
  all(names(df)==names(pj))
  #sapply(1:length(df),function(x){length(which(is.na(df[,x])))})
  #sapply(1:length(pj),function(x){length(which(is.na(pj[,x])))})
  
  {
    # #Calculate NumPeriods
    # TMP<-chr(seq(2005,2015,5))
    # av<-data.frame(matrix(NA,!!1273,3));names(av)<-c('PtID','NumPeriods_AN','NumPeriods_PL')
    # av$PtID<-chr(1:nrow(av))
    # am<-vector('list',2);names(am)<-c('AN','PL')
    # for (x in c('AN','PL')){
    #   vrnm<-paste0('OffshrAq_',x,'_projht_')
    #   am[[x]]<-readRDS(paste0('F:/',scope,'_data/Code files/Datasets/',vrnm,'_',aqrs,'.rds'))
    #   
    #   am[[x]]$Period[am[[x]]$Period=='2003']<-'2005'
    #   am[[x]]$Val[am[[x]]$Val>0]<-1;am[[x]]$Val<-chr(am[[x]]$Val)
    #   
    #   for (i in unique(am[[x]]$PtID)){
    #     #av[av$PtID==i,paste0('NumPeriods_',x)]<-length(which(am[[x]]$Val[am[[x]]$PtID==i & is.element(am[[x]]$Period,TMP)]=='1'))
    #     df[df$PtID==i,paste0('NumPeriods_',x)]<-length(which(am[[x]]$Val[am[[x]]$PtID==i & is.element(am[[x]]$Period,TMP)]=='1'))
    #   }
    # }
    # #table(df$NumPeriods_AN)
    # #table(df$NumPeriods_PL)
    
    # #Scale water quality with baseline data for projections (multiplicative; old)
    # am<-rbind(df,pj)
    # # For visual checking
    # # kk<-'chl';xx<-data.frame(matrix(NA,length(seq(2015,2100,5)),6));row.names(xx)<-seq(2015,2100,5)
    # # l<-1;for (jj in seq(2015,2100,5)){xx[l,]<-summary(am[am$Period==jj,kk]);l<-l+1};xx
    # cc<-c('thetao','chl','so','vhm0','ph')
    # sclfct<-colMeans(df[df$Period=='2020',cc])/colMeans(pj[pj$Period=='2030',cc]-(pj[pj$Period=='2040',cc]-pj[pj$Period=='2030',cc]))
    # for (i in 1:4){pj[,cc[i]]<-pj[,cc[i]]*sclfct[i]}
    # for (cc in c('chl','thetao','vhm0')){pj[pj$Period=='2025',cc]<-rowMeans(cbind(df[df$Period=='2020',cc],pj[pj$Period=='2030',cc]))}
    
    #Not used, just for format
    # df$OffshrAqBuffer___AN_D<-factor(df$OffshrAqBuffer___AN_D,levels=0:1)
    # df$OffshrAqBuffer___PL_D<-factor(df$OffshrAqBuffer___PL_D,levels=0:1)
    # pj$OffshrAqBuffer___AN_D<-factor(pj$OffshrAqBuffer___AN_D,levels=0:1)
    # pj$OffshrAqBuffer___PL_D<-factor(pj$OffshrAqBuffer___PL_D,levels=0:1)
    
    # #Adjust sample extent
    # df<-df[df$y>=(-60) & df$y<=80,]
    # pj<-pj[pj$y>=(-60) & pj$y<=80,]
    # 
    # #PtID adjustment
    # Already done
  }#Old
  
  #Rescale projections to match baseline and recalculate 2020 difference
  am<-rbind(df,pj)
  {
    #library(ggplot2)
    #am[am$PtID==50000,]
    #ggplot(am[is.element(am$Period,chr(seq(2015,2030,5))),],aes(x=Period,y=Population_D))+geom_boxplot()
    
    # k<-'chl'
    # par(mfrow=c(2,2))
    # plot(lt(am[am$Period=='2015',k]),lt(am[am$Period=='2020',k]),pch=16,col=rgb(0,0,1,0.05))
    # plot(lt(am[am$Period=='2020',k]),lt(am[am$Period=='2025',k]),pch=16,col=rgb(0,1,1,0.05))
    # plot(lt(am[am$Period=='2025',k]),lt(am[am$Period=='2030',k]),pch=16,col=rgb(1,0,1,0.05))
    # plot(am[am$Period=='2015',k],am[am$Period=='2020',k],pch=16,col=rgb(0,0,1,0.05))
    # points(am[am$Period=='2020',k],am[am$Period=='2025',k],pch=16,col=rgb(0,1,1,0.05))
    # points(am[am$Period=='2025',k],am[am$Period=='2030',k],pch=16,col=rgb(1,0,1,0.05))
  }#Old
  for (i in c('Population','GDPC')){
  
    k1<-am[am$Period=='2025',paste0(i)]-am[am$Period=='2020',paste0(i)] #Intercept
    k2<-rowMeans(cbind(am[am$Period=='2015',paste0(i,'_D')],am[am$Period=='2025',paste0(i,'_D')])) #Slope
    
    # print(i)
    # print(summary(k1+k2))
    
    #am[am$PtID==60000,]
    for (j in chr(seq(2025,2050,5))){am[am$Period==j,paste0(i)]<-am[am$Period==j,paste0(i)]-k1+k2*1}
    am[am$Period=='2020',paste0(i,'_D')]<-k2
  }
  for (i in c("LUC_Crp","LUC_Frs","LUC_Urb","LUC_Wtr","thetao","so","chl","ph","sws")){
    
    k1<-am[am$Period=='2030',paste0(i)]-am[am$Period=='2020',paste0(i)] #Intercept
    k2<-rowMeans(cbind(am[am$Period=='2015',paste0(i,'_D')],am[am$Period=='2030',paste0(i,'_D')]))
    
    # print(i)
    # print(summary(k1+k2))
    
    #am[am$PtID==40000,]
    for (j in chr(seq(2030,2050,5))){am[am$Period==j,paste0(i)]<-am[am$Period==j,paste0(i)]-k1+k2*2}
    am[am$Period=='2020',paste0(i,'_D')]<-k2
    am[am$Period=='2025',paste0(i)]<-am[am$Period=='2020',paste0(i)]+k2
    am[am$Period=='2025',paste0(i,'_D')]<-k2
    
  }
  
  #Recalculate 2020 differences
  if (scope=='China'){
    ac1x<-c('PriSect','SecSect','TerSect')
  }
  if (scope=='Global'){
    ac1x<-c("TradeNorm_Exports","TradeNorm_Imports","TradeNorm_Consumption","TradeNorm_Production")
  }
  for (i in c("Fsh",
              ac1x,
              "ProtAreas")){
    am[am$Period=='2020',paste0(i,'_D')]<-am[am$Period=='2025',paste0(i)]-am[am$Period=='2020',paste0(i)]
  }
  
  #Recalculate GDPC
  am$GDPC<-am$GDPC/am$Population
  am$GDPC[!is.finite(am$GDPC)]<-0 #Omitting infinites length(am$GDPC[is.infinite(am$GDPC)])
  TMP<-chr(seq(2020,2045,5))
  for (yr in TMP){
    am$GDPC_D[am$Period==chr(yr)]<-am$GDPC[am$Period==chr(nmr(yr)+5)]-am$GDPC[am$Period==chr(nmr(yr))]
  }
  
  #list(df,pj)
  am
  
} #Includes recalculation of GDP
ApplyRSAssumptions<-function(am){
  for (spc in c('AN','PL')){
    ww<-am$PtID[am[,'Period']=='2015' & am[,paste0('OffshrAq___',spc)]=='1' & am[,paste0('OffshrAq___',spc,'_F')]=='0'] #Points where there are farms in 2015 but not 2020
    #print(table(ww))
    #!is.unsorted(ww)
    av<-am[am[,'Period']=='2015' & is.element(am$PtID,ww),paste0('OffshrAq_',spc)] #Area at 2015
    #Update 2015 data
    am[am[,'Period']=='2015' & is.element(am$PtID,ww),paste0('OffshrAq___',spc,'_F')]<-'1'
    am[am[,'Period']=='2015' & is.element(am$PtID,ww),paste0('OffshrAq_',spc,'_F')]<-av
    #Update 2020 data
    am[am[,'Period']=='2020' & is.element(am$PtID,ww),paste0('OffshrAq___',spc)]<-'1'
    am[am[,'Period']=='2020' & is.element(am$PtID,ww),paste0('OffshrAq_',spc)]<-av
  }
  am
}
GenValGLMM<-function(sm,mm,distr,lnkfx,posintercept){
  
  am0<-mm
  am0[,"(Intercept)"]<-1
  #am0$`I(OffshrAq_AN^2)`<-am0$OffshrAq_AN^2
  #am0$`I(OffshrAq_PL^2)`<-am0$OffshrAq_PL^2
  
  cf<-summary(sm)$coefficients
  #Conditional model (fixed effects)
  am1<-data.frame(cf$cond)
  #if (posintercept==T & am1['(Intercept)','Estimate']<0){am1['(Intercept)','Estimate']<-0}
  if (nrow(am1)>0){
    cls<-row.names(am1)#[!is.element(row.names(am),c("LUC1","LUC_DNoCh"))]
    smc1<-as.matrix(am0[,cls])%*%as.matrix(mapply(function(x,y){rnorm(1,x,y)},x=am1[cls,'Estimate'],y=0)) #am1[cls,'Std..Error']
  }else{
    smc1<-as.matrix(1)
  }
  #Dispersion model (fixed effects)
  am4<-data.frame(cf$disp)
  if (nrow(am4)>0){
    cls<-row.names(am4)
    smd<-as.matrix(am0[,cls])%*%as.matrix(mapply(function(x,y){rnorm(1,x,y)},x=am4[cls,'Estimate'],y=0)) #am4[cls,'Std..Error']
  }else{
    smd<-as.matrix(1)
  }
  
  #Expected value
  if (lnkfx=='inverse'){mu<-1/(c(smc1))}
  if (lnkfx=='log'){mu<-exp(c(smc1))}
  if (lnkfx=='identity'){mu<-smc1}
  # if (lnkfx=='sqrt'){mu<-(smc1)^2}
  # if (lnkfx=='cloglog'){mu<-exp(log(-log(1-smc1)))}
  p<-family_params(sm) #Model parameter
  if (nrow(am4)>0){
    phi<-exp(smd)[,1] #Dispersion; Log-link
  }else{
    phi<-sigma(sm) #Dispersion, ~1
    if (distr=='Gamma'){phi<-sigma(sm)^2}
  }
  
  if (distr=='Gaussian'){vl<-mapply(function(x,y){rnorm(1,x,y)},x=mu,y=phi)}
  if (distr=='Gamma'){
    vl<-mapply(function(x,z){rgamma(1,shape=x,scale=z)},x=1/phi,z=mu*phi)
    #https://github.com/glmmTMB/glmmTMB/issues/990 for squared term in gamma
    #Old x=mu^2/phi^2,z=phi^2/mu #phi is the scale parameter, not the shape (as documentation says), assuming equation is correct
    #https://github.com/glmmTMB/glmmTMB/issues/990
  }
  if (distr=='Tweedie'){
    if (posintercept==F){vl<-mapply(function(x,y,z){mgcv::rTweedie(x,y,z)},x=mu,y=p,z=phi)}
    if (posintercept==T){print(c(length(which(mu<0)),length(mu)));vl<-mapply(function(x,y,z){mgcv::rTweedie(x,y,z)},x=sapply(mu,function(x){max(0,x)}),y=p,z=phi)}
  }
  
  # vv<-predict(sm,am0,type="conditional",se.fit=T,re.form=NULL) #se,fit gives error value compared with expected value (deterministic)
  # plot(vv$fit,am0[,paste0('OffshrAq_',spc,'_F')]);lines(0:5,0:5)
  
  vl
}

TrnsfrmData<-function(am){
  #k<-'VesselFrq';par(mfrow=c(3,1));hist(am[,k]);hist(sqrt(am[,k]));hist(lt(am[,k]))
  av<-c('OffshrAq_AN','OffshrAq_PL','OffshrAq_CB',
        'OffshrAq_AN_F','OffshrAq_PL_F','OffshrAq_CB_F',
        'Population','GDPC','HDI',#'VesselFrq','InfrstrFrq','TourismValue','Bth','ShrLength',
        'Population_D','GDPC_D',#'HDI_D','VesselFrq_D','InfrstrFrq_D','TourismValue_D',
        'Fsh','Fsh_D')#'PriSect','SecSect','TerSect','PriSect_D','SecSect_D','TerSect_D','ProtAreas'
  for (i in names(am)[is.element(names(am),av)]){
    am[,i]<-lt(am[,i])
  }
  # am$OffshrAq_AN<-round(am$OffshrAq_AN*1000,0) #COMPoisson(link="loglambda")
  # am$OffshrAq_AN[am$OffshrAq_AN==0]<-1 #Gamma(link=log)
  am
}
BackTrnsfrmData<-function(am){
  #k<-'VesselFrq';par(mfrow=c(3,1));hist(am[,k]);hist(sqrt(am[,k]));hist(lt(am[,k]))
  av<-c('OffshrAq_AN','OffshrAq_PL','OffshrAq_CB',
        'OffshrAq_AN_D','OffshrAq_PL_D','OffshrAq_CB_D',
        'OffshrAq_AN_F','OffshrAq_PL_F','OffshrAq_CB_F',
        'Population','HDI','GDPC','VesselFrq','InfrstrFrq','TourismValue','Bth','ShrLength',
        'Population_D','HDI_D','GDPC_D','VesselFrq_D','InfrstrFrq_D','TourismValue_D',
        'Fsh','PriSect','SecSect','TerSect','ProtAreas',
        'Fsh_D','PriSect_D','SecSect_D','TerSect_D')
  for (i in names(am)[is.element(names(am),av)]){
    am[,i]<-et(am[,i])
  }
  # am$OffshrAq_AN<-round(am$OffshrAq_AN*1000,0) #COMPoisson(link="loglambda")
  # am$OffshrAq_AN[am$OffshrAq_AN==0]<-1 #Gamma(link=log)
  am
}
NormData<-function(am,type,ac,save=F,scope){
  
  #an<-sapply(1:ncol(am),function(x){class(data.frame(am)[,x])});an<-an=='numeric'
  
  if (save==T){
    np<-data.frame(Var=chr(NA),Par1=rep(nmr(NA),length(am)),Par2=rep(nmr(NA),length(am)))
    for (i in 1:length(am)){
      
      if (is.element(names(am)[i],c('PtID','x','y','Period')) | class(am[,i])=='factor'){next()}
      if (type=='sd'){
        np[i,1]<-names(am)[i]
        np[i,2]<-mean(am[,i])
        np[i,3]<-sd(am[,i])
        am[,i]<-(am[,i]-np[i,2])/np[i,3]
      }
      if (type=='mm'){
        np[i,1]<-names(am)[i]
        np[i,2]<-min(am[,i]) #,na.rm=T
        np[i,3]<-max(am[,i])
        am[,i]<-(am[,i]-np[i,2])/(np[i,3]-np[i,2])
      }
      np<-np[!is.na(np$Var),]
      saveRDS(np,paste0("F:/",scope,"_data/Code files/Modeling/",ac,".rds"))
    }
    return(np)
  }
  
  if (save==F){
    np<-readRDS(paste0("F:/",scope,"_data/Code files/Modeling/",ac,".rds"))
    row.names(np)<-np$Var
    for (i in row.names(np)){
      if (type=='sd'){am[,i]<-(am[,i]-np[i,2])/np[i,3]}
      if (type=='mm'){am[,i]<-(am[,i]-np[i,2])/(np[i,3]-np[i,2])}
    }
    return(am)
  }
  
}
BackNormData<-function(am,type,ac,scope){
  
  # np<-readRDS(paste0("F:/China_data/Code files/Modeling/",ac,".rds"))
  # row.names(np)<-np$Var
  # for (i in np$Var){
  #   if (type=='sd'){am[,i]<-am[,i]*np[i,3]+np[i,2]}
  #   if (type=='mm'){am[,i]<-am[,i]*(np[i,3]-np[i,2])+np[i,2]}
  # }
  
  np<-readRDS(paste0("F:/",scope,"_data/Code files/Modeling/",ac,".rds"))
  row.names(np)<-np$Var
  for (i in row.names(np)){
    if (type=='sd'){am[,i]<-am[,i]*np[i,3]+np[i,2]}
    if (type=='mm'){am[,i]<-am[,i]*(np[i,3]-np[i,2])+np[i,2]}
  }
  
  am
  
}
TrnsNormVect<-function(av,ss,ac,scope){
  
  np<-readRDS(paste0("F:/",scope,"_data/Code files/Modeling/",ac,".rds"))
  row.names(np)<-np$Var
  
  av<-(lt(av)-np[ss,2])/(np[ss,3]-np[ss,2])
  
  av
}
BackTrnsNormVect<-function(av,ss,ac,scope){
  
  np<-readRDS(paste0("F:/",scope,"_data/Code files/Modeling/",ac,".rds"))
  row.names(np)<-np$Var
  
  av<-et(av*(np[ss,3]-np[ss,2])+np[ss,2])
  
  av
}

FitModel2<-function(rfave,rfdsp,rfzin,am,famlnk,mdllbl,spc,mdltp,optmth,aa){
  
  # glmmTMB
  # #?family_glmmTMB
  
  #Data preparation
  #https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
  if (!is.null(aa)){ #now as AA
    #Spatial groups
    #Based on custom distance matrix (different correlation coefficient for different groups)
    {
      # am$group1 <- NA
      # for (i in 1:length(aa)){
      #   am$group1[is.element(am$PtID,aa[[i]])]<-i
      # }
      # am$group1 <- factor(am$group1,sort(unique(am$group1)))
      # #table(am$group1);length(which(is.na(am$group1)))
    }
    #All points together(same correlation coefficient for different groups)
    {
      #am$group1 <- factor(rep(1,nrow(am))) #Constant across periods
      am$group1 <- factor(am$Period,sort(unique(am$Period))) #Different for each period
    }
    am$pos <- numFactor(am$x,am$y)
    #Temporal groups
    #Transform to unit time step
    am$Period<-chr(am$Period)
    av<-sort(unique(am$Period));for (i in 1:length(av)){am$Period[am$Period==av[i]]<-i}
    am$Period<-factor(am$Period,levels=sort(unique(am$Period)))
    am$Prd<-numFactor(am$Period) #For checking if somehow it changes results (probably not)
    am$group2 <- factor(am$PtID,levels=sort(unique(am$PtID))) #Set of time series
    #table(am$group2);length(which(is.na(am$group2)))
  }
  
  library(optimx)
  #checkallsolvers()
  
  if (optmth=='bobyqa'){ctrl<-glmmTMBControl(optimizer=optimr,optArgs=list(method='bobyqa',maxfun=1000),rank_check="adjust")}
  
  if (optmth=='nlminb'){ctrl<-glmmTMBControl(optimizer=nlminb,optArgs=list(iter.max=1000),parallel=20,rank_check="adjust")} #Default optimizer
  if (optmth=='BFGS'){ctrl<-glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS",maxit=1000),parallel=20,rank_check="adjust")}
  if (optmth=='L-BFGS-B'){ctrl<-glmmTMBControl(optimizer=optim,optArgs=list(method="L-BFGS-B",maxit=1000),parallel=20,rank_check="adjust")}
  if (optmth=='Nelder-Mead'){ctrl<-glmmTMBControl(optimizer=optim,optArgs=list(method="Nelder-Mead",maxit=5000),parallel=20,rank_check="adjust")}
  if (optmth=='SANN'){ctrl<-glmmTMBControl(optimizer=optim,optArgs=list(method="SANN",maxit=1000),parallel=20,rank_check="adjust")}
  
  sm<-tryCatch(glmmTMB(formula=rfave,
                       family=famlnk,
                       ziformula=rfzin,
                       dispformula=rfdsp,
                       data=am,
                       REML=F,
                       verbose=F,
                       control=ctrl),
               error=function(er){'Error'})
  #AIC(sm);summary(sm);plot(sm)
  
  # system.time(sm<-tryCatch(glmmTMB(formula=rfave,
  #                      family=famlnk,
  #                      ziformula=rfzin,
  #                      dispformula=rfdsp,
  #                      data=am,
  #                      REML=F,
  #                      verbose=F,
  #                      control=ctrl),
  #              error=function(er){'Error'}))
  
  if (sm[1]!='Error'){
    saveRDS(sm,paste0("F:/China_data/Code files/Modeling/",mdllbl,"_",spc,"_",mdltp,'_',optmth,".rds"))
    #print(c(rfave,'Saved'))
    return(sm)
  }else{
    print(sm)
  }
  
} #glmmTMB

GetSummaryTable<-function(am){
  uu<-data.frame(matrix(NA,length(names(am)),7+1));row.names(uu)<-names(am);colnames(uu)<-c('Var','Min','1stQu','2ndQu','3rdQu','Max','NmNA')
  uu$Var<-names(am)
  #sapply(am,class)
  for (i in names(am)){
    if (is.element(i,c('PtID','x','y'))){next()}
    if (is.element(i,c('Period'))){print(table(am[,i]));next()}
    if (class(am[,i])=='factor'){
      av<-table(am[,i])
      uu[i,c('Min')]<-av[1]
      uu[i,c('Max')]<-av[2]
      next()
    }
    av<-quantile(am[,i],seq(0,1,0.25),na.rm=T)
    uu[i,c('Min','1stQu','2ndQu','3rdQu','Max')]<-formatC(av,2,format='f',flag='0')
    av<-length(which(is.na(am[,i])))
    uu[i,c('NmNA')]<-formatC(av,0,format='f',flag='0')
  }
  uu
}

pdp<-function(sm,am,vr,fx){
  #Based on paper "pdp: An R Package for Constructing". Verified with examples from documentation.
  k<-500;if (nrow(am)>k){am<-am[sample(1:nrow(am),k,F),]}
  y<-sapply(am[,vr],function(x){
    d<-am;d[,vr]<-x;d[,paste0(vr,'2')]<-x^2
    fx(predict(sm,d,type='response',se.fit=F))
  })
  data.frame(am[,vr],y)
}


GetSpatStr<-function(am0,dco){
  
  # if (scope=='China'){ddm<-readRDS(paste0("F:/",scope,"_data/Code files/Boundaries/SmplPtsDist_10km_v2.rds"))}
  # if (scope=='Global'){ddm<-readRDS(paste0("F:/",scope,"_data/Code files/Administrative boundaries/SmplPtsDist_10km.rds"))}
  spt<-vect(am0,geom=c('x','y'),crs='epsg:4326')
  ddm<-distance(spt,spt,unit="m") #m
  ddm<-as.matrix(ddm);colnames(ddm)<-1:ncol(ddm);row.names(ddm)<-1:nrow(ddm)
  ddd<-ddm/1000 #To km
  #https://cran.r-project.org/web/packages/SDPDmod/vignettes/spatial_matrices.html
  W<-SDPDmod::InvDistMat(ddd,distCutOff=dco,powr=1) #50 used before
  W<-SDPDmod::eignor(W) #W<-SDPDmod::rownor(W)
  #table(rowSums(W))
  library(igraph)
  TT<-W #[1:50,1:50] #For testing
  TT[TT>0]<-1 #Assumes all points below distance threshold above are neighbours
  topology <- graph_from_adjacency_matrix(TT,'undirected')
  g3 <- simplify(topology)
  #https://stackoverflow.com/questions/23686729/how-to-identify-fully-connected-node-clusters-with-igraph
  #aa<-max_cliques(g3) #This is not right, points repeated in different groups (?)
  #https://stackoverflow.com/questions/55442939/count-number-of-disconnected-sub-networks
  AA<-components(g3)
  #plot(g3,mark.shape=AA$membership)
  print(c(AA$no,nrow(am0)))
  
  am0$SptStr<-factor(AA$membership,levels=1:AA$no)
  ll<-vect(am0,geom=c('x','y'),crs='epsg:4326')
  plot(ll,'SptStr')
  
  list(am0,ll)
}
FeatSel<-function(spc,rfave,amt,tgt,rfp,j,caseweights,classweights){
  
  library(Boruta)
  ac<-strsplit(chr(rfave)[3],' \\+ ',perl=T)[[1]]
  ac<-gsub('\n','',ac,fixed=T)
  ac<-gsub(' ','',ac,fixed=T)
  rfmvs<-tryCatch(
    Boruta(amt[,ac],amt[,tgt],
           pValue = 0.01,maxRuns = 100,
           doTrace = 0,holdHistory = T,
           #getImp = getImpLegacyRfZ,ntree=rfp[1,1],mtry=rfp[1,2],nodesize=rfp[1,3]) #Uses randomForest implementation
           getImp = getImpRfZ,num.trees=rfp[j,1],mtry=rfp[j,2],min.node.size=rfp[j,3],
           probability=T,case.weights=caseweights,class.weights=classweights), #Uses ranger implementation  
    error=function(e){c()}
  )
  
  if (length(rfmvs)==0){return(NULL)}
  
  #if (roughfix==T){rfmvs<-TentativeRoughFix(rfmvs)}
  av1<-rfmvs$finalDecision
  
  av2<-data.frame(rfmvs$ImpHistory)
  av2<-data.frame(tidyr::pivot_longer(av2,names(av2),names_to='Var',values_to='Val'))
  av2$Var<-factor(av2$Var,levels=c(
    'OffshrAq_AN','OffshrAq_PL',
    'OffshrAqBuffer_AN','OffshrAqBuffer_PL',
    'OffshrAq___AN','OffshrAq___PL',
    'OffshrAqBuffer___AN','OffshrAqBuffer___PL',
    'Population','GDPC','Fsh','Trade_ImpRat','Trade_ExpRat','PriSect','SecSect','TerSect',
    'Population_D','GDPC_D','Trade_ImpRat_D','Trade_ExpRat_D','Fsh_D','PriSect_D','SecSect_D','TerSect_D',
    'TradeNonNorm_Consumption','TradeNonNorm_Exports','TradeNonNorm_Imports','TradeNonNorm_Production',
    'TradeNorm_Consumption','TradeNorm_Exports','TradeNorm_Imports','TradeNorm_Production',
    'TradeNonNorm_Consumption_D','TradeNonNorm_Exports_D','TradeNonNorm_Imports_D','TradeNonNorm_Production_D',
    'TradeNorm_Consumption_D','TradeNorm_Exports_D','TradeNorm_Imports_D','TradeNorm_Production_D',
    'thetao','so','chl','ph','sws',
    'thetao_D','so_D','chl_D','ph_D','sws_D',
    'LUC_Crp','LUC_Frs','LUC_Urb','LUC_Wtr',
    'LUC_Crp_D','LUC_Frs_D','LUC_Urb_D','LUC_Wtr_D',
    'Bth','ShrLength','ProtAreas',
    'Subregion',
    'shadowMax','shadowMean','shadowMin'
  ))
  av2$Spc<-spc
  
  list(av1,av2)
}
GetWeights<-function(amt,tgt){
  
  #Observation weights (continuous)
  if (class(amt[,tgt])=='numeric'){
    ww<-cut(amt[,tgt],seq(0,1,0.1),right=F);ww[is.na(ww)]<-levels(ww)[length(levels(ww))]
    tt<-table(ww);tt[length(tt)]<-tt[length(tt)]+length(which(is.na(ww)))
    tt<-length(tt)/tt;tt[!is.finite(tt)]<-0
    tt<-data.frame(Cat=names(tt),Val=nmr(tt))
    W<-dplyr::left_join(data.frame(Cat=ww),tt,'Cat')$Val
  }
  
  #Observation weights (discrete)
  if (class(amt[,tgt])=='factor'){
    # #https://datascience.stackexchange.com/questions/44755/why-doesnt-class-weight-resolve-the-imbalanced-classification-problem?noredirect=1&lq=1
    W<-table(amt[,tgt]);W<-sum(W)/W
    #if (scope=='Global'){W<-table(amt[,tgt]);W<-rev(W/sum(W));names(W)<-rev(names(W))}
  }
  
  W
}
GetAccMetrics<-function(sm,amv,tgt){
  
  vv<-predict(sm,amv)$predictions
  p1<-vv;p2<-amv[,tgt]
  
  if (class(amv[,tgt])=='factor'){
    xx<-cbind(p1,p2)
    cf<-matrix(0,2,2)
    #Rows: Predicted, Columns: Truth #Ordered when assigning levels to factors
    for (z in 1:nrow(xx)){cf[xx[z,1],xx[z,2]]<-cf[xx[z,1],xx[z,2]]+1}
    cv<-sum(diag(cf))/sum(cf)
  }
  
  if (class(amv[,tgt])=='numeric'){
    cv<-sum((vv-amv[,tgt])^2)/nrow(amv)
    cf<-NULL
  }
  
  list(cv,cf)
  
}
#
AddVarCategories<-function(am,scope){
  
  if (scope=='Global'){
    #x<-strsplit(chr(formula(eval(rfave0)))[3],' + ',fixed=T)[[1]]
    #paste(x,collapse="','")
    x<-c('OffshrAq___AN','OffshrAq___PL','OffshrAqBuffer___AN','OffshrAqBuffer___PL',
         'OffshrAq_AN','OffshrAq_PL','OffshrAqBuffer_AN','OffshrAqBuffer_PL',
         'Population','GDPC','Fsh','TradeNorm_Exports','TradeNorm_Imports',
         'Population_D','GDPC_D','Fsh_D','TradeNorm_Exports_D','TradeNorm_Imports_D',
         'thetao','so','chl','ph','sws','thetao_D','so_D','chl_D','ph_D','sws_D',
         'LUC_Crp','LUC_Frs','LUC_Urb','LUC_Wtr','LUC_Crp_D','LUC_Frs_D','LUC_Urb_D','LUC_Wtr_D','Bth','ShrLength',
         'ProtAreas','Subregion')
    y<-c(rep('Industrial',4),
         rep('Industrial',4),
         rep('Socioeconomic',2),rep('Industrial',1),rep('Socioeconomic',2),
         rep('Socioeconomic',2),rep('Industrial',1),rep('Socioeconomic',2),
         rep('Environmental',10),
         rep('Geographical',10),
         rep('Institutional',2))
    ad<-data.frame(Var=x,Cat=y)
    ad$Cat<-factor(ad$Cat,levels=c('Socioeconomic','Industrial','Environmental','Geographical','Institutional'))
    return(left_join(am,ad,join_by(Var)))
  }
  
  if (scope=='China'){
    #x<-strsplit(chr(formula(eval(rfave0)))[3],' + ',fixed=T)[[1]]
    #paste(x,collapse="','")
    x<-c('OffshrAq___AN','OffshrAq___PL','OffshrAqBuffer___AN','OffshrAqBuffer___PL',
         'OffshrAq_AN','OffshrAq_PL','OffshrAqBuffer_AN','OffshrAqBuffer_PL',
         'Population','GDPC','Fsh',"PriSect","SecSect","TerSect",
         'Population_D','GDPC_D','Fsh_D',"PriSect_D","SecSect_D","TerSect_D",
         'thetao','so','chl','ph','sws','thetao_D','so_D','chl_D','ph_D','sws_D',
         'LUC_Crp','LUC_Frs','LUC_Urb','LUC_Wtr','LUC_Crp_D','LUC_Frs_D','LUC_Urb_D','LUC_Wtr_D','Bth','ShrLength',
         'ProtAreas','Subregion')
    y<-c(rep('Industrial',4),
         rep('Industrial',4),
         rep('Socioeconomic',2),rep('Industrial',1),rep('Socioeconomic',3),
         rep('Socioeconomic',2),rep('Industrial',1),rep('Socioeconomic',3),
         rep('Environmental',10),
         rep('Geographical',10),
         rep('Institutional',2))
    ad<-data.frame(Var=x,Cat=y)
    ad$Cat<-factor(ad$Cat,levels=c('Socioeconomic','Industrial','Environmental','Geographical','Institutional'))
    return(left_join(am,ad,join_by(Var)))
  }
  
}
ChangeVarNames<-function(am,scope){
  names(am)[which(names(am)=='Var')]<-'VarOld'
  
  if (scope=='Global'){
    x<-c('OffshrAq___AN','OffshrAq___PL','OffshrAqBuffer___AN','OffshrAqBuffer___PL',
         'OffshrAq_AN','OffshrAq_PL','OffshrAqBuffer_AN','OffshrAqBuffer_PL',
         'Population','GDPC','Fsh','TradeNorm_Exports','TradeNorm_Imports',
         'Population_D','GDPC_D','Fsh_D','TradeNorm_Exports_D','TradeNorm_Imports_D',
         'thetao','so','chl','ph','sws','thetao_D','so_D','chl_D','ph_D','sws_D',
         'LUC_Crp','LUC_Frs','LUC_Urb','LUC_Wtr','LUC_Crp_D','LUC_Frs_D','LUC_Urb_D','LUC_Wtr_D','Bth','ShrLength',
         'ProtAreas','Subregion')
    y<-c('Fish farm presence (at point)','Seaweed farm presence (at point)','Fish farm presence (at buffer)','Seaweed farm presence (at buffer)',
         'Fish farm area (at point)','Seaweed farm area (at point)','Fish farm area (at buffer)','Seaweed farm area (at buffer)',
         'Population','GDP per capita','Wild catch','Aquatic exports','Aquatic imports',
         'Population (diff)','GDP per capita (diff)','Wild catch (diff)','Aquatic exports (diff)','Aquatic imports (diff)',
         'Sea temperature','Salinity','Chlorophyll-A','pH','Seawater velocity',
         'Sea temperature (diff)','Salinity (diff)','Chlorophyll-A (diff)','pH (diff)','Seawater velocity (diff)',
         'Cropland proportion','Forest proportion','Urban proportion','Water proportion',
         'Cropland proportion (diff)','Forest proportion (diff)','Urban proportion (diff)','Water proportion (diff)',
         'Elevation/depth','Coastline complexity',
         'Protected area','Subregion')
    lvls<-c(
      'Population','GDP per capita','Aquatic exports','Aquatic imports',
      'Population (diff)','GDP per capita (diff)','Aquatic exports (diff)','Aquatic imports (diff)',
      'Fish farm presence (at point)','Seaweed farm presence (at point)','Fish farm presence (at buffer)','Seaweed farm presence (at buffer)',
      'Fish farm area (at point)','Seaweed farm area (at point)','Fish farm area (at buffer)','Seaweed farm area (at buffer)',
      'Wild catch','Wild catch (diff)',
      'Sea temperature','Salinity','Chlorophyll-A','pH','Seawater velocity',
      'Sea temperature (diff)','Salinity (diff)','Chlorophyll-A (diff)','pH (diff)','Seawater velocity (diff)',
      'Cropland proportion','Forest proportion','Urban proportion','Water proportion',
      'Cropland proportion (diff)','Forest proportion (diff)','Urban proportion (diff)','Water proportion (diff)',
      'Elevation/depth','Coastline complexity',
      'Protected area','Subregion')
    ad<-data.frame(VarOld=x,Var=y)
    ad$Var<-factor(ad$Var,levels=lvls)
    return(left_join(am,ad,join_by(VarOld)))
  }
    
  if (scope=='China'){
    x<-c('OffshrAq___AN','OffshrAq___PL','OffshrAqBuffer___AN','OffshrAqBuffer___PL',
         'OffshrAq_AN','OffshrAq_PL','OffshrAqBuffer_AN','OffshrAqBuffer_PL',
         'Population','GDPC','Fsh',"PriSect","SecSect","TerSect",
         'Population_D','GDPC_D','Fsh_D',"PriSect_D","SecSect_D","TerSect_D",
         'thetao','so','chl','ph','sws','thetao_D','so_D','chl_D','ph_D','sws_D',
         'LUC_Crp','LUC_Frs','LUC_Urb','LUC_Wtr','LUC_Crp_D','LUC_Frs_D','LUC_Urb_D','LUC_Wtr_D','Bth','ShrLength',
         'ProtAreas','Subregion')
    y<-c('Fish farm presence (at point)','Seaweed farm presence (at point)','Fish farm presence (at buffer)','Seaweed farm presence (at buffer)',
         'Fish farm area (at point)','Seaweed farm area (at point)','Fish farm area (at buffer)','Seaweed farm area (at buffer)',
         'Population','GDP per capita','Wild catch','Primary sector','Secondary sector','Tertiary sector',
         'Population (diff)','GDP per capita (diff)','Wild catch (diff)','Primary sector (diff)','Secondary sector (diff)','Tertiary sector (diff)',
         'Sea temperature','Salinity','Chlorophyll-A','pH','Seawater velocity',
         'Sea temperature (diff)','Salinity (diff)','Chlorophyll-A (diff)','pH (diff)','Seawater velocity (diff)',
         'Cropland proportion','Forest proportion','Urban proportion','Water proportion',
         'Cropland proportion (diff)','Forest proportion (diff)','Urban proportion (diff)','Water proportion (diff)',
         'Elevation/depth','Coastline complexity',
         'Protected area','Subregion')
    lvls<-c(
      'Population','GDP per capita','Primary sector','Secondary sector','Tertiary sector',
      'Population (diff)','GDP per capita (diff)','Primary sector (diff)','Secondary sector (diff)','Tertiary sector (diff)',
      'Fish farm presence (at point)','Seaweed farm presence (at point)','Fish farm presence (at buffer)','Seaweed farm presence (at buffer)',
      'Fish farm area (at point)','Seaweed farm area (at point)','Fish farm area (at buffer)','Seaweed farm area (at buffer)',
      'Wild catch','Wild catch (diff)',
      'Sea temperature','Salinity','Chlorophyll-A','pH','Seawater velocity',
      'Sea temperature (diff)','Salinity (diff)','Chlorophyll-A (diff)','pH (diff)','Seawater velocity (diff)',
      'Cropland proportion','Forest proportion','Urban proportion','Water proportion',
      'Cropland proportion (diff)','Forest proportion (diff)','Urban proportion (diff)','Water proportion (diff)',
      'Elevation/depth','Coastline complexity',
      'Protected area','Subregion')
    ad<-data.frame(VarOld=x,Var=y)
    ad$Var<-factor(ad$Var,levels=lvls)
    return(left_join(am,ad,join_by(VarOld)))
  }
  
  
}
#LogLoss <- function(pred, actual){-mean(actual * log(pred) + (1 - actual) * log(1 - pred))}
#LogSkillScore <- function(pred, actual){-mean(actual * log(pred) + (1 - actual) * log(1 - pred))}
BrierScore <- function(x,y){mean((x-y)^2)}
ProbScore <- function(x,y,R,mtrc){
  
  #https://stats.stackexchange.com/questions/403544/how-to-compute-the-brier-score-for-more-than-two-classes
  #https://mclust-org.github.io/mclust/reference/BrierScore.html
  am<-data.frame(fact=chr(y))
  am$fact<-factor(am$fact,levels=R)
  am<-model.matrix(formula("~0+fact"),am)
  
  if (mtrc=='Brier'){
    return(sum(
      (x-am)^2
      )/(2*nrow(am)))
  }
  if (mtrc=='LogLoss'){
    return(-sum(
      am*log(x)#+(1-am)*log(1-x)
      )/(nrow(am)))
  }
  
  # R<-1:5
  # y<-c(5,5,5,2,5,3,1,2,1,1)
  # x<-t(matrix(c(
  #   0.14, 0.38, 0.4 , 0.04, 0.05,
  #   0.55, 0.05, 0.34, 0.04, 0.01,
  #   0.3 , 0.35, 0.18, 0.09, 0.08,
  #   0.23, 0.22, 0.04, 0.05, 0.46,
  #   0.  , 0.15, 0.47, 0.28, 0.09,
  #   0.23, 0.13, 0.34, 0.27, 0.03,
  #   0.32, 0.06, 0.59, 0.02, 0.01,
  #   0.01, 0.19, 0.01, 0.03, 0.75,
  #   0.27, 0.38, 0.03, 0.12, 0.2 ,
  #   0.17, 0.45, 0.11, 0.25, 0.01
  # ),5,10))
  # 1.0068899999999998
  
  #MLmetrics::MultiLogLoss()
  
  #Focal loss
}
SkillScore <- function(x,y,R,mtrc){
  # av<-vector('numeric',10)
  # for (i in 1:length(av)){
  #   #av[i]<-1-mean((x-y)^2)/mean((sample(0:1,length(y),T,table(y)/sum(table(y)))-y)^2) #Compared with naive classification based on proportions
  #   av[i]<-1-mean((x-y)^2)/mean((0-y)^2) #Compared with naive classification based on always zero
  # }
  # mean(av)
  
  av<-vector('numeric',10)
  for (i in 1:length(av)){
    BS<-ProbScore(x,y,R,mtrc)
    BSrf<-ProbScore(x,sample(R,length(y),T,table(y)/sum(table(y))),R,mtrc)
    av[i]<-1-BS/BSrf
  }
  mean(av)
}

#### S1: Model fitting ####

library(terra)
library(dplyr)

scope<-'Global'
if (scope=='Global'){
  NN<-2500;dco<-50 #Sample and block distance model 1
  NN2<-250;dco2<-50 #Sample and block distance model 2
  {
    # ac1x<-c('Trade_ImpRat','Trade_ExpRat')
    # ac1<-'Trade_ImpRat + Trade_ExpRat +';ac2<-'Trade_ImpRat_D + Trade_ExpRat_D +'
    ac1x<-c(
      #"TradeNonNorm_Consumption","TradeNonNorm_Exports","TradeNonNorm_Imports","TradeNonNorm_Production"
      "TradeNorm_Exports","TradeNorm_Imports" #,"TradeNorm_Consumption","TradeNorm_Production"
    )
    ac2x<-c(
      #"TradeNonNorm_Consumption","TradeNonNorm_Exports","TradeNonNorm_Imports","TradeNonNorm_Production"
      "TradeNorm_Exports_D","TradeNorm_Imports_D" #,"TradeNorm_Consumption","TradeNorm_Production"
    )
    ac1<-c(
      #"TradeNonNorm_Consumption+TradeNonNorm_Exports+TradeNonNorm_Imports+TradeNonNorm_Production+"
      "TradeNorm_Exports+TradeNorm_Imports+" #TradeNorm_Consumption+TradeNorm_Production+"
    )
    ac2<-c(
      #"TradeNonNorm_Consumption_D+TradeNonNorm_Exports_D+TradeNonNorm_Imports_D+TradeNonNorm_Production_D+"
      "TradeNorm_Exports_D+TradeNorm_Imports_D+" #TradeNorm_Consumption_D+TradeNorm_Production_D+
    ) 
  }
  ar0<-rast(ext(c(-180,180,-90,90)),resolution=0.1,crs='epsg:4326')
  N<-10;N2<-10 #CV runs for model 1 and 2
}

#Prepare dataset for random forest
{
  set.seed(1000)
  #Get dataset
  am0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_','base','_','0','.rds'))
  ##Adjust sample period
  if (scope=='Global'){am0<-am0[am0$Period==2015,]}
  ##Add corrections
  if (scope=='Global'){am0<-ApplyRSAssumptions(am0)}
  ##
  amm<-vector('list',2);names(amm)<-c('AN','PL')
  np<-vector('list',2);names(amm)<-c('AN','PL')
  for (spc in c('AN','PL')){
    
    am<-am0
    
    ##Adjust sample size
    gs00<-10;spt<-vect(paste0('F:/',scope,'_data/Code files/Administrative boundaries/PopulationPointsBufferv2_',gs00,'km.shp'))
    av<-sample(spt$PtID,NN,F,spt$AvArCell)
    am<-am[is.element(am$PtID,av),]
    
    am<-TrnsfrmData(am)
    
    av<-GetSpatStr(am,dco)
    amm[[spc]]<-av[[1]]
    
    #Block sampling
    av<-sapply(unique(amm[[spc]]$SptStr),function(x){sample(which(amm[[spc]]$SptStr==x),1,F)})
    amm[[spc]]<-amm[[spc]][av,]
    
    np[[spc]]<-NormData(amm[[spc]],'mm',paste0('MinMax',spc),save=T,scope)
    amm[[spc]]<-NormData(amm[[spc]],'mm',paste0('MinMax',spc),save=F,scope)
  }
  saveRDS(amm,paste0('F:/',scope,'_data/Code files/Modeling/Model_df_TransNorm_Mod.rds'))
}
#Farm change (RF)
{
  rfp<-expand.grid(c(500,1000,1500,2000),c(15,10,5),c(5,10,20,25)) #num.trees,mtry,min.node.size 
  row.names(rfp)<-apply(rfp,1,paste0,collapse='-')
  mdltgt<-'___'
  if (mdltgt=='___'){ac0<-'FarmPresence';ii<-'factor'}
  if (mdltgt=='_'){ac0<-'FarmArea';ii<-'numeric'}
  UseWeights<-F
  VarSel<-F
  rfave0<-str2lang(
    "paste0('OffshrAq',mdltgt,spc,'_F ~','
        OffshrAq',mdltgt,'AN + OffshrAq',mdltgt,'PL +
        OffshrAqBuffer',mdltgt,'AN + OffshrAqBuffer',mdltgt,'PL +
        Population + GDPC + Fsh +',ac1,'
        Population_D + GDPC_D + Fsh_D +',ac2,'
        LUC_Crp + LUC_Frs + LUC_Urb + LUC_Wtr +
        thetao + so + chl + ph + sws +
        Bth + ShrLength + ProtAreas+
        Subregion')")
  if (scope=='Global'){Prob<-'Prob';N<-1}
  VarSelItr<-1 #Iterations of Boruta
  VarSelVal<-c('Confirmed','Tentative') #Values accepted #'Tentative'
  VarSelCtt<-1 #Variable selection iteration cutoff
  
  set.seed(2000)
  av<-vector('list',2);names(av)<-c('AN','PL')
  for (i in 1:length(av)){
    av[[i]]<-vector('list',nrow(rfp));names(av[[i]])<-row.names(rfp)
    for (j in 1:length(av[[i]])){
      av[[i]][[j]]<-vector('list',N);names(av[[i]][[j]])<-1:N
    }
  } #spc, rfp, cv
  sm<-av;cm<-av;cv<-av;cv1<-av;cv2<-av;cv3<-av;cv4<-av;vl<-av;vi<-av;cf<-av;aplt<-av #;p1<-av;p2<-av
  
  if (scope=='Global'){Clbr<-'isotonic'}
  window_size=0.1;step_size=0.05
  #Do cross-validation
  for (spc in c('AN','PL')){
    
    am0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_TransNorm_Mod.rds'))[[spc]]
    
    #For detecting growth (does not affect sample size, done here)
    if (scope=='Global'){
      av<-am0[,paste0('OffshrAq_',spc,'_F')]-am0[,paste0('OffshrAq_',spc)]>0
      am0[,paste0('OffshrAq___',spc,'_F')][av]<-'1'
      am0[,paste0('OffshrAq___',spc,'_F')][!av]<-'0'
    }
    
    set.seed(2000)
    print(spc)
    
    #Get folds
    rr<-data.frame(SptStr=sort(unique(am0$SptStr)),Fold=sample(1:N,length(sort(unique(am0$SptStr))),T)) #For ten-fold block cross-validation
    
    #Target variable
    tgt<-paste0('OffshrAq',mdltgt,spc,'_F')
      
    for (j in 1:nrow(rfp)){
        
      for (i in 1:N){
        
        #Get training and validation sets
        amt<-am0[!is.element(am0$SptStr,rr$SptStr[rr$Fold==i]),]
        amv<-am0[is.element(am0$SptStr,rr$SptStr[rr$Fold==i]),]
        if (N==1){amt<-am0;amv<-NULL}
        
        #Observation weights
        W<-GetWeights(amt,tgt)
        if (UseWeights==F){W<-NULL}
        rfave<-eval(rfave0)
        #rfave<-gsub(paste0('OffshrAq___',spc,' +'),'',rfave,fixed=T)
        rfave<-formula(rfave)
        numvar<-Inf
        
        if (ii=='factor'){classweights<-W}else{classweights<-NULL};if (ii=='numeric'){caseweights<-W}else{caseweights<-NULL};if (UseWeights==F){classweights<-NULL;caseweights<-NULL}
        
        #Feature selection
        if (VarSel==T){
          
          ftv1<-vector('list',VarSelItr);ftv2<-vector('list',VarSelItr)
          for (x in 1:VarSelItr){
            av<-tryCatch(FeatSel(spc,rfave,amt,tgt,rfp,j,caseweights,classweights),
                         error=function(e){c()})
            if (is.null(av)){next()}
            ftv1[[x]]<-av[[1]];ftv2[[x]]<-av[[2]]
            #View(ftv1b)
          }
          ftv1b<-data.frame(bind_rows(ftv1))
          ftv1b<-unlist(sapply(ftv1b,function(x){length(which(is.element(x,VarSelVal)))}))
          selvar<-ftv1b[ftv1b>=VarSelCtt]
          print(selvar)
          if (is.null(selvar)){next()}
          
          av<-ifelse(is.element(names(ftv1b),names(selvar)),'Confirmed','Rejected');names(av)<-names(ftv1b)
          vl[[spc]][[j]][[i]]<-av
          vi[[spc]][[j]][[i]]<-bind_rows(ftv2)
          numvar<-length(selvar)
          rfave<-paste0(tgt,'~',paste0(names(selvar),collapse='+'))
        }
        
        #Probability forest
        if (Prob=='Prob'){
          library(ranger);options(ranger.num.threads=8)
          if (rfp[j,2]>numvar){next()}
          if (ii=='factor'){classweights<-W}else{classweights<-NULL};if (ii=='numeric'){caseweights<-W}else{caseweights<-NULL};if (UseWeights==F){classweights<-NULL;caseweights<-NULL}
          sm[[spc]][[j]][[i]]<-tryCatch(
            ranger(rfave,amt,num.trees=rfp[j,1],mtry=rfp[j,2],min.node.size=rfp[j,3],probability=T,case.weights=caseweights,class.weights=classweights,keep.inbag=T), #importance='impurity_corrected'
            error=function(e){c()})
          if (length(sm[[spc]][[j]][[i]])==0){next()}
          
          ###Get accuracy metrics
          cv[[spc]][[j]][[i]]<-sm[[spc]][[j]][[i]]$prediction.error
          
          ###Calibrate probability values
          bnpr<-0.025 #Binning parameter
          prd<-vector('list',length(levels(amt[,tgt])));names(prd)<-levels(amt[,tgt])
          cma<-vector('list',length(levels(amt[,tgt])));names(cma)<-levels(amt[,tgt])
          
          #Prepare data for calibration
          #par(mfrow=c(length(levels(amt[,tgt])),1))
          for (pp in 1:length(prd)){
            
            #Get average OOB predictions
            ibcs <- sm[[spc]][[j]][[i]]$inbag.counts
            ibcs <- do.call(cbind, ibcs)
            preds <- predict(sm[[spc]][[j]][[i]],amt,predict.all=T)$predictions[,pp,]
            preds[which(ibcs > 0)] <- NA
            preds <- rowMeans(preds,na.rm=T)
            #all.equal(sm[[spc]][[j]][[i]]$predictions[,pp],preds)
            prd[[pp]]<-preds
            
          }
          #par(mfrow=c(1,1))
          ac<-paste0('P',1:length(levels(amt[,tgt])))
          am<-data.frame(prd,amt[,tgt]);names(am)<-c(ac,'class')
          
          library(probably)
          if (Clbr=='isotonic'){ #Fit binary calibration model
            cmg<-tryCatch(
              cal_estimate_isotonic(am,class,estimate=all_of(ac)),
              error=function(e){NULL}
            )
          }
          if (Clbr=='beta'){ #Fit binary calibration model
            cmg<-tryCatch(
              cal_estimate_beta(am,class,estimate=all_of(ac)),
              error=function(e){NULL}
            )
          }
          if (Clbr=='multinomial'){ #Fit multiple calibration model
            cmg<-tryCatch(
              cal_estimate_multinomial(am,class,all_of(ac),smooth=T),
              error=function(e){NULL}
            )
            if (is.null(cmg)){
              cmg<-cal_estimate_multinomial(am,class,all_of(ac),smooth=F)
            }
          }
          cg<-cal_apply(am,cmg)
          cm[[spc]][[j]][[i]]<-cmg #Save main model
          if (all(amt[,tgt]==am$class)==F){stop()}
          aplt[[spc]][[j]][[i]]<-cal_plot_windowed(cg,truth=class,estimate=all_of(ac),window_size=window_size,step_size=step_size,include_rug=F,include_ribbon=T)
          
          ###Get accuracy metrics (assuming conditional mean from calibration, no dispersion)
          NNNN<-10
          aaa<-vector('list',NNNN)
          aslm<-sample(1:NNNN,nrow(am),T)
          for (m in 1:NNNN){
            amdl<-cal_estimate_isotonic(am[aslm!=m,],class,estimate=all_of(ac))
            aaa[[m]]<-cal_apply(am[aslm==m,],amdl)
          }
          aaa<-bind_rows(aaa)
          av<-length(which(!is.finite(aaa[,'P1'])));if (av>0){print(c('NaN: ',av))}
          av<-which(is.finite(aaa[,'P1'])) #Select finite rows for calculation, just omit NaN
          if (length(levels(amt[,tgt]))==2){xx<-aaa[av,paste0('P',1:2)]}
          if (length(levels(amt[,tgt]))>2){xx<-aaa[av,paste0('P',1:3)]}
          cv1[[spc]][[j]][[i]]<-ProbScore(xx,amt[av,tgt],levels(amt[,tgt]),'Brier')
          cv3[[spc]][[j]][[i]]<-SkillScore(xx,amt[av,tgt],levels(amt[,tgt]),'Brier')
          cv2[[spc]][[j]][[i]]<-ProbScore(xx,amt[av,tgt],levels(amt[,tgt]),'LogLoss')
          cv4[[spc]][[j]][[i]]<-SkillScore(xx,amt[av,tgt],levels(amt[,tgt]),'LogLoss')
          MLmetrics::MultiLogLoss(data.frame(`0`=xx[,1],`1`=xx[,2]),amt[av,tgt])
          
        }
        
      }
      
      print(row.names(rfp)[j])
      #print(selvar)
      print(sm[[spc]][[j]][[i]]$prediction.error)
      
      
    }
    
    saveRDS(sm,paste0('F:/',scope,'_data/Code files/Modeling/sm_',ac0,'_',Prob,'.rds'))
    saveRDS(cm,paste0('F:/',scope,'_data/Code files/Modeling/cm_',ac0,'_',Prob,'.rds'))
    saveRDS(aplt,paste0('F:/',scope,'_data/Code files/Modeling/aplt_',ac0,'_',Prob,'.rds'))
    saveRDS(cv,paste0('F:/',scope,'_data/Code files/Modeling/cv_',ac0,'_',Prob,'.rds'))
    saveRDS(cv1,paste0('F:/',scope,'_data/Code files/Modeling/cv1_',ac0,'_',Prob,'.rds'))
    saveRDS(cv2,paste0('F:/',scope,'_data/Code files/Modeling/cv2_',ac0,'_',Prob,'.rds'))
    saveRDS(cv3,paste0('F:/',scope,'_data/Code files/Modeling/cv3_',ac0,'_',Prob,'.rds'))
    saveRDS(cv4,paste0('F:/',scope,'_data/Code files/Modeling/cv4_',ac0,'_',Prob,'.rds'))
    saveRDS(vl,paste0('F:/',scope,'_data/Code files/Modeling/vl_',ac0,'_',Prob,'.rds'))
    saveRDS(vi,paste0('F:/',scope,'_data/Code files/Modeling/vi_',ac0,'_',Prob,'.rds'))
    
  }
  
  #Check parameter results
  if (scope=='Global'){k<-1} #Brier
  {
    cv<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/cv',k,'_',ac0,'_',Prob,'.rds'))
    
    MI<-vector('numeric',2);names(MI)<-c('AN','PL')
    MV<-vector('character',2);names(MV)<-c('AN','PL')
    av<-expand.grid(Spc=c('AN','PL'),Par=1:nrow(rfp),Val=nmr(NA))
    for (spc in c('AN','PL')){
      for (j in 1:nrow(rfp)){
        av$Val[av$Spc==spc & av$Par==j]<-mean(unlist(cv[[spc]][[j]]))
      }
    }
    for (spc in c('AN','PL')){
      MI[spc]<-av$Par[av$Spc==spc][which.min(av$Val[av$Spc==spc])]
      if (k==3){MI[spc]<-av$Par[av$Spc==spc][which.max(av$Val[av$Spc==spc])]}
    }
    for (spc in c('AN','PL')){
      MV[spc]<-row.names(rfp)[MI[spc]]
    }
    
    print(MI)
    print(MV)
    print(c(av$Val[av$Spc=='AN'][MI['AN']],av$Val[av$Spc=='PL'][MI['PL']]))
    if (k==3){print(c(av$Val[av$Spc=='AN'][which.max(av$Val[av$Spc=='AN'])],av$Val[av$Spc=='PL'][which.max(av$Val[av$Spc=='PL'])]))}
      
    saveRDS(MV,paste0('F:/',scope,'_data/Code files/Modeling/SelectedCV_',ac0,'_',Prob,'.rds'))
  }
  #Error
  {
    MV<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/SelectedCV_',ac0,'_',Prob,'.rds'))
    MI<-vector('numeric',2);names(MI)<-c('AN','PL');for (spc in c('AN','PL')){MI[[spc]]<-which(row.names(rfp)==MV[spc])}
    cv1<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/cv1_',ac0,'_',Prob,'.rds'))
    cv2<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/cv2_',ac0,'_',Prob,'.rds'))
    cv3<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/cv3_',ac0,'_',Prob,'.rds'))
    cv4<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/cv4_',ac0,'_',Prob,'.rds'))
    dat_text <- data.frame(
      # x=rep(0.6,4),
      # y=rep(0,4),
      label=c(paste0('Brier Score: ',round(cv1[['AN']][[MI['AN']]][[1]],4)),
              paste0('Brier Skill: ',round(cv3[['AN']][[MI['AN']]][[1]]*100,1),'%'),
              paste0('LogLoss Score: ',round(cv2[['AN']][[MI['AN']]][[1]],4)),
              paste0('LogLoss Skill: ',round(cv4[['AN']][[MI['AN']]][[1]]*100,1),'%'),
              paste0('Brier Score: ',round(cv1[['PL']][[MI['PL']]][[1]],4)),
              paste0('Brier Skill: ',round(cv3[['PL']][[MI['PL']]][[1]]*100,1),'%'),
              paste0('LogLoss Score: ',round(cv2[['PL']][[MI['PL']]][[1]],4)),
              paste0('LogLoss Skill: ',round(cv4[['PL']][[MI['PL']]][[1]]*100,1),'%')
      ),
      Spc=c('Fish','Fish','Seaweed','Seaweed')
    )
    print(dat_text)
  }
  
  library(ranger)
  library(mgcv)
  library(glmnet)
  library(probably)
  
  #Reliability diagram
  {
    sm0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/sm_',ac0,'_',Prob,'.rds'))
    cm0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/cm_',ac0,'_',Prob,'.rds'))
    aplt0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/aplt_',ac0,'_',Prob,'.rds'))
    #aplt0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/_',ac0,'_',Prob,'.rds'))
    MV<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/SelectedCV_',ac0,'_',Prob,'.rds'))
    MI<-vector('numeric',2);names(MI)<-c('AN','PL');for (spc in c('AN','PL')){MI[[spc]]<-which(row.names(rfp)==MV[spc])}
    # cv1<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/cv1_',ac0,'_',Prob,'.rds'))
    # cv3<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/cv3_',ac0,'_',Prob,'.rds'))
    
    #Figure
    library(ggplot2)
    library(RColorBrewer) #RColorBrewer::display.brewer.all() 
    library(cowplot)
    library(extrafont)
    loadfonts(device="win",quiet=T)
    #font_import(paths="C:/WINDOWS/FONTS",pattern='Times')
    #fonts()
    
    prv<-vector('list',2);names(prv)<-c('AN','PL')
    for (spc in names(prv)){
      
      am0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_TransNorm_Mod.rds'))[[spc]]
      #For detecting growth (does no affect sample size, done here)
      if (scope=='Global'){
        av<-am0[,paste0('OffshrAq_',spc,'_F')]-am0[,paste0('OffshrAq_',spc)]>0
        am0[,paste0('OffshrAq___',spc,'_F')][av]<-'1'
        am0[,paste0('OffshrAq___',spc,'_F')][!av]<-'0'
      }
      #For detecting change type (does no affect sample size, done here)
      if (scope=='China'){
        # av<-mapply(function(x,y){paste0(x,y,collapse='-')},x=am0[,paste0('OffshrAq___',spc)],y=am0[,paste0('OffshrAq___',spc,'_F')]) #table(av)
        # am0[,paste0('OffshrAq___',spc,'_F')]<-fct(av)
        av<-am0[,paste0('OffshrAq_',spc,'_F')]-am0[,paste0('OffshrAq_',spc)]
        am0[,paste0('OffshrAq___',spc,'_F')]<-chr(am0[,paste0('OffshrAq___',spc,'_F')])
        am0[,paste0('OffshrAq___',spc,'_F')][av>0]<-'1'
        am0[,paste0('OffshrAq___',spc,'_F')][av==0]<-'0'
        am0[,paste0('OffshrAq___',spc,'_F')][av<0]<-'-1'
        am0[,paste0('OffshrAq___',spc,'_F')]<-fct(am0[,paste0('OffshrAq___',spc,'_F')])
      }
      
      set.seed(2000)
      
      #Get folds
      amt<-am0
      
      #Target variable
      tgt<-paste0('OffshrAq',mdltgt,spc,'_F')
      
      ac<-1:length(levels(amt[,tgt]))
      preds1<-data.frame(predict(sm0[[spc]][[MI[[spc]]]][[1]],amt)$predictions);names(preds1)<-paste0('P',ac)
      preds2<-cal_apply(preds1,cm0[[spc]][[MI[[spc]]]][[1]])
      prv[[spc]]<-data.frame(preds1,PrdsCal=preds2,Spc=spc);names(prv[[spc]])<-c(paste0('PrdUnc',ac),paste0('PrdCal',ac),'Spc')
      
    }
    prv<-bind_rows(prv)
    if (scope=='Global'){
      prv<-bind_rows(
        data.frame(PrdUnc=prv[,1],PrdCal=prv[,1+2],Spc=prv[,'Spc'],class='0'),
        data.frame(PrdUnc=prv[,2],PrdCal=prv[,2+2],Spc=prv[,'Spc'],class='1')
      )
    }
    #
    spc<-'AN'
    aa<-aplt0[[spc]][[MI[spc]]][[1]]
    spc<-'PL'
    bb<-aplt0[[spc]][[MI[spc]]][[1]]
    am<-rbind(data.frame(aa$data,Spc='Fish'),data.frame(bb$data,Spc='Seaweed'))
    #Add other class for two-class case
    if (scope=='China'){
      am$Lbl<-ifelse(am$class=='0','Null\ngrowth',ifelse(am$class=='-1','Strictly\nnegative\ngrowth','Strictly\npositive\ngrowth'))
      am$Lbl<-factor(am$Lbl,levels=c('Strictly\nnegative\ngrowth','Null\ngrowth','Strictly\npositive\ngrowth'))
      
      prv$Spc<-ifelse(prv$Spc=='AN','Fish','Seaweed')
      prv$Lbl<-ifelse(prv$class=='0','Null\ngrowth',ifelse(prv$class=='-1','Strictly\nnegative\ngrowth','Strictly\npositive\ngrowth'))
      prv$Lbl<-factor(prv$Lbl,levels=c('Strictly\nnegative\ngrowth','Null\ngrowth','Strictly\npositive\ngrowth'))
    }
    if (scope=='Global'){
      am$class<-'1'
      am2<-am;am2$class<-'0'
      am2$event_rate<-1-am2$event_rate
      am2$predicted_midpoint<-1-am2$predicted_midpoint
      am2$lower<-1-am2$lower
      am2$upper<-1-am2$upper
      am<-rbind(am,am2)
      am$Lbl<-ifelse(am$class=='0','Null\ngrowth','Strictly\npositive\ngrowth')
      am$Lbl<-factor(am$Lbl,levels=c('Null\ngrowth','Strictly\npositive\ngrowth'))
      
      prv$Spc<-ifelse(prv$Spc=='AN','Fish','Seaweed')
      prv$Lbl<-ifelse(prv$class=='0','Null\ngrowth','Strictly\npositive\ngrowth')
      prv$Lbl<-factor(prv$Lbl,levels=c('Null\ngrowth','Strictly\npositive\ngrowth'))
    }
    
    hght<-15/2*length(levels(amt[,tgt]))
    aplt<-ggplot(am,aes(x=predicted_midpoint,y=event_rate,color=Spc,fill=Spc))+
      geom_point(pch=16,size=2,alpha=1)+ #position=ggforce::position_jitternormal(0.025,0.025),
      geom_line(linewidth=1.2,linetype=1)+
      geom_abline(intercept=0,slope=1,linewidth=1,linetype=3)+
      #geom_point(aes(x=PrdUnc,y=PrdCal,color=Spc,fill=Spc),data=prv,pch=16,size=2,alpha=0.1,position=ggforce::position_jitternormal(0.025,0.025),)+
      geom_rug(aes(x=PrdUnc,y=PrdCal,color=Spc,alpha=1),data=prv,sides="b")+
      #General
      theme(text=element_text(size=16,family="Times New Roman",color='black'))+
      theme(panel.background = element_rect(fill="White",colour="black"),
            panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
      theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
            panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
      #Legend
      theme(legend.position="none",legend.title=element_blank(),
            legend.text=element_text(size=16,color='black'),
            legend.key = element_rect(color=NA,fill=NA),
            #legend.key.height=unit(0.75,"cm"),
            #legend.key.size=unit(1,"cm"),
            #legend.key.width=unit(1,"cm"),
            legend.spacing.y=unit(0.4,'cm'))+
      #guides(guide=guide_legend(reverse=T))+
      #Axis
      #Margins from top, clockwise
      xlab(expression("Predicted probability"))+
      ylab(expression("Actual probability"))+
      theme(axis.title.x=element_text(size=16,margin=ggplot2::margin(t=10,r=0,b=0,l=0)),
            axis.title.y=element_text(size=16,margin=ggplot2::margin(t=0,r=20,b=0,l=0)))+
      scale_x_continuous(breaks=seq(0,1,0.2),minor_breaks=seq(0,1,0.2),expand=expansion(add=c(0.1,0.1),mult=c(0,0)))+ #,position='bottom',
      scale_y_continuous(breaks=seq(0,1,0.2),minor_breaks=seq(0,1,0.2),expand=expansion(add=c(0.1,0.1),mult=c(0,0)))+ #,position='bottom',
      #scale_x_discrete(limits=rev)+
      theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+ 
      theme(axis.text.x=element_text(size=14,color="black"))+
      theme(axis.text.y=element_text(size=14,color="black"))+
      theme(axis.ticks.length=unit(.1,"cm"))+
      #Others
      scale_color_manual(values=brewer.pal(9,'Set1')[c(5,2)])+
      #geom_text(data=dat_text,aes(x=x,y=y,label=label,hjust=0),color='black',size=5,family="Times New Roman")+
      facet_grid(rows=vars(Lbl),cols=vars(Spc),scales="free",space='free')+
      theme(strip.background=element_rect(fill="NA"),
            strip.text.x=element_text(size=16,color="black",angle=0),
            strip.text.y=element_text(size=16,color="black",angle=0))
    pp<-aplt #plot_grid(plotlist=aplt,nrow=1,ncol=1)
    #aplt+theme(legend.position="top") #Check legend
    # ggsave(paste0("F:/Global_data/Code files/Figures/Jiji.jpg"),aplt,
    #        device='jpg',width=20,height=16,units=c("cm"),dpi=1000,limitsize=F)
    jpeg(paste0("F:/",scope,"_data/Code files/Figures/ReliabilityDiagram.jpg"),units="cm",width=30,height=hght,res=800);plot(pp);dev.off()
    jpeg(paste0("F:/",scope,"_data/Code files/Figures/ReliabilityDiagram_Legend.jpg"),units="cm",width=30,height=hght,res=800);plot(aplt+theme(legend.position="top"));dev.off()
    
  }
  
  #Get selected variables
  varcutoff<-6
  if (Prob=='Prob'){varcutoff<-1}
  if (VarSel==T){
    vl<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/vl_',ac0,'_',Prob,'.rds'))
    vrblnm<-vector('list',2);names(vrblnm)<-c('AN','PL')
    for (spc in c('AN','PL')){
      av<-as.data.frame(vl[[spc]][[MI[[spc]]]])
      #vl<-apply(vl,1,function(x){length(which(x=='Confirmed'))})
      av<-apply(av,1,function(x){length(which(x!='Rejected'))})
      print(av)
      vrblnm[[spc]]<-names(av)[av>=varcutoff]
    }
    
    print(vrblnm)
    
    saveRDS(vrblnm,paste0('F:/',scope,'_data/Code files/Modeling/SelectedVariables_',ac0,'_',Prob,'.rds'))
  }
  #Get variable importance (including selection cut-off)
  if (VarSel==T){
    vi<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/vi_',ac0,'_',Prob,'.rds'))
    vrblnm<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/SelectedVariables_',ac0,'_',Prob,'.rds'))
    vrblim<-vector('list',2);names(vrblim)<-c('AN','PL')
    for (spc in c('AN','PL')){
      vrblim[[spc]]<-dplyr::bind_rows(vi[[spc]][[MI[[spc]]]]) #as.data.frame(vi[[spc]][[MI[[spc]]]])
      vrblim[[spc]]$Imp<-'No'
      vrblim[[spc]]$Imp[is.element(vrblim[[spc]]$Var,vrblnm[[spc]])]<-'Yes'
    }
    am<-dplyr::bind_rows(vrblim);names(am)<-c('Var','Val','Spc','Imp')
    
    library(ggplot2)
    library(RColorBrewer) #RColorBrewer::display.brewer.all() 
    library(cowplot)
    library(extrafont)
    loadfonts(device="win",quiet=T)
    #font_import(paths="C:/WINDOWS/FONTS",pattern='Times')
    #fonts()
    
    df<-am
    au<-expand.grid(Spc=c('AN','PL'),Var=unique(df$Var))
    for (p in c('AN','PL')){for (q in unique(au$Var)){
      for (r in seq(0,1,0.05)){
        av<-df[df$Spc==p & df$Var==q,'Val']
        au[au$Spc==p & au$Var==q,paste0('q',r)]<-quantile(av[is.finite(av)],r)
      }
      au[au$Spc==p & au$Var==q,'Mean']<-mean(av[is.finite(av)])
    }}
    au$Spc<-factor(chr(au$Spc),levels=c('PL','AN'))
    av<-data.frame(Spc=chr(df[,'Spc']),Var=chr(df[,'Var']),Imp=chr(df[,'Imp']));av<-av[!duplicated(av),]
    au$Imp<-left_join(data.frame(Spc=chr(au[,c('Spc')]),Var=chr(au[,c('Var')])),av,join_by(Var,Spc))[,3]
    au$Imp<-factor(chr(au$Imp),levels=c('Yes','No'))
    au$Var<-factor(chr(au$Var),levels=rev(levels(au$Var)))
    amfi<-au
    
    aplt<-ggplot(au[au$Imp=='Yes',])+
      # geom_boxplot(aes(xlower=q0.25,xmiddle=q0.5,xupper=q0.75,xmin=q0,xmax=q1,y=Var,color=Spc),stat="identity",staplewidth=0.8,linewidth=1.2)+
      #geom_boxplot(aes(lower=q0.25,middle=q0.5,upper=q0.75,min=q0,max=q1,x=Var,color=Spc),stat="identity",position=position_dodge(width=1.2))+
      geom_point(aes(y=q0.5,x=Var,color=Spc),size=4,shape=16,stroke=0,alpha=1.0,position=position_dodge(width=1.4))+
      geom_boxplot(aes(lower=q0.25,middle=q0.5,upper=q0.75,min=q0.05,max=q0.95,x=Var,color=Spc),stat="identity",staplewidth=0,width=0,linewidth=0.8,alpha=1.0,position=position_dodge(width=1.4))+
      coord_flip()+
      #General
      theme(text=element_text(size=16,family="Times New Roman",color='black'))+
      theme(panel.background = element_rect(fill="White",colour="black"),
            panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
      theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
            panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
      #Legend
      theme(legend.position="right",legend.title=element_blank(),
            legend.text=element_text(size=16,color='black'),
            legend.key = element_rect(color=NA,fill=NA),
            #legend.key.height=unit(0.75,"cm"),
            #legend.key.size=unit(1,"cm"),
            #legend.key.width=unit(1,"cm"),
            legend.spacing.y=unit(0.4,'cm'))+
      #guides(guide=guide_legend(reverse=T))+
      #Axis
      #Margins from top, clockwise
      xlab(expression("Variable"))+
      ylab(expression("Z-Score"))+
      theme(axis.title.x=element_text(size=16,margin=ggplot2::margin(t=10,r=0,b=0,l=0)),
            axis.title.y=element_text(size=16,margin=ggplot2::margin(t=0,r=20,b=0,l=0)))+
      #scale_x_continuous(limits=c(-5,35),breaks=seq(-5,35,5),minor_breaks=seq(-5,35,2.5),expand=expansion(add=c(0,0),mult=c(0,0)))+ #,position='bottom',
      #scale_x_discrete(limits=rev)+
      scale_y_log10()+
      theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+ 
      theme(axis.text.x=element_text(size=14,color="black"))+
      theme(axis.text.y=element_text(size=14,color="black"))+
      #theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
      theme(axis.ticks.length=unit(.1,"cm"))+
      #Others
      scale_color_manual(values=brewer.pal(9,'Set1')[c(2,5)])+
      scale_fill_discrete(breaks=c("r","q"))+
      #facet_grid(rows=vars(Imp),scales="free",space='free')+
      theme(strip.background=element_rect(fill="NA"),
            strip.text.y=element_text(size=12,color="black",angle=0))
    pp<-aplt
    jpeg(paste0("E:/Global_data/Figures/FeatImp.jpg"),units="cm",width=20,height=20,res=800);pp;dev.off()
    
    # an<-which(sapply(am0,class)!='factor')
    # par(mfrow=c(1,1))
    # corrplot::corrplot(as.matrix(cor(am0[,an[4:length(an)]],method='spearman')),method="number",is.corr=T) #method='color'
    
  }
  #Only implemented for N=1, variable selection not part of fitting procedure
  if (VarSel==F){
    
    MV<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/SelectedCV_',ac0,'_',Prob,'.rds'))
    MI<-vector('numeric',2);names(MI)<-c('AN','PL');for (spc in c('AN','PL')){MI[[spc]]<-which(row.names(rfp)==MV[spc])}
    
    vl<-vector('list',2);names(vl)<-c('AN','PL')
    vi<-vector('list',2);names(vi)<-c('AN','PL')
    
    for (spc in c('AN','PL')){
      
      am0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_TransNorm_Mod.rds'))[[spc]]
      
      #For detecting growth (does no affect sample size, done here)
      if (scope=='Global'){
        av<-am0[,paste0('OffshrAq_',spc,'_F')]-am0[,paste0('OffshrAq_',spc)]>0
        am0[,paste0('OffshrAq___',spc,'_F')][av]<-'1'
        am0[,paste0('OffshrAq___',spc,'_F')][!av]<-'0'
      }
      #For detecting change type (does no affect sample size, done here)
      if (scope=='China'){
        # av<-mapply(function(x,y){paste0(x,y,collapse='-')},x=am0[,paste0('OffshrAq___',spc)],y=am0[,paste0('OffshrAq___',spc,'_F')]) #table(av)
        # am0[,paste0('OffshrAq___',spc,'_F')]<-fct(av)
        av<-am0[,paste0('OffshrAq_',spc,'_F')]-am0[,paste0('OffshrAq_',spc)]
        am0[,paste0('OffshrAq___',spc,'_F')]<-chr(am0[,paste0('OffshrAq___',spc,'_F')])
        am0[,paste0('OffshrAq___',spc,'_F')][av>0]<-'1'
        am0[,paste0('OffshrAq___',spc,'_F')][av==0]<-'0'
        am0[,paste0('OffshrAq___',spc,'_F')][av<0]<-'-1'
        am0[,paste0('OffshrAq___',spc,'_F')]<-fct(am0[,paste0('OffshrAq___',spc,'_F')])
      }
      
      set.seed(2000)
      print(spc)
      
      #Get folds
      rr<-data.frame(SptStr=sort(unique(am0$SptStr)),Fold=sample(1:N,length(sort(unique(am0$SptStr))),T)) #For ten-fold block cross-validation
      
      #Target variable
      tgt<-paste0('OffshrAq',mdltgt,spc,'_F')
      
      j<-MI[[spc]]
      
      i<-1
      
      #Get training and validation sets
      amt<-am0[!is.element(am0$SptStr,rr$SptStr[rr$Fold==i]),]
      amv<-am0[is.element(am0$SptStr,rr$SptStr[rr$Fold==i]),]
      if (N==1){amt<-am0;amv<-NULL}
      
      #Observation weights
      W<-GetWeights(amt,tgt)
      if (UseWeights==F){W<-NULL}
      rfave<-eval(rfave0)
      #rfave<-gsub(paste0('OffshrAq___',spc,' +'),'',rfave,fixed=T)
      rfave<-formula(rfave)
      numvar<-Inf
      
      if (ii=='factor'){classweights<-W}else{classweights<-NULL};if (ii=='numeric'){caseweights<-W}else{caseweights<-NULL};if (UseWeights==F){classweights<-NULL;caseweights<-NULL}
      
      #Feature selection
      #https://datascience.stackexchange.com/questions/56672/unimportant-features-impact-on-models-performance
      ftv1<-vector('list',VarSelItr);ftv2<-vector('list',VarSelItr)
      for (x in 1:VarSelItr){
        av<-tryCatch(FeatSel(spc,rfave,amt,tgt,rfp,j,caseweights,classweights),
                     error=function(e){c()})
        if (is.null(av)){next()}
        ftv1[[x]]<-av[[1]];ftv2[[x]]<-av[[2]]
        #View(ftv1b)
      }
      ftv1b<-data.frame(bind_rows(ftv1))
      ftv1b<-unlist(sapply(ftv1b,function(x){length(which(is.element(x,VarSelVal)))}))
      selvar<-ftv1b[ftv1b>=VarSelCtt]
      print(selvar)
      if (is.null(selvar)){next()}
      
      av<-ifelse(is.element(names(ftv1b),names(selvar)),'Confirmed','Rejected');names(av)<-names(ftv1b)
      vl[[spc]]<-av
      vi[[spc]]<-bind_rows(ftv2)
      numvar<-length(selvar)
      rfave<-paste0(tgt,'~',paste0(names(selvar),collapse='+'))
      
    }
    
    saveRDS(vl,paste0('F:/',scope,'_data/Code files/Modeling/vl_',ac0,'_',Prob,'.rds'))
    saveRDS(vi,paste0('F:/',scope,'_data/Code files/Modeling/vi_',ac0,'_',Prob,'.rds'))
    
  } #Data
  if (VarSel==F){
    vi<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/vi_',ac0,'_',Prob,'.rds'))
    ac<-strsplit(chr(formula(eval(rfave0)))[3],' \\+ ',perl=T)[[1]]
    ac<-gsub('\n','',ac,fixed=T)
    va<-gsub(' ','',ac,fixed=T)
    vrblnm<-va#readRDS(paste0('F:/',scope,'_data/Code files/Modeling/SelectedVariables_',ac0,'_',Prob,'.rds'))
    vrblim<-vector('list',2);names(vrblim)<-c('AN','PL')
    for (spc in c('AN','PL')){
      vrblim[[spc]]<-dplyr::bind_rows(vi[[spc]]) #as.data.frame(vi[[spc]][[MI[[spc]]]])
      vrblim[[spc]]$Imp<-'Yes'
    }
    am<-dplyr::bind_rows(vrblim);names(am)<-c('Var','Val','Spc','Imp')
    
    library(ggplot2)
    library(RColorBrewer) #RColorBrewer::display.brewer.all() 
    library(colorspace)
    library(cowplot)
    library(extrafont)
    loadfonts(device="win",quiet=T)
    #font_import(paths="C:/WINDOWS/FONTS",pattern='Times')
    #fonts()
    
    df<-am
    au<-expand.grid(Spc=c('AN','PL'),Var=unique(df$Var))
    for (p in c('AN','PL')){for (q in unique(au$Var)){
      for (r in seq(0,1,0.05)){
        av<-df[df$Spc==p & df$Var==q,'Val']
        au[au$Spc==p & au$Var==q,paste0('q',r)]<-quantile(av[is.finite(av)],r)
      }
      au[au$Spc==p & au$Var==q,'Mean']<-mean(av[is.finite(av)])
    }}
    au$Spc<-factor(chr(au$Spc),levels=c('PL','AN'))
    av<-data.frame(Spc=chr(df[,'Spc']),Var=chr(df[,'Var']),Imp=chr(df[,'Imp']));av<-av[!duplicated(av),]
    au$Imp<-left_join(data.frame(Spc=chr(au[,c('Spc')]),Var=chr(au[,c('Var')])),av,join_by(Var,Spc))[,3]
    au$Imp<-factor(chr(au$Imp),levels=c('Yes','No'))
    au$Var<-factor(chr(au$Var),levels=rev(levels(au$Var)))
    #ao<-data.frame(tidyr::pivot_longer(au,names(au)[3:7],names_to='Stt',values_to='Val'))
    # ao<-vector('list',1000)
    # for (n in 1:nrow(au)){
    #   a<-au$q0.25[n]
    #   b<-au$q0.75[n]
    #   ao[[n]]<-data.frame(Spc=au$Spc[n],Var=au$Var[n],Val=seq(a,b,(b-a)/1000),Imp=au$Imp[n])
    #   n<-n+1
    # }
    # ao<-bind_rows(ao)
    # ggplot()+
    #   geom_point(data=au[au$Imp=='Yes',],aes(y=q0.5,x=Var,color=Spc),size=4,shape=16,stroke=0,alpha=1.0,position=position_dodge2(width=1))+
    #   geom_point(data=ao[ao$Imp=='Yes',],aes(y=Val,x=Var,color=Spc),size=1,shape=16,stroke=0,alpha=1.0,position=position_dodge(width=1))+
    au<-AddVarCategories(au,scope)
    au$Spc<-factor(ifelse(au$Spc=='AN','Fish','Seaweed'),levels=c('Seaweed','Fish'))
    au$Lbl<-NA
    au$Lbl[au$Spc=='Fish']<-au$q0.5[au$Spc=='Fish']/max(au$q0.5[au$Spc=='Fish'])*100
    au$Lbl[au$Spc=='Seaweed']<-au$q0.5[au$Spc=='Seaweed']/max(au$q0.5[au$Spc=='Seaweed'])*100
    au$Lbl[au$Lbl<0]<-0
    au<-ChangeVarNames(au,scope)
    amfi<-au
    
    # #http://colorspace.r-forge.r-project.org/articles/manipulation_utilities.html
    # clr0<-RColorBrewer::brewer.pal(8,'Dark2')[c(1,2,3)]
    # clr1<-colorspace::lighten(clr0,0.15,"absolute")
    # clr2<-colorspace::lighten(clr1,0.15,"absolute")
    # scales::show_col(c(clr0,clr1,clr2))
    # #clr<-c(clr0,clr1,clr2)
    # clr<-c(t(matrix(c(clr2,clr1,clr0),3,3)))
    #scales::show_col(clr0)
    #library(colorspace)
    #hcl_palettes()
    #divergingx_palettes(n = 7, plot = TRUE)
    #scale_color_continuous_sequential(palette="Oranges")
    #https://stackoverflow.com/questions/78926052/make-color-blocks-with-text-using-ggplot2
    
    #Additional changes to details
    av<-chr(au$Var);lvl<-levels(au$Var)
    av<-gsub('(diff)','(diff.)',av,fixed=T);lvl<-gsub('(diff)','(diff.)',lvl,fixed=T)
    av<-gsub('Sea temperature','Seawater temperature',av,fixed=T);lvl<-gsub('Sea temperature','Seawater temperature',lvl,fixed=T)
    av<-gsub('Chlorophyll-A','Chlorophyll-a',av,fixed=T);lvl<-gsub('Chlorophyll-A','Chlorophyll-a',lvl,fixed=T)
    av<-gsub('Elevation/depth','Sea depth',av,fixed=T);lvl<-gsub('Elevation/depth','Sea depth',lvl,fixed=T)
    lvl<-lvl[c(1:19,22,20,21,23:40)]
    au$Var<-factor(av,levels=lvl)
    
    aplt<-ggplot(au[!is.na(au$Cat),])+
      geom_boxplot(aes(xlower=q0.25,xmiddle=q0.5,xupper=q0.75,xmin=q0,xmax=q1,y=Var,color=Spc),stat="identity",staplewidth=0.8,linewidth=1)+
      #geom_boxplot(aes(lower=q0.25,middle=q0.5,upper=q0.75,min=q0,max=q1,x=Var,color=Spc),stat="identity",position=position_dodge(width=1.2))+
      #geom_point(aes(y=0,x=Var,color=round(Lbl,0)),size=10,shape=15,stroke=0,alpha=1.0)+
      #geom_text(aes(y=0,x=Var,label=round(Lbl,0),family="Times New Roman"),color=ifelse(round(au$Lbl[!is.na(au$Cat)],0)>=50,'white','black'))+
      #geom_boxplot(aes(lower=q0.25,middle=q0.5,upper=q0.75,min=q0.05,max=q0.95,x=Var,color=Spc),stat="identity",staplewidth=0,width=0,linewidth=0.8,alpha=1.0,position=position_dodge(width=1.4))+
      #coord_flip()+
      #General
      theme(text=element_text(size=14,family="Times New Roman",color='black'))+
      theme(panel.background = element_rect(fill="White",colour="black"),
            panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
      theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
            panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
      #Legend
      theme(legend.position="bottom",legend.title=element_blank(),
            legend.text=element_text(size=14,color='black'),
            legend.key = element_rect(color=NA,fill=NA),
            #legend.key.height=unit(0.75,"cm"),
            #legend.key.size=unit(1,"cm"),
            #legend.key.width=unit(1,"cm"),
            legend.spacing.y=unit(0.4,'cm'))+
      #guides(guide=guide_legend(reverse=T))+
      #Axis
      #Margins from top, clockwise
      xlab(expression("Z-score"))+
      ylab(expression("Covariate"))+
      theme(axis.title.x=element_text(size=14,margin=ggplot2::margin(t=10,r=0,b=0,l=0)),
            axis.title.y=element_text(size=14,margin=ggplot2::margin(t=0,r=20,b=0,l=0)))+
      #scale_x_continuous(limits=c(-5,35),breaks=seq(-5,35,5),minor_breaks=seq(-5,35,2.5),expand=expansion(add=c(0,0),mult=c(0,0)))+ #,position='bottom',
      scale_y_discrete(limits=rev)+
      #scale_y_log10()+
      theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+ 
      theme(axis.text.x=element_text(size=12,color="black"))+
      theme(axis.text.y=element_text(size=12,color="black"))+
      theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
      #Others
      scale_color_manual(values=brewer.pal(9,'Set1')[c(2,5)])+
      guides(color=guide_legend(reverse=T))+
      facet_grid(rows=vars(Cat),scales="free",space='free')+ #cols=vars(Spc)
      theme(strip.background=element_rect(fill="NA"),
            strip.text.x=element_text(size=14,color="black",angle=45),
            strip.text.y=element_text(size=14,color="black",angle=0)) #element_text(size=12,color="black",angle=0)
    pp<-aplt
    if (scope=='China'){hght<-22}
    if (scope=='Global'){hght<-26}
    jpeg(paste0("F:/",scope,"_data/Code files/Figures/FeatureImportance_v2.jpg"),units="cm",width=26,height=hght,res=800);plot(pp);dev.off()
  
  } #Plot
  
  #Get partial dependency plots
  {
    set.seed(4000)
    
    if (scope=='China'){fx<-mean;lbl<-'Mean';R<-c('-1','0','1')}
    if (scope=='Global'){fx<-mean;lbl<-'Mean';R<-c('0','1')}
    
    #Data
    {
      library(dplyr)
      sm0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/sm_',ac0,'_',Prob,'.rds'))
      cm0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/cm_',ac0,'_',Prob,'.rds'))
      MV<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/SelectedCV_',ac0,'_',Prob,'.rds'))
      MI<-vector('numeric',2);names(MI)<-c('AN','PL');for (spc in c('AN','PL')){MI[[spc]]<-which(row.names(rfp)==MV[spc])}
      #vrblnm<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/SelectedVariables_',ac0,'_',Prob,'.rds'))
      library(ranger)
      library(probably)
      af<-function(sm,am,vr,fx,cm,R){
        
        #Based on paper "pdp: An R Package for Constructing". Verified with examples from documentation.
        if (!is.null(cm)){
          y<-sapply(am[,vr],function(x){
            d<-am;d[,vr]<-x
            #a<-fx(predict(sm,d)$predictions)
            #b<-predict(cm,data.frame(x=a),type='response')
            a<-predict(sm,d)$predictions
            b<-data.frame(a);names(b)<-paste0('P',1:length(R))
            b<-cal_apply(b,cm)
            b<-apply(b,2,fx)
            b
          })
          dd<-data.frame(am[,vr],t(y));names(dd)[1]<-vr
          return(dd)
        }
        
        # if (is.null(cm)){
        #   y<-sapply(am[,vr],function(x){
        #     d<-am;d[,vr]<-x
        #     predict(sm,d)$predictions
        #   })
        #   return(data.frame(am[,vr],y))
        # }
        
      }
      afsmpl<-function(sm,df,vr,fx,cm,R){
        y<-sapply(df[,vr],function(x){
          d<-df;d[,vr]<-x
          #a<-fx(predict(sm,d)$predictions)
          #b<-predict(cm,data.frame(x=a),type='response')
          a<-predict(sm,d)$predictions
          b<-data.frame(a);names(b)<-paste0('P',1:length(R))
          b<-cal_apply(b,cm)
          b
        })
        z<-data.frame(y)
        z<-lapply(1:length(z),function(x){data.frame(Val=df[x,vr],z[[x]])})
        z<-bind_rows(z)
        dd<-data.frame(vr,z);names(dd)[1]<-'Var'
        return(dd)
      }
      ad<-vector('list',2*40)
      adsmpl<-vector('list',2*40)
      n<-1;for (spc in c('AN','PL')){
        
        df<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_TransNorm_Mod.rds'))[[spc]]
        ac<-strsplit(chr(formula(eval(rfave0)))[3],' \\+ ',perl=T)[[1]]
        ac<-gsub('\n','',ac,fixed=T)
        va<-gsub(' ','',ac,fixed=T)
        
        for (vr in va){
          
          print(c(spc,vr))
          ad[[n]]<-af(sm0[[spc]][[MI[[spc]]]][[1]],df,vr,fx,cm0[[spc]][[MI[[spc]]]][[1]],R)
          #adsmpl[[n]]<-afsmpl(sm0[[spc]][[MI[[spc]]]][[1]],df,vr,fx,cm0[[spc]][[MI[[spc]]]][[1]],R)
          # ggplot(adsmpl[[n]],aes(x=Val,y=P2))+
          #   stat_density_2d(aes(fill=log1p(..density..)))
          
          ii<-class(df[,vr])
          if (ii=='numeric'){ad[[n]]<-cbind(spc,vr,ad[[n]][,1],NA,ad[[n]][,2:length(ad[[n]])])}
          if (ii=='factor'){ad[[n]]<-cbind(spc,vr,NA,ad[[n]][,1],ad[[n]][,2:length(ad[[n]])])}
          ad[[n]]<-data.frame(ad[[n]])
          names(ad[[n]])[1:4]<-c('Spc','Var','Val_Nmr','Val_Chr')
          ad[[n]]$Val_Nmr<-nmr(ad[[n]]$Val_Nmr)
          ad[[n]]$Val_Chr<-chr(ad[[n]]$Val_Chr)
          #ad[[n]]$Prob<-nmr(ad[[n]]$Prob)
          ad[[n]]$Class<-ii
          names(ad)[n]<-vr
          n<-n+1
          
        }
      }
      dd<-bind_rows(ad)
      dd<-AddVarCategories(dd,scope)
      dd$Mtrc<-lbl
      saveRDS(dd,paste0("E:/",scope,"_data/Figures/PDP_FarmPresence_",lbl,".rds"))
    }
    
    #PDP
    {
      library(ggplot2)
      library(RColorBrewer) #RColorBrewer::display.brewer.all() 
      library(cowplot)
      library(extrafont)
      loadfonts(device="win",quiet=T)
      #font_import(paths="C:/WINDOWS/FONTS",pattern='Times')
      #fonts()
      
      #dd<-readRDS("E:/Global_data/Figures/PDP_FarmPresence.rds")
      dd<-readRDS(paste0("E:/",scope,"_data/Figures/PDP_FarmPresence_",lbl,".rds"))
      dd<-ChangeVarNames(dd,scope)
      dd$Rng<-NA;dd$Step<-NA
      #spc<-'AN' #Arbitrary
      # ac<-strsplit(chr(formula(eval(rfave0)))[3],' \\+ ',perl=T)[[1]]
      # ac<-gsub('\n','',ac,fixed=T)
      # nn<-gsub(' ','',ac,fixed=T)
      # nn<-unique(nn)
      nn<-unique(dd$Var)
      dd$Lbl<-chr(dd$Var)
      dd$Lbl[dd$Lbl=='Fish farm presence (at point)']<-'Fish farm presence\n(at point)'
      dd$Lbl[dd$Lbl=='Seaweed farm presence (at point)']<-'Seaweed farm presence\n(at point)'
      dd$Lbl[dd$Lbl=='Fish farm presence (at buffer)']<-'Fish farm presence\n(at buffer)'
      dd$Lbl[dd$Lbl=='Seaweed farm presence (at buffer)']<-'Seaweed farm presence\n(at buffer)'
      #Additional changes to details
      dd$Lbl<-gsub('(diff)','(diff.)',dd$Lbl,fixed=T)
      dd$Lbl<-gsub('Sea temperature','Seawater temperature',dd$Lbl,fixed=T)
      dd$Lbl<-gsub('Chlorophyll-A','Chlorophyll-a',dd$Lbl,fixed=T)
      dd$Lbl<-gsub('Elevation/depth','Sea depth',dd$Lbl,fixed=T)
      
      #
      dd0<-dd
      
      usecutoff<-F
      if (scope=='Global'){PP0<-'P2'}
      if (scope=='China'){PP0<-c('P1','P3')}
      for (PP in PP0){
        
        #PDP 1
        dd<-dd0
        vthr<-1e-3
        if (scope=='Global'){
          #Arbitrary cut-off
          #sv<-amfi[amfi$Lbl>40,c('Spc','Var','q0.5','Lbl')]
          #sv<-amfi[amfi$Lbl>=50,c('Spc','Var','q0.5','Lbl')]
          
          #Cut-off based on limit of shadow feature importance
          ctf1<-amfi$q1[amfi$VarOld=="shadowMax" & amfi$Spc=='Fish']#/max(amfi$q0.5)*100 #Cut-off
          ctf2<-amfi$q1[amfi$VarOld=="shadowMax" & amfi$Spc=='Seaweed']#/max(amfi$q0.5)*100 #Cut-off
          ctfvar<-'q0'
          if (usecutoff==F){ctf1<--99;ctf2<--99}
          sv1<-amfi[amfi$Spc=='Fish' & amfi[,ctfvar]>ctf1,c('Spc','Var','q0.5','Lbl')]
          sv2<-amfi[amfi$Spc=='Seaweed' & amfi[,ctfvar]>ctf2,c('Spc','Var','q0.5','Lbl')]
          sv<-rbind(sv1,sv2)
          
          print(length(unique(sv$Var)))
          
        } #Based on feature importance
        if (scope=='China'){
          #Arbitrary cut-off
          #sv<-amfi[amfi$Lbl>20,c('Spc','Var','q0.5','Lbl')]
          #sv<-amfi[amfi$Lbl>=50,c('Spc','Var','q0.5','Lbl')]
          
          #Cut-off based on limit of shadow feature importance
          ctf1<-amfi$q0.5[amfi$VarOld=="shadowMax" & amfi$Spc=='Fish']#/max(amfi$q0.5)*100 #Cut-off
          ctf2<-amfi$q0.5[amfi$VarOld=="shadowMax" & amfi$Spc=='Seaweed']#/max(amfi$q0.5)*100 #Cut-off
          ctfvar<-'q0'
          if (usecutoff==F){ctf1<--99;ctf2<--99}
          sv1<-amfi[amfi$Spc=='Fish' & amfi[,ctfvar]>ctf1,c('Spc','Var','q0.5','Lbl')]
          sv2<-amfi[amfi$Spc=='Seaweed' & amfi[,ctfvar]>ctf2,c('Spc','Var','q0.5','Lbl')]
          sv<-rbind(sv1,sv2)
          
          print(length(unique(sv$Var)))
          
        } #Based on feature importance
        #Smooth curve (same as below)
        movavrg<-function(am,stp=0.1,wnd=0.1,fx=mean){
          data.frame(x=seq(0,1,stp),y=sapply(seq(0,1,stp),function(x){fx(am[am[,1]>=x-wnd & am[,1]<=x+wnd,2])}))
        }
        uu<-vector('list',2);names(uu)<-c('AN','PL')
        for (spc in c('AN','PL')){
          amdl<-vector('list',length(nn));names(amdl)<-nn
          for (l in nn){
            amdl[[l]]<-data.frame(Val_Nmr=dd[dd$Var==l & dd$Spc==spc,c('Val_Nmr')],
                                  Val_Chr=dd[dd$Var==l & dd$Spc==spc,c('Val_Chr')],
                                  y=dd[dd$Var==l & dd$Spc==spc,PP])
            if (!any(dd$Class[dd$Var==l]=='factor')){
              ##
              amdl[[l]]<-movavrg(amdl[[l]][,c(1,3)])
              #Missing values linearly interpolated
              av<-is.finite(amdl[[l]][,'y'])
              amdl[[l]][!av,]<-data.frame(approx(x=amdl[[l]][av,'x'],y=amdl[[l]][av,'y'],xout=amdl[[l]][!av,'x']))
              ##
              # av<-av[order(av$Val_Nmr),];av<-av[,2] #Already sorted
              amdl[[l]]$AvSlope<-c(0,amdl[[l]]$y[2:length(amdl[[l]]$y)]-amdl[[l]]$y[1:(length(amdl[[l]]$y)-1)])
            }else{
              ##
              amdl[[l]]<-amdl[[l]][!duplicated(amdl[[l]][,'Val_Chr']),]
              amdl[[l]]<-data.frame(x=0,y=amdl[[l]][,'y'])
              ##
              amdl[[l]]$AvSlope<-0
            }
            amdl[[l]]$Var<-l
            amdl[[l]]$Spc<-spc
            amdl[[l]]$Mean<-mean(amdl[[l]]$y)
            amdl[[l]]$Variance<-var(amdl[[l]]$y)
            amdl[[l]]$Cat<-unique(dd[dd$Var==l & dd$Spc==spc,c('Cat')])
          }
          uu[[spc]]<-bind_rows(amdl)
        }
        amdl<-bind_rows(uu);names(amdl)[1:2]<-c('Val',PP)
        amdl0<-amdl
        aplt<-vector('list',length(nn));names(aplt)<-nn
        varvAN<-vector('numeric',length(nn));names(varvAN)<-nn
        varvPL<-vector('numeric',length(nn));names(varvPL)<-nn
        ballsize<-c(0.5,5)
        balllimits<-c(0.01,3.5)
        ballbreaks=seq(1,3.5,0.5)
        if (scope=='Global'){
          aa<-data.frame(Subregion=c('Australia and New Zealand','Central Asia','Eastern Asia','Eastern Europe','Latin America and the Caribbean','Melanesia','Micronesia','Northern Africa','Northern America','Northern Europe','Other','Polynesia','South-eastern Asia','Southern Asia','Southern Europe','Sub-Saharan Africa','Western Asia','Western Europe'),
                         Val=c('ANZ','C-AS','E-AS','E-EU','LAC','MELA','MICR','N-AF','N-AM','N-EU','Other','POLY','SE-AS','S-AS','S-EU','SS-AF','W-AS','W-EU'))
        }
        if (scope=='China'){
          #unique(dd$Val_Chr[dd$Var=='Subregion'])
          aa<-data.frame(Subregion=unique(dd$Val_Chr[dd$Var=='Subregion']),
                         Val=chr(1:length(unique(dd$Val_Chr[dd$Var=='Subregion']))))
        }
        k<-1;for (i in nn){
          
          if (unique(dd$Class[dd$Var==i])=='numeric'){
            
            av<-diff(range(dd[dd$Var==i,PP]))
            av<-range(dd[dd$Var==i,PP])[2] #Range of plot calculated from zero
            dd$Rng[dd$Var==i]<-av 
            dd$Step[dd$Var==i]<-0.01
            dd$Lnwdth[dd$Var==i]<-0.05
            
            #Smooth curve and segments
            amdl<-amdl0[amdl0$Var==i,]
            ffAN<-data.frame(approx(x=amdl[amdl$Spc=='AN','Val'],y=amdl[amdl$Spc=='AN',PP],xout=dd[dd$Var==i,'Val_Nmr']),Spc='AN')
            ffPL<-data.frame(approx(x=amdl[amdl$Spc=='PL','Val'],y=amdl[amdl$Spc=='PL',PP],xout=dd[dd$Var==i,'Val_Nmr']),Spc='PL')
            ff<-rbind(ffAN,ffPL)
            
            #Variable filtering
            varvAN[k]<-mean(amdl$Variance[amdl$Spc=='AN'])
            varvPL[k]<-mean(amdl$Variance[amdl$Spc=='PL'])
            if (!is.element(i,sv$Var[sv$Spc=='Fish'])){
              ff<-ff[-which(ff$Spc=='AN'),]
              dd<-dd[-which(dd$Var==i & dd$Spc=='AN'),]
              amdl<-amdl[-which(amdl$Var==i & amdl$Spc=='AN'),]
            }
            if (!is.element(i,sv$Var[sv$Spc=='Seaweed'])){
              ff<-ff[-which(ff$Spc=='PL'),]
              dd<-dd[-which(dd$Var==i & dd$Spc=='PL'),]
              amdl<-amdl[-which(amdl$Var==i & amdl$Spc=='PL'),]
            }
            # if (varvAN[k]<vthr){
            #   ff<-ff[-which(ff$Spc=='AN'),]
            #   dd<-dd[-which(dd$Var==i & dd$Spc=='AN'),]
            #   amdl<-amdl[-which(amdl$Var==i & amdl$Spc=='AN'),]
            # }
            # if (varvPL[k]<vthr){
            #   ff<-ff[-which(ff$Spc=='PL'),]
            #   dd<-dd[-which(dd$Var==i & dd$Spc=='PL'),]
            #   amdl<-amdl[-which(amdl$Var==i & amdl$Spc=='PL'),]
            # }
            if (nrow(amdl)==0){k<-k+1;next()}
            
            clr<-unique(amdl$Spc[amdl$Var==i])
            if (length(clr)==1){clr<-ifelse(clr=='AN',5,2)}
            if (length(clr)==2){clr<-c(5,2)}
            aplt[[k]]<-ggplot(mapping=aes(x=Val_Nmr,y=!!sym(PP),color=Spc,fill=Spc))+
              #geom_rug(aes(y=2,linewidth=Lnwdth,color=Spc),data=dd[dd$Var==i,],sides="b",outside=F)+#coord_cartesian(clip="off")+
              geom_segment(aes(color=Spc,xend=Val_Nmr,
                               # y=ifelse(Spc=='AN',(-0.1*Rng-Step+0.1*Rng/2),(-0.1*Rng-Step+0.0*Rng/2)),
                               # yend=ifelse(Spc=='AN',(-0.1*Rng-Step+0.2*Rng/2),(-0.1*Rng-Step+0.1*Rng/2))),
                               y=ifelse(Spc=='AN',-Rng*(Step+Lnwdth),-Rng*2*(Step+Lnwdth)),
                               yend=ifelse(Spc=='AN',-Rng*(Step),-Rng*(2*Step+Lnwdth))),
                           data=dd[dd$Var==i,])+
              # geom_point(aes(x=x,y=y,color=Spc,fill=Spc),data=ff,size=1.8,shape='\u007c',alpha=0.8,
              #            position=ggforce::position_jitternormal(0.05,0.01))+
              geom_step(data=dd[dd$Var==i,],linetype=1,linewidth=2,alpha=0.2)+
              geom_line(aes(x=Val,y=!!sym(PP),color=Spc),data=amdl,linetype=1,linewidth=1.2,alpha=1.)+
              #General
              theme(text=element_text(size=14,family="Times New Roman",color='black'))+
              theme(panel.background = element_rect(fill="White",colour="black"),
                    panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
              # theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
              #       panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
              #Legend
              theme(legend.position="none",legend.title=element_blank(),
                    legend.text=element_text(size=14,color='black'),
                    legend.key = element_rect(color=NA,fill=NA),
                    #legend.key.height=unit(0.75,"cm"),
                    #legend.key.size=unit(1,"cm"),
                    #legend.key.width=unit(1,"cm"),
                    legend.spacing.y=unit(0.4,'cm'),
                    plot.title=element_text(size=14,hjust=0.5))+
              labs(title=unique(dd$Lbl[dd$Var==i]))+
              #Axis
              #Margins from top, clockwise
              xlab(expression(" "))+
              ylab(expression(" "))+
              # ylab(paste0(
              #   paste0(round(diff(range(dd[dd$Spc=='AN' & dd$Var==i,PP]))*100,2),'%'),'\n',
              #   paste0(round(diff(range(dd[dd$Spc=='PL' & dd$Var==i,PP]))*100,2),'%')
              # ))+
              theme(axis.title.x=element_text(size=14,margin=ggplot2::margin(t=0,r=0,b=0,l=0)), #t=10
                    axis.title.y=element_text(size=14,margin=ggplot2::margin(t=0,r=15,b=0,l=0)))+ #r=20
              scale_x_continuous(breaks=seq(0,1,0.5),minor_breaks=seq(0,1,0.25))+ #,position='bottom',limits=rev #expand=expansion(add=c(0.05,0.05),mult=c(0,0))
              #scale_y_continuous(expand=expansion(mult=c(0.1,0.1)))+
              #coord_trans(y="log1p")+
              theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+ 
              theme(axis.text.x=element_text(size=12,color="black"))+
              theme(axis.text.y=element_text(size=12,color="black"))+
              theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+ #
              theme(axis.ticks.length=unit(.1,"cm"))+
              #Others
              scale_color_manual(values=brewer.pal(9,'Set1')[clr])+
              scale_fill_manual(values=brewer.pal(9,'Set1')[clr])
            
            k<-k+1
          }
          
          if (unique(dd$Class[dd$Var==i])=='factor'){
            
            au<-expand.grid(Spc=c('AN','PL'),Val_Chr=sort(unique(dd[dd$Var==i,'Val_Chr'])),lngth=nmr(NA),
                            q0=nmr(NA),q0.25=nmr(NA),q0.5=nmr(NA),q0.75=nmr(NA),q1=nmr(NA),Min=nmr(NA),Max=nmr(NA),Variance=nmr(NA))
            for (p in unique(au$Spc)){
              for (s in unique(au$Val_Chr)){
                for (r in seq(0,1,0.25)){
                  av<-dd[dd$Spc==p & dd$Var==i & dd$Val_Chr==s,PP]
                  au[au$Spc==p & au$Val_Chr==s,paste0('q',r)]<-quantile(av[is.finite(av)],r)
                  au[au$Spc==p & au$Val_Chr==s,'Min']<-min(av[is.finite(av)])
                  au[au$Spc==p & au$Val_Chr==s,'Max']<-max(av[is.finite(av)])
                  au[au$Spc==p & au$Val_Chr==s,'lngth']<-length(av[is.finite(av)])
                }
              }
              av<-dd[dd$Spc==p & dd$Var==i,PP]
              au[au$Spc==p,'Variance']<-var(unique(av)) #Variance across all levels
            }
            print(levels(au$Val_Chr))
            aulgnd<-au
            
            varvAN[k]<-mean(au$Variance[au$Spc=='AN'])
            varvPL[k]<-mean(au$Variance[au$Spc=='PL'])
            if (!is.element(i,sv$Var[sv$Spc=='Fish'])){au<-au[-which(au$Spc=='AN'),]}
            if (!is.element(i,sv$Var[sv$Spc=='Seaweed'])){au<-au[-which(au$Spc=='PL'),]}
            # if (varvAN[k]<vthr){au<-au[-which(au$Spc=='AN'),]}
            # if (varvPL[k]<vthr){au<-au[-which(au$Spc=='PL'),]}
            if (nrow(au)==0){k<-k+1;next()}
            
            clr<-unique(au$Spc)
            if (length(clr)==1){clr<-ifelse(clr=='AN',5,2)}
            if (length(clr)==2){clr<-c(5,2)}
            aplt[[k]]<-ggplot()+
              #geom_boxplot(aes(lower=q0.25,middle=q0.5,upper=q0.75,min=q0,max=q1,x=Val_Chr,color=Spc),stat="identity",staplewidth=0.8,linewidth=1.2,position=position_dodge2())+
              #geom_point(aes(x=Val_Chr,y=!!sym(PP),color=Spc),data=dd[dd$Var==i,],size=1,shape=16,alpha=0.25,
              #position=ggforce::position_jitternormal(0.05,0.05))+
              #position=position_jitterdodge(0.25,0.25,1))+
              #geom_text(aes(y=0.1,x=Val_Chr,label=paste0('n = ',lngth),family='Times New Roman'),data=au,position=position_dodge2(width=1))+
              geom_point(aes(y=q0.5,x=Val_Chr,color=Spc,size=log10(lngth)),data=au,position=position_dodge(width=1))+
              #General
              theme(text=element_text(size=16,family="Times New Roman",color='black'))+
              theme(panel.background = element_rect(fill="White",colour="black"),
                    panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
              # theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
              #       panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
              #Legend
              theme(legend.position="none",legend.title=element_blank(),
                    #legend.position.inside = c(0,0.04),
                    legend.background=element_blank(),
                    legend.text=element_text(size=14,color='black'),
                    legend.key = element_rect(color=NA,fill=NA),
                    #legend.key.height=unit(0.75,"cm"),
                    #legend.key.size=unit(1,"cm"),
                    #legend.key.width=unit(1,"cm"),
                    legend.spacing.y=unit(0.4,'cm'),
                    plot.title=element_text(size=14,hjust=0.5))+
              labs(title=unique(dd$Lbl[dd$Var==i]))+
              #Axis
              #Margins from top, clockwise
              xlab(expression(" "))+
              ylab(expression(" "))+
              # ylab(paste0(
              #   paste0(round(diff(range(dd[dd$Spc=='AN' & dd$Var==i,PP]))*100,2),'%'),'\n',
              #   paste0(round(diff(range(dd[dd$Spc=='PL' & dd$Var==i,PP]))*100,2),'%')
              # ))+
              theme(axis.title.x=element_text(size=14,margin=ggplot2::margin(t=0,r=0,b=0,l=0)), #t=10
                    axis.title.y=element_text(size=14,margin=ggplot2::margin(t=0,r=15,b=0,l=0)))+ #r=20
              #scale_x_continuous(breaks=seq(0,1,0.5),minor_breaks=seq(0,1,0.25),expand=expansion(add=c(0.05,0.05),mult=c(0,0)))+ #,position='bottom',limits=rev
              #coord_trans(y="log1p")+
              scale_x_discrete(labels=c('No','Yes'))+
              scale_y_continuous(expand=expansion(mult=c(0.1,0.1)))+
              theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+ 
              theme(axis.text.x=element_text(size=12,color="black"))+
              theme(axis.text.y=element_text(size=12,color="black"))+
              theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
              theme(axis.ticks.length=unit(.1,"cm"))+
              #Others
              scale_color_manual(values=brewer.pal(9,'Set1')[c(clr,clr)])+
              scale_fill_manual(values=brewer.pal(9,'Set1')[c(clr,clr)])+
              scale_size_continuous(limits=balllimits,range=ballsize,breaks=ballbreaks)
            
            if (i=='Subregion'){
              aplt[[k]]<-ggplot()+
                #geom_boxplot(aes(lower=q0.25,middle=q0.5,upper=q0.75,min=q0,max=q1,x=Val_Chr,color=Spc),stat="identity",staplewidth=0.8,linewidth=1.2,position=position_dodge2())+
                #geom_point(aes(x=Val_Chr,y=!!sym(PP),color=Spc),data=dd[dd$Var==i,],size=1,shape=16,alpha=0.25,
                #position=ggforce::position_jitternormal(0.05,0.05))+
                #position=position_jitterdodge(0.25,0.25,1))+
                #geom_text(aes(y=q0.5+rep(c(0.005,0.003),9),x=Val_Chr,label=paste0('n = ',lngth),family='Times New Roman',size=10,angle=45),data=au,position=position_dodge2(width=1))+
                geom_point(aes(y=q0.5,x=Val_Chr,color=Spc,size=log10(lngth)),data=au,position=position_dodge(width=1))+
                #General
                theme(text=element_text(size=16,family="Times New Roman",color='black'))+
                theme(panel.background = element_rect(fill="White",colour="black"),
                      panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
                # theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
                #       panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
                #Legend
                theme(legend.position="none",legend.title=element_blank(),
                      legend.position.inside = c(0.5,0.75),
                      legend.background=element_blank(),
                      legend.text=element_text(size=14,color='black'),
                      legend.key = element_rect(color=NA,fill=NA),
                      #legend.key.height=unit(0.75,"cm"),
                      #legend.key.size=unit(1,"cm"),
                      #legend.key.width=unit(1,"cm"),
                      legend.spacing.y=unit(0.4,'cm'),
                      plot.title=element_text(size=14,hjust=0.5))+
                guides(size=guide_legend(ncol=4))+
                labs(title=unique(dd$Lbl[dd$Var==i]))+
                #Axis
                #Margins from top, clockwise
                xlab(expression(" "))+
                ylab(expression(" "))+
                # ylab(paste0(
                #   paste0(round(diff(range(dd[dd$Spc=='AN' & dd$Var==i,PP]))*100,2),'%'),'\n',
                #   paste0(round(diff(range(dd[dd$Spc=='PL' & dd$Var==i,PP]))*100,2),'%')
                # ))+
                theme(axis.title.x=element_text(size=14,margin=ggplot2::margin(t=0,r=0,b=0,l=0)), #t=10
                      axis.title.y=element_text(size=14,margin=ggplot2::margin(t=0,r=15,b=0,l=0)))+ #r=20
                #scale_x_continuous(breaks=seq(0,1,0.5),minor_breaks=seq(0,1,0.25),expand=expansion(add=c(0.05,0.05),mult=c(0,0)))+ #,position='bottom',limits=rev
                #coord_trans(y="log1p")+
                scale_x_discrete(labels=left_join(data.frame(Subregion=levels(au$Val_Chr)),aa,join_by(Subregion))$Val)+
                scale_y_continuous(expand=expansion(mult=c(0.1,0.1)))+
                theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+ 
                theme(axis.text.x=element_text(size=8,color="black",angle=45,margin=margin(t=10)))+
                theme(axis.text.y=element_text(size=12,color="black"))+
                theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
                theme(axis.ticks.length=unit(.1,"cm"))+
                #Others
                scale_color_manual(values=brewer.pal(9,'Set1')[c(clr,clr)])+
                scale_fill_manual(values=brewer.pal(9,'Set1')[c(clr,clr)])+
                scale_size_continuous(limits=balllimits,range=ballsize,breaks=ballbreaks)
            }
            
            k<-k+1
          }
          
        }
        #Legend
        {
          #https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
          av<-aplt[['Population']]+
            theme(legend.position="right",legend.title=element_blank(),
                  legend.text=element_text(size=14,color='black'),
                  legend.key = element_rect(color=NA,fill=NA),
                  #legend.key.height=unit(0.75,"cm"),
                  #legend.key.size=unit(1,"cm"),
                  legend.key.width=unit(1.2,"cm"),
                  legend.spacing.x=unit(0.4,'cm'),
                  legend.spacing.y=unit(0.4,'cm'))+
            scale_color_hue(labels=c("Fish", "Seaweed"))
          apltlgnd1<-cowplot::get_plot_component(av,'guide-box-right',return_all=T)
          av<-aplt[['Subregion']]+
            theme(legend.position="right",legend.title=element_blank(),
                  legend.text=element_text(size=14,color='black'),
                  legend.key = element_rect(color=NA,fill=NA),
                  #legend.key.height=unit(0.75,"cm"),
                  #legend.key.size=unit(1,"cm"),
                  legend.key.width=unit(1.2,"cm"),
                  legend.spacing.x=unit(0.4,'cm'),
                  legend.spacing.y=unit(0.4,'cm'))+
            scale_color_hue(labels=c("Fish", "Seaweed"))
          apltlgnd2<-cowplot::get_plot_component(av,'guide-box-right',return_all=T)
        }
        #sv<-unique(c(names(varvAN)[varvAN>vthr],names(varvPL)[varvPL>vthr]))  #Based on variance of probability
        #aplt[[i]]+theme(legend.position="top") #Check legend
        #pp<-plot_grid(plotlist=aplt[sv],nrow=3,ncol=4)
        # if (scope=='Global'){nrowp<-4+1;ncolp<-7;wdth<-30}
        # if (scope=='China'){nrowp<-3;ncolp<-6;wdth<-30}
        ##pp<-plot_grid(plotlist=aplt[chr(unique(sv$Var))],nrow=nrowp,ncol=ncolp)
        {
          ## pp1<-plot_grid(plotlist=aplt[chr(unique(sv$Var))][1:21],nrow=3,ncol=7)
          ## pp2<-plot_grid(plotlist=aplt[chr(unique(sv$Var))][22:26],nrow=1,ncol=5)
          ## pp3<-plot_grid(pp2,aplt[['Subregion']],nrow=1,ncol=2,rel_widths=c(2.5,1))
          ## pp4<-plot_grid(pp1,pp3,nrow=2,ncol=1,rel_heights=c(3,1))
          ## ppl<-plot_grid(apltlgnd1,apltlgnd2,nrow=1,ncol=2)
          ## pp5<-plot_grid(pp4,NULL,ppl,nrow=3,ncol=1,rel_heights=c(10,0.5,1))
          ## pp<-pp5
        }
        #Ordered
        vv<-c("Population","GDP per capita","Aquatic exports","Aquatic imports",
              "Population (diff)","GDP per capita (diff)","Aquatic exports (diff)","Aquatic imports (diff)",
              "Fish farm presence (at point)","Seaweed farm presence (at point)","Fish farm presence (at buffer)","Seaweed farm presence (at buffer)",
              "Wild catch","Wild catch (diff)",
              "Sea temperature","pH","Salinity","Chlorophyll-A","Seawater velocity",
              "Cropland proportion","Forest proportion","Urban proportion","Water proportion",
              "Elevation/depth","Coastline complexity",
              "Protected area","Subregion")
        aplt<-aplt[vv]
        pp1<-plot_grid(plotlist=aplt[1:24],nrow=6,ncol=4)
        pp2<-plot_grid(plotlist=aplt[25:26],nrow=1,ncol=2)
        pp3<-plot_grid(pp2,aplt[['Subregion']],nrow=1,ncol=2,rel_widths=c(1,1))
        pp4<-plot_grid(pp1,pp3,nrow=2,ncol=1,rel_heights=c(6,1))
        ppl<-plot_grid(apltlgnd1,apltlgnd2,nrow=1,ncol=2)
        pp5<-plot_grid(pp4,NULL,ppl,nrow=3,ncol=1,rel_heights=c(2.5*7,0.5,1))
        pp<-pp5
        #Not ordered
        # pp1<-plot_grid(plotlist=aplt[chr(unique(sv$Var))][1:24],nrow=6,ncol=4)
        # pp2<-plot_grid(plotlist=aplt[chr(unique(sv$Var))][25:26],nrow=1,ncol=2)
        # pp3<-plot_grid(pp2,aplt[['Subregion']],nrow=1,ncol=2,rel_widths=c(1,1))
        # pp4<-plot_grid(pp1,pp3,nrow=2,ncol=1,rel_heights=c(6,1))
        # ppl<-plot_grid(apltlgnd1,apltlgnd2,nrow=1,ncol=2)
        # pp5<-plot_grid(pp4,NULL,ppl,nrow=3,ncol=1,rel_heights=c(2.5*7,0.5,1))
        # pp<-pp5
        ## jpeg(paste0("F:/",scope,"_data/Code files/Figures/PDP_",PP,".jpg"),units="cm",width=wdth,height=15,res=800);plot(pp);dev.off()
        jpeg(paste0("F:/",scope,"_data/Code files/Figures/PDP_",PP,".jpg"),units="cm",width=25,height=35,res=800);plot(pp);dev.off() #width=35,height=25
        #pp<-apltlgnd
        #jpeg(paste0("F:/",scope,"_data/Code files/Figures/PDP_",PP,"_Legend.jpg"),units="cm",width=wdth,height=15,res=800);plot(pp);dev.off()
        
        
        
        #PDP 2  (not rechecked; used for superplot below it seems)
        dd<-dd0
        #Smooth curve
        movavrg<-function(am,stp=0.1,wnd=0.1,fx=mean){
          data.frame(x=seq(0,1,stp),y=sapply(seq(0,1,stp),function(x){fx(am[am[,1]>=x-wnd & am[,1]<=x+wnd,2])}))
        }
        uu<-vector('list',2);names(uu)<-c('AN','PL')
        for (spc in c('AN','PL')){
          amdl<-vector('list',length(nn));names(amdl)<-nn
          for (l in nn){
            amdl[[l]]<-data.frame(Val_Nmr=dd[dd$Var==l & dd$Spc==spc,c('Val_Nmr')],
                                  Val_Chr=dd[dd$Var==l & dd$Spc==spc,c('Val_Chr')],
                                  y=dd[dd$Var==l & dd$Spc==spc,PP])
            if (!any(dd$Class[dd$Var==l]=='factor')){
              ##
              amdl[[l]]<-movavrg(amdl[[l]][,c(1,3)])
              #Missing values linearly interpolated
              av<-is.finite(amdl[[l]][,'y'])
              amdl[[l]][!av,]<-data.frame(approx(x=amdl[[l]][av,'x'],y=amdl[[l]][av,'y'],xout=amdl[[l]][!av,'x']))
              ##
              # av<-av[order(av$Val_Nmr),];av<-av[,2] #Already sorted
              amdl[[l]]$AvSlope<-c(0,amdl[[l]]$y[2:length(amdl[[l]]$y)]-amdl[[l]]$y[1:(length(amdl[[l]]$y)-1)])
            }else{
              ##
              amdl[[l]]<-amdl[[l]][!duplicated(amdl[[l]][,'Val_Chr']),]
              amdl[[l]]<-data.frame(x=0,y=amdl[[l]][,'y'])
              ##
              amdl[[l]]$AvSlope<-0
            }
            amdl[[l]]$Var<-l
            amdl[[l]]$Spc<-spc
            amdl[[l]]$Mean<-mean(amdl[[l]]$y)
            amdl[[l]]$Variance<-var(amdl[[l]]$y)
            amdl[[l]]$Cat<-unique(dd[dd$Var==l & dd$Spc==spc,c('Cat')])
          }
          uu[[spc]]<-bind_rows(amdl)
        }
        amdl<-bind_rows(uu);names(amdl)[1:2]<-c('Val',PP)
        
        #Summary statistics
        mm<-amdl
        au<-expand.grid(Spc=c('AN','PL'),Var=unique(mm$Var))
        for (p in c('AN','PL')){
          for (q in unique(mm$Var)){
            an1<-au$Spc==p & au$Var==q
            an2<-mm$Spc==p & mm$Var==q
            for (r in seq(0,1,0.25)){au[an1,paste0('q',r)]<-quantile(mm[an2,PP],r)}
            au[an1,'Mean']<-mean(mm[an2,PP])
            au[an1,'AvSlope']<-mean(mm[an2,'AvSlope'])
            au[an1,'Cat']<-unique(mm[an2,'Cat'])
          }
        }
        au$Sgn<-sign(au$AvSlope)
        au$Sgn<-ifelse(au$Sgn==0,23,ifelse(au$Sgn==1,24,25))
        au$Spc<-factor(chr(au$Spc),levels=c('AN','PL'))
        au$Lbl<-chr(au$Var)
        # au$Lbl[au$Lbl=='Fish farm presence (at point)']<-'Fish farm presence\n(at point)'
        # au$Lbl[au$Lbl=='Seaweed farm presence (at point)']<-'Seaweed farm presence\n(at point)'
        # au$Lbl[au$Lbl=='Fish farm presence (at buffer)']<-'Fish farm presence\n(at buffer)'
        # au$Lbl[au$Lbl=='Seaweed farm presence (at buffer)']<-'Seaweed farm presence\n(at buffer)'
        apdp<-au
        #
        if (F){
          
          aplt<-ggplot(au)+
            # geom_boxplot(aes(xlower=q0.25,xmiddle=q0.5,xupper=q0.75,xmin=q0,xmax=q1,y=Var,color=Spc),stat="identity",staplewidth=0.8,linewidth=1.2)+
            #geom_boxplot(aes(lower=q0.25,middle=q0.5,upper=q0.75,min=q0,max=q1,x=Var,color=Spc),stat="identity",position=position_dodge(width=1.2))+
            geom_point(aes(y=q0.5,x=Lbl,color=Spc,fill=Spc,shape=I(Sgn)),size=4,stroke=0,alpha=1.0,position=position_dodge(width=1.))+
            geom_point(aes(y=Mean,x=Lbl,color=Spc,fill=Spc,shape=I(Sgn)),size=4,stroke=0,alpha=1.0,position=position_dodge(width=1.))+
            geom_boxplot(aes(lower=q0.25,middle=q0.5,upper=q0.75,min=q0,max=q1,x=Lbl,color=Spc),stat="identity",staplewidth=0,width=0,linewidth=0.8,alpha=1.0,position=position_dodge(width=1.))+
            #coord_flip()+
            #General
            theme(text=element_text(size=14,family="Times New Roman",color='black'))+
            theme(panel.background = element_rect(fill="White",colour="black"),
                  panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
            theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
                  panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
            #Legend
            theme(legend.position="right",legend.title=element_blank(),
                  legend.text=element_text(size=14,color='black'),
                  legend.key = element_rect(color=NA,fill=NA),
                  #legend.key.height=unit(0.75,"cm"),
                  #legend.key.size=unit(1,"cm"),
                  #legend.key.width=unit(1,"cm"),
                  legend.spacing.y=unit(0.4,'cm'))+
            #guides(guide=guide_legend(reverse=T))+
            #Axis
            #Margins from top, clockwise
            xlab(expression(" "))+
            ylab(expression(" "))+
            theme(axis.title.x=element_text(size=14,margin=ggplot2::margin(t=10,r=0,b=0,l=0)),
                  axis.title.y=element_text(size=14,margin=ggplot2::margin(t=0,r=20,b=0,l=0)))+
            scale_x_discrete(position='top')+
            scale_y_continuous(limits=c(0,0.6),breaks=seq(0,0.6,0.1),minor_breaks=seq(0,0.6,0.1),expand=expansion(add=c(0,0),mult=c(0,0)))+
            #scale_y_log10()+
            theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+ 
            theme(axis.text.x=element_text(size=12,color="black",angle=45,hjust=0))+
            theme(axis.text.y=element_text(size=12,color="black"))+
            #theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
            theme(axis.ticks.length=unit(.1,"cm"))+
            #Others
            scale_color_manual(values=brewer.pal(9,'Set1')[c(5,2)])+
            scale_fill_manual(values=brewer.pal(9,'Set1')[c(5,2)])+
            facet_grid(cols=vars(Cat),scales="free",space='free')+
            theme(strip.background=element_blank(),
                  strip.text.x=element_blank())
          pp<-aplt
          jpeg(paste0("F:/",scope,"_data/Code files/Figures/PDP2_",PP,"_Legend.jpg"),units="cm",width=wdth,height=15,res=800);plot(pp);dev.off()
        }
        
      }
     
    }
    
  }
  
  #Get and save final model from above
  {
    library(ranger)
    sm0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/sm_',ac0,'_',Prob,'.rds'))
    cm0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/cm_',ac0,'_',Prob,'.rds'))
    av<-vector('list',2);names(av)<-c('AN','PL')
    sm<-av;cm<-av
    #vrblnm<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/SelectedVariables_',ac0,'_',Prob,'.rds'))
    MV<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/SelectedCV_',ac0,'_',Prob,'.rds'))
    for (spc in c('AN','PL')){
      
      sm[[spc]]<-sm0[[spc]][[MI[[spc]]]][[1]]
      cm[[spc]]<-cm0[[spc]][[MI[[spc]]]][[1]]
      
    }
    saveRDS(sm,paste0('F:/',scope,'_data/Code files/Modeling/RF_PR_','FarmPresence','_',Prob,'.rds'))
    saveRDS(cm,paste0('F:/',scope,'_data/Code files/Modeling/RF_CL_','FarmPresence','_',Prob,'.rds'))
  }
  
}

#Prepare dataset for glm
{
  set.seed(1000)
  #Get dataset
  am0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_','base','_','0','.rds'))
  ##Adjust sample period
  if (scope=='Global'){am0<-am0[am0$Period==2015,]}
  ##Add corrections
  if (scope=='Global'){am0<-ApplyRSAssumptions(am0)}
  ##
  amm<-vector('list',4);names(amm)<-c('AN+','AN-','PL+','PL-')
  np<-vector('list',4);names(amm)<-c('AN+','AN-','PL+','PL-')
  for (spc in c('AN','PL')){
    for (chng in c('-','+')){
      
      am<-am0
      
      ##Adjust sample scope
      if (scope=='Global'){
        av<-am[,paste0('OffshrAq_',spc,'_F')]-am[,paste0('OffshrAq_',spc)]>0
        am[,paste0('OffshrAq___',spc,'_F')][av]<-'1';am[,paste0('OffshrAq___',spc,'_F')][!av]<-'0'
        am<-am[am[,paste0('OffshrAq___',spc,'_F')]=='1',]
      }
      
      #Model assumption
      if (scope=='Global' & chng=='-'){next()}
      
      ##Adjust sample size
      gs00<-10;spt<-vect(paste0('F:/',scope,'_data/Code files/Administrative boundaries/PopulationPointsBufferv2_',gs00,'km.shp'))
      av<-sample(spt$PtID[is.element(spt$PtID,am$PtID)],min(c(NN2,length(unique(spt$PtID[is.element(spt$PtID,am$PtID)])))),F,spt$AvArCell[is.element(spt$PtID,am$PtID)])
      am<-am[is.element(am$PtID,av),]
      
      am<-TrnsfrmData(am)
      
      av<-GetSpatStr(am,dco2)
      #writeVector(av[[2]],paste0('F:/',scope,'_data/Code files/Modeling/Model_df_FarmArea_TransNorm.shp'),overwrite=T)
      amm[[paste0(spc,chng)]]<-av[[1]]
      
      #Block sampling
      av<-sapply(unique(amm[[paste0(spc,chng)]]$SptStr),function(x){sample(which(amm[[paste0(spc,chng)]]$SptStr==x),1,F)})
      amm[[paste0(spc,chng)]]<-amm[[paste0(spc,chng)]][av,]
      
      np[[paste0(spc,chng)]]<-NormData(amm[[paste0(spc,chng)]],'mm',paste0('FarmArea_MinMax',paste0(spc,chng)),save=T,scope)
      amm[[paste0(spc,chng)]]<-NormData(amm[[paste0(spc,chng)]],'mm',paste0('FarmArea_MinMax',paste0(spc,chng)),save=F,scope)
    }
  }
  saveRDS(amm,paste0('F:/',scope,'_data/Code files/Modeling/Model_df_FarmArea_TransNorm_Mod.rds'))
}
#Farm area (GLM)
{
  
  ac0<-'FarmArea'
  ii<-'numeric'
  UseWeights<-F
  VarSel<-F
  
  set.seed(4000)
  
  if (scope=='Global'){type<-'glmnet'}
  
  ppp<-vector('numeric',4);names(ppp)<-c('AN-','AN+','PL-','PL+') #For setting max change bounds, and renormalization of OffShrAq_AN_F and OffShrAq_PL_F
  #pvct<-0.1
  av<-vector('list',4);names(av)<-c('AN-','AN+','PL-','PL+')
  xt<-av;sm0<-av;yt<-av;em0<-av;aplt<-av
  #Growth instead of area
  for (spc in c('AN','PL')){
    for (chng in c('-','+')){
      
      #Model assumption
      if (scope=='Global' & chng=='-'){next()}
      
      #Prepare data
      am0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_FarmArea_TransNorm_Mod.rds'))[[paste0(spc,chng)]]
      print(dim(am0))
      #For use with GAM and RE (non-normalized)
      am0[,paste0('OffshrAq_',spc,'_F')]<-
        lt(
          BackTrnsNormVect(am0[,paste0('OffshrAq_',spc,'_F')],paste0('OffshrAq_',spc,'_F'),paste0('FarmArea_MinMax',paste0(spc,chng)),scope)-
            BackTrnsNormVect(am0[,paste0('OffshrAq_',spc)],paste0('OffshrAq_',spc),paste0('FarmArea_MinMax',paste0(spc,chng)),scope)
        )
      #plot(am0[,paste0('OffshrAq_',spc)],am0[,paste0('OffshrAq_',spc,'_F')])
      #am0[,paste0('OffshrAq_',spc,'_F')][am0[,paste0('OffshrAq_',spc,'_F')]<=0]<-0 #Not needed anymore
      ppp[[paste0(spc,chng)]]<-max(abs(am0[,paste0('OffshrAq_',spc,'_F')]))
      am0[,paste0('OffshrAq_',spc,'_F')]<-abs(am0[,paste0('OffshrAq_',spc,'_F')])/ppp[[paste0(spc,chng)]] #For use with GAM but not RE (renormalized); omitting min due to low value and time cost.
      
      #Select (approximately) independent observations (cross-validation not trivial); done above
      #am0<-am0[!duplicated(am0$SptStr),]
      # av<-sapply(unique(am0$SptStr),function(x){sample(which(am0$SptStr==x),1,F)})
      # am0<-am0[av,]
      
      #Select variables
      if (F){
        #https://stats.stackexchange.com/questions/577/is-there-any-reason-to-prefer-the-aic-or-bic-over-the-other
        ## Variable range across all years
        ac<-c('Population','GDPC','Fsh',
              ac1x,
              'LUC_Wtr','LUC_Crp','LUC_Frs','LUC_Urb',
              'thetao','ph','so','chl','sws',
              'Bth','ShrLength','ProtAreas')
        xxx<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/Model_df_','base','_','0','.rds'))
        corrplot::corrplot(as.matrix(cor(xxx[,ac],method='spearman')),method="number",is.corr=T) #method='color'
        #pairs(am[,c(paste0('OffshrAq_',spc),ac,paste0('OffshrAq_',spc,'_F'))])
        acc<-c('Population','Fsh','LUC_Wtr','ph','sws','ShrLength','ProtAreas') #AN; Cutoff of 0.5
        acc<-c('Population','Fsh','LUC_Wtr','ph','chl','sws','ShrLength') #PL; Cutoff of 0.5
      }
      #acc<-c('Population','Fsh','thetao','ph','chl','sws','Bth','ShrLength') #Bth better fit than LUC_Wtr #Old
      #acc<-c('Population','Fsh','thetao','chl','sws','Bth','ShrLength') #ph narrow prediction range #Old
      if (scope=='Global'){acc<-c('Population','Fsh','LUC_Wtr','ph','sws','ShrLength')} #Cutoff of 0.5 #ProtAreas and Subregion assumed low relevance.
      if (scope=='China'){acc<-c('Population','Fsh','PriSect','TerSect','thetao','ph','sws','ProtAreas')} #Cutoff of 0.5 #ProtAreas and Subregion assumed low relevance.
      #Manual selection below
      
      #am<-am0[am0[,paste0('OffshrAq_',spc,'_F')]>0,] #Change #
      am<-am0 #Expansion #
      
      #Build model
      # rfave<-eval(rfave0)
      # rfave<-formula(rfave)
      W<-GetWeights(am,paste0('OffshrAq_',spc))
      if (UseWeights==F){W<-NULL}
      
      if (type=='glmnet'){
        ac<-c('Population','GDPC','Fsh',
              ac1x,
              'LUC_Wtr','LUC_Crp','LUC_Frs','LUC_Urb',
              'thetao','ph','so','chl','sws',
              'Bth','ShrLength','ProtAreas')
        an<-c(paste0('OffshrAq_',spc),ac)
        
        
        
        # ac<-c('OffshrAq',mdltgt,spc,,
        # OffshrAq_'AN',OffshrAq_'PL'
        # OffshrAqBuffer',mdltgt,'AN + OffshrAqBuffer',mdltgt,'PL +
        # Population,GDPC,Fsh,ac1x,
        # Population_D,GDPC_D,Fsh_D,'ac2x',
        # LUC_Crp,LUC_Frs,LUC_Urb,LUC_Wtr,
        # thetao,so,chl,ph,sws,
        # Bth,ShrLength,ProtAreas,
        # Subregion')
        # an<-c(paste0('OffshrAq_',spc),ac)
        
        
        
        # glmnet (only relevant without dispersion model)
        # #https://www.geeksforgeeks.org/ridge-regression-vs-lasso-regression/
        # #Lasso (L1): Focus on subset of key predictors (alpha=1)
        # #Ridge (L2): Distribute importance among all predictors (alpha=0)
        library(glmnet)
        #an<-strsplit(chr(rfave),'\\+|~',perl=T);an<-an[[3]];an<-gsub(' ','',an,fixed=T)
        #uu<-glmnet(as.matrix(am[,an]),am[,paste0('OffshrAq_',spc,'_F')],family=Gamma(link='log'),alpha=0,maxit=2e5)
        #uu<-glmnet(as.matrix(am[,an]),am[,paste0('OffshrAq_',spc,'_F')],family=Gamma(link='identity'),alpha=0,maxit=2e5)
        vv<-cv.glmnet(as.matrix(am[,an]),am[,paste0('OffshrAq_',spc,'_F')],family=gaussian(link='identity'),alpha=0,maxit=2e5,standardize=F)
        #vv<-cv.glmnet(as.matrix(am[,an]),am[,paste0('OffshrAq_',spc,'_F')],family=Gamma(link='identity'),alpha=0,maxit=2e5,standardize=F)
        vv$MSE<-mean((predict(vv,as.matrix(am[,an]),s='lambda.min',type='response')-am[,paste0('OffshrAq_',spc,'_F')])^2)
        vv$Error<-mean(abs(predict(vv,as.matrix(am[,an]),s='lambda.min',type='response')-am[,paste0('OffshrAq_',spc,'_F')]))
        vv$Lambda<-vv$lambda.min
        print(vv)
        #plot(glmnet(as.matrix(am[,an]),am[,paste0('OffshrAq_',spc,'_F')],family=gaussian(link='identity'),alpha=0,maxit=2e5))
        sm<-vv
        #coef(uu)[,ncol(coef(uu))]
        #coef(vv,s='lambda.min')
        
        plot((predict(vv,as.matrix(am[,an]),s='lambda.min',type='response')),am[,paste0('OffshrAq_',spc,'_F')]);lines(0:4,0:4)
        print(sum((predict(vv,as.matrix(am[,an]),s='lambda.min',type='response')-am[,paste0('OffshrAq_',spc,'_F')])^2)/length(am[,paste0('OffshrAq_',spc,'_F')]))
        
        #plot(summary(sm)$coefficients$cond[,1],coef(uu)[,ncol(coef(uu))]);lines(-100:100,-100:100)
      }
      
      sm0[[paste0(spc,chng)]]<-sm
      
    }
  }
  saveRDS(sm0,paste0('F:/',scope,'_data/Code files/Modeling/GLMM_AR_','FarmArea','_X','.rds'))
  saveRDS(ppp,paste0('F:/',scope,'_data/Code files/Modeling/PPP.rds')) #Only transformed, not normalized
  
  #Plots
  {
    library(glmnet)
    sm0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/GLMM_AR_','FarmArea','_X','.rds'))
    uu<-vector('list',length(sm0));names(uu)<-names(sm0)
    for (i in 1:length(uu)){
      #Model assumption
      if (scope=='Global' & is.element(i,c(1,3))){next()}
      
      cfs<-coef(sm0[[i]],s='lambda.min')
      uu[[i]]<-data.frame(Var=cfs@Dimnames[[1]],Coef=cfs@x)
      #uu[[i]]<-uu[[i]][-1,] #Intercept
      uu[[i]]$q0.5<-abs(uu[[i]]$Coef)
      #uu[[i]]$Lbl<-uu[[i]]$q0.5/max(uu[[i]]$q0.5)*100
      uu[[i]]$Lbl<-uu[[i]]$q0.5/max(uu[[i]]$q0.5[uu[[i]]$Var!='(Intercept)'])*100
      uu[[i]]$Model<-gsub('[A-Z]*','',names(sm0)[i],perl=T)
      uu[[i]]$Spc<-gsub('\\+|\\-','',names(sm0)[i],perl=T)
      uu[[i]]<-AddVarCategories(uu[[i]],scope)
      uu[[i]]<-ChangeVarNames(uu[[i]],scope)
    }
    uu<-bind_rows(uu)
    uu$Spc<-ifelse(uu$Spc=='AN','Fish','Seaweed')
    ##
    vv<-amfi[,-which(names(amfi)=='Imp')]
    vv$Coef<-NA
    vv$Model<-'+-'
    au<-rbind(vv[,names(uu)],uu)
    au$Mdl<-apply(au[,c('Spc','Model')],1,function(x){paste0(x,collapse='')})
    au$Mdl<-gsub('+-',' area change',au$Mdl,fixed=T)
    au$Mdl<-gsub('-',' area decrease',au$Mdl,fixed=T)
    au$Mdl<-gsub('+',' area increase',au$Mdl,fixed=T)
    #lvl<-c("Fish area change","Seaweed area change","Fish area decrease","Seaweed area decrease","Fish area increase","Seaweed area increase")
    lvl<-c("Fish area change","Fish area decrease","Fish area increase","Seaweed area change","Seaweed area decrease","Seaweed area increase")
    au$Mdl<-factor(au$Mdl,levels=lvl[is.element(lvl,unique(au$Mdl))])
    lvl<-levels(au$Var)
    lvl<-gsub('farm presence|farm area','farm presence\\/area',lvl,perl=T)
    lvl<-lvl[!duplicated(lvl)]
    avv<-chr(au$Var)
    avv<-gsub('farm presence|farm area','farm presence\\/area',au$Var,perl=T)
    au$Var<-factor(avv,levels=lvl)
    
    library(ggplot2)
    library(RColorBrewer) #RColorBrewer::display.brewer.all() 
    library(colorspace)
    library(cowplot)
    library(extrafont)
    loadfonts(device="win",quiet=T)
    #font_import(paths="C:/WINDOWS/FONTS",pattern='Times')
    #fonts()
    
    ##
    uu<-apdp
    uu$Mdl<-'Change type'
    uu$Coef<-uu$Mean
    uu$Lbl<-factor(uu$Lbl,levels=levels(uu$Var))
    uu$Sgn2<-factor(ifelse(uu$Coef==0,'0',ifelse(uu$Coef>0,'1','-1')))
    ##
    am<-au[au$Model!='+-',]
    lvl<-levels(am$Var)
    lvl[lvl=='Fish farm presence/area (at point)']<-'Fish farm area (at point)';lvl[lvl=='Seaweed farm presence/area (at point)']<-'Seaweed farm area (at point)'
    av<-chr(am$Var)
    av[is.na(av)]<-'Intercept'
    av[av=='Fish farm presence/area (at point)']<-'Fish farm area (at point)';av[av=='Seaweed farm presence/area (at point)']<-'Seaweed farm area (at point)'
    am$Var<-factor(av,levels=c('Intercept',lvl))
    vv<-am
    vv$Lbl<-vv$Var
    vv$Sgn<-ifelse(vv$Coef==0,23,ifelse(vv$Coef>0,24,25)) #vv$Sgn<-21
    vv$Sgn2<-factor(ifelse(vv$Coef==0,'0',ifelse(vv$Coef>0,'1','-1')))
    vv$Coef<-abs(vv$Coef)
    vv$q0<-NA
    vv$q1<-NA
    vv$Lbl<-factor(vv$Lbl,levels=levels(vv$Var))
    #
    lvl<-levels(vv$Cat)
    avv<-chr(vv$Cat);avv[is.na(avv)]<-' '
    vv$Cat<-factor(avv,levels=c(' ',lvl))
    #
    ##
    
    print(chr(unique(uu$Lbl[uu$q1>=0.1])))
    #ac<-'black'
    ac<-'grey90'
    uu$SpcSgn<-factor(apply(uu[,c('Spc','Sgn')],1,function(x){paste0(ifelse(x[1]=='AN','Fish','Seaweed'),' ',ifelse(x[2]==23,'No sign',ifelse(x[2]==24,'Positive','Negative')))}))
    lgnlbl<-levels(uu$SpcSgn)
    aplt<-ggplot(uu)+
      geom_point(aes(y=Coef,x=Lbl,color=SpcSgn,fill=SpcSgn,shape=SpcSgn),size=4,stroke=0,alpha=1.0,position=position_dodge(width=1.))+
      geom_segment(aes(y=q0,yend=q1,x=Lbl,color=SpcSgn),stat="identity",linewidth=0.8,alpha=1.0,position=position_dodge(width=1.))+
      #General
      theme(text=element_text(size=14,family="Times New Roman",color='black'))+
      theme(panel.background = element_rect(fill="White",colour="black"),
            panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
      theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
            panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
      #Legend
      theme(legend.position="none",legend.title=element_blank(),
            legend.text=element_text(size=14,color='black'),
            legend.key = element_rect(color=NA,fill=NA),
            #legend.key.height=unit(0.75,"cm"),
            #legend.key.size=unit(1,"cm"),
            #legend.key.width=unit(1,"cm"),
            legend.spacing.y=unit(0.4,'cm'))+
      #guides(guide=guide_legend(reverse=T))+
      #Axis
      #Margins from top, clockwise
      xlab(expression(""))+
      ylab(expression(""))+
      theme(axis.title.x=element_text(size=14,margin=ggplot2::margin(t=10,r=0,b=0,l=0)),
            axis.title.y=element_text(size=14,margin=ggplot2::margin(t=0,r=20,b=0,l=0)))+
      scale_x_discrete(position='top')+
      #scale_y_continuous(limits=c(0,0.6),breaks=seq(0,0.6,0.1),minor_breaks=seq(0,0.6,0.1),expand=expansion(add=c(0,0),mult=c(0,0)))+
      #scale_y_log10()+
      theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+ 
      theme(axis.text.x=element_text(size=12,color=ac,angle=45,hjust=0))+
      theme(axis.text.y=element_text(size=12,color="black"))+
      #theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
      theme(axis.ticks.length=unit(.1,"cm"))+
      #Others
      scale_color_manual(values=brewer.pal(9,'Set1')[c(5,5,5,2,2,2)],labels=lgnlbl)+
      scale_fill_manual(values=brewer.pal(9,'Set1')[c(5,5,5,2,2,2)],labels=lgnlbl)+
      scale_shape_manual(values=c(25,23,24,25,23,24),labels=lgnlbl)+
      facet_grid(cols=vars(Cat),scales="free",space='free')+ #rows=vars(Mdl),
      theme(strip.background=element_blank(),
            strip.text.x=element_blank(),
            strip.text.y=element_blank())
    pp<-aplt
    jpeg(paste0("F:/",scope,"_data/Code files/Figures/Model_1_",ac,".jpg"),units="cm",width=30,height=15,res=800);plot(pp);dev.off()
    
    av<-aplt+
      theme(legend.position="right",legend.title=element_blank(),
            legend.text=element_text(size=14,color='black'),
            legend.key = element_rect(color=NA,fill=NA),
            #legend.key.height=unit(0.75,"cm"),
            #legend.key.size=unit(1,"cm"),
            legend.key.width=unit(1.2,"cm"),
            legend.spacing.x=unit(0.4,'cm'),
            legend.spacing.y=unit(0.4,'cm'))
    pp<-cowplot::get_plot_component(av,'guide-box-right',return_all=T);plot(pp)
    jpeg(paste0("F:/",scope,"_data/Code files/Figures/Model_1_Legend.jpg"),units="cm",width=10,height=10,res=800);plot(pp);dev.off()
    
    
    
    print(chr(unique(vv$Lbl[vv$Coef>=0.1])))
    #ac<-'black'
    ac<-'grey90'
    vv$SpcSgn<-factor(apply(vv[,c('Spc','Sgn2')],1,function(x){paste0(x[1],' ',ifelse(x[2]==-1,'decrease','increase'))}))
    lgnlbl<-c('Fish decrease','Fish increase','Seaweed decrease','Seaweed increase')
    aplt<-ggplot(vv)+
      geom_point(aes(y=Coef,x=Lbl,color=SpcSgn,fill=SpcSgn,shape=SpcSgn),size=4,stroke=0,alpha=1.0,position=position_dodge(width=1.))+
      #geom_boxplot(aes(lower=q0,middle=Coef,upper=q1,min=q0,max=q1,x=Lbl,color=Spc),stat="identity",staplewidth=0,width=0,linewidth=0.8,alpha=1.0,position=position_dodge(width=1.))+
      #General
      theme(text=element_text(size=14,family="Times New Roman",color='black'))+
      theme(panel.background = element_rect(fill="White",colour="black"),
            panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
      theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
            panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
      #Legend
      theme(legend.position="none",legend.title=element_blank(),
            legend.text=element_text(size=14,color='black'),
            legend.key = element_rect(color=NA,fill=NA),
            #legend.key.height=unit(0.75,"cm"),
            #legend.key.size=unit(1,"cm"),
            #legend.key.width=unit(1,"cm"),
            legend.spacing.y=unit(0.4,'cm'))+
      guides(guide=guide_legend(reverse=T))+
      #Axis
      #Margins from top, clockwise
      xlab(expression(""))+
      ylab(expression(""))+
      theme(axis.title.x=element_text(size=14,margin=ggplot2::margin(t=10,r=0,b=0,l=0)),
            axis.title.y=element_text(size=14,margin=ggplot2::margin(t=0,r=20,b=0,l=0)))+
      scale_x_discrete(position='top')+
      #scale_y_continuous(limits=c(0,0.6),breaks=seq(0,0.6,0.1),minor_breaks=seq(0,0.6,0.1),expand=expansion(add=c(0,0),mult=c(0,0)))+
      #scale_y_log10()+
      theme(plot.margin=unit(c(0.8,1.6,0.4,0.4),"cm"))+ 
      theme(axis.text.x=element_text(size=12,color=ac,angle=45,hjust=0))+
      theme(axis.text.y=element_text(size=12,color="black"))+
      #theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
      theme(axis.ticks.length=unit(.1,"cm"))+
      #Others
      scale_color_manual(values=brewer.pal(9,'Set1')[c(5,5,2,2)],labels=lgnlbl)+
      scale_fill_manual(values=brewer.pal(9,'Set1')[c(5,5,2,2)],labels=lgnlbl)+
      scale_shape_manual(values=c(24,25,24,25),labels=lgnlbl)+
      #guide=guide_legend(override.aes=list(color='black',fill="black"))
      facet_grid(cols=vars(Cat),scales="free",space='free')+ #rows=vars(Mdl),
      theme(strip.background=element_blank(),
            strip.text.x=element_blank(),
            strip.text.y=element_blank())
    pp<-aplt
    jpeg(paste0("F:/",scope,"_data/Code files/Figures/Model_2_",ac,".jpg"),units="cm",width=30,height=15,res=800);plot(pp);dev.off()
    
    av<-aplt+
      theme(legend.position="right",legend.title=element_blank(),
            legend.text=element_text(size=14,color='black'),
            legend.key = element_rect(color=NA,fill=NA),
            #legend.key.height=unit(0.75,"cm"),
            #legend.key.size=unit(1,"cm"),
            legend.key.width=unit(1.2,"cm"),
            legend.spacing.x=unit(0.4,'cm'),
            legend.spacing.y=unit(0.4,'cm'))
    pp<-cowplot::get_plot_component(av,'guide-box-right',return_all=T);plot(pp)
    jpeg(paste0("F:/",scope,"_data/Code files/Figures/Model_2_Legend.jpg"),units="cm",width=10,height=10,res=800);plot(pp);dev.off()
    
  }
  
  #More plots
  {
    if (F){
      am1<-au[!is.na(au$Cat),] #Boruta values and GlmnetCoefs
      names(am1)[is.element(names(am1),c('Lbl'))]<-'Val'
      am1$SuperMdl<-'Variable\nimportance'
      am1$Mdl<-chr(am1$Mdl)
      am1$Mdl[am1$Mdl=='Fish area change']<-1.1
      am1$Mdl[am1$Mdl=='Seaweed area change']<-1
      am1$Mdl[am1$Mdl=='Fish area increase']<-1.3
      am1$Mdl[am1$Mdl=='Seaweed area increase']<-1.2
      am1$Mdl<-nmr(am1$Mdl)
      av1<-chr(am1$Var);av2<-levels(am1$Var)
      av1<-gsub('/area','',av1,fixed=T);av2<-gsub('/area','',av2,fixed=T)
      am1$Var<-factor(av1,levels=av2)
      am2<-uu #RF curves
      am2$Mdl<-NA
      am2$SuperMdl<-'Partial\ndependence'
      am3<-vv #glmnet coefficients
      am3$Mdl<-NA
      am3$SuperMdl<-'Penalized\ncoefficients'
      av1<-chr(am3$Var);av2<-levels(am3$Var)
      av1<-gsub('farm area','farm presence/area',av1,fixed=T);av2<-gsub('farm area','farm presence/area',av2,fixed=T)
      am3$Var<-factor(av1,levels=av2)
      ##
      amm<-bind_rows(am1,am2,am3)
      amm$SuperMdl<-factor(amm$SuperMdl,levels=c('Variable\nimportance','Partial\ndependence','Penalized\ncoefficients'))
      saveRDS(amm,paste0("F:/",scope,"_data/Code files/Figures/FIPDP_data.rds"))
    }
    
    library(ggplot2)
    library(RColorBrewer) #RColorBrewer::display.brewer.all() 
    library(colorspace)
    library(cowplot)
    library(extrafont)
    loadfonts(device="win",quiet=T)
    #font_import(paths="C:/WINDOWS/FONTS",pattern='Times')
    #fonts()
    library(ggnewscale)
    library(ggh4x)
    
    amm<-readRDS(paste0("F:/",scope,"_data/Code files/Figures/FIPDP_data.rds"))
    #
    av<-chr(amm$SuperMdl);al<-levels(amm$SuperMdl)
    av[av=='Variable\nimportance']<-'Relative\nimportance\n(%)';al[al=='Variable\nimportance']<-'Relative\nimportance\n(%)'
    av[av=='Partial\ndependence']<-'Partial\ndependence\n(probability)';al[al=='Partial\ndependence']<-'Partial\ndependence\n(probability)'
    amm$SuperMdl<-factor(av,levels=al)
    #
    #Additional changes to details
    av<-chr(amm$Var);lvl<-levels(amm$Var)
    av<-gsub('(diff)','(diff.)',av,fixed=T);lvl<-gsub('(diff)','(diff.)',lvl,fixed=T)
    av<-gsub('Sea temperature','Seawater temperature',av,fixed=T);lvl<-gsub('Sea temperature','Seawater temperature',lvl,fixed=T)
    av<-gsub('Chlorophyll-A','Chlorophyll-a',av,fixed=T);lvl<-gsub('Chlorophyll-A','Chlorophyll-a',lvl,fixed=T)
    #av<-gsub(' (','\n(',av,fixed=T);lvl<-gsub(' (','\n(',lvl,fixed=T)
    #lvl<-lvl[c(1:19,22,20,21,23:40)]
    amm$Var<-factor(av,levels=lvl)
    ###amm$Spc<-factor(amm$Spc,levels=c('Seaweed','Fish','PL','AN'))
    #Order
    amm$Ord<-NA
    ov<-data.frame(amm[!is.na(amm$Mdl) & amm$Mdl<=1.1,] %>% group_by(Var,SuperMdl) %>% summarise(Order = mean(Val), n = n()))
    amm$Ord[!is.na(amm$Mdl) & amm$Mdl<=1.1]<-left_join(amm[!is.na(amm$Mdl) & amm$Mdl<=1.1,],ov,join_by(Var,SuperMdl))$Order
    #
    print(chr(unique(amm$Lbl[amm$q1>=0.1 & amm$SuperMdl=='Partial\ndependence\n(probability)'])))
    #'Viridis'
    #clr<-'Oranges'
    #clr<-'Blues'
    ac<-'black'
    #ac<-'grey80'
    lgnlbl1<-c(rep('',3),'Probability range','Probability range','Probability range')
    lgnlbl2<-c(rep('',3),'Negative correlation','Average effect','Positive correlation')
    adp1<-amm[amm$SuperMdl=='Relative\nimportance\n(%)' & is.element(amm$Mdl,c(1,1.1)),]
    adp2<-amm[amm$SuperMdl=='Partial\ndependence\n(probability)',]
    adp2$SpcSgn<-factor(apply(adp2[,c('Spc','Sgn')],1,function(x){paste0(ifelse(x[1]=='AN','Fish','Seaweed'),' ',ifelse(x[2]==23,'No sign',ifelse(x[2]==24,'Positive','Negative')))}))
    adp2$SpcSgn<-factor(chr(adp2$SpcSgn),levels=c('Seaweed Negative','Seaweed No sign','Seaweed Positive','Fish Negative','Fish No sign','Fish Positive'))
    #Cutoff
    ctv<-unique(amm$Var[!is.na(amm$Ord) & amm$Ord>=40])
    adp1<-adp1[is.element(adp1$Var,ctv),]
    adp2<-adp2[is.element(adp2$Var,ctv),]
    #
    aplt<-ggplot()+
      geom_point(aes(x=-Mdl,y=reorder(Var,-Ord,mean),color=Val,group=Spc),data=adp1,size=10,shape=15,stroke=0,alpha=1.0)+
      geom_text(aes(x=-Mdl,y=reorder(Var,-Ord,mean),label=round(Val,0),family="Times New Roman"),data=adp1,color=ifelse(round(adp1$Val[!is.na(adp1$Cat)],0)>=60,'white','black'))+
      scale_color_continuous_sequential(palette=clr,guide=guide_colourbar(title='Relative importance',title.position="top",direction="vertical"))+
      new_scale("color")+
      geom_point(aes(x=Coef,y=reorder(Var,-Ord,mean),color=SpcSgn,fill=SpcSgn,shape=SpcSgn),data=adp2,size=4,stroke=0,alpha=1.0,position=position_dodge(width=1.))+
      geom_segment(aes(x=q0,xend=q1,y=reorder(Var,-Ord,mean),color=SpcSgn),data=adp2,stat="identity",linewidth=0.8,alpha=1.0,position=position_dodge(width=1.))+
      #
      #General
      theme(text=element_text(size=16,family="Times New Roman",color='black'))+
      theme(panel.background = element_rect(fill="White",colour="black"),
            panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
      theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
            panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
      #Legend
      theme(legend.position="right",legend.title=element_text(size=12,color='black',margin=margin(t=0,r=0,b=10,l=0)),
            legend.text=element_text(size=12,color='black'),
            legend.key = element_rect(color=NA,fill=NA),
            #legend.key.height=unit(0.75,"cm"),
            #legend.key.size=unit(1,"cm"),
            #legend.key.width=unit(1,"cm"),
            legend.spacing.y=unit(0.4,'cm'),
            plot.title=element_text(size=14,face='bold',color='black',margin=margin(t=0,r=0,b=10,l=0)))+
      guides(guide=guide_legend(reverse=F),
             color=guide_legend(title='Covariate effect',ncol=2), #order=1
             fill=guide_legend(title='Covariate effect',ncol=2), #order=1
             shape=guide_legend(title='Covariate effect',ncol=2,
                                override.aes=list(color=brewer.pal(9,'Set1')[c(5,5,5,2,2,2)],
                                                  fill=brewer.pal(9,'Set1')[c(5,5,5,2,2,2)])))+
      labs(title="Farm growth probability")+
      #Axis
      #Margins from top, clockwise
      xlab(expression(" "))+
      ylab(expression(" "))+
      theme(axis.title.x=element_text(size=12,margin=ggplot2::margin(t=10,r=0,b=0,l=0)),
            axis.title.y=element_text(size=12,margin=ggplot2::margin(t=0,r=20,b=0,l=0)))+
      #scale_x_continuous(limits=c(-5,35),breaks=seq(-5,35,5),minor_breaks=seq(-5,35,2.5),expand=expansion(add=c(0,0),mult=c(0,0)))+ #,position='bottom',
      #scale_x_discrete(limits=rev)+
      scale_y_discrete(limits=rev,position="left")+
      #scale_y_log10()+
      theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+
      theme(axis.text.x=element_text(size=12,color=ac,angle=0))+
      theme(axis.text.y=element_text(size=12,color="black"))+
      #theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
      theme(axis.ticks.length=unit(0,"cm"))+
      #Others
      scale_color_manual(values=brewer.pal(9,'Set1')[c(2,2,2,5,5,5)],labels=lgnlbl1)+
      scale_fill_manual(values=brewer.pal(9,'Set1')[c(2,2,2,5,5,5)],labels=lgnlbl1)+
      scale_shape_manual(values=c(25,23,24,25,23,24),labels=lgnlbl2)+
      facet_grid(cols=vars(SuperMdl),rows=vars(Cat),scales="free",space='free')+
      theme(strip.background=element_rect(fill="NA"),
            strip.text.x=element_text(size=12,color="black",angle=0),
            strip.text.y=element_text(size=12,color="black",angle=0))+
      facetted_pos_scales(x=list(
        SuperMdl=='Relative\nimportance\n(%)'~scale_x_continuous(breaks=c(-1.1,-1.),minor_breaks=c(-1.1,-1.),expand=expansion(add=c(0.05,0.05),mult=c(0,0)),labels=c('Fish','Seaweed')),
        SuperMdl=='Partial\ndependence\n(probability)'~scale_x_continuous(limits=c(0,0.6),breaks=seq(0,0.6,0.2),minor_breaks=seq(0,0.6,0.2),expand=expansion(add=c(0.05,0.05),mult=c(0,0)))
      ))+
      force_panelsizes(cols=c(1,1.6))
    pp<-aplt
    jpeg(paste0("F:/",scope,"_data/Code files/Figures/FIPDP_Changev2_",clr,'_',ac,".jpg"),units="cm",width=24,height=15,res=800);plot(pp);dev.off()
    
    amm<-readRDS(paste0("F:/",scope,"_data/Code files/Figures/FIPDP_data.rds"))
    #
    av<-chr(amm$SuperMdl);al<-levels(amm$SuperMdl)
    av[av=='Variable\nimportance']<-'Relative\nimportance\n(%)';al[al=='Variable\nimportance']<-'Relative\nimportance\n(%)'
    av[av=='Penalized\ncoefficients']<-'Penalized\ncoefficients\n(normalized area)';al[al=='Penalized\ncoefficients']<-'Penalized\ncoefficients\n(normalized area)'
    amm$SuperMdl<-factor(av,levels=al)
    #
    av<-chr(amm$Var);al<-levels(amm$Var)
    av<-gsub('presence/area','area',av,fixed=T);al<-gsub('presence/area','area',al,fixed=T)
    av<-gsub('presence','area',av,fixed=T);al<-gsub('presence','area',al,fixed=T)
    av<-gsub('(diff)','(diff.)',av,fixed=T);al<-gsub('(diff)','(diff.)',al,fixed=T)
    av<-gsub('Sea temperature','Seawater temperature',av,fixed=T);al<-gsub('Sea temperature','Seawater temperature',al,fixed=T)
    av<-gsub('Chlorophyll-A','Chlorophyll-a',av,fixed=T);al<-gsub('Chlorophyll-A','Chlorophyll-a',al,fixed=T)
    #av<-gsub(' (','\n(',av,fixed=T);al<-gsub(' (','\n(',al,fixed=T)
    amm$Var<-factor(av,levels=al[!duplicated(al)])
    #
    amm$Sgn[amm$Var=='Intercept']<-23
    #
    #Order
    amm$Ord<-NA
    ov<-data.frame(amm[!is.na(amm$Mdl) & amm$Mdl>1.1,] %>% group_by(Var,SuperMdl) %>% summarise(Order = mean(Val), n = n()))
    amm$Ord[!is.na(amm$Mdl) & amm$Mdl>1.1]<-left_join(amm[!is.na(amm$Mdl) & amm$Mdl>1.1,],ov,join_by(Var,SuperMdl))$Order
    #
    #print(chr(unique(amm$Var[amm$Coef>=0.1 & amm$SuperMdl=='Penalized\ncoefficients'])))
    #clr<-'Oranges'
    #clr<-'Blues'
    ac<-'black'
    #ac<-'grey80'
    adp1<-amm[amm$SuperMdl=='Relative\nimportance\n(%)' & is.element(amm$Mdl,c(1.2,1.3)),]
    adp1$Mdl<-adp1$Mdl-0.2
    adp2<-amm[amm$SuperMdl=='Penalized\ncoefficients\n(normalized area)',]
    av1<-chr(adp1$Var);av2<-levels(adp1$Var)
    av1<-gsub('farm presence','farm area',av1,fixed=T);av2<-gsub('farm presence','farm area',av2,fixed=T)
    adp1$Var<-factor(av1,levels=av2)
    adp2$SpcSgn<-factor(apply(adp2[,c('Spc','Sgn')],1,function(x){paste0(ifelse(x[1]=='Fish','Fish','Seaweed'),' ',ifelse(x[2]==23,'No sign',ifelse(x[2]==24,'Positive','Negative')))}))
    lgnlbl<-levels(adp2$SpcSgn)
    adp2$SpcSgn<-factor(chr(adp2$SpcSgn),levels=c('Seaweed Negative','Seaweed No sign','Seaweed Positive','Fish Negative','Fish No sign','Fish Positive'))
    #Cutoff
    ctv<-unique(amm$Var[!is.na(amm$Ord) & amm$Ord>=15]) #20
    adp1<-adp1[is.element(adp1$Var,ctv),]
    adp2<-adp2#[is.element(adp2$Var,ctv),] #No cutoff
    #
    aplt<-ggplot()+
      geom_point(aes(x=-Mdl,y=reorder(Var,-Ord,mean),color=Val,group=Mdl),data=adp1,size=10,shape=15,stroke=0,alpha=1.0)+
      geom_text(aes(x=-Mdl,y=reorder(Var,-Ord,mean),label=round(Val,0),family="Times New Roman"),data=adp1,color=ifelse(round(adp1$Val[!is.na(adp1$Cat)],0)>=50,'white','black'))+
      scale_color_continuous_sequential(palette=clr)+
      new_scale("color")+
      geom_point(aes(x=Coef,y=reorder(Var,-Ord,mean),color=SpcSgn,fill=SpcSgn,shape=SpcSgn),data=adp2,size=4,stroke=0,alpha=1.0,position=position_dodge(width=1.))+
      #
      #General
      theme(text=element_text(size=16,family="Times New Roman",color='black'))+
      theme(panel.background = element_rect(fill="White",colour="black"),
            panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
      theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
            panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
      #Legend
      theme(legend.position="right",legend.title=element_blank(),
            legend.text=element_text(size=12,color='black'),
            legend.key = element_rect(color=NA,fill=NA),
            #legend.key.height=unit(0.75,"cm"),
            #legend.key.size=unit(1,"cm"),
            #legend.key.width=unit(1,"cm"),
            legend.spacing.y=unit(0.4,'cm'),
            plot.title=element_text(size=14,face='bold',color='black',margin=margin(t=0,r=0,b=10,l=0)))+
      #guides(guide=guide_legend(reverse=T))+
      labs(title="Farm growth extent")+
      #Axis
      #Margins from top, clockwise
      xlab(expression(" "))+
      ylab(expression(" "))+
      theme(axis.title.x=element_text(size=12,margin=ggplot2::margin(t=10,r=0,b=0,l=0)),
            axis.title.y=element_text(size=12,margin=ggplot2::margin(t=0,r=20,b=0,l=0)))+
      #scale_x_continuous(limits=c(-5,35),breaks=seq(-5,35,5),minor_breaks=seq(-5,35,2.5),expand=expansion(add=c(0,0),mult=c(0,0)))+ #,position='bottom',
      #scale_x_discrete(limits=rev)+
      scale_y_discrete(limits=rev,position="left")+
      #scale_y_log10()+
      theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+
      theme(axis.text.x=element_text(size=12,color=ac,angle=0))+
      theme(axis.text.y=element_text(size=12,color="black"))+
      #theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
      theme(axis.ticks.length=unit(0,"cm"))+
      #Others
      scale_color_manual(values=brewer.pal(9,'Set1')[c(2,2,5,5)],labels=lgnlbl)+
      scale_fill_manual(values=brewer.pal(9,'Set1')[c(2,2,5,5)],labels=lgnlbl)+
      scale_shape_manual(values=c(25,24,25,24),labels=lgnlbl)+
      facet_grid(cols=vars(SuperMdl),rows=vars(Cat),scales="free",space='free')+
      theme(strip.background=element_rect(fill="NA"),
            strip.text.x=element_text(size=12,color="black",angle=0),
            strip.text.y=element_text(size=12,color="black",angle=0))+
      facetted_pos_scales(x=list(
        SuperMdl=="Relative\nimportance\n(%)"~scale_x_continuous(breaks=c(-1.1,-1.),minor_breaks=c(-1.1,-1.),expand=expansion(add=c(0.05,0.05),mult=c(0,0)),labels=c('Fish','Seaweed')),
        SuperMdl=="Penalized\ncoefficients\n(normalized area)"~scale_x_continuous(limits=c(0,0.6),breaks=seq(0,0.6,0.2),minor_breaks=seq(0,0.6,0.2),expand=expansion(add=c(0.05,0.05),mult=c(0,0)))
      ))+
      force_panelsizes(cols=c(1,1.6))
    pp<-aplt
    jpeg(paste0("F:/",scope,"_data/Code files/Figures/FIPDP_Areav2_",clr,'_',ac,".jpg"),units="cm",width=22,height=15,res=800);plot(pp);dev.off()
    
    
    
    #Suppl. Fig. 11
    av<-adp2$Sgn=='23' & adp2$Sgn2=='-1'
    adp2$Sgn[av]<-25;adp2$SpcSgn[av]<-'Fish Negative'
    av<-adp2$Sgn=='23' & adp2$Sgn2=='1'
    adp2$Sgn[av]<-24;adp2$SpcSgn[av]<-'Seaweed Positive'
    lgnlbl<-c('Seaweed Negative','Seaweed Positive','Fish Negative','Fish Positive')
    aplt<-ggplot()+
      geom_point(aes(x=Coef,y=reorder(Var,-Ord,mean),color=SpcSgn,fill=SpcSgn,shape=SpcSgn),data=adp2,size=4,stroke=0,alpha=1.0,position=position_dodge(width=1.))+
      #General
      theme(text=element_text(size=14,family="Times New Roman",color='black'))+
      theme(panel.background = element_rect(fill="White",colour="black"),
            panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
      theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
            panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
      #Legend
      theme(legend.position="bottom",legend.title=element_blank(),
            legend.text=element_text(size=14,color='black'),
            legend.key = element_rect(color=NA,fill=NA),
            #legend.key.height=unit(0.75,"cm"),
            #legend.key.size=unit(1,"cm"),
            #legend.key.width=unit(1,"cm"),
            legend.spacing.y=unit(0.4,'cm'),
            plot.title=element_text(size=14,face='bold',color='black',margin=margin(t=0,r=0,b=10,l=0)))+
      #guides(guide=guide_legend(reverse=T))+
      #labs(title="Farm growth extent")+
      #Axis
      #Margins from top, clockwise
      xlab(expression("Penalized coefficient"))+
      ylab(expression("Covariate"))+
      theme(axis.title.x=element_text(size=14,margin=ggplot2::margin(t=10,r=0,b=0,l=0)),
            axis.title.y=element_text(size=14,margin=ggplot2::margin(t=0,r=20,b=0,l=0)))+
      #scale_x_continuous(limits=c(-5,35),breaks=seq(-5,35,5),minor_breaks=seq(-5,35,2.5),expand=expansion(add=c(0,0),mult=c(0,0)))+ #,position='bottom',
      #scale_x_discrete(limits=rev)+
      scale_y_discrete(limits=rev,position="left")+
      #scale_y_log10()+
      theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+
      theme(axis.text.x=element_text(size=12,color=ac,angle=0))+
      theme(axis.text.y=element_text(size=12,color="black"))+
      #theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
      theme(axis.ticks.length=unit(0,"cm"))+
      #Others
      scale_color_manual(values=brewer.pal(9,'Set1')[c(2,2,5,5)],labels=lgnlbl)+
      scale_fill_manual(values=brewer.pal(9,'Set1')[c(2,2,5,5)],labels=lgnlbl)+
      scale_shape_manual(values=c(25,24,25,24),labels=lgnlbl)+
      guides(color=guide_legend(ncol=2),fill=guide_legend(ncol=2),shape=guide_legend(ncol=2))+
      facet_grid(cols=vars(SuperMdl),rows=vars(Cat),scales="free",space='free')+
      theme(strip.background=element_rect(fill="NA"),
            strip.text.x=element_blank(),
            strip.text.y=element_text(size=14,color="black",angle=0))
    pp<-aplt
    jpeg(paste0("F:/",scope,"_data/Code files/Figures/PenalizedCoefficients_v2.jpg"),units="cm",width=24,height=24,res=800);plot(pp);dev.off()
    
  }
  
  
  
  #Errors
  sm0<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/GLMM_AR_','FarmArea','_X','.rds'))
  for (i in 1:length(sm0)){
    if (scope=='Global' & is.element(i,c(1,3))){next()}
    print(c(names(sm0)[i],round(sm0[[i]]$MSE,4)))
    #print(c(names(sm0)[i],sm0[[i]]$Error)) #Mean absolute error
  }
  
}

#Prepare dataset for difference model
{
  #Wild fisheries
  {
    library(dplyr)
    
    am1<-read.csv("F:/Global_data/FAO/FI_Trade_2023.1.0/TRADE_VALUE.csv")
    am2<-read.csv("F:/Global_data/FAO/FI_Trade_2023.1.0/TRADE_QUANTITY.csv")
    am3<-read.table("F:/Global_data/FAO/FI_Trade_2023.1.0/CL_FI_COUNTRY_GROUPS.csv",sep=',',header=T) #Countries and regions
    am4<-read.table("F:/Global_data/FAO/GlobalProduction_2023.1.1/CL_FI_SPECIES_GROUPS.csv",sep=',',header=T) #Species groups
    am5<-read.csv("F:/Global_data/FAO/FI_Trade_2023.1.0/CL_FI_COMMODITY_ISSCFC.csv",sep=',',header=T) #Species groups
    
    am6<-read.table("F:/Global_data/FAO/GlobalProduction_2023.1.1/Global_production_quantity.csv",sep=',',header=T)
    am6<-left_join(am6,am3[,c('UN_Code','Name_En')],join_by(COUNTRY.UN_CODE==UN_Code))
    am6<-left_join(am6,am4[,c('X3A_Code','Yearbook_Group_En','Major_Group')],join_by(SPECIES.ALPHA_3_CODE==X3A_Code))
    #table(am6$Yearbook_Group_En)
    #table(am6$MEASURE)
    am6<-am6[am6$PRODUCTION_SOURCE_DET.CODE=='CAPTURE',] #Only wild fisheries
    am6<-am6[am6$MEASURE=='Q_tlw',] #Excluding species measured in units (usually mammals and some amphibians)
    ############################am6<-am6[is.element(am6$PERIOD,c('2015','2019')),] #Only periods of interest
    am6$Yearbook_Group_En[!is.element(am6$Yearbook_Group_En,c("Fish, crustaceans and molluscs, etc.","Other aq. animals & products"))]<-'PL'
    am6$Yearbook_Group_En[is.element(am6$Yearbook_Group_En,c("Fish, crustaceans and molluscs, etc.","Other aq. animals & products"))]<-'AN'
    names(am6)[names(am6)=='PERIOD']<-'Period'
    names(am6)[names(am6)=='Name_En']<-'Region'
    names(am6)[names(am6)=='Yearbook_Group_En']<-'ProductCategory'
    names(am6)[names(am6)=='VALUE']<-'Val'
    
    am0<-HarmFAOtoSHP(am6)
    ############################am0<-data.frame(am0[,c('Period','Region','ProductCategory','Val')] %>% group_by(Period,Region,ProductCategory) %>% summarise(Val=sum(Val))) #n = n()
    am0<-data.frame(am0[,c('Period','Region','ProductCategory','Val')] %>% group_by(Period,Region) %>% summarise(Val=sum(Val))) #n = n()
    names(am0)[3]<-'Fsh'
    amFs<-am0
    
  }
  #Aquaculture
  {
    library(dplyr)
    
    am1<-read.table('F:/Global_data/FAO/Aquaculture_2023.1.1/Aquaculture_Quantity.csv',sep=',',header=T) #Main data
    am2<-read.table("F:/Global_data/FAO/Aquaculture_2023.1.1/CL_FI_COUNTRY_GROUPS.csv",sep=',',header=T) #Countries and regions
    am3<-read.table("F:/Global_data/FAO/Aquaculture_2023.1.1/CL_FI_SPECIES_GROUPS.csv",sep=',',header=T) #Species groups
    am4<-read.table("F:/Global_data/FAO/Aquaculture_2023.1.1/CL_FI_SYMBOL.csv",sep=',',header=T) #Source and processing
    
    #Preprocessing
    am0<-am1
    #am0<-am0[am0$ENVIRONMENT.ALPHA_2_CODE=='MA' & is.element(am0$PERIOD,c(2015,2019,2020)),]############
    print(unique(am0$MEASURE))
    ad<-am3[,c('X3A_Code','Yearbook_Group_En')];ad<-ad[!duplicated(ad),] #'Major_Group'
    which(!is.element(unique(am0$SPECIES.ALPHA_3_CODE),unique(ad$X3A_Code)))
    for (i in unique(am0$SPECIES.ALPHA_3_CODE)){ 
      am0$SPECIES.ALPHA_3_CODE[am0$SPECIES.ALPHA_3_CODE==i]<-ad$Yearbook_Group_En[ad$X3A_Code==i]
    }
    am0$SPECIES.ALPHA_3_CODE[is.element(am0$SPECIES.ALPHA_3_CODE,c('Fish, crustaceans and molluscs, etc.','Other aq. animals & products'))]<-'AN'
    am0$SPECIES.ALPHA_3_CODE[is.element(am0$SPECIES.ALPHA_3_CODE,c('Aquatic plants'))]<-'PL'
    am0<-am0[,-c(3)];colnames(am0)<-c('Region','ProductGroup','WaterEnv','Unit','Period','Prod_FAO','Prod_Nts');am0<-am0[,c(1,5,3,2,6,4,7)]
    am0$Prod_Nts[am0$Prod_Nts=='']<-'F' #Official value
    ############################am0<- am0 %>% group_by(Region,ProductGroup,WaterEnv,Unit,Period) %>% summarise(Prod_FAO=sum(Prod_FAO),Prod_Nts=paste(Prod_Nts,collapse='')) #n = n()
    am0<- am0 %>% group_by(Region,WaterEnv,Period) %>% summarise(Prod_FAO=sum(Prod_FAO),Prod_Nts=paste(Prod_Nts,collapse='')) #n = n()
    av<-left_join(am0,data.frame(Region=am2$UN_Code,Val=am2$Name_En),'Region');am0$Region<-av$Val
    
    am0<-HarmFAOtoSHP(am0)
    am0<-am0[am0$WaterEnv=='MA',]
    am0<-data.frame(am0)[,c('Region','Period','Prod_FAO')]
    names(am0)[3]<-'Aq'
    amAq<-am0
  }
  #Trade
  {
    library(tidyr)
    library(dplyr)
    
    am<-read.csv("D:/Global_data/Fisheries and aquaculture/OECD.TAD.ATM,DSD_AGR@DF_OUTLOOK_2024_2033,1.1+all.csv")
    av<-unique(am$Commodity);av<-av[grepl('Fish|fish',av,ignore.case=T,perl=T)];av1<-av[-1]
    av<-unique(am$Measure);av2<-av[is.element(av,c("Production","Exports","Imports","Consumption"))]
    am<-am[is.element(am$Commodity,av1) & is.element(am$Measure,av2),]
    
    av<-data.frame(am %>% group_by(TIME_PERIOD,Reference.area,Measure) %>% summarise(Val=sum(OBS_VALUE)))
    av<-data.frame(pivot_wider(av,id_cols=names(av)[1:2],names_from='Measure',values_from='Val'))
    av<-data.frame(av %>% mutate(CnsRat=Consumption/Production,ExpRat=Exports/Production,ImpRat=Imports/Production))
    am<-av
    names(am)[1:2]<-c('Period','Region')
    
    amTr<-am
  }
  #WB
  {
    am<-read.table("F:/Global_data/Various WorldBank/WDI_CSV/WDICSV.csv",sep=',',header=T)
    CstReg<-sort(unique(readRDS('F:/Global_data/Code files/Datasets/FAO_MAqProd.rds')$Region))
    #sort(unique(am$Country.Name))[!is.element(sort(unique(am$Country.Name)),CstReg)]
    v1<-CstReg[!is.element(CstReg,sort(unique(am$Country.Name)))]
    v2<-unique(am$Country.Name)
    # an1<-sapply(v1,function(x){u<-adist(x,v2);v2[which(u==min(u))]})
    # an2<-sapply(v1,function(x){u<-adist(x,v2);min(u)})
    # View(cbind(v1,an1,an2))
    va1<-v1[c( 1,  2,  6,  7,  8, 11, 12, 15, 16, 17, 18, 21, 22, 25, 26)]
    va2<-v2[c(63,136,141,153,154,178,189,231,232,234,182,241,248,256,260)] #k<-'China';which(grepl(k,v2,perl=T));v2[which(grepl(k,v2,fixed=T))]
    # View(cbind(v1[va1],v2[va2]))
    # View(cbind(va1,va2))
    for (i in 1:length(va1)){
      am$Country.Name[am$Country.Name==va2[i]]<-va1[i]
    }
    #FAO data excludes Macao, Antilles, and other small states
    #Taiwan missing data in WDI
    
    #Search for indicators
    x<-unique(am$Indicator.Name);k<-'Trade';x[(grepl(k,x,perl=T))]
    
    ac<-c(
      #Sociodemographic
      'Population, total',
      #' #'Population, female', 
      #' 'Population, female (% of total population)',
      #' 'Population ages 15-64 (% of total population)',
      #' 'Educational attainment, at least completed upper secondary, population 25+, total (%) (cumulative)',
      #' #'Labor force with intermediate education (% of total working-age population with intermediate education)'
      #' #'Households and NPISHs Final consumption expenditure (current US$)',
      #' #'Final consumption expenditure (current US$)',
      #' #Technoeconomic
      'GDP (current US$)'
      #' #'Labor force, total',
      #' 'Employment in agriculture (% of total employment) (modeled ILO estimate)',
      #' 'Employment in industry (% of total employment) (modeled ILO estimate)',
      #' 'Employment in services (% of total employment) (modeled ILO estimate)',
      #' #Institutional
      #' 'Ease of doing business score (0 = lowest performance to 100 = best performance)',
      #' # 'Start-up procedures to register a business (number)',
      #' # 'Time required to start a business (days)',
      #' # 'Time to export, documentary compliance (hours)',
      #' # 'Domestic credit to private sector (% of GDP)',
      #' # 'Lending interest rate (%)',
      #' # Others
      #' 'Firms using banks to finance investment (% of firms)',
      #' #Too general
      #' #'Agriculture, forestry, and fishing, value added (current US$)',
      #' #'Travel services (% of service exports, BoP)',
      #' #'International tourism, receipts (current US$)',
      #' #'Transport services (% of commercial service imports)',
      #' #'Transport services (% of service exports, BoP)',
      #' #'Service exports (BoP, current US$)',
      #' #'Air transport, freight (million ton-km)',
      #' 'Container port traffic (TEU: 20 foot equivalent units)',
      #' #'Air transport, registered carrier departures worldwide',
      #' #'Food imports (% of merchandise imports)',
      #' #'Arable land (% of land area)',
      #' #'Land area (sq. km)',
      #' #'Rural land area (sq. km)',
      #' #'Urban land area (sq. km)',
      #' #'Surface area (sq. km)',
      #' #'Agricultural land (sq. km)'
      #' #'Merchandise imports (current US$)'
      #' #'Merchandise exports (current US$)',
      #' #'Agricultural raw materials exports (% of merchandise exports)',
      #' 'Marine protected areas (% of territorial waters)'
    )
    am<-am[is.element(am$Indicator.Name,ac),]
    #am<-HarmWBtoSHP(am) #See below
    {
      avR<-vect("F:/Global_data/Administrative boundaries/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp")
      avR<-SetSHPLabels(avR)
      amR<-read.csv("F:/Global_data/Various WorldBank/HDR23-24_Composite_indices_complete_time_series.csv")
      
      #x<-unique(avc$Territory);k<-'Bahamas';which(grepl(k,x,perl=T,ignore.case=T));x[(grepl(k,x,perl=T,ignore.case=T))]
      va1<-unique(amR$country)[!is.element(unique(amR$country),unique(avR$Territory))][1:28]
      va2<-unique(avR$Territory)[c(267,97,68,162,188,235,106,259,249,NA,42,30,115,250,272,154,48,251,263,NA,227,81,12,127,77,93,181,24)]
      #View(cbind(va1,chr(va2)))
      
      for (i in 1:length(va1)){am$Country.Name[am$Country.Name==va1[i]]<-va2[i]}
    }
    names(am)[c(1,3)]<-c('Region','Var')
    
    am<-am[!is.na(am$Region),]
    am<-data.frame(pivot_longer(am,cols=names(am)[5:length(am)],names_to='Period',values_to='Val'))
    am$Period<-gsub('X','',am$Period,fixed=T)
    am$Period<-nmr(am$Period)
    am$Var[am$Var=="GDP (current US$)"]<-'GDP'
    am$Var[am$Var=="Population, total"]<-"Pop"
    am<-am[,c('Region','Period','Var','Val')]
    am<-pivot_wider(am,id_cols=names(am)[1:2],names_from='Var',values_from='Val')
    
    amWB<-data.frame(am)
  }
  #Check country names (based on FAO fisheries)
  {
    r0<-unique(amFs$Region)
    k<-amAq;unique(k$Region[which(!is.element(k$Region,r0))])
    k<-amTr;unique(k$Region[which(!is.element(k$Region,r0))])
    amTr$Region[amTr$Region=="China (People’s Republic of)"]<-"China"
    amTr$Region[amTr$Region=="Korea"]<-"South Korea"
    amTr$Region[amTr$Region=="Türkiye"]<-"Turkey"
    amTr$Region[amTr$Region=="Viet Nam"]<-"Vietnam"
    k<-amWB;unique(k$Region[which(!is.element(k$Region,r0))])
    k<-am0;k$Region[which(!is.element(k$Region,r0))]
    #r0[grep('Cu',r0,ignore.case=T,fixed=T)]
    ad<-data.frame(t(matrix(c(
      "Channel Islands",NA,
      "Congo, Dem. Rep.","Democratic Republic of the Congo",
      "Congo, Rep.","Republic of the Congo",
      "Cote d'Ivoire","Ivory Coast",
      "Curacao","Curaçao",
      "Egypt, Arab Rep.","Egypt",
      "Eswatini",NA,
      "Faroe Islands","Faeroe",
      "Gambia, The","Gambia",
      "Gibraltar",NA, #"Gibraltar_Spain"
      "China, Hong Kong SAR",NA,
      "Iran (Islamic Rep. of)","Iran",
      "Isle of Man",NA,
      "Korea, Dem. People's Rep","North Korea",
      "Korea, Republic of","South Korea",
      "Kosovo",NA,
      "Kyrgyz Republic","Kyrgyzstan",
      "Lao PDR","Laos",
      "Macao SAR, China",NA,
      "Micronesia (Fed. States)","Micronesia",
      "Serbia and Montenegro","Montenegro",
      "Netherlands (Kingdom of the)","Netherlands",
      "Sint Maarten (Dutch part)","Sint-Maarten",
      "Slovak Republic","Slovakia",
      "St. Martin (French part)","Collectivity of Saint Martin",
      "Saint Vincent/Grenadines","Saint Vincent and the Grenadines",
      "Tanzania, United Rep. of","Tanzania",
      "Türkiye","Turkey",
      "United States of America","United States",
      "Venezuela (Boliv Rep of)","Venezuela",
      "Virgin Islands (U.S.)","United States Virgin Islands",
      "Yemen, Rep.","Yemen"),2,32)))
    for (i in 1:nrow(ad)){amWB$Region[amWB$Region==ad$X1[i]]<-ad$X2[i]}
  }
  #Calculate observations for larger areas
  GetCntrGroups<-function(am,p,gg){
    
    avc0<-vect(paste0('F:/Global_data/Code files/Administrative boundaries/RegionBufferv2_',gs00,'km.shp'))
    cnt0<-centroids(avc0)
    avx<-avc0
    
    #UN countries
    UNL<-c('Austria','Belgium','Bulgaria','Croatia','Cyprus','Czech Republic','Denmark','Estonia','Finland','France',
           'Germany','Greece','Hungary','Ireland','Italy','Latvia','Lithuania','Luxembourg','Malta','Netherlands',
           'Poland','Portugal','Romania','Slovakia','Slovenia','Spain','Sweden')
    library(terra)
    library(sf)
    av<-vect(t(matrix(c(-11.415347026814477,37.16233965175856,
                        10.531579008688489,39.25839438548637,
                        14.230499127031685,34.69639290619643,
                        26.930124866676657,33.710014207971575,
                        31.12223433413228,45.42326124939169,
                        25.573854156617486,66.63040326122602,
                        -12.155131050483114,56.150129592586964),2,7)),crs='epsg:4326')
    av1<-vect(st_cast(st_combine(st_as_sf(av)),'POLYGON'))
    av<-vect(t(matrix(c(-18.68988992622276,36.1759609535337,
                        7.819037588570147,39.38169172276447,
                        14.230499127031685,33.710014207971575,
                        32.231910369635244,32.230446160634294,
                        35.06774912703169,26.312173971285176,
                        46.041212144783174,12.25627752158104,
                        65.2755967601678,-21.034003543507737,
                        32.97169439330388,-37.43254940149589,
                        -1.3049653700097394,-37.30925206421779,
                        -7.469832233915078,-17.211786087886395,
                        -29.170163594861826,15.831900302646162,
                        -20.78594465995058,33.21682485885918),2,12)),crs='epsg:4326')
    av2<-vect(st_cast(st_combine(st_as_sf(av)),'POLYGON'))
    # plot(avc0,ext=ext(-40,80,-90,90))
    # plot(av1,col=rgb(1,1,0,0.5),add=T)
    # plot(av2,col=rgb(1,0,1,0.5),add=T)
    
    ll<-vector('list',3)
    ll[[1]]<-avx$Region[is.related(cnt0,av1,"within") & is.element(avx$Region,UNL)]
    ll[[2]]<-avx$Region[is.related(cnt0,av1,"within") & !is.element(avx$Region,UNL)]
    ll[[3]]<-avx$Region[is.related(cnt0,av2,"within") & !is.element(avx$Region,UNL)]
    ss<-c("European Union","Europe","Africa")
    
    am<-am[,c("Period","Region",p)]
    #am$ArType<-'Nrm'
    for (j in gg){
      an<-is.element(am$Region,ll[[j]])
      if (!any(an)==T){print(j);next()}
      av<-am[an,] #Get countries
      aa<-data.frame(av %>% group_by(Period) %>% summarise(Var=sum(!!sym(p),na.rm=T))) #Summarize and create new observation (remove NAs as approximation, NAs assumed to be low values)
      aa$Region<-ss[j]
      names(aa)<-c("Period",p,"Region");aa<-aa[,names(am)[1:3]]
      ##
      am<-rbind(am[!an,],aa)
      ##
      #aa$ArType<-'Agg'
      #am<-rbind(am,aa)
      ##
    }
    am$Var<-p
    #names(am)<-c('Period','Region','Val','ArType','Var')
    names(am)<-c('Period','Region','Val','Var')
    
    am
    
  }  #gg "European Union","Europe","Africa"
  #unique(amTr$Region)[which(!is.element(unique(amTr$Region),unique(pp$Region)))]
  {
    am1<-GetCntrGroups(amWB,'Pop',2:3)
    am2<-GetCntrGroups(amWB,'GDP',2:3)
    am3<-GetCntrGroups(amFs,'Fsh',1:3)
    am4<-GetCntrGroups(amAq,'Aq',1:3)
    
    am<-rbind(am1,am2);am<-pivot_wider(am[!is.na(am$Region),],id_cols=names(am)[c(1:2)],names_from='Var',values_from='Val')
    amWB<-data.frame(am)
    am<-am3;am<-pivot_wider(am[!is.na(am$Region),],id_cols=names(am)[c(1:2)],names_from='Var',values_from='Val')
    amFs<-data.frame(am)
    am<-am4
    #av<-which(!is.na(am$Region) & am$Region=='Montenegro' & am$Val<100);am<-am[-av,]
    am<-pivot_wider(am[!is.na(am$Region),],id_cols=names(am)[c(1:2)],names_from='Var',values_from='Val')
    amAq<-data.frame(am)
  }
  #Get statistics dataset
  {
    #Combine data
    am0<-amTr
    am0<-left_join(am0,amFs,join_by(Region,Period))
    am0<-left_join(am0,amAq,join_by(Region,Period))
    am0<-left_join(am0,amWB[!is.na(amWB$Region),],join_by(Region,Period))
    #unique(amTr$Region)[which(!is.element(unique(amTr$Region),unique(am0$Region)))]
    am0<-am0[,c('Period','Region','Pop','GDP','Fsh','Consumption','Exports','Imports','Production','Aq')]
    
    #Add difference variables
    am<-am0#[,!grepl('Aq',names(am0),fixed=T)]
    for (i in names(am)[3:length(am)]){
      am[,paste0(i,'_D')]<-NA
      for (rg in unique(am$Region)){
        for (yr in unique(am$Period)){
          av<-am[am$Region==rg & am$Period==yr+5,paste0(i)]-am[am$Region==rg & am$Period==yr,paste0(i)]
          if (length(av)>0){am[am$Region==rg & am$Period==yr,paste0(i,'_D')]<-av}
        }
      }
    }
    #av<-complete.cases(am);am<-am[av,] #table(av)
    am$Region<-fct(am$Region)
    am$Period<-fct(am$Period)
    am0<-am
    
  }
  saveRDS(am0,"F:/Global_data/Code files/Various/TradeStats.rds")
  
  
  
  af<-function(qq,sb,fx,astr,scn=NULL){
    if (is.null(scn)){scn<-'0'}
    #print(length(astr))
    if (length(astr)==1){ar<-rast(astr)}
    if (length(astr)>1){
      ar<-vector('list',length(astr))
      for (l in 1:length(ar)){ar[[l]]<-rast(astr[l])}
      ar<-app(rast(ar),'sum',na.rm=T)
    }
    ar<-resample(ar,sb,'bilinear');names(ar)<-'Val'
    av<-zonal(ar,sb,fx,na.rm=T,as.raster=F)
    av$Period<-yr;av$Var<-k;av$Scn<-nmr(scn)
    qq<-rows_update(qq,av,c('Region','Period','Var','Scn'),unmatched="ignore")
    qq
  } #,k
  #Calculate values from rasters
  {
    sb<-rast(paste0("F:/Global_data/Code files/Administrative boundaries/MarZones_Global_","0.1",".tif"))
    qq1<-expand.grid(Period=c(2015,2019),Region=unique(sb$Region)[,1],Scn=nmr(c(0)),Var=c('Pop','GDP','Fsh','Aq',paste0(c('Pop','GDP','Fsh','Aq'),'_D')),Val=nmr(NA))
    qq2<-expand.grid(Period=seq(2025,2050,5),Region=unique(sb$Region)[,1],Scn=nmr(c(2,4,5)),Var=c('Pop','GDP','Fsh','Aq',paste0(c('Pop','GDP','Fsh','Aq'),'_D')),Val=nmr(NA))
    qq<-rbind(qq1,qq2)
    
    k<-'Pop'
    for (yr in c(2015,2019)){qq<-af(qq,sb,'sum',paste0('F:/Global_data/Population and others/ppp_',yr,'_1km_Aggregated.tif'))}
    for (yr in seq(2025,2050,5)){for (SCN in c(2,4,5)){qq<-af(qq,sb,'sum',paste0('D:/China_data/Population, Weather, and Projections/Population/Global/SSP/SSP',SCN,'_',yr,'.tif'),SCN)}}
    
    k<-'GDP'
    for (yr in c(2015,2019)){qq<-af(qq,sb,'sum',paste0('F:/Global_data/Code files/GDP/GDPC_Global_',yr,'.tif'))}
    for (yr in seq(2025,2050,5)){for (SCN in c(2,4,5)){qq<-af(qq,sb,'sum',paste0('D:/Global_data/Code files/GDPC/GDPC_Global_SSP',SCN,'_',yr,'.tif'),SCN)}}
    
    k<-'Fsh'
    for (yr in c(2015,2019)){qq<-af(qq,sb,'mean',c(paste0('F:/Global_data/Code files/Various/WldFshNonNorm_',yr,'_AN_10km.tif'),paste0('F:/Global_data/Code files/Various/WldFshNonNorm_',yr,'_PL_10km.tif')))}
    for (yr in seq(2025,2050,5)){for (SCN in c(2,4,5)){qq<-af(qq,sb,'mean',paste0("D:/Global_data/Code files/Wild fisheries/ProjectionsNonNorm_",SCN,"_",yr,".tif"),SCN)}}
    
    k<-'Aq'
    for (yr in c(2015,2019)){qq<-af(qq,sb,'sum',c(paste0('F:/Global_data/Code files/Offshore aquaculture/OffshrAq_Area_',yr,'_AN_0.1.tif'),paste0('F:/Global_data/Code files/Offshore aquaculture/OffshrAq_Area_',yr,'_PL_0.1.tif')))}
    
    QQ<-qq
  }
  #Recalculate observations for larger areas (prediction data)
  {
    {
      # ashp<-vect("F:/Global_data/Administrative boundaries/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp")
      # ashp<-SetSHPLabels(ashp)
      # avc0<-ashp
      # cnt0<-centroids(avc0)
      # avx<-avc0
      # 
      # #UN countries
      # UNL<-c('Austria','Belgium','Bulgaria','Croatia','Cyprus','Czech Republic','Denmark','Estonia','Finland','France',
      #        'Germany','Greece','Hungary','Ireland','Italy','Latvia','Lithuania','Luxembourg','Malta','Netherlands',
      #        'Poland','Portugal','Romania','Slovakia','Slovenia','Spain','Sweden')
      # library(terra)
      # library(sf)
      # av<-vect(t(matrix(c(-11.415347026814477,37.16233965175856,
      #                     10.531579008688489,39.25839438548637,
      #                     14.230499127031685,34.69639290619643,
      #                     26.930124866676657,33.710014207971575,
      #                     31.12223433413228,45.42326124939169,
      #                     25.573854156617486,66.63040326122602,
      #                     -12.155131050483114,56.150129592586964),2,7)),crs='epsg:4326')
      # av1<-vect(st_cast(st_combine(st_as_sf(av)),'POLYGON'))
      # av<-vect(t(matrix(c(-18.68988992622276,36.1759609535337,
      #                     7.819037588570147,39.38169172276447,
      #                     14.230499127031685,33.710014207971575,
      #                     32.231910369635244,32.230446160634294,
      #                     35.06774912703169,26.312173971285176,
      #                     46.041212144783174,12.25627752158104,
      #                     65.2755967601678,-21.034003543507737,
      #                     32.97169439330388,-37.43254940149589,
      #                     -1.3049653700097394,-37.30925206421779,
      #                     -7.469832233915078,-17.211786087886395,
      #                     -29.170163594861826,15.831900302646162,
      #                     -20.78594465995058,33.21682485885918),2,12)),crs='epsg:4326')
      # av2<-vect(st_cast(st_combine(st_as_sf(av)),'POLYGON'))
      # # plot(avc0,ext=ext(-40,80,-90,90))
      # # plot(av1,col=rgb(1,1,0,0.5),add=T)
      # # plot(av2,col=rgb(1,0,1,0.5),add=T)
      # 
      # ll<-vector('list',3)
      # ll[[1]]<-avx$Region[is.related(cnt0,av1,"within") & is.element(avx$Region,UNL)]
      # ll[[2]]<-avx$Region[is.related(cnt0,av1,"within") & !is.element(avx$Region,UNL)]
      # ll[[3]]<-avx$Region[is.related(cnt0,av2,"within") & !is.element(avx$Region,UNL)]
      # ss<-c("European Union","Europe","Africa")
      # 
      # pl1<-avx[is.element(avx$Region,ll[[1]]),];pl1$Region<-ss[1]
      # pl2<-avx[is.element(avx$Region,ll[[2]]),];pl2$Region<-ss[2]
      # pl3<-avx[is.element(avx$Region,ll[[3]]),];pl3$Region<-ss[3]
      # 
      # pl<-vect(list(pl1,pl2,pl3)) #Combine
      # pl0<-aggregate(pl,'Region') #Dissolve
      # sb<-rasterize(pl0,ar0,'Region') #plot(sb)
      
    } #Get polygons for large areas
    
    k<-'Pop'
    # for (yr in c(2015,2019)){qq<-af(qq,sb,'sum',paste0('F:/Global_data/Population and others/ppp_',yr,'_1km_Aggregated.tif'))}
    # for (yr in seq(2025,2050,5)){for (SCN in c(2,4,5)){qq<-af(qq,sb,'sum',paste0('D:/China_data/Population, Weather, and Projections/Population/Global/SSP/SSP',SCN,'_',yr,'.tif'),SCN)}}
    am<-vector('list',2);names(am)<-c(2015,2019);for (yr in chr(c(2015,2019))){
      am[[yr]]<-data.frame(pivot_wider(qq[qq$Period==nmr(yr),],id_cols=names(qq)[1:3],names_from='Var',values_from='Val'))
      am[[yr]]<-GetCntrGroups(am[[yr]],k,c(1:3));am[[yr]]$Scn<-0
    };am1<-bind_rows(am)
    am<-vector('list',3);names(am)<-c(2,4,5);for (i in 1:length(am)){am[[i]]<-vector('list',6);names(am[[i]])<-seq(2025,2050,5)}
    for (yr in chr(seq(2025,2050,5))){for (SCN in chr(c(2,4,5))){
      am[[SCN]][[yr]]<-data.frame(pivot_wider(qq[qq$Period==nmr(yr) & qq$Scn==SCN,],id_cols=names(qq)[1:3],names_from='Var',values_from='Val'))
      am[[SCN]][[yr]]<-GetCntrGroups(am[[SCN]][[yr]],k,c(1:3));am[[SCN]][[yr]]$Scn<-nmr(SCN)
    }}
    for (SCN in chr(c(2,4,5))){am[[SCN]]<-bind_rows(am[[SCN]])};am2<-bind_rows(am)
    amV1<-bind_rows(am1,am2)
    #av<-which(qq$Var==k);qq<-qq[-av,];qq<-bind_rows(qq,am[,names(qq)])
    
    k<-'GDP'
    # for (yr in c(2015,2019)){qq<-af(qq,sb,'sum',paste0('F:/Global_data/Code files/GDP/GDPC_Global_',yr,'.tif'))}
    # for (yr in seq(2025,2050,5)){for (SCN in c(2,4,5)){qq<-af(qq,sb,'sum',paste0('D:/Global_data/Code files/GDPC/GDPC_Global_SSP',SCN,'_',yr,'.tif'),SCN)}}
    am<-vector('list',2);names(am)<-c(2015,2019);for (yr in chr(c(2015,2019))){
      am[[yr]]<-data.frame(pivot_wider(qq[qq$Period==nmr(yr),],id_cols=names(qq)[1:3],names_from='Var',values_from='Val'))
      am[[yr]]<-GetCntrGroups(am[[yr]],k,c(1:3));am[[yr]]$Scn<-0
    };am1<-bind_rows(am)
    am<-vector('list',3);names(am)<-c(2,4,5);for (i in 1:length(am)){am[[i]]<-vector('list',6);names(am[[i]])<-seq(2025,2050,5)}
    for (yr in chr(seq(2025,2050,5))){for (SCN in chr(c(2,4,5))){
      am[[SCN]][[yr]]<-data.frame(pivot_wider(qq[qq$Period==nmr(yr) & qq$Scn==SCN,],id_cols=names(qq)[1:3],names_from='Var',values_from='Val'))
      am[[SCN]][[yr]]<-GetCntrGroups(am[[SCN]][[yr]],k,c(1:3));am[[SCN]][[yr]]$Scn<-nmr(SCN)
    }}
    for (SCN in chr(c(2,4,5))){am[[SCN]]<-bind_rows(am[[SCN]])};am2<-bind_rows(am)
    amV2<-bind_rows(am1,am2)
    #av<-which(qq$Var==k);qq<-qq[-av,];qq<-bind_rows(qq,am[,names(qq)])
    
    k<-'Fsh'
    # for (yr in c(2015,2019)){qq<-af(qq,sb,'mean',c(paste0('F:/Global_data/Code files/Various/WldFshNonNorm_',yr,'_AN_10km.tif'),paste0('F:/Global_data/Code files/Various/WldFshNonNorm_',yr,'_PL_10km.tif')))}
    # for (yr in seq(2025,2050,5)){for (SCN in c(2,4,5)){qq<-af(qq,sb,'mean',paste0("D:/Global_data/Code files/Wild fisheries/ProjectionsNonNorm_",SCN,"_",yr,".tif"),SCN)}}
    am<-vector('list',2);names(am)<-c(2015,2019);for (yr in chr(c(2015,2019))){
      am[[yr]]<-data.frame(pivot_wider(qq[qq$Period==nmr(yr),],id_cols=names(qq)[1:3],names_from='Var',values_from='Val'))
      am[[yr]]<-GetCntrGroups(am[[yr]],k,c(1:3));am[[yr]]$Scn<-0
    };am1<-bind_rows(am)
    am<-vector('list',3);names(am)<-c(2,4,5);for (i in 1:length(am)){am[[i]]<-vector('list',6);names(am[[i]])<-seq(2025,2050,5)}
    for (yr in chr(seq(2025,2050,5))){for (SCN in chr(c(2,4,5))){
      am[[SCN]][[yr]]<-data.frame(pivot_wider(qq[qq$Period==nmr(yr) & qq$Scn==SCN,],id_cols=names(qq)[1:3],names_from='Var',values_from='Val'))
      am[[SCN]][[yr]]<-GetCntrGroups(am[[SCN]][[yr]],k,c(1:3));am[[SCN]][[yr]]$Scn<-nmr(SCN)
    }}
    for (SCN in chr(c(2,4,5))){am[[SCN]]<-bind_rows(am[[SCN]])};am2<-bind_rows(am)
    amV3<-bind_rows(am1,am2)
    #av<-which(qq$Var==k);qq<-qq[-av,];qq<-bind_rows(qq,am[,names(qq)])
    
    k<-'Aq'
    #for (yr in c(2015,2019)){qq<-af(qq,sb,'sum',c(paste0('F:/Global_data/Code files/Offshore aquaculture/OffshrAq_Area_',yr,'_AN_0.1.tif'),paste0('F:/Global_data/Code files/Offshore aquaculture/OffshrAq_Area_',yr,'_PL_0.1.tif')))}
    am<-vector('list',2);names(am)<-c(2015,2019);for (yr in chr(c(2015,2019))){
      am[[yr]]<-data.frame(pivot_wider(qq[qq$Period==nmr(yr),],id_cols=names(qq)[1:3],names_from='Var',values_from='Val'))
      am[[yr]]<-GetCntrGroups(am[[yr]],k,c(1:3));am[[yr]]$Scn<-0
    };am1<-bind_rows(am)
    amV4<-am1
    #av<-which(qq$Var==k);qq<-qq[-av,];qq<-bind_rows(qq,am[,names(qq)])
    
    qq<-bind_rows(amV1,amV2,amV3,amV4)
    
  }
  saveRDS(QQ,"F:/Global_data/Code files/Various/TradeRst_Nrm.rds")
  saveRDS(qq,"F:/Global_data/Code files/Various/TradeRst_Agg.rds")
}
#Growth model (DF)
{
  library(tidyr)
  library(dplyr)
  am0<-readRDS("F:/Global_data/Code files/Various/TradeStats.rds")
  #View(am0[!complete.cases(am0),])
  av<-complete.cases(am0[,c('Period','Region','Pop','GDP','Fsh','Consumption','Exports','Imports','Production')]);am0<-am0[av,] #table(av)
  mm<-data.frame(pivot_longer(am0,names(am0)[3:length(am0)],names_to='Var',values_to='Val'))
  qq<-readRDS("F:/Global_data/Code files/Various/TradeRst_Agg.rds")
  #Scale statistics dataset (to baseline years of raster data)
  {
    # am1<-am0
    # am2<-qq[is.element(qq$Scn,c(0)),]
    # ##Difference variables for raster data
    # pp<-am2
    # pp<-data.frame(pivot_wider(pp,id_cols=names(pp)[1:3],names_from='Var',values_from='Val'))
    # pp$Period[pp$Period==2019]<-2020
    # for (i in c('Pop','GDP','Fsh','Aq')){
    #   for (rg in unique(pp$Region)){
    #     for (yr in unique(pp$Period)){
    #       av<-pp[pp$Region==rg & pp$Period==yr+5,paste0(i)]-pp[pp$Region==rg & pp$Period==yr,paste0(i)]
    #       if (length(av)>0){pp[pp$Region==rg & pp$Period==yr,paste0(i,'_D')]<-av}
    #     }
    #   }
    # }
    # am2<-pp
    
    # Old
    am<-am0
    yr<-2015 #Just one year for simplicity. 2019 also not directly available for 5-year differences.
    par(mfrow=c(3,2))
    for (k in c('Pop','GDP','Fsh')){
      ad1<-qq[qq$Period==yr & qq$Var==k,] #Raster values
      ad2<-mm[mm$Period==yr & mm$Var==k,];ad2$Period<-nmr(chr(ad2$Period)) #Statistic values
      ad<-left_join(ad1,ad2,join_by(Period,Region,Var))
      plot((ad$Val.x),(ad$Val.y),main=k)
      plot(lt(ad$Val.x),lt(ad$Val.y),main=k)
      
      rm<-lm(Val.x~Val.y,ad)
      kk<-summary(rm)$coefficients['Val.y','Estimate']
      print(kk)
      am[,k]<-am[,k]*kk
    }
    par(mfrow=c(1,1))
    
    # for (i in c('Pop','GDP','Fsh')){
    #   
    #   for (rg in unique(am2$Region)){
    #     
    #     k1<-am1[am1$Period==2010 & am1$Region==rg,paste0(i)]-am2[am2$Period==2015 & am2$Region==rg,paste0(i)] #Intercept
    #     k2<-mean(c(am1[am1$Period==2005 & am1$Region==rg,paste0(i,'_D')],am2[am2$Period==2015 & am2$Region==rg,paste0(i,'_D')])) #Slope
    #     
    #     if (length(k1)>0 & length(k2)>0){
    #       for (j in unique(am1$Period)){am1[am1$Period==j & am1$Region==rg,paste0(i)]<-am1[am1$Period==j & am1$Region==rg,paste0(i)]-k1+k2*1}
    #       am1[am1$Period==2010,paste0(i,'_D')]<-k2
    #     }
    #     
    #   }
    #   
    # }
    
    xx<-am
  } #xx
  #Scale raster dataset (to baseline years of raster data)
  {
    ee<-vector('list',3);names(ee)<-c(2,4,5)
    for (scn in c(2,4,5)){
      
      pp<-qq[is.element(qq$Scn,c(0,scn)),]
      pp<-pp[,c("Period","Region","Scn","Val","Var")]
      ##Difference variables for raster data
      pp<-data.frame(pivot_wider(pp,id_cols=names(pp)[1:3],names_from='Var',values_from='Val'))
      pp$Period[pp$Period==2019]<-2020
      for (i in c('Pop','GDP','Fsh','Aq')){
        for (rg in unique(pp$Region)){
          for (yr in unique(pp$Period)){
            av<-pp[pp$Region==rg & pp$Period==yr+5,paste0(i)]-pp[pp$Region==rg & pp$Period==yr,paste0(i)]
            if (length(av)>0){pp[pp$Region==rg & pp$Period==yr,paste0(i,'_D')]<-av}
          }
        }
      }
      ##
      
      #Scale
      for (i in c('Pop','GDP','Fsh')){
        
        k1<-pp[pp$Period==2025,paste0(i)]-pp[pp$Period==2020,paste0(i)] #Intercept
        k2<-rowMeans(cbind(pp[pp$Period==2015,paste0(i,'_D')],pp[pp$Period==2025,paste0(i,'_D')])) #Slope
        
        for (j in chr(seq(2025,2050,5))){pp[pp$Period==j,paste0(i)]<-pp[pp$Period==j,paste0(i)]-k1+k2*1}
        pp[pp$Period==2020,paste0(i,'_D')]<-k2
      }
      
      ee[[chr(scn)]]<-pp
      
    }
    
  } #ee
  #Build trade model
  {
    am<-xx
    am$Period<-nmr(chr(am$Period))
    av<-data.frame(pivot_wider(qq[qq$Var=='Aq',],id_cols=names(qq)[c(1,2,5)],names_from='Var',values_from='Val'));names(av)[4]<-'OA'
    am<-left_join(am,av[,c('Period','Region','OA')],join_by(Period,Region))
    ap<-ee
    
    set.seed(1111)
    
    #Transformations
    for (i in c('Pop','GDP','Fsh')){
      am[,i]<-lt(am[,i])
      for (scn in c(2,4,5)){ap[[chr(scn)]][,i]<-lt(ap[[chr(scn)]][,i])}
    }
    #for (i in paste0(c('Pop','GDP','Fsh','Aq'),'_D')){am[,i]<-lt(am[,i])}
    
    #Predict values
    #famlnk<-gaussian(link='identity') #family=tw(link='log')
    library(mgcv)
    #famlnk<-tw(link='log') #Exxagerates production
    famlnk<-gaussian(link='identity')
    for (vr in c("Production","Imports","Exports","Consumption")){
      
      #rfave<-formula(Imports~Pop+GDP+Fsh+Aq+Pop_D+GDP_D+Fsh_D+Aq_D+(1|Region)+(1|Period)) #mat(0+pos|group1) + ar1(0+Period|group2)
      #library(lme4)
      if (F){
        # rfave<-paste0(vr,'~',"
        #                 #s(Pop)+s(GDP)+s(Fsh)+
        #                 s(Pop,Region,bs='re')+s(GDP,Region,bs='re')+s(Fsh,Region,bs='re')+
        #                 #s(Pop_D)+s(GDP_D)+s(Fsh_D)+
        #                 s(Region,bs='re')+
        #                 #s(Period,Region,bs='re')
        #                 s(Period)
        #                 ") #mat(0+pos|group1) + ar1(0+Period|group2)
        # sm<-glmer(rfave,family=famlnk,data=am)
      }
      #library(DHARMa)
      if (F){
        rm<-lmer(rfave,am,REML=T)
        summary(rm)
        plot(rm)
        so<-simulateResiduals(rm,500) #recalculateResiduals
        plot(so)
      }
      #library(glmmTMB)
      if (F){
        famlnk<-gaussian(link='identity')
        rfzin<-formula(~0)
        rfdsp<-formula(~1)
        ctrl<-glmmTMBControl(optimizer=nlminb,optArgs=list(iter.max=1000),parallel=20,rank_check="adjust") #Default optimizer
        sm<-glmmTMB(formula=rfave,
                    family=famlnk,
                    ziformula=rfzin,
                    dispformula=rfdsp,
                    data=am,REML=F,verbose=F,control=ctrl)
        AIC(sm);BIC(sm);summary(sm);plot(sm)
      }
      
      library(mgcv)
      if (T){
        
        # library(mgcv)
        #G (can be used)
        # sm<-gam(Val~s(Period,k=3,bs="tp")+s(Period,PtID,k=3,bs="re"),
        #         family=gaussian('identity'),data=am,method="REML",unconditional=T)
        #GS (a bit slow)
        # sm<-gam(Val~s(Period,k=3,bs="tp",m=2)+s(Period,PtID,k=3,bs="fs",m=2),
        #         family=gaussian('identity'),data=am,method="REML",unconditional=T)
        #GI (crashes with many groups)
        # sm<-gam(Val~s(Period,k=3,m=2,bs="tp")+s(Period,by=PtID,k=3,m=1,bs="tp")+s(PtID,k=3,bs="re"),
        #         family=gaussian('identity'),data=am,method="REML",unconditional=T)
        #S (decent fit)
        # sm<-gam(Val~s(Period,PtID,k=3,m=2,bs="fs"),
        #         family=gaussian('identity'),data=am,method="REML",unconditional=T)
        #I (slow)
        # sm<-gam(Val~s(Period,by=PtID,k=3,m=1,bs="tp")+s(PtID,k=1,bs="re"),
        #         family=gaussian('identity'),data=am,method="REML",unconditional=T)
        
        #https://stats.stackexchange.com/questions/552880/by-group-random-effect-gam
        #Old (not terrible but improved by below)
        rfave<-paste0(vr,'~',"
                      #s(Pop)+s(GDP)+s(Fsh)+
                      s(Pop,Region,bs='re')+s(GDP,Region,bs='re')+s(Fsh,Region,bs='re')+
                      #s(Pop_D)+s(GDP_D)+s(Fsh_D)+
                      s(Region,bs='re')+
                      #s(Period,Region,bs='re')
                      s(Period)
                      ") #mat(0+pos|group1) + ar1(0+Period|group2)
        
        #Rev202604(better fit)
        # rfave<-paste0(vr,'~',"
        #                 #s(Pop)+s(GDP)+s(Fsh)+
        #                 s(Pop,Region,bs='re')+s(GDP,Region,bs='re')+s(Fsh,Region,bs='re')+
        #                 #s(Pop_D)+s(GDP_D)+s(Fsh_D)+
        #                 s(Region,bs='re')+
        #                 s(Period,Region,bs='re')
        #                 ") #mat(0+pos|group1) + ar1(0+Period|group2)
        rfave<-formula(rfave)
        sm<-gam(rfave,family=famlnk,data=am,unconditional=T)
        plot(sm,pages=1,residuals=T,main=vr,pch=16)  ## show partial residuals
        #k.check(sm)
        #par(mfrow=c(2,2));gam.check(sm);lines(0:10000,0:10000);par(mfrow=c(1,1))
        print(AIC(sm))#;BIC(sm)
        print(summary(sm))
        
        #Plot
        if (F){
          ll<-am
          nn<-quantile(am[,vr],c(0.25,0.50,0.75))
          # ll$Cat<-chr(NA)
          # for (i in 1:length(nn)){ll$Cat[ll[,vr]<nn[i] & is.na(ll$Cat)]<-chr(i)};ll$Cat[is.na(ll$Cat)]<-chr(i+1)
          # ll$Cat<-fct(ll$Cat)
          ll$Prediction<-predict(sm,type='response',se.fit=T)$fit
          library(ggplot2)
          library(cowplot)
          aplt<-vector('list',length(sort(unique(ll$Region))));names(aplt)<-sort(unique(ll$Region))
          for (i in sort(unique(ll$Region))){
            aplt[[i]]<-ggplot()+
              geom_line(aes(x=Period,y=!!sym(vr),color=Region),data=ll[ll$Region==i,])+
              geom_point(aes(x=Period,y=Prediction,color=Region),data=ll[ll$Region==i,])+
              labs(title=i)+
              theme(legend.position='none')
            #ggh4x::facet_grid2(rows=vars(Cat),scales='free')
          }
          cowplot::plot_grid(plotlist=aplt,nrow=4)
        }
      }
      
      for (scn in c(2,4,5)){
        #ap[[chr(scn)]][,vr]<-NA
        ap[[chr(scn)]][ap[[chr(scn)]]$Period>=2020,vr]<-predict(sm,ap[[chr(scn)]][ap[[chr(scn)]]$Period>=2020,],type="response",se.fit=F)
      }
      
    } #"OA" Satellite data for OA not good enough
    
    #u<-'China'
    #u<-'Africa' #China European Union Turkey United States
    #u<-'Brazil'
    #(am[am$Region==u,])
    #head(ap[[1]][ap[[1]]$Region==u,])
    
    aa<-am;aa$Scn<-0;aa<-aa[,names(ap[[1]])]
    aa<-data.frame(pivot_longer(aa,names(aa)[4:length(aa)],names_to='Var',values_to='Val'))
    bb<-ap
    for (scn in c(2,4,5)){
      bb[[chr(scn)]]<-bb[[chr(scn)]][bb[[chr(scn)]]$Period>=2020,]
      bb[[chr(scn)]]<-data.frame(pivot_longer(bb[[chr(scn)]],names(bb[[chr(scn)]])[4:length(bb[[chr(scn)]])],names_to='Var',values_to='Val'))
    }
    cc<-bind_rows(aa,bind_rows(bb))
    
    #Check results
    if (F){
      library(ggplot2)
      for (k in c("Production","Imports","Exports","Consumption")){
        print(
          ggplot(cc[cc$Var==k & !is.na(cc$Val),],aes(x=Period,y=Val,color=Region,fill=Region))+
            geom_line()+
            geom_point()+
            theme(legend.position='bottom')#+
          #ggrepel::geom_label_repel(aes(label=Region),nudge_x=1,na.rm=T)+
          #coord_trans(y="log1p")
          #scale_y_continuous(limits=c(0,200000))
        )
      }
      
      for (l in 1:length(unique(cc$Region))){
        #l<-1
        print(
          ggplot(cc[cc$Var==k & !is.na(cc$Val) & cc$Region==unique(cc$Region)[l],],aes(x=Period,y=Val,color=Scn,fill=Scn))+
            geom_line()+
            geom_point()+
            coord_trans(y="log1p")+
            labs(title=unique(cc$Region)[l]) 
        )
        Sys.sleep(2)
      }
      
    }
    
  }
  #Rasterize trade datasets
  {
    #Baseline
    xx<-readRDS("F:/Global_data/Code files/Various/TradeStats.rds")
    #No need to normalize trade variables to raster dataset. other variables not used here.
    xx<-xx[is.element(xx$Period,c(2015,2019)),]
    xx$Scn<-0
    
    #Projections
    dd<-bind_rows(bb)
    dd<-dd[dd$Period>=2025,]
    dd<-data.frame(pivot_wider(dd,id_cols=names(dd)[1:3],names_from='Var',values_from='Val'))
    
    nrml<-T;ac<-'Norm'
    #nrml<-F;ac<-'NonNorm'
    {
      #Use subregion classification
      #https://unstats.un.org/unsd/methodology/m49
      # am0<-openxlsx::read.xlsx("F:/Global_data/Administrative boundaries/UNSD — Methodology.xlsx")
      # avc<-vect("F:/Global_data/Administrative boundaries/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp")
      # avc<-SetSHPLabels(avc)
      # avc$ISOTer<-apply(data.frame(avc)[,c('ISO_TER1','ISO_TER2','ISO_TER3')],1,function(x){paste0(unique(x[!is.na(x)]),collapse='_')})
      # uu<-left_join(data.frame(avc[,c('Territory','ISO_TER1')]),am0[,c('ISO-alpha3.Code','Region.Name','Sub-region.Name','Country.or.Area')],join_by(ISO_TER1==`ISO-alpha3.Code`))
      # avc$SubRegion<-left_join(data.frame(avc[,c('Territory','ISO_TER1')]),am0[,c('ISO-alpha3.Code','Region.Name','Sub-region.Name','Country.or.Area')],join_by(ISO_TER1==`ISO-alpha3.Code`))[,4]
      # xx<-left_join(data.frame(avc0[,c('Region')]),am0[,c('Region.Name','Sub-region.Name','Country.or.Area')],join_by(Region==Country.or.Area),keep=T)
      
      #UN countries
      UNL<-c('Austria','Belgium','Bulgaria','Croatia','Cyprus','Czech Republic','Denmark','Estonia','Finland','France',
             'Germany','Greece','Hungary','Ireland','Italy','Latvia','Lithuania','Luxembourg','Malta','Netherlands',
             'Poland','Portugal','Romania','Slovakia','Slovenia','Spain','Sweden')
      library(terra)
      library(sf)
      av<-vect(t(matrix(c(-11.415347026814477,37.16233965175856,
                          10.531579008688489,39.25839438548637,
                          14.230499127031685,34.69639290619643,
                          26.930124866676657,33.710014207971575,
                          31.12223433413228,45.42326124939169,
                          25.573854156617486,66.63040326122602,
                          -12.155131050483114,56.150129592586964),2,7)),crs='epsg:4326')
      av1<-vect(st_cast(st_combine(st_as_sf(av)),'POLYGON'))
      av<-vect(t(matrix(c(-18.68988992622276,36.1759609535337,
                          7.819037588570147,39.38169172276447,
                          14.230499127031685,33.710014207971575,
                          32.231910369635244,32.230446160634294,
                          35.06774912703169,26.312173971285176,
                          46.041212144783174,12.25627752158104,
                          65.2755967601678,-21.034003543507737,
                          32.97169439330388,-37.43254940149589,
                          -1.3049653700097394,-37.30925206421779,
                          -7.469832233915078,-17.211786087886395,
                          -29.170163594861826,15.831900302646162,
                          -20.78594465995058,33.21682485885918),2,12)),crs='epsg:4326')
      av2<-vect(st_cast(st_combine(st_as_sf(av)),'POLYGON'))
      # plot(avc0,ext=ext(-40,80,-90,90))
      # plot(av1,col=rgb(1,1,0,0.5),add=T)
      # plot(av2,col=rgb(1,0,1,0.5),add=T)
      
    }#Corrections to selected continent
    for (am0 in list(xx,dd)){
      gs00<-10
      #HarmFAOtoSHP(am)
      #av<-cbind(am$Region,am0$Region);View(av[av[,1]!=av[,2] | is.na(av[,2]),]) #Names transformed
      avc0<-vect(paste0('F:/Global_data/Code files/Administrative boundaries/RegionBufferv2_',gs00,'km.shp'))
      avc0$Area<-expanse(avc0,'km')
      cnt0<-centroids(avc0)
      for (i in unique(am0$Period)){
        for (k in c('Production','Imports','Exports','Consumption')){
          for (scn in unique(am0$Scn)){
            print(c(i,k))
            avx<-avc0
            am<-am0[am0$Period==i & am0$Scn==scn,]
            #am$Region[which(!is.element(am$Region,avx$Region))]
            avx$Value<-data.frame(left_join(data.frame(avx[,"Region"]),am[,c("Region",k)],join_by(Region)))[,k]
            if (nrml==T){avx$Value<-avx$Value/avx$Area} #Normalize by area
            
            #Corrections to selected continent
            av<-is.na(avx$Value) & is.related(cnt0,av1,"within") & is.element(avx$Region,UNL)
            avx$Value[av]<-am0[am0$Period==i & am0$Scn==scn & am0$Region=="European Union",k]
            if (nrml==T){avx$Value[av]<-avx$Value[av]/sum(avx$Area[av])}
            av<-is.na(avx$Value) & is.related(cnt0,av1,"within") & !is.element(avx$Region,UNL)
            avx$Value[av]<-am0[am0$Period==i & am0$Scn==scn & am0$Region=="Europe",k]
            if (nrml==T){avx$Value[av]<-avx$Value[av]/sum(avx$Area[av])}
            av<-is.na(avx$Value) & is.related(cnt0,av2,"within") & !is.element(avx$Region,UNL)
            avx$Value[av]<-am0[am0$Period==i & am0$Scn==scn & am0$Region=="Africa",k]
            if (nrml==T){avx$Value[av]<-avx$Value[av]/sum(avx$Area[av])}
            
            rasterize(avx,ar0,'Value',filename=paste0('F:/Global_data/Code files/Trade demand consumption/Trade_',ac,'_',scn,'_',i,'_',k,'_',gs00,'km.tif'),overwrite=T)
          }
        }
      }
    }
    
  }
  #Get different values for offshore aquaculture classes (AN and PL)
  {
    am<-data.frame(matrix(NA,4,3));names(am)<-c('Period','Var','Val')
    n<-1;for (spc in c('AN','PL')){
      for (yr in c(2015,2019)){
        av<-rast(paste0('F:/Global_data/Code files/Offshore aquaculture/OffshrAq_Area_',yr,'_',spc,'_0.1.tif'))
        am[n,]<-c(yr,spc,sum(values(av),na.rm=T));n<-n+1
      }  
    }
    am$Period<-nmr(am$Period)
    am$Val<-nmr(am$Val)
    am$Dff<-NA
    am$Dff[1]<-am$Val[2]-am$Val[1]
    am$Dff[3]<-am$Val[4]-am$Val[3]
    oa<-am[am$Period==2015,]
  }
  #Add trade information
  {
    #Baseline
    xj<-readRDS("F:/Global_data/Code files/Various/TradeStats.rds")
    #No need to normalize trade variables to raster dataset. other variables not used here.
    xj<-xj[is.element(xj$Period,c(2015,2019)),]
    xj<-xj[,c('Period','Region','Production','Imports','Exports','Consumption')]
    xj<-data.frame(pivot_longer(xj,names(xj)[3:length(xj)],names_to='Var',values_to='Val'))
    xj<-xj[is.element(xj$Region,unique(bb[[1]]$Region)),]
    xj<-data.frame(xj %>% group_by(Period,Var) %>% summarise(Val=sum(Val,na.rm=T)))
    xj$Period<-nmr(chr(xj$Period))
    xj$Period[xj$Period==2019]<-2020
    xj$Scn<-0
    
    #Projections
    dd<-bind_rows(bb)
    dd<-dd[dd$Period>=2025,]
    dd<-dd[is.element(dd$Var,c('Production','Imports','Exports','Consumption')),]
    dd<-data.frame(dd %>% group_by(Period,Var,Scn) %>% summarise(Val=sum(Val[Val>0],na.rm=T)))
    
    jj<-bind_rows(xj,dd[,names(xj)])
  }
  #Build growth model
  {
    #Historical, raster
    am<-ee[[1]][ee[[1]]$Scn==0,]
    am<-data.frame(pivot_longer(am,names(am)[4:7],names_to='Var',values_to='Val'))[,c('Period','Region','Var','Val')]
    am<-data.frame(am%>%group_by(Period,Var)%>%summarise(Val=sum(Val,na.rm=T)))
    am<-am[is.element(am$Period,c(2015,2020)),]
    am[1:(nrow(am)/2),'Dff']<-(am[(nrow(am)/2+1):nrow(am),'Val']-am[1:(nrow(am)/2),'Val'])#/am[1:(nrow(am)/2),'Val']*100
    am<-am[is.element(am$Period,c(2015)),]
    #
    ae<-jj[is.element(jj$Scn,c(0,i)) & is.element(jj$Period,c(2015,2020)),];ae$Dff<-NA
    ae$Dff[is.element(ae$Period,2015)]<-ae$Val[is.element(ae$Period,2020)]-ae$Val[is.element(ae$Period,2015)]
    #
    am0<-bind_rows(am,oa,ae[ae$Period==2015,names(oa)])
    
    #Projected, raster
    am<-bind_rows(ee);am<-am[ee[[1]]$Scn!=0,];am<-bind_rows(ee[[1]][ee[[1]]$Scn==0,],am)
    am<-data.frame(pivot_longer(am,names(am)[4:7],names_to='Var',values_to='Val'))[,c('Period','Region','Scn','Var','Val')]
    am<-data.frame(am%>%group_by(Period,Scn,Var)%>%summarise(Val=sum(Val,na.rm=T)))
    am2<-vector('list',3);names(am2)<-c(2,4,5)
    for (i in c(2,4,5)){
      ad<-am[is.element(am$Scn,c(0,i)),];ad$Dff<-NA
      ad$Dff[is.element(ad$Period,seq(2015,2045,5))]<-(ad$Val[is.element(ad$Period,seq(2020,2050,5))]-
                                                         ad$Val[is.element(ad$Period,seq(2015,2045,5))])#/ad$Val[is.element(ad$Period,seq(2015,2045,5))]*100
      #View(cbind(ad[is.element(ad$Period,seq(2020,2050,5)),],ad[is.element(ad$Period,seq(2015,2045,5)),]))
      ad$Dff[is.element(ad$Period,2050)]<-NA
      
      ae<-jj[is.element(jj$Scn,c(0,i)),];ae$Dff<-NA
      ae$Dff[is.element(ae$Period,seq(2015,2045,5))]<-(ae$Val[is.element(ae$Period,seq(2020,2050,5))]-
                                                         ae$Val[is.element(ae$Period,seq(2015,2045,5))])#/ae$Val[is.element(ae$Period,seq(2015,2045,5))]*100
      #View(cbind(ae[is.element(ae$Period,seq(2020,2050,5)),],ae[is.element(ae$Period,seq(2015,2045,5)),]))
      ae$Dff[is.element(ae$Period,2050)]<-NA
      
      am2[[chr(i)]]<-bind_rows(ad[complete.cases(ad),],ae[complete.cases(ae),])
    }
    
    ww<-vector('list',3);names(ww)<-c(2,4,5);for (i in 1:length(ww)){ww[[i]]<-vector('list',2);names(ww[[i]])<-c('AN','PL')}
    for (spc in c('AN','PL')){
      K<-data.frame(Var=c('GDP','Pop','Production','Imports','Exports','Consumption'));K$K<-NA;for (i in 1:nrow(K)){K$K[i]<-am0$Dff[am0$Var==spc]/am0$Dff[am0$Var==K$Var[i]]}
      #Fsh removed, not proportional to growth.
      for (i in c(2,4,5)){
        am<-left_join(am2[[chr(i)]],K,join_by(Var))
        am$Chng<-am$Dff*am$K
        ww[[chr(i)]][[spc]]<-data.frame(am%>%group_by(Period)%>%summarise(Mean=mean(Chng,na.rm=T),StdE=sd(Chng,na.rm=T),Median=median(Chng,na.rm=T),Min=min(Chng,na.rm=T),Max=max(Chng,na.rm=T)))
        #am[order(am$Period,am$Var),]
        #print(ww)
      }
      
    }
    saveRDS(ww,"E:/Global_data/Code files/Modeling/OA_GrowthModel.rds")
    
  }
}
#Plots
{
  if (scope=='Global'){TMP<-c(2015,2019);gs<-0.1}
  if (scope=='China'){TMP<-c(2005,2010,2015,2019);gs<-0.1}
  
  #Projections
  ww0<-readRDS(paste0('E:/',scope,'_data/Code files/Modeling/OA_GrowthModel.rds')) #Global model
  for (i in 1:length(ww0)){
    for (j in 1:length(ww0[[i]])){
      #ww0[[i]][[j]]<-ww0[[i]][[j]][ww0[[i]][[j]]$Period>=2020,]
      ww0[[i]][[j]]$Scn<-names(ww0)[i]
      ww0[[i]][[j]]$Spc<-names(ww0[[i]])[j]
      #ww0[[i]][[j]]$Val<-Reduce('+',ww0[[i]][[j]][,'Mean'],accumulate=T)
    }
    ww0[[i]]<-data.frame(bind_rows(ww0[[i]]))
  }
  ww0<-data.frame(bind_rows(ww0))#[,names(am0)]
  ww0$Mdl<-gsub('AN|PL','',ww0$Spc,perl=T)
  ww0$Mdl[ww0$Mdl=='']<-'Overall'
  ww0$Spc<-gsub('IncN|DecN|IncX|DecX','',ww0$Spc,perl=T)
  ww0$Spc<-ifelse(ww0$Spc=='AN','Fish','Seaweed')
  
  #Figure
  library(ggplot2)
  library(RColorBrewer) #RColorBrewer::display.brewer.all() 
  library(cowplot)
  library(extrafont)
  loadfonts(device="win",quiet=T)
  #font_import(paths="C:/WINDOWS/FONTS",pattern='Times')
  #fonts()
  
  ichazo<-'Median'
  #http://colorspace.r-forge.r-project.org/articles/manipulation_utilities.html
  clr<-RColorBrewer::brewer.pal(8,'Dark2')[c(1,2,3)]
  lbl<-mapply(function(x,y){paste0(x,'-',y)},x=seq(2015,2045,5),y=seq(2020,2050,5))
  aplt<-ggplot(ww0)+
    geom_line(aes(x=Period,y=!!sym(ichazo),color=Scn),linewidth=1.6)+
    geom_ribbon(aes(x=Period,ymin=Min,ymax=Max,fill=Scn),alpha=0.2)+
    geom_hline(yintercept=0,linetype=2,linewidth=1,color='grey')+
    #General
    theme(text=element_text(size=16,family="Times New Roman",color='black'))+
    theme(panel.background = element_rect(fill="White",colour="black"),
          panel.border=element_rect(colour="black",fill=NA,linewidth=0.5))+
    theme(panel.grid.major=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5),
          panel.grid.minor=element_line(colour=brewer.pal(9,'Greys')[2],linewidth=0.5))+
    #Legend
    theme(legend.position="right",legend.title=element_blank(),
          legend.text=element_text(size=16,color='black'),
          legend.key = element_rect(color=NA,fill=NA),
          #legend.key.height=unit(0.75,"cm"),
          #legend.key.size=unit(1,"cm"),
          #legend.key.width=unit(1,"cm"),
          legend.spacing.y=unit(0.4,'cm'),
          plot.title=element_text(size=24,hjust=0.5))+
    #Axis
    #Margins from top, clockwise
    xlab(expression("Period"))+
    ylab(expression("Farm area change (km"^2*")"))+
    theme(axis.title.x=element_text(size=16,margin=ggplot2::margin(t=10,r=0,b=0,l=0)),
          axis.title.y=element_text(size=16,margin=ggplot2::margin(t=0,r=20,b=0,l=0)))+
    scale_x_continuous(labels=lbl,
      breaks=seq(min(TMP),2045,5),minor_breaks=seq(min(TMP),2045,5),expand=expansion(add=c(2,2),mult=c(0,0)))+
    #scale_y_discrete(labels=c('No','Yes'))+
    theme(plot.margin=unit(c(0.8,0.8,0.4,0.4),"cm"))+ 
    theme(axis.text.x=element_text(size=14,color="black",angle=45,hjust=1))+
    theme(axis.text.y=element_text(size=14,color="black"))+
    #theme(axis.ticks.x=element_line(color='black'),axis.ticks.y=element_line(color='black'))+
    theme(axis.ticks.length=unit(.1,"cm"))+
    facet_grid(cols=vars(Spc),rows=vars(Mdl),scales='free')+
    theme(strip.background=element_rect(fill="NA"),
          strip.text.y=element_blank(), #text(size=12,color="black",angle=0)
          strip.text.x=element_text(size=16,color="black",angle=0))+
    scale_color_manual(values=clr,labels=paste0('SSP',c(2,4,5)))+
    scale_fill_manual(values=clr,labels=paste0('SSP',c(2,4,5)))
  pp<-aplt
  jpeg(paste0("F:/",scope,"_data/Code files/Figures/OverallGrowthModel.jpg"),units="cm",width=25,height=15,res=800);plot(pp);dev.off()

}

#### S2: Scenario simulation ####

library(terra)
library(parallel)
library(dplyr)

N<-50 #Number of iterations
scope<-'Global'
tgt<-vector('list',2);names(tgt)<-c('AN','PL');for (i in 1:length(tgt)){tgt[[i]]<-paste0('OffshrAq_',names(tgt)[i])}
gs<-0.2 #Grid resolution
gs00<-10 #Coastal buffer size (km)
spt<-vect(paste0('F:/',scope,'_data/Code files/Administrative boundaries/PopulationPointsBufferv2_',gs00,'km.shp'))
ar0<-rast(ext(c(-180,180,-90,90)),resolution=gs,crs='epsg:4326')
bf1<-20;bf2<-50
MD<-'Prb'
arr<-rasterize(spt,ar0,'PtID')
for (S in list(01:25,26:50)){ #Set of iterations
for (md in MD){ #Type of prediction
for (scn in c('2','4','5')){ #Set of scenarios
for (yr in chr(c(2030,2040,2050))){
    
    ## Get data of initial conditions and projections
    am<-PrepareProjData(scn,scope)
    
    #Prepare data
    am<-TrnsfrmData(am)
    
    #Factors to characters
    ac<-sapply(1:length(am),function(x){class(am[,x])=='factor'})
    for (j in which(ac==T)){am[,j]<-chr(am[,j])}
    am$OffshrAqBuffer___AN_D<-nmr(am$OffshrAqBuffer___AN_D)
    am$OffshrAqBuffer___PL_D<-nmr(am$OffshrAqBuffer___PL_D)
    
    #Get baseline
    if (scope=='Global'){am<-ApplyRSAssumptions(am)}
    for (spc in c('AN','PL')){
      av<-arr
      an<-which(values(arr)>0) #Index of valid cells (within coastal buffer)
      values(av)[an]<-am[am[,'Period']=='2020',paste0('OffshrAq_',spc)]
      print(sum(et(am[am[,'Period']=='2020',paste0('OffshrAq_',spc)])))
      writeRaster(av,paste0('F:/',scope,'_data/Code files/Modeling/SSP/FarmBaseline_',spc,'_',gs,'.tif'),overwrite=T)
    }
    
    #Get iterations
    cl<-makeCluster(4,type='PSOCK')
    clusterExport(cl,c("scope","yr","am","md","scn","md")) #"N" "tgt"
    CCC<-parLapply(cl,S,function(nnn)
    {
      set.seed(nnn)
      
      print(nnn)
      
      if (scope=='Global'){R<-c('0','1')}
      if (scope=='China'){R<-c('-1','0','1')}
      
      #Preparation
      ppp<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/PPP.rds')) #Max farm area values in 2020 for AN and PL, also used in normalizing separate models
      ww0<-readRDS(paste0('E:/',scope,'_data/Code files/Modeling/OA_GrowthModel.rds')) #Global model
      GetProbArea<-function(scope,scn,spc,vr0,yi,amp,ww0,rr0,mm0,xxx,yyy){
        
        ac<-'Median'
        
        #Get change
        if (scope=='Global'){
          #Overall change
          ww<-ww0[[scn]][[spc]][ww0[[scn]][[spc]]$Period==(nmr(yi)-5),c('Period','Min','Max',ac)]
          rr<-EnvStats::rtri(1,max(c(0,ww[1,'Min'])),ww[1,'Max'],ww[1,'Median'])
          asgn<-sign(1)
        }
        if (scope=='China'){
          
          #Proportion growth
          vr<-'Dec'
          ww<-ww0[[scn]][[paste0(spc,vr,'N')]][ww0[[scn]][[paste0(spc,vr,'N')]]$Period==(nmr(yi)-5),c('Period','Min','Max',ac)]
          rr1<-unlist(abs(round(ww[ac],0)))
          ww<-ww0[[scn]][[paste0(spc,vr,'X')]][ww0[[scn]][[paste0(spc,vr,'N')]]$Period==(nmr(yi)-5),c('Period','Min','Max',ac)]
          rr2<-unlist((ww[ac]))
          
          vr<-'Inc'
          ww<-ww0[[scn]][[paste0(spc,vr,'N')]][ww0[[scn]][[paste0(spc,vr,'N')]]$Period==(nmr(yi)-5),c('Period','Min','Max',ac)]
          rr3<-unlist(abs(round(ww[ac],0)))
          ww<-ww0[[scn]][[paste0(spc,vr,'X')]][ww0[[scn]][[paste0(spc,vr,'N')]]$Period==(nmr(yi)-5),c('Period','Min','Max',ac)]
          rr4<-unlist((ww[ac]))
          
          #Inc+Dec=rr0
          #Inc/(-Dec)=P
          P<-(rr3*rr4)/(-rr1*rr2)
          IncAr<-rr0*(-P)/(1-P)
          DecAr<-rr0/(1-P)
          
          if (vr0=='Dec'){rr<-abs(DecAr);asgn<-sign(DecAr)}
          if (vr0=='Inc'){rr<-abs(IncAr);asgn<-sign(IncAr)}
        }
        
        #Calculate accumulated area
        XX<-data.frame(PtID=amp$PtID,Prev=amp[,paste0('OffshrAq___',spc)],Prob=yyy,Growth=et(c(xxx)))
        av<-XX$Prob>0;table(av);XX<-XX[sample(XX$PtID[av],length(XX$PtID[av]),F,XX$Prob[av]),]
        XX$AccArea<-Reduce('+',XX$Growth,accumulate=T)
        if (!any(XX$AccArea<=rr)){
          avl<-XX$PtID[1]
          sv<-rep(0,length(amp$PtID))
          kk<-0
        }else{
          avn<-XX$PtID[1:max(which(XX$AccArea<=rr))]
          avl<-XX$PtID[max(which(XX$AccArea<=rr))+1]
          sv<-ifelse(is.element(amp$PtID,c(avn,avl)),1,0)
          kk<-sum(XX$Growth[1:max(which(XX$AccArea<=rr))])
          if (all(sv==1)){print(c(scn,spc,vr0,yi,' Not enough change'))}
        }
        
        #Calculate final area
        xxx[sv==0]<-0 #Growth only in areas selected by probability
        xxx[which(amp$PtID==avl)]<-lt(rr-kk) #Adjust growth in last point
        xxx<-xxx*asgn #Asign sign
        
        list(xxx,rr*asgn)
        
        #Scaling (problem, probably exceeds maximum growth in existing areas)
        # KK<-rr/sum(et(mm2));print(KK)
        # if (KK>=1){
        #   mm2<-et(mm2)*KK
        #   ##
        #   mm2[mm1==0]<-0 #Growth only in areas with actual farms
        # }
        
        
      }
      {
        library(terra)
        library(ranger)
        library(glmmTMB)
        library(glmnet)
        library(mgcv)
        library(probably)
        library(dplyr)
        lt<-function(x){sign(x)*log1p(abs(x))} #log(abs(x)+1)
        et<-function(x){sign(x)*expm1(abs(x))} #exp(abs(x))-1
        nmr<-function(x){as.numeric(x)}
        chr<-function(x){as.character(x)}
        fct<-function(x){factor(x,levels=sort(unique(x)))}
        mdl<-function(x){
          ux<-unique(x)
          tab<-tabulate(match(x, ux))
          ux[tab==max(tab)]
        }
        ws<-function(x,p){
          av<-quantile(x,p);x[x<av]<-av
          av<-quantile(x,1-p);x[x>av]<-av
          x
        }
        GetOffshrAqBufferData<-function(xx,spt,sptb,spc,fld,yi,scope,ar0){
          
          tmp<-yi
          
          if (fld=='_'){fx<-function(x){x[is.na(x)]<-0;sum(x)}}
          #if (fld=='___'){fx<-function(x){any(x[is.na(x)]=='1')}}
          if (fld=='___'){fx<-function(x){x[is.na(x)]<-'0';chr(max(nmr(x)))}} 
          ars<-rasterize(xx,ar0,field=paste0('OffshrAq',fld,spc))
          am<-terra::extract(ars,sptb,fx,touches=T) #table(am$last)
          #if (fld=='___'){am[,2][!is.finite(am[,2])]<-0} #Using #max(x,na.rm=T) above
          colnames(am)<-c('PtID','Val')
          am$Var<-paste0('OffshrAqBuffer',fld,spc)
          am$Period<-nmr(yi)
          av<-data.frame(geom(spt))
          am$x<-av$x;am$y<-av$y
          am<-am[,c('Period','x','y','PtID','Var','Val')]
          
          am
          
        }
        #
        TrnsfrmData<-function(am){
          #k<-'VesselFrq';par(mfrow=c(3,1));hist(am[,k]);hist(sqrt(am[,k]));hist(lt(am[,k]))
          av<-c('OffshrAq_AN','OffshrAq_PL','OffshrAq_CB',
                'OffshrAq_AN_F','OffshrAq_PL_F','OffshrAq_CB_F',
                'Population','GDPC','HDI',#'VesselFrq','InfrstrFrq','TourismValue','Bth','ShrLength',
                'Population_D','GDPC_D',#'HDI_D','VesselFrq_D','InfrstrFrq_D','TourismValue_D',
                'Fsh','Fsh_D')#'PriSect','SecSect','TerSect','PriSect_D','SecSect_D','TerSect_D','ProtAreas'
          for (i in names(am)[is.element(names(am),av)]){
            am[,i]<-lt(am[,i])
          }
          # am$OffshrAq_AN<-round(am$OffshrAq_AN*1000,0) #COMPoisson(link="loglambda")
          # am$OffshrAq_AN[am$OffshrAq_AN==0]<-1 #Gamma(link=log)
          am
        }
        BackTrnsfrmData<-function(am){
          #k<-'VesselFrq';par(mfrow=c(3,1));hist(am[,k]);hist(sqrt(am[,k]));hist(lt(am[,k]))
          av<-c('OffshrAq_AN','OffshrAq_PL','OffshrAq_CB',
                'OffshrAq_AN_D','OffshrAq_PL_D','OffshrAq_CB_D',
                'OffshrAq_AN_F','OffshrAq_PL_F','OffshrAq_CB_F',
                'Population','HDI','GDPC','VesselFrq','InfrstrFrq','TourismValue','Bth','ShrLength',
                'Population_D','HDI_D','GDPC_D','VesselFrq_D','InfrstrFrq_D','TourismValue_D',
                'Fsh','PriSect','SecSect','TerSect','ProtAreas',
                'Fsh_D','PriSect_D','SecSect_D','TerSect_D')
          for (i in names(am)[is.element(names(am),av)]){
            am[,i]<-et(am[,i])
          }
          # am$OffshrAq_AN<-round(am$OffshrAq_AN*1000,0) #COMPoisson(link="loglambda")
          # am$OffshrAq_AN[am$OffshrAq_AN==0]<-1 #Gamma(link=log)
          am
        }
        NormData<-function(am,type,ac,save=F,scope){
          
          #an<-sapply(1:ncol(am),function(x){class(data.frame(am)[,x])});an<-an=='numeric'
          
          if (save==T){
            np<-data.frame(Var=chr(NA),Par1=rep(nmr(NA),length(am)),Par2=rep(nmr(NA),length(am)))
            for (i in 1:length(am)){
              
              if (is.element(names(am)[i],c('PtID','x','y','Period')) | class(am[,i])=='factor'){next()}
              if (type=='sd'){
                np[i,1]<-names(am)[i]
                np[i,2]<-mean(am[,i])
                np[i,3]<-sd(am[,i])
                am[,i]<-(am[,i]-np[i,2])/np[i,3]
              }
              if (type=='mm'){
                np[i,1]<-names(am)[i]
                np[i,2]<-min(am[,i]) #,na.rm=T
                np[i,3]<-max(am[,i])
                am[,i]<-(am[,i]-np[i,2])/(np[i,3]-np[i,2])
              }
              np<-np[!is.na(np$Var),]
              saveRDS(np,paste0("F:/",scope,"_data/Code files/Modeling/",ac,".rds"))
            }
            return(np)
          }
          
          if (save==F){
            np<-readRDS(paste0("F:/",scope,"_data/Code files/Modeling/",ac,".rds"))
            row.names(np)<-np$Var
            for (i in row.names(np)){
              if (type=='sd'){am[,i]<-(am[,i]-np[i,2])/np[i,3]}
              if (type=='mm'){am[,i]<-(am[,i]-np[i,2])/(np[i,3]-np[i,2])}
            }
            return(am)
          }
          
        }
        BackNormData<-function(am,type,ac,scope){
          
          # np<-readRDS(paste0("F:/China_data/Code files/Modeling/",ac,".rds"))
          # row.names(np)<-np$Var
          # for (i in np$Var){
          #   if (type=='sd'){am[,i]<-am[,i]*np[i,3]+np[i,2]}
          #   if (type=='mm'){am[,i]<-am[,i]*(np[i,3]-np[i,2])+np[i,2]}
          # }
          
          np<-readRDS(paste0("F:/",scope,"_data/Code files/Modeling/",ac,".rds"))
          row.names(np)<-np$Var
          for (i in row.names(np)){
            if (type=='sd'){am[,i]<-am[,i]*np[i,3]+np[i,2]}
            if (type=='mm'){am[,i]<-am[,i]*(np[i,3]-np[i,2])+np[i,2]}
          }
          
          am
          
        }
        TrnsNormVect<-function(av,ss,ac,scope){
          
          np<-readRDS(paste0("F:/",scope,"_data/Code files/Modeling/",ac,".rds"))
          row.names(np)<-np$Var
          
          av<-(lt(av)-np[ss,2])/(np[ss,3]-np[ss,2])
          
          av
        }
        BackTrnsNormVect<-function(av,ss,ac,scope){
          
          np<-readRDS(paste0("F:/",scope,"_data/Code files/Modeling/",ac,".rds"))
          row.names(np)<-np$Var
          
          av<-et(av*(np[ss,3]-np[ss,2])+np[ss,2])
          
          av
        }
        #
        
        # Estimation of mean and error
        GenValGLMM<-function(sm,mm,distr,lnkfx,posintercept){
          
          am0<-mm
          am0[,"(Intercept)"]<-1
          #am0$`I(OffshrAq_AN^2)`<-am0$OffshrAq_AN^2
          #am0$`I(OffshrAq_PL^2)`<-am0$OffshrAq_PL^2
          
          cf<-summary(sm)$coefficients
          #Conditional model (fixed effects)
          am1<-data.frame(cf$cond)
          #if (posintercept==T & am1['(Intercept)','Estimate']<0){am1['(Intercept)','Estimate']<-0}
          if (nrow(am1)>0){
            cls<-row.names(am1)#[!is.element(row.names(am),c("LUC1","LUC_DNoCh"))]
            smc1<-as.matrix(am0[,cls])%*%as.matrix(mapply(function(x,y){rnorm(1,x,y)},x=am1[cls,'Estimate'],y=0)) #am1[cls,'Std..Error']
          }else{
            smc1<-as.matrix(1)
          }
          #Dispersion model (fixed effects)
          am4<-data.frame(cf$disp)
          if (nrow(am4)>0){
            cls<-row.names(am4)
            smd<-as.matrix(am0[,cls])%*%as.matrix(mapply(function(x,y){rnorm(1,x,y)},x=am4[cls,'Estimate'],y=0)) #am4[cls,'Std..Error']
          }else{
            smd<-as.matrix(1)
          }
          
          #Expected value
          if (lnkfx=='inverse'){mu<-1/(c(smc1))}
          if (lnkfx=='log'){mu<-exp(c(smc1))}
          if (lnkfx=='identity'){mu<-smc1}
          # if (lnkfx=='sqrt'){mu<-(smc1)^2}
          # if (lnkfx=='cloglog'){mu<-exp(log(-log(1-smc1)))}
          p<-family_params(sm) #Model parameter
          if (nrow(am4)>0){
            phi<-exp(smd)[,1] #Dispersion; Log-link
          }else{
            phi<-sigma(sm) #Dispersion, ~1
            if (distr=='Gamma'){phi<-sigma(sm)^2}
          }
          
          if (distr=='Gaussian'){vl<-mapply(function(x,y){rnorm(1,x,y)},x=mu,y=phi)}
          if (distr=='Gamma'){
            vl<-mapply(function(x,z){rgamma(1,shape=x,scale=z)},x=1/phi,z=mu*phi)
            #https://github.com/glmmTMB/glmmTMB/issues/990 for squared term in gamma
            #Old x=mu^2/phi^2,z=phi^2/mu #phi is the scale parameter, not the shape (as documentation says), assuming equation is correct
            #https://github.com/glmmTMB/glmmTMB/issues/990
          }
          if (distr=='Tweedie'){
            if (posintercept==F){vl<-mapply(function(x,y,z){mgcv::rTweedie(x,y,z)},x=mu,y=p,z=phi)}
            if (posintercept==T){
              print(c(length(which(mu<0)),length(mu)));vl<-mapply(function(x,y,z){mgcv::rTweedie(x,y,z)},x=sapply(mu,function(x){max(0,x)}),y=p,z=phi)
            }
          }
          
          # vv<-predict(sm,am0,type="conditional",se.fit=T,re.form=NULL) #se,fit gives error value compared with expected value (deterministic)
          # plot(vv$fit,am0[,paste0('OffshrAq_',spc,'_F')]);lines(0:5,0:5)
          
          vl
        }
        GenValGLMNET<-function(sm,mm,scope){
          if (scope=='China'){ac1x<-c("PriSect","SecSect","TerSect")}
          if (scope=='Global'){ac1x<-c("TradeNorm_Exports","TradeNorm_Imports")}
          ac<-c('Population','GDPC','Fsh',
                ac1x,
                'LUC_Wtr','LUC_Crp','LUC_Frs','LUC_Urb',
                'thetao','ph','so','chl','sws',
                'Bth','ShrLength','ProtAreas')
          ac<-c(paste0('OffshrAq_',spc),ac)
          predict(sm,as.matrix(mm[,ac]),s='lambda.min',type='response')
        }
        GenValRF<-function(sm,mm,cf){
          
          library(ranger)
          vv<-predict(sm,mm)$predictions
          
          #Use cross-validation iterations
          # av1<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/CNFM_',spc,'_00.rds'))
          # av2<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/CNFM_',spc,'_11.rds'))
          # ff<-sapply(vv,function(x){
          #   if (x=='0'){av<-av1}
          #   if (x=='1'){av<-av2}
          #   y<-sample(av,1,F)
          #   if (x=='0'){p<-c(y,1-y)}
          #   if (x=='1'){p<-c(1-y,y)}
          #   return(sample(chr(0:1),1,F,p))
          # })
          
          #Use cross-validation summary 
          av0<-cf/rowSums(cf)
          ff<-sapply(vv,function(x){
            return(sample(colnames(cf),1,F,av0[chr(x),]))
          })
          #plot(nmr(vv)+rnorm(length(ff),0,0.1),nmr(ff)+rnorm(length(ff),0,0.1))
          
          ff
          
        }
        smP<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/RF_PR_','FarmPresence','_','Prob','.rds'))
        smC<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/RF_CL_','FarmPresence','_','Prob','.rds'))
        sm4X<-readRDS(paste0('F:/',scope,'_data/Code files/Modeling/GLMM_AR_','FarmArea','_X','.rds'))
        
        ## Get points, grid, and area constraints
        if (scope=='China'){
          gs<-0.1 #Grid resolution
          gs00<-10 #km
          spt<-vect(paste0('F:/',scope,'_data/Code files/Administrative boundaries/PopulationPointsBufferv2_',gs00,'km.shp'))
          ar0<-rast(ext(c(106,127,15,43)),resolution=gs,crs='epsg:4326')
          bf1<-10;bf2<-20
        }
        if (scope=='Global'){
          gs<-0.2 #Grid resolution
          gs00<-10 #km
          spt<-vect(paste0('F:/',scope,'_data/Code files/Administrative boundaries/PopulationPointsBufferv2_',gs00,'km.shp'))
          ar0<-rast(ext(c(-180,180,-90,90)),resolution=gs,crs='epsg:4326')
          bf1<-20;bf2<-50
          
        }
        #sptb<-vect(sapply(1:length(spt),function(x){buffer(spt[x,],bf2*1000)-buffer(spt[x,],bf1*1000)}))
        #writeVector(sptb,paste0('F:/',scope,'_data/Code files/Modeling/PopBuffer_Crl_',bf2,'-',bf1,'.rds'))
        sptb<-vect(paste0('F:/',scope,'_data/Code files/Modeling/PopBuffer_Crl_',bf2,'-',bf1,'.rds')) #Needed to calculate buffer variables in subsequent iterations
        arr<-rasterize(spt,ar0,'PtID')
        an<-which(values(arr)>0) #Index of valid cells (within buffer)
        #k<-10000;plot(arr,ext=ext(buffer(sptb[k,],2500)));plot(sptb[k,],col=rgb(0,0,0,0.5),add=T)
        AvAr<-spt$AvArCell
      }
      
      #Generate simulation
      ap<-vector('list',2);names(ap)<-c('AN','PL');for (j in 1:length(ap)){
        ap[[j]]<-vector('list',length(seq(2025,yr,5)));names(ap[[j]])<-seq(2025,yr,5)
      }
      ar<-vector('list',2);names(ar)<-c('AN','PL')
      tx<-vector('list',2);names(tx)<-c('AN','PL');for (j in 1:length(tx)){
        tx[[j]]<-vector('list',length(seq(2025,yr,5)));names(tx[[j]])<-seq(2025,yr,5)
      }
      ta<-vector('list',2);names(ta)<-c('AN','PL');for (j in 1:length(ta)){
        ta[[j]]<-vector('numeric',length(seq(2025,yr,5)));names(ta[[j]])<-seq(2025,yr,5)
      }
      tc<-vector('list',2);names(tc)<-c('AN','PL');for (j in 1:length(tc)){
        tc[[j]]<-vector('numeric',length(seq(2025,yr,5)));names(tc[[j]])<-seq(2025,yr,5)
      }
      wv<-vector('list',2);names(wv)<-c('AN','PL');for (j in 1:length(wv)){
        wv[[j]]<-vector('list',length(seq(2025,yr,5)));names(wv[[j]])<-seq(2025,yr,5)
      }
      `wv-`<-wv;`wv+`<-wv
      #par(mfrow=c(4,4)) #For printing results
      for (yi in chr(seq(2025,yr,5))){
        
        #Generate predictions
        for (spc in c('AN','PL')){
          
          amp<-am[am$Period==chr(nmr(yi)-5),]
          #sapply(1:length(amp),function(x){length(which(is.na(amp[,x])))})
          am1<-NormData(amp,'mm',paste0('MinMax',spc),save=F,scope)
          am2<-NormData(amp,'mm',paste0('FarmArea_MinMax',spc),save=F,scope);am2<-cbind(am2,data.frame(Spc=spc))
          if (scope=='China'){
            `am2-`<-NormData(amp,'mm',paste0('FarmArea_MinMax',paste0(spc,'-')),save=F,scope);`am2-`<-cbind(`am2-`,data.frame(Spc=spc))
            `am2+`<-NormData(amp,'mm',paste0('FarmArea_MinMax',paste0(spc,'+')),save=F,scope);`am2+`<-cbind(`am2+`,data.frame(Spc=spc))
          }
          
          #Predict growth change type
          if (md=='Prb'){
            mmP<-predict(smP[[spc]],am1,type='response')$predictions #Prediction
            mmP<-data.frame(mmP);names(mmP)<-paste0('P',1:length(mmP))
            mmP<-cal_apply(mmP,smC[[spc]]) #Calibration
            names(mmP)<-R
            mm0<-mmP #Growth probabilities
            mm1<-sapply(1:nrow(mmP),function(x){sample(R,1,F,mmP[x,])})
          }
          if (md=='Cls'){
            mmP<-NULL
            mm0<-NULL
            mm1<-nmr(GenValRF(sm0[[spc]],am1,t(sm0[[spc]]$confusion.matrix)))
          }
          #print(table(mm1))
          
          #Predict change amount
          if (scope=='Global'){
            
            amm<-am2
            spcsgn<-paste0(spc,'+')
            mm2<-c(GenValGLMNET(sm4X[[spcsgn]],amm,scope))+rnorm(nrow(amm),0,sm4X[[spcsgn]]$Error)
            #
            #xx0<-min(abs(lt(et(am[am$Period=='2020',paste0('OffshrAq_',spc)])-et(am[am$Period=='2015',paste0('OffshrAq_',spc)]))))
            mm2[mm2>1]<-1 #Change cannot surpass historical rates
            mm2[mm2<1e-6]<-1e-6 #Minimum change assumed, needed due to normal distribution
            #
            xx1<-lt(BackTrnsNormVect(amm[,paste0('OffshrAq_',spc)],paste0('OffshrAq_',spc),paste0('FarmArea_MinMax',spcsgn),scope)) #Previous farm area
            xx2<-max(am[am$Period=='2020',paste0('OffshrAq_',spc)]) #Maximum farm area observed (2020)
            mm2[xx1>xx2]<-0 #Changes in areas larger than farm area in 2020 are zero
            `mm2+A`<-mm2*ppp[spcsgn]
            
          }
          if (scope=='China'){
            
            amm<-`am2-`
            spcsgn<-paste0(spc,'-')
            mm2<-c(GenValGLMNET(sm4X[[spcsgn]],amm,scope))+rnorm(nrow(amm),0,sm4X[[spcsgn]]$Error)
            #
            #xx0<-abs(min(lt(et(am[am$Period=='2020',paste0('OffshrAq_',spc)])-et(am[am$Period=='2015',paste0('OffshrAq_',spc)]))))
            mm2[mm2>1]<-1 #Change cannot surpass historical rates
            mm2[mm2<1e-6]<-1e-6 #Needed due to normal distribution
            #
            xx1<-BackTrnsNormVect(amm[,paste0('OffshrAq_',spc)],paste0('OffshrAq_',spc),paste0('FarmArea_MinMax',spcsgn),scope) #Previous farm area
            xx2<-et(mm2*ppp[spcsgn]) #Reduction in farm area
            mm2[xx1-xx2<0]<-lt(xx1[xx1-xx2<0])/ppp[spcsgn] #Reduction cannot be larger than current area
            `mm2-A`<-mm2*ppp[spcsgn]
            
            amm<-`am2+`
            spcsgn<-paste0(spc,'+')
            mm2<-c(GenValGLMNET(sm4X[[spcsgn]],amm,scope))+rnorm(nrow(amm),0,sm4X[[spcsgn]]$Error)
            #
            #xx0<-min(abs(lt(et(am[am$Period=='2020',paste0('OffshrAq_',spc)])-et(am[am$Period=='2015',paste0('OffshrAq_',spc)]))))
            mm2[mm2>1]<-1 #Change cannot surpass historical rates
            mm2[mm2<1e-6]<-1e-6 #Needed due to normal distribution
            #
            xx1<-lt(BackTrnsNormVect(amm[,paste0('OffshrAq_',spc)],paste0('OffshrAq_',spc),paste0('FarmArea_MinMax',spcsgn),scope)) #Previous farm area
            xx2<-max(am[am$Period=='2020',paste0('OffshrAq_',spc)]) #Maximum farm area observed (2020)
            mm2[xx1>xx2]<-0 #Changes in areas larger than farm area in 2020 are zero
            `mm2+A`<-mm2*ppp[spcsgn]
            
          }
          #print(summary(`mm2-A`))
          #print(summary(`mm2+A`))
          
          #Adjust change to overall growth
          if (scope=='Global'){
            ac<-'Median'
            ww<-ww0[[scn]][[spc]][ww0[[scn]][[spc]]$Period==(nmr(yi)-5),c('Period','Min','Max',ac)]
            rr0<-EnvStats::rtri(1,ww[1,'Min'],ww[1,'Max'],ww[1,'Median']) #Put out, same for growth and contraction
            wv[[spc]][[yi]]<-rr0
            
            av<-GetProbArea(scope,scn,spc,paste0(spc),yi,amp,ww0,rr0,mm0,mm2,mm0[,'1'])
            mm2<-av[[1]]
            `wv+`[[spc]][[yi]]<-av[[2]]
          }
          if (scope=='China'){
            
            #Overall change
            ac<-'Median'
            ww<-ww0[[scn]][[spc]][ww0[[scn]][[spc]]$Period==(nmr(yi)-5),c('Period','Min','Max',ac)]
            rr0<-EnvStats::rtri(1,ww[1,'Min'],ww[1,'Max'],ww[1,'Median']) #Put out, same for growth and contraction
            wv[[spc]][[yi]]<-rr0
            
            av<-GetProbArea(scope,scn,spc,'Dec',yi,amp,ww0,rr0,mm0,`mm2-A`,mm0[,'-1'])
            `mm2-`<-av[[1]]
            `wv-`[[spc]][[yi]]<-av[[2]]
            av<-GetProbArea(scope,scn,spc,'Inc',yi,amp,ww0,rr0,mm0,`mm2+A`,mm0[,'1'])
            `mm2+`<-av[[1]]
            `wv+`[[spc]][[yi]]<-av[[2]]
          }
          #print(summary(`mm2-`))
          #print(summary(`mm2+`))
          
          #Predict final area
          area0<-amp[,paste0('OffshrAq_',spc)] #Previous farm area, transformed
          par(mfrow=c(2,1))
          if (scope=='Global'){
            plot(area0,mm2,main=paste0('Initial area and expansion ',spc,' ',yi))
            mm3<-lt(et(area0)+et(mm2))
          } #Old (global)
          if (scope=='China'){
            plot(area0,et(ifelse(mm1=='-1',`mm2-`,0)),main=paste0('Initial area and expansion ',spc,' ',yi),ylim=c(-5,5))
            points(area0,et(ifelse(mm1=='1',`mm2+`,0)))
            mm3<-lt(et(area0)-et(ifelse(mm1=='-1',`mm2-`,0))+et(ifelse(mm1=='1',`mm2+`,0))) #Final farm area
          }
          mm3[mm3<0]<-0 #Approximation needed due to overall growth adjustment (only relevnt for simultaneous expansion contractino simulation).
          print(c(yi,length(which(!is.finite(mm3)))))
          #sum(et(area0));sum(et(ifelse(mm1=='-1',`mm2-`,0)));sum(et(ifelse(mm1=='1',`mm2+`,0)))
          plot(area0,mm3,main=paste0('Initial and final area ',spc,' ',yi))
          par(mfrow=c(1,1))
          #print(summary(mm3))
          
          #Farm area
          am[am$Period==chr(nmr(yi)-5),paste0('OffshrAq_',spc,'_F')]<-mm3 #Only transform, normalization done at each step
          am[am$Period==chr(nmr(yi)),paste0('OffshrAq_',spc)]<-mm3 #Only transform, normalization done at each step
          #Farm absence/presence
          am[am$Period==chr(nmr(yi)-5),paste0('OffshrAq___',spc,'_F')]<-ifelse(mm3>0,'1','0') #mm1
          am[am$Period==chr(nmr(yi)),paste0('OffshrAq___',spc)]<-ifelse(mm3>0,'1','0') #mm1 #No need to normalize and transform
          
          ap[[spc]][[yi]]<-mm0 #For saving probability map of growth change
          
          #Check change patterns
          # print(c(yi,spc))
          # print(round(table(apply(am[am$Period==chr(nmr(yi)-5),c(paste0('OffshrAq___',spc),paste0('OffshrAq___',spc,'_F'))],1,function(x){paste0(x,collapse='-')}))/nrow(am[am$Period==chr(nmr(yi)-5),])*100,2))
        }
        
        #Check area within cell bounds
        if (T){
          bb<-AvAr #Area available (km2)
          av1<-am[am$Period==chr(nmr(yi)),'OffshrAq_AN'] #table(is.finite(av1))
          av2<-am[am$Period==chr(nmr(yi)),'OffshrAq_PL'] #table(is.finite(av2))
          if (length(which(!is.finite(av1)))>0){print(c('AN finite',length(which(!is.finite(av1)))))}
          if (length(which(!is.finite(av2)))>0){print(c('PL finite',length(which(!is.finite(av2)))))}
          av1[!is.finite(av1)]<-max(av1[is.finite(av1)])
          av2[!is.finite(av2)]<-max(av2[is.finite(av2)])
          av1<-et(av1)
          av2<-et(av2)
          av<-av1+av2 #Total area occupied
          ad<-av-bb;ad[ad<0]<-0 #Excess area
          if (any(ad>0)){ #am$PtID[am$Period==chr(nmr(yi)) & !is.finite(qq)];am[am$PtID==44958,]
            print(c(yi,'ad',length(which(ad>0))))
            #Weighted reduction (proportional to area already occupied)
            #qq<-av!=0 & ad>0
            qq<-ad>0
            # av1[qq]<-av1[qq]-ad[qq]*av1[qq]/av[qq]
            # av2[qq]<-av2[qq]-ad[qq]*av2[qq]/av[qq]
            av1[qq]<-av1[qq]/(av1[qq]+av2[qq])*bb[qq]
            av2[qq]<-av2[qq]/(av1[qq]+av2[qq])*bb[qq]
            #Update values
            am[am$Period==chr(nmr(yi)-5),'OffshrAq_AN_F'][qq]<-lt(av1[qq])
            am[am$Period==chr(nmr(yi)-5),'OffshrAq_PL_F'][qq]<-lt(av2[qq])
            am[am$Period==chr(nmr(yi)),'OffshrAq_AN'][qq]<-lt(av1[qq])
            am[am$Period==chr(nmr(yi)),'OffshrAq_PL'][qq]<-lt(av2[qq])
          }
        }
        
        #Calculate buffer and lag variables
        for (spc in c('AN','PL')){
          
          #Calculate buffer variables
          for (fld in c('___','_')){
            xx<-vect(am[am$Period==chr(nmr(yi)),],geom=c('x','y'),crs='epsg:4326',keepgeom=T)
            if (fld=='_'){xx[,paste0('OffshrAq',fld,spc)]<-et(data.frame(xx[,paste0('OffshrAq_',spc)]))}
            av<-GetOffshrAqBufferData(xx,spt,sptb,spc,fld,yi,scope,ar0)
            #av[,'Val'][!is.finite(av[,'Val'])]<-0
            if (fld=='_'){
              am[am$Period==chr(nmr(yi)),paste0('OffshrAqBuffer_',spc)]<-lt(av[,'Val'])
            }
            if (fld=='___'){
              am[am$Period==chr(nmr(yi)),paste0('OffshrAqBuffer___',spc)]<-av$Val
            }
          }
          
        }
        
        next()
        
      }
      
      #Get cells, area, and density
      for (spc in c('AN','PL')){
        for (yi in chr(seq(2025,yr,5))){
          
          #tx[[spc]][[yi]]<-am[am$Period==chr(nmr(yi)),]
          tx[[spc]][[yi]]<-sum(nmr(am[am$Period==chr(nmr(yi)),paste0('OffshrAq___',spc)]))
          
          av<-am[am$Period==chr(nmr(yi)),paste0('OffshrAq_',spc)]
          av<-ws(av,0.001) #Does not affect lower bound for data with many zeroes.
          ta[[spc]][yi]<-sum(et(av))
          #tc[[spc]][yi]<-mean(nmr(am[am$Period==chr(nmr(yi)),paste0('OffshrAq___',spc)])/AvAr) #Not used
        }
      }
      
      #plot(c(data.frame(`wv-`)[1,]+data.frame(`wv+`)[1,]),data.frame(wv)[1,]);lines(0:200,0:200)
      #plot(unlist(tx[[1]]),unlist(ta[[1]]))
      #plot(unlist(tx[[2]]),unlist(ta[[2]]))
      
      #Get probability map
      for (spc in c('AN','PL')){
        av<-arr
        values(av)[an]<-ap[[spc]][[chr(nmr(yr))]][,'1']
        ap[[spc]]<-as.data.frame(av,xy=T)
      }

      #Get outcome map
      for (spc in c('AN','PL')){
        av<-arr
        values(av)[an]<-am[am$Period==chr(nmr(yr)),paste0('OffshrAq_',spc)]
        ar[[spc]]<-as.data.frame(av,xy=T)
      }
      
      XX<-list(PROB=ap,RAST=ar,TA=ta,TX=tx,WW=wv,DT=am[,c('PtID','x','y','Period','Subregion',grep('OffshrAq_',names(am),fixed=T,value=T))])
    })
    stopCluster(cl)
    saveRDS(CCC,paste0('E:/',scope,'_data/Code files/Modeling/SSP/RawData202512_',scn,'_',yr,'_',N,'_',paste0(c(S[1],S[length(S)]),collapse='-'),'_',md,'.rds'))
    
    gc()
    
}
}}}

#Organize data
library(tidyr)
library(dplyr)
AvAr<-spt$AvArCell
tgt<-vector('list',2);names(tgt)<-c('AN','PL');for (i in 1:length(tgt)){tgt[[i]]<-paste0('OffshrAq_',names(tgt)[i])}
for (md in MD){
  for (scn in c('2','4','5')){
    for (yr in chr(c(2030,2040,2050))){
      S<-1:25
      XX1<-readRDS(paste0('E:/',scope,'_data/Code files/Modeling/SSP/RawData202512_',scn,'_',yr,'_',N,'_',paste0(c(S[1],S[length(S)]),collapse='-'),'_',md,'.rds'))
      S<-26:50
      XX2<-readRDS(paste0('E:/',scope,'_data/Code files/Modeling/SSP/RawData202512_',scn,'_',yr,'_',N,'_',paste0(c(S[1],S[length(S)]),collapse='-'),'_',md,'.rds'))
      XX<-c(XX1,XX2)
      
      #Arrange total area
      ta<-expand.grid(Period=chr(seq(2005,yr,5)),It=1:N,Spc=c('AN','PL'),Val=nmr(NA))
      for (i in 1:N){
        
        am<-XX[[i]][['DT']]
        
        fx<-function(x){sum(et(x))}
        am<-data.frame(am %>% group_by(Period) %>% summarise(AN=fx(OffshrAq_AN),PL=fx(OffshrAq_PL)))
        am<-data.frame(pivot_longer(am,names(am)[2:3],names_to='Spc',values_to='Val'));am$It<-i
        #ta<-dplyr::rows_update(ta,am[!is.na(am$Val) & !is.element(am$Period,c('2015','2020')),],c('Period','It','Spc'))
        ta<-dplyr::rows_update(ta,am[!is.na(am$Val),],c('Period','It','Spc'))
      }
      data.frame(pivot_wider(ta[ta$Spc=='AN',],id_cols=names(ta)[1],names_from=It,values_from=Val))
      data.frame(pivot_wider(ta[ta$Spc=='PL',],id_cols=names(ta)[1],names_from=It,values_from=Val))
      
      tx<-expand.grid(Period=chr(seq(2025,yr,5)),It=1:N,Spc=c('AN','PL'),Val=nmr(NA))
      for (spc in c('AN','PL')){
        for (i in 1:N){
          av<-XX[[i]][['TX']][[spc]]
          av<-unlist(av)
          am<-data.frame(Period=names(av),It=i,Spc=spc,Val=av)
          tx<-dplyr::rows_update(tx,am,c('Period','It','Spc'))
        }
      }
      data.frame(pivot_wider(tx[tx$Spc=='AN',],id_cols=names(tx)[1],names_from=It,values_from=Val))
      data.frame(pivot_wider(tx[tx$Spc=='PL',],id_cols=names(tx)[1],names_from=It,values_from=Val))
      
      #Arrange probability map
      ap<-vector('list',2);names(ap)<-c('AN','PL');for (j in 1:length(ap)){ap[[j]]<-vector('list',N)}
      for (spc in c('AN','PL')){
        for (i in 1:N){
          # av<-vect(XX[[i]][['PROB']][[spc]],geom=c('x','y'),crs='epsg:4326',keepgeom=T)
          # av<-av[,c('x','y','Prob1')]
          # ap[[spc]][[i]]<-rasterize(av,arr,'Prob1')
          av<-XX[[i]][['PROB']][[spc]]
          av<-av[,c('x','y','last')]
          av$Spc<-spc
          av$It<-i
          av<-av[,c('x','y','Spc','It','last')]
          names(av)<-c('x','y','Spc','It','Val')
          ap[[spc]][[i]]<-av
        }
        ap[[spc]]<-dplyr::bind_rows(ap[[spc]])
      }
      
      #Arrange outcome map
      ar<-vector('list',2);names(ar)<-c('AN','PL');for (j in 1:length(ar)){ar[[j]]<-vector('list',N)}
      for (spc in c('AN','PL')){
        for (i in 1:N){
          # av<-vect(XX[[i]][['RAST']][[spc]],geom=c('x','y'),crs='epsg:4326',keepgeom=T)
          # names(av)<-c('x','y',tgt[[spc]])
          # ar[[spc]][[i]]<-rasterize(av,arr,tgt[[spc]])
          av<-XX[[i]][['RAST']][[spc]]
          names(av)<-c('x','y',tgt[[spc]])
          av$Spc<-spc
          av$It<-i
          av<-av[,c('x','y','Spc','It',tgt[[spc]])]
          names(av)<-c('x','y','Spc','It','Val')
          ar[[spc]][[i]]<-av
        }
        ar[[spc]]<-dplyr::bind_rows(ar[[spc]])
      }
      
      #Get simulation metrics
      #saveRDS(ta,paste0('F:/',scope,'_data/Code files/Modeling/SSP/TA_',scn,'_',yr,'.rds'))
      for (spc in c('AN','PL')){
        
        cc<-data.frame(ap[[spc]] %>% group_by(x,y,Spc) %>% summarise(
          Mdn=median(Val),Q0=quantile(Val,0),Q5=quantile(Val,0.05),Q10=quantile(Val,0.10),Q25=quantile(Val,0.25),Q50=quantile(Val,0.50),Q75=quantile(Val,0.75),Q90=quantile(Val,0.90),Q95=quantile(Val,0.95),Q100=quantile(Val,1),
          IQR=quantile(Val,0.75)-quantile(Val,0.25),IDR=quantile(Val,0.90)-quantile(Val,0.10),
          AppPrb=1-ecdf(Val)(0))) #(1-min(which(sort(Val)>0))/length(Val)
        arp<-rast(cc[cc$Spc==spc,c('x','y','Mdn','Q0','Q5','Q10','Q25','Q50','Q75','Q90','Q95','Q100','IQR','IDR','AppPrb')],type='xyz',crs='epsg:4326')
        writeRaster(arp,paste0('E:/',scope,'_data/Code files/Modeling/SSP/Scn_Prob_',scn,'_',yr,'_',spc,'_',md,'.tif'),overwrite=T)
        
        cc<-data.frame(ar[[spc]] %>% group_by(x,y,Spc) %>% summarise(
          Mdn=median(Val),Q0=quantile(Val,0),Q5=quantile(Val,0.05),Q10=quantile(Val,0.10),Q25=quantile(Val,0.25),Q50=quantile(Val,0.50),Q75=quantile(Val,0.75),Q90=quantile(Val,0.90),Q95=quantile(Val,0.95),Q100=quantile(Val,1),
          IQR=quantile(Val,0.75)-quantile(Val,0.25),IDR=quantile(Val,0.90)-quantile(Val,0.10),
          AppPrb=1-ecdf(Val)(0))) #(1-min(which(sort(Val)>0))/length(Val)
        arp<-rast(cc[cc$Spc==spc,c('x','y','Mdn','Q0','Q5','Q10','Q25','Q50','Q75','Q90','Q95','Q100','IQR','IDR','AppPrb')],type='xyz',crs='epsg:4326')
        writeRaster(arp,paste0('E:/',scope,'_data/Code files/Modeling/SSP/Scn_Outc_',scn,'_',yr,'_',spc,'_',md,'.tif'),overwrite=T)
        
      }
      library(dplyr)
      bb<-data.frame(tx %>% group_by(Period,Spc) %>% summarise(Mdn=median(Val),Q0=quantile(Val,0),Q5=quantile(Val,0.05),Q10=quantile(Val,0.10),Q25=quantile(Val,0.25),Q50=quantile(Val,0.50),Q75=quantile(Val,0.75),Q90=quantile(Val,0.90),Q95=quantile(Val,0.95),Q100=quantile(Val,1),
                                                               IQR=quantile(Val,0.75)-quantile(Val,0.25),IDR=quantile(Val,0.90)-quantile(Val,0.10),Q100=quantile(Val,1)))
      aa<-data.frame(ta[!is.na(ta$Val),] %>% group_by(Period,Spc) %>% summarise(Mdn=median(Val),Q0=quantile(Val,0),Q5=quantile(Val,0.05),Q10=quantile(Val,0.10),Q25=quantile(Val,0.25),Q50=quantile(Val,0.50),Q75=quantile(Val,0.75),Q90=quantile(Val,0.90),Q95=quantile(Val,0.95),Q100=quantile(Val,1),
                                                               IQR=quantile(Val,0.75)-quantile(Val,0.25),IDR=quantile(Val,0.90)-quantile(Val,0.10),Q100=quantile(Val,1)))
      openxlsx::write.xlsx(bb,paste0('E:/',scope,'_data/Code files/Modeling/SSP/FarmExpansion_',scn,'_',yr,'_',md,'.xlsx'))
      openxlsx::write.xlsx(aa,paste0('E:/',scope,'_data/Code files/Modeling/SSP/FarmArea_',scn,'_',yr,'_',md,'.xlsx'))
      
    }
  }
}

