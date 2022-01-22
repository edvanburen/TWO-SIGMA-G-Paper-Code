# boxplot_functions.R

# This script provides functions which can be used to create boxplots of results
# as we did in the TWO-SIGMA-G paper. Included are functions to:

# 1. Plot type-I error for both summary statistic based and raw data based simulations
# 2. Plot power at various test and reference set configurations for both
# summary statistic based and raw data based simulations
power<-function(x){mean(x<.05,na.rm=T)}
create_boxplot_null-function(num_list,ts=30,rs=100,meth="corr",axis.cex=2.5,re=FALSE,ylim=c(0,.3)){
  par(mar=c(5, 5.3, 4, 2) + 0.1)
  #a<-lapply(ls(pattern="p.val_corr_twosigmag_test80_ref20_adjust_ind"),FUN=get)
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test_",ts,"_ref_",rs,"_sim",x)}
  a1<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_refset2_test_",ts,"_ref_",rs,"_sim",x)}
  a1c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_adjust_ind_refset2_test_",ts,"_ref_",rs,"_sim",x)}
  a2<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_adjust_refset2_test_",ts,"_ref_",rs,"_sim",x)}
  a2c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_adjust_refset2_test_",ts,"_ref_",rs,"_sim",x)}
  a2m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test_",ts,"_ref_",rs,"_sim",x)}
  a2g<-lapply(num_list,temp)
  
  getlist<-function(x){c(sapply(unlist(x),get))}
  
  b1<-lapply(a1,FUN=getlist)
  b2<-lapply(a2,FUN=getlist)
  
  b1c<-lapply(a1c,FUN=getlist)
  b2c<-lapply(a2c,FUN=getlist)
  
  b2m<-lapply(a2m,FUN=getlist)
  b2g<-lapply(a2g,FUN=getlist)
  
  c1<-lapply(b1,FUN=power)
  c2<-lapply(b2,FUN=power)
  
  c1c<-lapply(b1c,FUN=power)
  c2c<-lapply(b2c,FUN=power)
  
  c2m<-lapply(b2m,FUN=power)
  c2g<-lapply(b2g,FUN=power)
  
  pow_all<-c(unlist(c1),unlist(c2),unlist(c1c),unlist(c2c),unlist(c2m),unlist(c2g))
  group_all<-c(rep("Unadjusted",length(a1)),rep("Adjusted",length(a2)),rep("Unadjusted",length(a1c)),rep("Adjusted",length(a2c)),rep("Adjusted",length(a2m)),rep("Unadjusted",length(a2g)))
  method_all<-c(rep("TWO-SIGMA-G",length(a1)),rep("TWO-SIGMA-G",length(a2)),rep("CAMERA",length(a1c)),rep("CAMERA",length(a2c)),rep("MAST",length(a2m)),rep("GSEA",length(a2g)))
  
  dat<-as.data.frame(cbind(pow_all,group_all,method_all), stringsAsFactors = FALSE)
  colnames(dat)<-c("Power","Group","Method")
  dat$Method<-factor(dat$Method)
  dat$Power <- as.numeric(dat$Power)
  save(dat,file="dat.RData")
  p1<-ggplot(data = dat,aes(x=Group,y=Power,group=interaction(Group,Method)))+geom_boxplot(aes(fill=Method))+scale_fill_manual(values=cbPalette)+
    scale_x_discrete(limits=c("Unadjusted","Adjusted"),name="p-values")+
    scale_y_continuous(limits=ylim,name="Set-Level Type-I Error")+geom_hline(yintercept=.05,linetype="dashed")
  
  if(is.null(re)){
    p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nGene-Level Random Effects\nIncorrectly Absent","\nTest Size = ",ts,", Ref. Size = ",rs))
  }else{
    if(re==TRUE){
      p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nGene-Level Random Effects Present","\nTest Size = ",ts,", Ref. Size = ",rs))
    }else{
      if(re==FALSE){
        p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nNo Gene-Level Random Effects","\nTest Size = ",ts,", Ref. Size = ",rs))
      }
      
    }
  }
  p3<-p2+theme(legend.text = element_text(size=rel(1.3)),legend.key.size =unit(1.2, "cm"),
               legend.title = element_blank(),legend.position = "bottom",legend.background = element_blank(),
               axis.text.y=element_text(size=rel(2)),axis.text.x=element_text(size=rel(1.8)),plot.title=element_text(size=rel(2),hjust=.5)
               ,axis.title=element_text(size=rel(2.5)))
  print(p3)
  if(is.null(re)){
    save(p3,file=paste0("t1e_boxplot_mast_gsea_remiss",min(unlist(num_list)),"-"
                        ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
  }else
  {
    if(re==TRUE){
      save(p3,file=paste0("t1e_boxplot_mast_gsea_re",min(unlist(num_list)),"-"
                          ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
    }else{
      if(re==FALSE){
        save(p3,file=paste0("t1e_boxplot_mast_gsea_nore",min(unlist(num_list)),"-"
                            ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
      }
      
    }
  }
  
}
create_boxplot_alt<-function(num_list,ts=30,rs=100,meth="corr",axis.cex=2.5,re=FALSE,ylim=c(0,1)){
  par(mar=c(6.2, 5, 4, 2) + 0.1)
  #a<-lapply(ls(pattern="p.val_corr_twosigmag_test80_ref20_adjust_ind"),FUN=get)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_mixed_ref50_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a1<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_mixed_ref50_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a1c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_mixed_ref50_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a1m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_mixed_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a1g<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_ref50_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_ref50_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_ref50_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2g<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test80_mixed_ref50_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a3<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_test80_mixed_ref50_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a3c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_test80_mixed_ref50_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a3m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test80_mixed_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a3g<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test80_ref50_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a4<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_test80_ref50_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a4c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_test80_ref50_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a4m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test80_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a4g<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test80_mixed_ref20_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a5<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_test80_mixed_ref20_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a5c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_test80_mixed_ref20_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a5m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test80_mixed_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a5g<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test80_ref20_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a6<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_test80_ref20_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a6c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_test80_ref20_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a6m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test80_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a6g<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test50_mixed_ref20_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a7<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_test50_mixed_ref20_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a7c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_test50_mixed_ref20_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a7m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test50_mixed_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a7g<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test50_ref20_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a8<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_test50_ref20_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a8c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_test50_ref20_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a8m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test50_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a8g<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test20_mixed_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a9<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_test20_mixed_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a9c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_test20_mixed_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a9m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test20_mixed_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a9g<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test20_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a10<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_camera_test20_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a10c<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_mast_test20_adjust_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a10m<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a10g<-lapply(num_list,temp)
  
  getlist<-function(x){c(sapply(unlist(x),get))}
  
  b1<-lapply(a1,FUN=getlist)
  b2<-lapply(a2,FUN=getlist)
  b3<-lapply(a3,FUN=getlist)
  b4<-lapply(a4,FUN=getlist)
  b5<-lapply(a5,FUN=getlist)
  b6<-lapply(a6,FUN=getlist)
  b7<-lapply(a7,FUN=getlist)
  b8<-lapply(a8,FUN=getlist)
  b9<-lapply(a9,FUN=getlist)
  b10<-lapply(a10,FUN=getlist)
  
  
  b1c<-lapply(a1c,FUN=getlist)
  b2c<-lapply(a2c,FUN=getlist)
  b3c<-lapply(a3c,FUN=getlist)
  b4c<-lapply(a4c,FUN=getlist)
  b5c<-lapply(a5c,FUN=getlist)
  b6c<-lapply(a6c,FUN=getlist)
  b7c<-lapply(a7c,FUN=getlist)
  b8c<-lapply(a8c,FUN=getlist)
  b9c<-lapply(a9c,FUN=getlist)
  b10c<-lapply(a10c,FUN=getlist)
  
  b1m<-lapply(a1m,FUN=getlist)
  b2m<-lapply(a2m,FUN=getlist)
  b3m<-lapply(a3m,FUN=getlist)
  b4m<-lapply(a4m,FUN=getlist)
  b5m<-lapply(a5m,FUN=getlist)
  b6m<-lapply(a6m,FUN=getlist)
  b7m<-lapply(a7m,FUN=getlist)
  b8m<-lapply(a8m,FUN=getlist)
  b9m<-lapply(a9m,FUN=getlist)
  b10m<-lapply(a10m,FUN=getlist)
  
  b1g<-lapply(a1g,FUN=getlist)
  b2g<-lapply(a2g,FUN=getlist)
  b3g<-lapply(a3g,FUN=getlist)
  b4g<-lapply(a4g,FUN=getlist)
  b5g<-lapply(a5g,FUN=getlist)
  b6g<-lapply(a6g,FUN=getlist)
  b7g<-lapply(a7g,FUN=getlist)
  b8g<-lapply(a8g,FUN=getlist)
  b9g<-lapply(a9g,FUN=getlist)
  b10g<-lapply(a10g,FUN=getlist)
  
  c1<-lapply(b1,FUN=power)
  c2<-lapply(b2,FUN=power)
  c3<-lapply(b3,FUN=power)
  c4<-lapply(b4,FUN=power)
  c5<-lapply(b5,FUN=power)
  c6<-lapply(b6,FUN=power)
  c7<-lapply(b7,FUN=power)
  c8<-lapply(b8,FUN=power)
  c9<-lapply(b9,FUN=power)
  c10<-lapply(b10,FUN=power)
  
  c1c<-lapply(b1c,FUN=power)
  c2c<-lapply(b2c,FUN=power)
  c3c<-lapply(b3c,FUN=power)
  c4c<-lapply(b4c,FUN=power)
  c5c<-lapply(b5c,FUN=power)
  c6c<-lapply(b6c,FUN=power)
  c7c<-lapply(b7c,FUN=power)
  c8c<-lapply(b8c,FUN=power)
  c9c<-lapply(b9c,FUN=power)
  c10c<-lapply(b10c,FUN=power)
  
  c1m<-lapply(b1m,FUN=power)
  c2m<-lapply(b2m,FUN=power)
  c3m<-lapply(b3m,FUN=power)
  c4m<-lapply(b4m,FUN=power)
  c5m<-lapply(b5m,FUN=power)
  c6m<-lapply(b6m,FUN=power)
  c7m<-lapply(b7m,FUN=power)
  c8m<-lapply(b8m,FUN=power)
  c9m<-lapply(b9m,FUN=power)
  c10m<-lapply(b10m,FUN=power)
  
  c1g<-lapply(b1g,FUN=power)
  c2g<-lapply(b2g,FUN=power)
  c3g<-lapply(b3g,FUN=power)
  c4g<-lapply(b4g,FUN=power)
  c5g<-lapply(b5g,FUN=power)
  c6g<-lapply(b6g,FUN=power)
  c7g<-lapply(b7g,FUN=power)
  c8g<-lapply(b8g,FUN=power)
  c9g<-lapply(b9g,FUN=power)
  c10g<-lapply(b10g,FUN=power)
  
  pow_all<-c(unlist(c1),unlist(c1c),unlist(c1m),unlist(c1g),
             unlist(c2),unlist(c2c),unlist(c2m),unlist(c2g),unlist(c3),unlist(c3c),unlist(c3m),unlist(c3g),unlist(c4),unlist(c4c),unlist(c4m),unlist(c4g)
             ,unlist(c5),unlist(c5c),unlist(c5m),unlist(c5g),unlist(c6),unlist(c6c),unlist(c6m),unlist(c6g),unlist(c7),unlist(c7c),unlist(c7m),unlist(c7g)
             ,unlist(c8),unlist(c8c),unlist(c8m),unlist(c8g),unlist(c9),unlist(c9c),unlist(c9m),unlist(c9g),unlist(c10),unlist(c10c),unlist(c10m),unlist(c10g))
  group_all<-c(rep("T100M,R50",length(a1)),rep("T100M,R50",length(a1c)),rep("T100M,R50",length(a1m)),rep("T100M,R50",length(a1g)),
               rep("T100,R50",length(a2)),rep("T100,R50",length(a2c)),rep("T100,R50",length(a2m)),rep("T100,R50",length(a2g)),
               rep("T80M,R50",length(a3)),rep("T80M,R50",length(a3c)),rep("T80M,R50",length(a3m)),rep("T80M,R50",length(a3g)),
               rep("T80,R50",length(a4)),rep("T80,R50",length(a4c)),rep("T80,R50",length(a4m)),rep("T80,R50",length(a4g)),
               rep("T80M,R20",length(a5)),rep("T80M,R20",length(a5c)),rep("T80M,R20",length(a5m)),rep("T80M,R20",length(a5g)),
               rep("T80,R20",length(a6)),rep("T80,R20",length(a6c)),rep("T80,R20",length(a6m)),rep("T80,R20",length(a6g)),
               rep("T50M,R20",length(a7)),rep("T50M,R20",length(a7c)),rep("T50M,R20",length(a7m)),rep("T50M,R20",length(a7g)),
               rep("T50,R20",length(a8)),rep("T50,R20",length(a8c)),rep("T50,R20",length(a8m)),rep("T50,R20",length(a8g)),
               rep("T20M,R0",length(a9)),rep("T20M,R0",length(a9c)),rep("T20M,R0",length(a9m)),rep("T20M,R0",length(a9g)),
               rep("T20,R0",length(a10)),rep("T20,R0",length(a10c)),rep("T20,R0",length(a10m)),rep("T20,R0",length(a10g)))
  
  method_all<-c(rep("TWO-SIGMA-G",length(a1)),rep("CAMERA",length(a1c)),rep("MAST",length(a1m))
                ,rep("GSEA",length(a1g)),
                rep("TWO-SIGMA-G",length(a2)),rep("CAMERA",length(a2c)),rep("MAST",length(a2m))
                ,rep("GSEA",length(a2g)),
                rep("TWO-SIGMA-G",length(a3)),rep("CAMERA",length(a3c)),rep("MAST",length(a3m))
                ,rep("GSEA",length(a3g)),
                rep("TWO-SIGMA-G",length(a4)),rep("CAMERA",length(a4c)),rep("MAST",length(a4m))
                ,rep("GSEA",length(a4g)),
                rep("TWO-SIGMA-G",length(a5)),rep("CAMERA",length(a5c)),rep("MAST",length(a5m))
                ,rep("GSEA",length(a5g)),
                rep("TWO-SIGMA-G",length(a6)),rep("CAMERA",length(a6c)),rep("MAST",length(a6m))
                ,rep("GSEA",length(a6g)),
                rep("TWO-SIGMA-G",length(a7)),rep("CAMERA",length(a7c)),rep("MAST",length(a7m))
                ,rep("GSEA",length(a7g)),
                rep("TWO-SIGMA-G",length(a8)),rep("CAMERA",length(a8c)),rep("MAST",length(a8m))
                ,rep("GSEA",length(a8g)),
                rep("TWO-SIGMA-G",length(a9)),rep("CAMERA",length(a9c)),rep("MAST",length(a9m))
                ,rep("GSEA",length(a9g)),
                rep("TWO-SIGMA-G",length(a10)),rep("CAMERA",length(a10c)),rep("MAST",length(a10m))
                ,rep("GSEA",length(a10g)))
  
  dat<-as.data.frame(cbind(pow_all,group_all,method_all), stringsAsFactors = FALSE)
  colnames(dat)<-c("Power","Group","Method")
  dat$Method<-factor(dat$Method)
  dat$Power <- as.numeric(dat$Power)
  save(dat,file="dat.RData")
  p1<-ggplot(data = dat,aes(x=Group,y=Power,group=interaction(Group,Method)))+geom_boxplot(aes(fill=Method))+scale_fill_manual(values=cbPalette)+
    scale_x_discrete(limits=c("T100M,R50","T100,R50","T80M,R50","T80,R50","T80M,R20","T80,R20","T50M,R20","T50,R20","T20M,R0","T20,R0"),name="Configuration",labels=c("T100M,\nR50","T100,\nR50","T80M,\nR50","T80,\nR50","T80M,\nR20","T80,\nR20","T50M,\nR20","T50,\nR20","T20M,\nR0","T20,\nR0"))+
    scale_y_continuous(limits=c(0,1),name="Set-Level Power")
  
  if(is.null(re)){
    p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nGene-Level Random Effects\nIncorrectly Absent","\nTest Size = ",ts,", Ref. Size = ",rs))
  }else{
    if(re==TRUE){
      p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nGene-Level Random Effects Present","\nTest Size = ",ts,", Ref. Size = ",rs))
    }else{
      if(re==FALSE){
        p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nNo Gene-Level Random Effects","\nTest Size = ",ts,", Ref. Size = ",rs))
      }
      
    }
  }
  p3<-p2+theme(legend.text = element_text(size=rel(1.3)),legend.title = element_blank(),legend.position = "bottom",legend.background = element_blank(),legend.key.size =unit(1.2, "cm"),
               axis.text.y=element_text(size=rel(2)),axis.text.x=element_text(size=rel(1.3)),plot.title=element_text(size=rel(2),hjust=.5)
               ,axis.title=element_text(size=rel(2)))+ guides(colour = guide_legend(nrow = 1))
  print(p3)
  if(is.null(re)){
    save(p3,file=paste0("pow_boxplot2_mast_gsea_mixed5_remiss",min(unlist(num_list)),"-"
                        ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
  }else
  {
    if(re==TRUE){
      save(p3,file=paste0("pow_boxplot2_mast_gsea_mixed5_re",min(unlist(num_list)),"-"
                          ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
    }else{
      if(re==FALSE){
        save(p3,file=paste0("pow_boxplot2_mast_gsea_mixed5_nore",min(unlist(num_list)),"-"
                            ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
      }
      
    }
  }
  
}
create_sum_stat_boxplot_alt<-function(num_list,ts=30,rs=100,meth="corr",axis.cex=2.5,re=FALSE,ylim=c(0,1)){
  par(mar=c(6.2, 5, 4, 2) + 0.1)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_ref50_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_idea_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2i<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_page_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2p<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test80_ref50_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a4<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_idea_test80_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a4i<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_page_test80_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a4p<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test80_ref20_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a5<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_idea_test80_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a5i<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_page_test80_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a5p<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test50_ref20_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a7<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_idea_test50_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a7i<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_page_test50_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a7p<-lapply(num_list,temp)
  
  
  getlist<-function(x){c(sapply(unlist(x),get))}
  # b1<-lapply(a1,FUN=getlist)
  b2<-lapply(a2,FUN=getlist)
  # b3<-lapply(a3,FUN=getlist)
  b4<-lapply(a4,FUN=getlist)
  b5<-lapply(a5,FUN=getlist)
  # b6<-lapply(a6,FUN=getlist)
  b7<-lapply(a7,FUN=getlist)
  # b8<-lapply(a8,FUN=getlist)
  
  # b1i<-lapply(a1i,FUN=getlist)
  b2i<-lapply(a2i,FUN=getlist)
  # b3i<-lapply(a3i,FUN=getlist)
  b4i<-lapply(a4i,FUN=getlist)
  b5i<-lapply(a5i,FUN=getlist)
  # b6i<-lapply(a6i,FUN=getlist)
  b7i<-lapply(a7i,FUN=getlist)
  # b8i<-lapply(a8i,FUN=getlist)
  
  # b1p<-lapply(a1p,FUN=getlist)
  b2p<-lapply(a2p,FUN=getlist)
  # b3p<-lapply(a3p,FUN=getlist)
  b4p<-lapply(a4p,FUN=getlist)
  b5p<-lapply(a5p,FUN=getlist)
  # b6p<-lapply(a6p,FUN=getlist)
  b7p<-lapply(a7p,FUN=getlist)
  # b8p<-lapply(a8p,FUN=getlist)
  
  # c1<-lapply(b1,FUN=power)
  c2<-lapply(b2,FUN=power)
  # c3<-lapply(b3,FUN=power)
  c4<-lapply(b4,FUN=power)
  c5<-lapply(b5,FUN=power)
  # c6<-lapply(b6,FUN=power)
  c7<-lapply(b7,FUN=power)
  # c8<-lapply(b8,FUN=power)
  
  # c1i<-lapply(b1i,FUN=power)
  c2i<-lapply(b2i,FUN=power)
  # c3i<-lapply(b3i,FUN=power)
  c4i<-lapply(b4i,FUN=power)
  c5i<-lapply(b5i,FUN=power)
  # c6i<-lapply(b6i,FUN=power)
  c7i<-lapply(b7i,FUN=power)
  # c8i<-lapply(b8i,FUN=power)
  
  # c1p<-lapply(b1p,FUN=power)
  c2p<-lapply(b2p,FUN=power)
  # c3p<-lapply(b3p,FUN=power)
  c4p<-lapply(b4p,FUN=power)
  c5p<-lapply(b5p,FUN=power)
  # c6p<-lapply(b6p,FUN=power)
  c7p<-lapply(b7p,FUN=power)
  # c8p<-lapply(b8p,FUN=power)
  
  pow_all<-c(unlist(c2),unlist(c2i),unlist(c2p),unlist(c5),unlist(c5i),unlist(c5p),unlist(c4),unlist(c4i),unlist(c4p),unlist(c7),unlist(c7i),unlist(c7p))
  #pow_all<-c(unlist(c2),unlist(c2i),unlist(c5),unlist(c5i),unlist(c4),unlist(c4i),unlist(c7),unlist(c7i))
  group_all<-c(rep("T100,R50",length(a2)),rep("T100,R50",length(a2i)),rep("T100,R50",length(a2p)),
               rep("T80,R20",length(a5)),rep("T80,R20",length(a5i)),rep("T80,R20",length(a5p)),
               rep("T80,R50",length(a4)),rep("T80,R50",length(a4i)),rep("T80,R50",length(a4p)),
               rep("T50,R20",length(a7)),rep("T50,R20",length(a7i)),rep("T50,R20",length(a7i)))
  method_all<-c(rep("TWO-SIGMA-G",length(a2)),rep("iDEA",length(a2i)),rep("PAGE",length(a2p)),
                rep("TWO-SIGMA-G",length(a5)),rep("iDEA",length(a5i)),rep("PAGE",length(a5p)),
                rep("TWO-SIGMA-G",length(a4)),rep("iDEA",length(a4i)),rep("PAGE",length(a4p)),
                rep("TWO-SIGMA-G",length(a7)),rep("iDEA",length(a7i)),rep("PAGE",length(a7p)))
  
  dat<-as.data.frame(cbind(pow_all,group_all,method_all), stringsAsFactors = FALSE)
  colnames(dat)<-c("Power","Group","Method")
  dat$Method<-factor(dat$Method)
  dat$Power <- as.numeric(dat$Power)
  save(dat,file="dat.RData")
  p1<-ggplot(data = dat,aes(x=Group,y=Power,group=interaction(Group,Method)))+geom_boxplot(aes(fill=Method))+scale_fill_manual(values=cbPalette)+
    scale_x_discrete(limits=c("T100,R50","T80,R20","T80,R50","T50,R20"),name="Configuration (\\% of Genes \nDE in Test/Reference Sets)",labels=c("T100,\nR50","T80,\nR20","T80,\nR50","T50,\nR20"))+scale_y_continuous(limits=c(0,1),name="Set-Level Power")
  if(is.null(re)){
    p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nGene-Level Random Effects\nIncorrectly Absent","\nTest Size = ",ts,", Ref. Size = ",rs))
  }else{
    if(re==TRUE){
      p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nGene-Level Random Effects Present","\nTest Size = ",ts,", Ref. Size = ",rs))
    }else{
      if(re==FALSE){
        p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nNo Gene-Level Random Effects","\nTest Size = ",ts,", Ref. Size = ",rs))
      }
      
    }
  }
  p3<-p2+theme(legend.text = element_text(size=rel(1.5)),legend.title = element_blank(),legend.position = "bottom",legend.background = element_blank(),legend.key.size =unit(1.5, "cm"),
               axis.text.y=element_text(size=rel(2)),axis.text.x=element_text(size=rel(1.8)),plot.title=element_text(size=rel(2),hjust=.5)
               ,axis.title=element_text(size=rel(2.5)))+ guides(colour = guide_legend(nrow = 1))
  print(p3)
  if(is.null(re)){
    save(p3,file=paste0("pow_boxplot2_idea_page_combnull4_remiss",min(unlist(num_list)),"-"
                        ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
  }else
  {
    if(re==TRUE){
      save(p3,file=paste0("pow_boxplot2_idea_page_combnull4_re",min(unlist(num_list)),"-"
                          ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
    }else{
      if(re==FALSE){
        save(p3,file=paste0("pow_boxplot2_idea_page_combnull4_nore",min(unlist(num_list)),"-"
                            ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
      }
      
    }
  }
  
}
create_boxplot_sum_stat_null<-function(num_list,ts=30,rs=100,meth="corr",axis.cex=2.5,re=FALSE,ylim=c(0,.2),nulldiff=2300,cbPalette=c("#999999", "#E69F00", "#009E73","#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")){
  par(mar=c(6.2, 5, 4, 2) + 0.1)
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test20_ref20_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a1<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_idea_test20_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a1i<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_page_test20_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a1p<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_fgsea_test20_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a1f<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test20_ref20_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a1g<-lapply(num_list,temp)
  
  temp<-function(x){paste0("p.val_",meth,"_twosigmag_test50_ref50_adjust_ind_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_idea_test50_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2i<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_page_test50_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2p<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_fgsea_test50_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2f<-lapply(num_list,temp)
  temp<-function(x){paste0("p.val_",meth,"_gsea_test50_ref50_testsize_",ts,"_refsize_",rs,"_sim",x)}
  a2g<-lapply(num_list,temp)
  
  getlist<-function(x){c(sapply(unlist(x),get))}
  
  #b0<-lapply(a0,FUN=getlist)
  b1<-lapply(a1,FUN=getlist)
  b2<-lapply(a2,FUN=getlist)
  
  #b0i<-lapply(a0i,FUN=getlist)
  b1i<-lapply(a1i,FUN=getlist)
  b2i<-lapply(a2i,FUN=getlist)
  
  #b0p<-lapply(a0p,FUN=getlist)
  b1p<-lapply(a1p,FUN=getlist)
  b2p<-lapply(a2p,FUN=getlist)
  
  #b0f<-lapply(a0f,FUN=getlist)
  b1f<-lapply(a1f,FUN=getlist)
  b2f<-lapply(a2f,FUN=getlist)
  
  #b0g<-lapply(a0g,FUN=getlist)
  b1g<-lapply(a1g,FUN=getlist)
  b2g<-lapply(a2g,FUN=getlist)
  
  #c0<-lapply(b0,FUN=power)
  c1<-lapply(b1,FUN=power)
  c2<-lapply(b2,FUN=power)
  
  #c0i<-lapply(b0i,FUN=power)
  c1i<-lapply(b1i,FUN=power)
  c2i<-lapply(b2i,FUN=power)
  
  #c0p<-lapply(b0p,FUN=power)
  c1p<-lapply(b1p,FUN=power)
  c2p<-lapply(b2p,FUN=power)
  
  #c0f<-lapply(b0f,FUN=power)
  c1f<-lapply(b1f,FUN=power)
  c2f<-lapply(b2f,FUN=power)
  
  #c0g<-lapply(b0g,FUN=power)
  c1g<-lapply(b1g,FUN=power)
  c2g<-lapply(b2g,FUN=power)
  
  pow_all<-c(#unlist(c0),unlist(c0i),unlist(c0p),unlist(c0f),unlist(c0g),
    unlist(c1),unlist(c1i),unlist(c1p),unlist(c1f),unlist(c1g)
    ,unlist(c2),unlist(c2i),unlist(c2p),unlist(c2f),unlist(c2g))
  group_all<-c(#rep("T0,R0",length(a0)),rep("T0,R0",length(a0i)),rep("T0,R0",length(a0p)),rep("T0,R0",length(a0f)),rep("T0,R0",length(a0g)),
    rep("T20,R20",length(a1)),rep("T20,R20",length(a1i)),rep("T20,R20",length(a1p)),rep("T20,R20",length(a1f)),rep("T20,R20",length(a1g)),
    rep("T50,R50",length(a2)),rep("T50,R50",length(a2i)),rep("T50,R50",length(a2p)),rep("T50,R50",length(a2f)),rep("T50,R50",length(a2g)))
  
  method_all<-c(#rep("TWO-SIGMA-G",length(a0)),rep("iDEA",length(a0i)),rep("PAGE",length(a0p)),rep("fGSEA",length(a0f)),rep("GSEA",length(a0g)),
    rep("TWO-SIGMA-G",length(a1)),rep("iDEA",length(a1i)),rep("PAGE",length(a1p)),rep("fGSEA",length(a1f)),rep("GSEA",length(a1g)),
    rep("TWO-SIGMA-G",length(a2)),rep("iDEA",length(a2i)),rep("PAGE",length(a2p)),rep("fGSEA",length(a2f)),rep("GSEA",length(a2g)))
  
  dat<-as.data.frame(cbind(pow_all,group_all,method_all), stringsAsFactors = FALSE)
  colnames(dat)<-c("Power","Group","Method")
  dat$Method<-factor(dat$Method)
  dat$Power <- as.numeric(dat$Power)
  save(dat,file="dat.RData")
  p1<-ggplot(data = dat,aes(x=Group,y=Power,group=interaction(Group,Method)))+geom_boxplot(aes(fill=Method))+scale_fill_manual(values=cbPalette)+
    #scale_x_discrete(limits=c("T0,R0","T20,R20","T50,R50"),name="Configuration (\\% of Genes DE in \nTest/Reference Sets)",labels=c("T0,\nR0","T20,\nR20","T50,\nR50"))+
    scale_x_discrete(limits=c("T20,R20","T50,R50"),name="Configuration (\\% of Genes DE in \nTest/Reference Sets)",labels=c("T20,\nR20","T50,\nR50"))+
    scale_y_continuous(limits=ylim,name="Set-Level Type-I Error")+geom_hline(yintercept=.05,linetype="dashed")
  
  if(is.null(re)){
    p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nGene-Level Random Effects\nIncorrectly Absent","\nTest Size = ",ts,", Ref. Size = ",rs))
  }else{
    if(re==TRUE){
      p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nGene-Level Random Effects Present","\nTest Size = ",ts,", Ref. Size = ",rs))
    }else{
      if(re==FALSE){
        p2<-p1+ggtitle(paste0("Genes Simulated with IGC\nNo Gene-Level Random Effects","\nTest Size = ",ts,", Ref. Size = ",rs))
      }
      
    }
  }
  p3<-p2+theme(legend.text = element_text(size=rel(1.5)),legend.title = element_blank(),legend.position = "bottom",legend.background = element_blank(),legend.key.size =unit(1, "cm"),
               axis.text.y=element_text(size=rel(2)),axis.text.x=element_text(size=rel(1.8)),plot.title=element_text(size=rel(2),hjust=.5)
               # Changed 10/28/21
               ,axis.title=element_text(size=rel(2.5)))+ guides(colour = guide_legend(nrow = 1))
  print(p3)
  if(is.null(re)){
    save(p3,file=paste0("pow_boxplot2_idea_page_gsea_altnull_remiss",min(unlist(num_list)),"-"
                        ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
  }else
  {
    if(re==TRUE){
      save(p3,file=paste0("pow_boxplot2_idea_page_gsea_altnull_re",min(unlist(num_list)),"-"
                          ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
    }else{
      if(re==FALSE){
        save(p3,file=paste0("pow_boxplot2_idea_page_gsea_altnull_nore",min(unlist(num_list)),"-"
                            ,max(unlist(num_list)),"rs",rs,"ts",ts,"meth",meth,".RData"))
      }
      
    }
  }
  
}