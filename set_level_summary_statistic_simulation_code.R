# set_simulation_code.R

# This script is used to accomplish several set-level tasks
# related to the TWO-SIGMA-G paper simulations comparing
# to fGSEA, PAGE, and iDEA

# 1. Constructs test and reference sets for each scenario
# by randomly sampling genes from matched null and alternative
# hypotheses (using the arbitrary "sim_number" to ensure a match)
# 2. Computes and saves the IGC estimate using the individual-level
# residual correlation matrix for twosigmag
# 3. Runs set-level inference for CAMERA, MAST, GSEA, and twosigmag.  
# Including all methods in the same script is not ideal for distributed computing
# but serves as a way of being especially sure that the same
# test and reference sets are being used for all methods
# 4. Saves output



# Arbitrary Number used to uniquely identify simulations
sim_number_alt=as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

test_size<-30
# reference set size is used at both 30 and 100 
ref_size<-30

library(twosigma) 
library(fgsea) # for fGSEA function
library(iDEA) # for iDEA function
library(PGSEA) # For PAGE function
library(Biobase)
library(ArrayTools)
library(readr)
setwd("~/sim/array2/pathway")
capture.output(system.time({
  if(sim_number_alt>=14300 & sim_number_alt<=14399){
    sim_number_null<-sim_number_alt-2300
  }
  if(sim_number_alt>=14500 & sim_number_alt<=14599){
    sim_number_null<-sim_number_alt-2300
  }
  if(sim_number_alt>=14800 & sim_number_alt<=14899){
    sim_number_null<-sim_number_alt-2100
  }
  
  # How many sims are run in this submission
  # As with other code distributed computing is used to run multiple times
  nsims=250
  print(paste("Sim ",sim_number_alt, ", nsims = ",nsims))
  print(paste("Sim Number Null",sim_number_null, ", nsims = ",nsims))
  
  # Load precomputed results for null and alternative hypotheses
  load(paste0("/results_sim",sim_number_null,"_pathway.RData"))
  load(paste0("/residuals_sim",sim_number_null,"_pathway.RData"))
  load(paste0("/LRparams_sim",sim_number_null,"_pathway.RData"))
  ngenes_null_indep=length(get(paste0("LR_stat_orig_null_all_sim",sim_number_null)))

  load(paste0("/results_sim",sim_number_alt,"_pathway.RData"))
  load(paste0("/residuals_sim",sim_number_alt,"_pathway.RData"))
  load(paste0("/LRparams_sim",sim_number_alt,"_pathway.RData"))

  ngenes_alt_indep=length(get(paste0("LR_stat_orig_null_all_sim",sim_number_alt)))

  nruns=get(paste0("nruns_sim",sim_number_alt))
  sum_stats_alt<-get(paste0("sum_stats_sim",sim_number_alt))
  sum_stats_null<-get(paste0("sum_stats_sim",sim_number_null))
  
  # Variance for adjusted Wilcoxon test
  var0<-test_size*ref_size*(ref_size+test_size+1)/12

  # Way to drop sim_numbers from object names for easier referencing
  LR_stat_orig_alt_all<-get(paste0("LR_stat_orig_null_all_sim",sim_number_alt))
  LR_stat_large_noise_alt_all<-get(paste0("LR_stat_large_noise_null_all_sim",sim_number_alt))
  rm(list=paste0("LR_stat_large_noise_null_all_sim",sim_number_alt))

  LR_stat_orig_null_all<-get(paste0("LR_stat_orig_null_all_sim",sim_number_null))
  LR_stat_large_noise_null_all<-get(paste0("LR_stat_large_noise_null_all_sim",sim_number_null))
  rm(list=paste0("LR_stat_large_noise_null_all_sim",sim_number_null))

  dat_orig_alt_all<-get(paste0("dat_orig_null_all_sim",sim_number_alt))
  rm(list=paste0("dat_orig_null_all_sim",sim_number_alt))

  dat_orig_null_all<-get(paste0("dat_orig_null_all_sim",sim_number_null))
  rm(list=paste0("dat_orig_null_all_sim",sim_number_null))

  residuals_orig_alt_all<-get(paste0("residuals_orig_null_all_sim",sim_number_alt))
  rm(list=paste0("residuals_orig_null_all_sim",sim_number_alt))

  residuals_orig_null_all<-get(paste0("residuals_orig_null_all_sim",sim_number_null))
  rm(list=paste0("residuals_orig_null_all_sim",sim_number_null))

  #precomputed correlations when the test set is completely pure and from one simulation
  cor_large_noise<-get(paste0("cor_large_noise_sim",sim_number_alt))
  cor_ind<-get(paste0("cor_ind_sim",sim_number_alt))
  cor_ind_dat<-get(paste0("cor_ind_dat_sim",sim_number_alt))
  cor_dat<-get(paste0("cor_dat_sim",sim_number_alt))

  index_large_noise_null_all<-get(paste0("index_large_noise_null_all_sim",sim_number_alt))

  sim.seed<-get(paste0("sim.seed_sim",sim_number_alt))
  nind<-get(paste0("nind_sim",sim_number_alt))
  ncellsper<-get(paste0("ncellsper_sim",sim_number_alt))

  #Get covariates used in TWO-SIGMA for comparisons to CAMERAs
  set.seed(sim.seed[1])
  id.levels<-1:nind
  t2d_ind<-rbinom(nind,1,p=.4)
  t2d_sim<-rep(t2d_ind,times=ncellsper)
  nind<-length(id.levels)
  id<-rep(id.levels,times=ncellsper)
  cdr_sim<-rbeta(sum(ncellsper),3,6) #Is this good distn
  age_sim_ind<-sample(c(20:60),size=nind,replace = TRUE)
  age_sim<-rep(age_sim_ind,times=ncellsper)

  Z<-cbind(1,scale(t2d_sim),scale(age_sim),scale(cdr_sim))
  colnames(Z)<-c("Intercept","t2d_sim","age_sim","cdr_sim")
  X<-cbind(1,scale(t2d_sim),scale(age_sim),scale(cdr_sim))
  colnames(X)<-c("Intercept","t2d_sim","age_sim","cdr_sim")

  design_camera<-X[,c("Intercept","t2d_sim",setdiff(all.vars(get(paste0("mean_form_sim",sim_number_alt))[[3]]),"id"))]

  p.val_corr_wilcox<-numeric(length=nsims)
  p.val_corr_wilcox_adjust<-numeric(length=nsims)

  p.val_corr_wilcox_ref50<-numeric(length=nsims)
  p.val_corr_wilcox_ref50_adjust<-numeric(length=nsims)

  p.val_corr_wilcox_test50_adjust<-numeric(length=nsims)
  p.val_corr_wilcox_test50_ref50<-numeric(length=nsims)
  p.val_corr_wilcox_test50_ref50_adjust<-numeric(length=nsims)
  p.val_corr_wilcox_test50_ref20<-numeric(length=nsims)
  p.val_corr_wilcox_test50_ref20_adjust<-numeric(length=nsims)

  p.val_corr_wilcox_test50<-numeric(length=nsims)
  p.val_corr_wilcox_test80<-numeric(length=nsims)
  p.val_corr_wilcox_test80_adjust<-numeric(length=nsims)
  p.val_corr_wilcox_test80_ref50<-numeric(length=nsims)
  p.val_corr_wilcox_test80_ref50_adjust<-numeric(length=nsims)
  p.val_corr_wilcox_test80_ref20<-numeric(length=nsims)
  p.val_corr_wilcox_test80_ref20_adjust<-numeric(length=nsims)

  p.val_corr_wilcox_test20<-numeric(length=nsims)

  p.val_corr_wilcox_test0<-numeric(length=nsims)
  
  p.val_corr_wilcox_test20_ref20<-numeric(length=nsims)

  p.val_corr_wilcox_test20_adjust<-numeric(length=nsims)
  p.val_corr_wilcox_test0_adjust<-numeric(length=nsims)

  p.val_corr_wilcox_test20_ref20_adjust<-numeric(length=nsims)

  adj_stat_indep<-numeric(length=nsims)
  adj_stat_corr<-numeric(length=nsims)

  one_signif_in_test_set_indep<-logical(length=nsims)
  one_signif_in_test_set_corr<-logical(length=nsims)
  max_in_test_set_corr<-logical(length=nsims)
  max_in_test_set_corr_ref50<-logical(length=nsims)
  max_in_test_set_indep<-logical(length=nsims)
  gene_selected<-numeric(length=nsims)

  wilcox_stat_corr<-numeric(length=nsims)
  wilcox_stat_indep<-numeric(length=nsims)
  wilcox_stat_corr_ref50<-numeric(length=nsims)
  wilcox_stat_indep_ref50<-numeric(length=nsims)

  wilcox_stat_indep_test20<-numeric(length=nsims)
  wilcox_stat_indep_test20_ref20<-numeric(length=nsims)
  wilcox_stat_indep_test50<-numeric(length=nsims)
  wilcox_stat_indep_test80<-numeric(length=nsims)
  wilcox_stat_indep_test50_ref50<-numeric(length=nsims)
  wilcox_stat_indep_test50_ref20<-numeric(length=nsims)
  wilcox_stat_indep_test80_ref50<-numeric(length=nsims)
  wilcox_stat_indep_test80_ref20<-numeric(length=nsims)

  wilcox_stat_corr_test20<-numeric(length=nsims)
  wilcox_stat_corr_test0<-numeric(length=nsims)
  wilcox_stat_corr_test20_ref20<-numeric(length=nsims)
  wilcox_stat_corr_test50<-numeric(length=nsims)
  wilcox_stat_corr_test80<-numeric(length=nsims)
  wilcox_stat_corr_test50_ref50<-numeric(length=nsims)
  wilcox_stat_corr_test50_ref20<-numeric(length=nsims)
  wilcox_stat_corr_test80_ref50<-numeric(length=nsims)
  wilcox_stat_corr_test80_ref20<-numeric(length=nsims)

  ref_set_corr_0<-vector('list',length=nsims)
  ref_set_corr_50<-vector('list',length=nsims)
  ref_set_corr_20<-vector('list',length=nsims)

  sum_stats_ref_corr_0<-vector('list',length=nsims)
  sum_stats_ref_corr_20<-vector('list',length=nsims)
  sum_stats_ref_corr_50<-vector('list',length=nsims)

  ref_set_indep_0<-vector('list',length=nsims)
  ref_set_indep_50<-vector('list',length=nsims)
  ref_set_indep_20<-vector('list',length=nsims)

  sum_stats_test_corr<-vector('list',length=nsims)
  sum_stats_test_corr_test80<-vector('list',length=nsims)
  sum_stats_test_corr_test50<-vector('list',length=nsims)
  sum_stats_test_corr_test20<-vector('list',length=nsims)
  sum_stats_test_corr_test0<-vector('list',length=nsims)

  test_set_indep<-vector('list',length=nsims)
  test_set_corr<-vector('list',length=nsims)

  test_set_indep_test20<-vector('list',length=nsims)
  test_set_indep_test50<-vector('list',length=nsims)
  test_set_indep_test80<-vector('list',length=nsims)

  test_set_corr_test0<-vector('list',length=nsims)
  test_set_corr_test20<-vector('list',length=nsims)
  test_set_corr_test50<-vector('list',length=nsims)
  test_set_corr_test80<-vector('list',length=nsims)

  rho_indep_ind<-numeric(length=nsims)
  rho_indep_test50_ind<-numeric(length=nsims)
  rho_indep_test20_ind<-numeric(length=nsims)
  rho_indep_test80_ind<-numeric(length=nsims)

  rho_corr_ind<-numeric(length=nsims)
  rho_corr_test50_ind<-numeric(length=nsims)
  rho_corr_test20_ind<-numeric(length=nsims)
  rho_corr_test0_ind<-numeric(length=nsims)
  rho_corr_test80_ind<-numeric(length=nsims)

  rho_corr_dat<-numeric(length=nsims)
  rho_corr_test50_dat<-numeric(length=nsims)
  rho_corr_test20_dat<-numeric(length=nsims)
  rho_corr_test0_dat<-numeric(length=nsims)
  rho_corr_test80_dat<-numeric(length=nsims)

  rho_indep_camera<-numeric(length=nsims)
  rho_indep_test50_camera<-numeric(length=nsims)
  rho_indep_test20_camera<-numeric(length=nsims)
  rho_indep_test80_camera<-numeric(length=nsims)

  rho_corr_camera<-numeric(length=nsims)
  rho_corr_test50_camera<-numeric(length=nsims)
  rho_corr_test20_camera<-numeric(length=nsims)
  rho_corr_test0_camera<-numeric(length=nsims)
  rho_corr_test80_camera<-numeric(length=nsims)

  p.val_corr_mast_adjust<-numeric(length=nsims)
  p.val_corr_mast_ref50_adjust<-numeric(length=nsims)
  p.val_corr_mast_test50_adjust<-numeric(length=nsims)
  p.val_corr_mast_test50_ref50_adjust<-numeric(length=nsims)
  p.val_corr_mast_test50_ref20_adjust<-numeric(length=nsims)
  p.val_corr_mast_test80_adjust<-numeric(length=nsims)
  p.val_corr_mast_test80_ref50_adjust<-numeric(length=nsims)
  p.val_corr_mast_test80_ref20_adjust<-numeric(length=nsims)
  p.val_corr_mast_test20_adjust<-numeric(length=nsims)
  p.val_corr_mast_test20_ref20_adjust<-numeric(length=nsims)

  rho_corr_mast<-matrix(NA,nrow=nsims,ncol=2)
  rho_corr_test50_mast<-matrix(NA,nrow=nsims,ncol=2)
  rho_corr_test20_mast<-matrix(NA,nrow=nsims,ncol=2)
  rho_corr_test80_mast<-matrix(NA,nrow=nsims,ncol=2)

  if(ncol(design_camera)==2){mast_form<-~t2d_sim}
  if(ncol(design_camera)==4){mast_form<-~t2d_sim+age_sim+cdr_sim}

  if(ncol(design_camera)==2){mean_form<-count~t2d_sim;zi_form<-~t2d_sim}
  if(ncol(design_camera)==4){mean_form<-count~t2d_sim+age_sim+cdr_sim;zi_form<-~t2d_sim+age_sim+cdr_sim}
  compute_summary_stats_ts<-function(dat){
    beta<-numeric(length=test_size+ref_size)
    beta_se2<-numeric(length=test_size+ref_size)

    beta_zi<-numeric(length=test_size+ref_size)
    beta_zi_se2<-numeric(length=test_size+ref_size)
    for(xx in 1:(test_size+ref_size)){
      temp_ts<-twosigma_custom(count=dat[xx,]
        ,mean_form = mean_form
        ,zi_form=zi_form,id = id)
      beta[xx]<-summary(temp_ts)$coefficients$cond['t2d_sim',1]
      beta_zi[xx]<-summary(temp_ts)$coefficients$zi['t2d_sim',1]
      beta_se2[xx]<-summary(temp_ts)$coefficients$cond['t2d_sim',2]^2
      beta_zi_se2[xx]<-summary(temp_ts)$coefficients$zi['t2d_sim',2]^2
      if(xx%%10==0){print(xx)}}

    sum_stats<-data.frame(cbind(beta,beta_se2))
    rownames(sum_stats)<-paste("Gene",1:(test_size+ref_size))

    sum_stats_zi<-data.frame(cbind(beta_zi,beta_zi_se2))
    rownames(sum_stats)<-paste("Gene",1:(test_size+ref_size))

    return(sum_stats)
  }
  compute_summary_stats_zingeR<-function(dat){
    if(ncol(design_camera)==2){
      n_covar<-1
      colData <- data.frame(t2d_sim=t2d_sim)
      design <- model.matrix(~ t2d_sim)
      dse <- DESeqDataSetFromMatrix(countData = dat, colData = colData, design = ~ t2d_sim)
      weights <- zingeR::zeroWeightsLS(counts = dat, design = design, maxit = 500, normalization = "DESeq2_poscounts", colData = colData, designFormula = ~t2d_sim, verbose = TRUE,designZI = design)
      assays(dse)[["weights"]] <- weights
      dse <- DESeq2::estimateSizeFactors(dse, type="poscounts")
      dse <- estimateDispersions(dse,fitType="mean")
      # wald test contrast wrong??
      dse <- nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=F, useT=TRUE, df=rowSums(weights)-(n_covar-1) )
      res <- results(dse,cooksCutoff=FALSE)
      sum_stats_zinger<-cbind(res$log2FoldChange,res$lfcSE^2)
    }
    if(ncol(design_camera)==4){
      n_covar<-3
      colData <- data.frame(age_sim=age_sim,cdr_sim=cdr_sim,t2d_sim=t2d_sim)
      design <- model.matrix(~age_sim+cdr_sim+ t2d_sim)
      #dse <- DESeqDataSetFromMatrix(countData = dat, colData = colData, design = ~ age_sim+cdr_sim+t2d_sim)
      dse <- DESeqDataSetFromMatrix(countData = dat, colData = colData, design = ~ age_sim+cdr_sim+t2d_sim)
      #weights <- zingeR::zeroWeightsLS(counts = dat, design = design, maxit = 500, normalization = "DESeq2_poscounts", colData = colData, designFormula = ~age_sim+cdr_sim+t2d_sim, verbose = TRUE,designZI = design)
      weights <- zingeR::zeroWeightsLS(counts = dat, design = design, maxit = 500, normalization = "DESeq2_poscounts", colData = colData, designFormula = ~age_sim+cdr_sim+t2d_sim, verbose = TRUE,designZI = design)
      assays(dse)[["weights"]] <- weights
      dse <- DESeq2::estimateSizeFactors(dse, type="poscounts")
      dse <- estimateDispersions(dse,fitType="mean")
      # wald test contrast wrong??
      dse <- nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=F, useT=TRUE, df=rowSums(weights)-(n_covar-1) )
      res <- results(dse,cooksCutoff=FALSE)
      sum_stats_zinger<-cbind(res$log2FoldChange,res$lfcSE^2)
      return(sum_stats_zinger)
    }

  }
  compute_summary_stats_MAST<-function(dat){
      coldat<-as.data.frame(cbind(t2d_sim,age_sim,cdr_sim))
      rownames(coldat)<-colnames(dat)
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(dat+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      b<-cbind(coef(mast_model_corr,"C")[,'t2d_sim'],(se.coef(mast_model_corr,"C")[,'t2d_sim'])^2)
      return(b)
  }

  p.val_corr_idea_ref50<-rep(NA,nsims)
  p.val_corr_idea_test50_ref50<-rep(NA,nsims)

  p.val_corr_idea_test80_ref50<-rep(NA,nsims)
  p.val_corr_idea_test80_ref20<-rep(NA,nsims)
  p.val_corr_idea_test50_ref20<-rep(NA,nsims)
  p.val_corr_idea_test20_ref20<-rep(NA,nsims)

  p.val_corr_idea<-rep(NA,nsims)
  p.val_corr_idea_test80<-rep(NA,nsims)
  p.val_corr_idea_test50<-rep(NA,nsims)
  p.val_corr_idea_test20<-rep(NA,nsims)
  p.val_corr_idea_test0<-rep(NA,nsims)
  
  p.val_corr_page_ref50<-rep(NA,nsims)
  p.val_corr_page_test50_ref50<-rep(NA,nsims)
  
  p.val_corr_page_test80_ref50<-rep(NA,nsims)
  p.val_corr_page_test80_ref20<-rep(NA,nsims)
  p.val_corr_page_test50_ref20<-rep(NA,nsims)
  p.val_corr_page_test20_ref20<-rep(NA,nsims)
  
  p.val_corr_page<-rep(NA,nsims)
  p.val_corr_page_test80<-rep(NA,nsims)
  p.val_corr_page_test50<-rep(NA,nsims)
  p.val_corr_page_test20<-rep(NA,nsims)
  p.val_corr_page_test0<-rep(NA,nsims)
  
  p.val_corr_fgsea_ref50<-rep(NA,nsims)
  p.val_corr_fgsea_test50_ref50<-rep(NA,nsims)
  
  p.val_corr_fgsea_test80_ref50<-rep(NA,nsims)
  p.val_corr_fgsea_test80_ref20<-rep(NA,nsims)
  p.val_corr_fgsea_test50_ref20<-rep(NA,nsims)
  p.val_corr_fgsea_test20_ref20<-rep(NA,nsims)
  
  p.val_corr_fgsea<-rep(NA,nsims)
  p.val_corr_fgsea_test80<-rep(NA,nsims)
  p.val_corr_fgsea_test50<-rep(NA,nsims)
  p.val_corr_fgsea_test20<-rep(NA,nsims)
  p.val_corr_fgsea_test0<-rep(NA,nsims)
  
  p.val_corr_gsea_ref50<-rep(NA,nsims)
  p.val_corr_gsea_test80_ref50<-rep(NA,nsims)
  p.val_corr_gsea_test80_ref20<-rep(NA,nsims)
  p.val_corr_gsea_test50_ref20<-rep(NA,nsims)
  
  p.val_corr_gsea_test20_ref20<-rep(NA,nsims)
  p.val_corr_gsea_test50_ref50<-rep(NA,nsims)
  p.val_corr_gsea_test0<-rep(NA,nsims)
  
  p.val_corr_gsea<-rep(NA,nsims)
  p.val_corr_gsea_test80<-rep(NA,nsims)
  p.val_corr_gsea_test50<-rep(NA,nsims)
  p.val_corr_gsea_test20<-rep(NA,nsims)
  


  time<-proc.time()
  for(i in 1:nsims){
    # j is the "master gene" and corresponding test set
    # No matter what scenario evaluating we want to use this in the test set
    j<-sample(1:ngenes_alt_indep,1)
    gene_selected[i]<-j

    # Randomly pick correlated genes to be in the test set for correlated case
    # will be all if test_size is 30 based on simulation setup
    index<-sample(1:nruns,test_size-1,replace=test_size>30)

    index_null_indep<-sample(setdiff(1:ngenes_null_indep,j),ref_size,replace=ref_size>(ngenes_alt_indep-1))
    index_alt_indep<-sample(setdiff(1:ngenes_alt_indep,j),ref_size,replace=ref_size>(ngenes_alt_indep-1))

    #perfectly competitive ref set for independent case
    ref_data_0_indep<-dat_orig_null_all[index_null_indep,]
    ref_set_indep_0[[i]]<-LR_stat_orig_null_all[index_null_indep]

    #polluted reference set for independent case
    # Half of genes in reference under alt and half under null
    index_alt_ref50<-sample(index_alt_indep,round(.5*ref_size,0),replace=ref_size>(ngenes_alt_indep-1))
    index_null_ref50<-sample(index_null_indep,ref_size-round(.5*ref_size,0),replace=ref_size>(ngenes_alt_indep-1))
    ref_data_50_indep<-rbind(dat_orig_alt_all[index_alt_ref50,],dat_orig_null_all[index_null_ref50,])
    ref_set_indep_50[[i]]<-c(LR_stat_orig_alt_all[index_alt_ref50],LR_stat_orig_null_all[index_null_ref50])

    # 20% of ref genes alt
    index_alt_ref20<-sample(index_alt_indep,round(.20*ref_size,0))
    index_null_ref20<-sample(index_null_indep,ref_size-round(.20*ref_size,0))
    ref_data_20_indep<-rbind(dat_orig_alt_all[index_alt_ref20,],dat_orig_null_all[index_null_ref20,])
    ref_set_indep_20[[i]]<-c(LR_stat_orig_alt_all[index_alt_ref20],LR_stat_orig_null_all[index_null_ref20])

    #Perfectly competitive ref set but composed of the genes created with added noise
    # below we are constructing only the reference set for now
    #First randomly sample a set to include--one gene per set ensures no correlation exists

    # This will give genes that could be used in a reference set
    # where the could is depending on breakdown of alt vs null genes in ref set
    # can just exclude all of these from test set
    # to be sure we are not including genes in both test and reference sets
    # Test sets will only ever have alt genes associated with "j"
    # so excluding j here and the index_corr_ref genes for any null genes
    # makes sure that we will not double include any genes in test and reference sets
    index_corr_ref<-sample(setdiff(1:ngenes_alt_indep,j),ref_size)

    #index_corr_ref0<-sample(index_corr_ref,ref_size)
    ref_data_corr_0<-matrix(nrow=ref_size,ncol=sum(ncellsper))
    temp<-numeric(length=ref_size)
    temp2<-matrix(NA,nrow=ref_size,ncol=2)
    for(d in 1:ref_size){
      g<-index_corr_ref[d]
      # From the randomly sampled set choose a random gene to use
      # This preserves the general structure from the "correlated" genes but ensures we do not have correlated genes being chosen
      h<-sample(1:nruns,1)
      load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
      temp[d]<-LR_stat_large_noise_null_all[[g]][h]
      temp2[d,]<-sum_stats_null[[g]][h,]
      ref_data_corr_0[d,]<-test_dat[h,]
    }
    ref_set_corr_0[[i]]<-temp
    sum_stats_ref_corr_0[[i]]<-temp2

    #Polluted Competitive Reference Set made of genes added with noise
    # would still like to preserve the independence in the reference set though
    # so be sure no correlation exists in reference set
    #index_corr_ref50<-sample(setdiff(1:ngenes_alt,j),ref_size)
    ref_data_50_corr<-matrix(nrow=ref_size,ncol=sum(ncellsper))
    temp<-numeric(length=ref_size)
    temp2<-matrix(NA,nrow=ref_size,ncol=2)
    for(d in 1:ref_size){
      g<-index_corr_ref[d]
      # For the simulated set, randomly choose a given gene to include
      # Again we want to get data that is simulated like the "correlated" genes
      # But ensure that we do not have correlation in the reference set
      h<-sample(1:nruns,1)
      if(d<=.5*ref_size){
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_alt_all[[g]][h]
        temp2[d,]<-sum_stats_alt[[g]][h,]
      }
      if(d>.5*ref_size){
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_null_all[[g]][h]
        temp2[d,]<-sum_stats_null[[g]][h,]
      }
      ref_data_50_corr[d,]<-test_dat[h,]
    }
    ref_set_corr_50[[i]]<-temp
    sum_stats_ref_corr_50[[i]]<-temp2
    #index_corr_ref20<-sample(setdiff(1:ngenes_alt,j),ref_size)
    ref_data_20_corr<-matrix(nrow=ref_size,ncol=sum(ncellsper))
    temp<-numeric(length=ref_size)
    temp2<-matrix(NA,nrow=ref_size,ncol=2)
    for(d in 1:ref_size){
      g<-index_corr_ref[d]
      # For the simulated set, randomly choose a given gene to include
      # Again we want to get data that is simulated like the "correlated" genes
      # But ensure that we do not have correlation in the reference set
      h<-sample(1:nruns,1)
      if(d<=.2*ref_size){
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_alt_all[[g]][h]
        temp2[d,]<-sum_stats_alt[[g]][h,]
      }
      if(d>.2*ref_size){
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_null_all[[g]][h]
        temp2[d,]<-sum_stats_null[[g]][h,]
      }
      ref_data_20_corr[d,]<-test_dat[h,]
    }
    ref_set_corr_20[[i]]<-temp
    sum_stats_ref_corr_20[[i]]<-temp2
    #Randomly pick test set for indep case from statstics from independent genes
    # Being sure to exclude original and genes that are in any reference set (index_alt)
    index_alt_test<-sample(setdiff(1:ngenes_alt_indep,c(j,index_alt_indep)),test_size-1)

    # First set of statistics combines the original gene's statistic with either independent gene statistics

    test_set_indep[[i]]<-LR_stat_orig_alt_all[c(j,index_alt_test)]
    #sum_stats_indep_test[[i]]<-sum_stats[]
    one_signif_in_test_set_indep[i]<-sum(test_set_indep[[i]]>5.99,na.rm = T)>1

    test_set_corr[[i]]<-c(LR_stat_orig_alt_all[j],LR_stat_large_noise_alt_all[[j]][index])
    sum_stats_test_corr[[i]]<-sum_stats_alt[[j]]
    one_signif_in_test_set_corr[i]<-sum(test_set_corr[[i]]>5.99,na.rm = T)>1

    max_in_test_set_corr[i]<-max(c(test_set_corr[[i]],ref_set_corr_0[[i]]))%in%test_set_corr[[i]]
    max_in_test_set_indep[i]<-max(c(test_set_indep[[i]],ref_set_corr_0[[i]]))%in%test_set_indep[[i]]

    test_data_indep_100<-dat_orig_alt_all[c(j,index_alt_test),]
    test_resid_indep<-residuals_orig_alt_all[c(j,index_alt_test),]

    #loading this object will have the name test_dat
    load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",j,".RData"))
    test_data_corr_100<-rbind(matrix(dat_orig_alt_all[j,],nrow=1),test_dat)
    cor_temp<-cor(t(test_data_corr_100))
    rho_corr_dat[i]<-mean(cor_temp[upper.tri(cor_temp)])

    #Test set that is not perfectly DE
    test_data_corr_test50<-matrix(nrow=test_size,ncol=sum(ncellsper))
    test_data_corr_test50[1,]<-dat_orig_alt_all[j,]

    sum_stats_test_corr_test50[[i]]<-matrix(NA,nrow=test_size,ncol=2)
    sum_stats_test_corr_test50[[i]][1,]<-sum_stats_alt[[j]][1,]

    test_resid_corr_test50<-matrix(nrow=test_size,ncol=sum(ncellsper))
    test_resid_corr_test50[1,]<-residuals_orig_alt_all[j,]

    # gives the null genes that will be used
    # Being sure to remove any gene that could have been used in a reference set
    #Here
    index4<-sample(setdiff(1:ngenes_null_indep,c(j,index_corr_ref)),test_size,replace=test_size>30)
    temp<-numeric(length=test_size)

    # Gives the correlated genes that will be used
    # sample all possibilities (1:nruns) but will use less than all of them
    # depending on desired composition of test set of null and alternative genes
    # always pull
    #Here
    index4a<-sample(1:nruns,test_size-1,replace=test_size>30)
    temp[1]<-LR_stat_orig_alt_all[j]
    for(d in 2:(test_size)){
      if(d<=.5*test_size){
        g<-j
        h<-index4a[d]
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_alt_all[[g]][h]
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_alt,"/residuals_large_noise_sim",sim_number_alt,"_",g,"_",h,".RData"))
        test_resid_corr_test50[d,]<-residuals_large_noise
        sum_stats_test_corr_test50[[i]][d,]<-sum_stats_alt[[g]][h,]
      }
      if(d>.5*test_size){# use null data
        g<-index4[d]
        h<-sample(1:nruns,1)
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_null_all[[g]][h]
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/residuals_large_noise_sim",sim_number_null,"_",g,"_",h,".RData"))
        test_resid_corr_test50[d,]<-residuals_large_noise
        sum_stats_test_corr_test50[[i]][d,]<-sum_stats_null[[g]][h,]
      }
      test_data_corr_test50[d,]<-test_dat[h,]
    }
    test_set_corr_test50[[i]]<-temp

    cor_temp<-cor(t(test_data_corr_test50))
    rho_corr_test50_dat[i]<-mean(cor_temp[upper.tri(cor_temp)])

    #Test set that is not perfectly DE
    test_data_corr_test80<-matrix(nrow=test_size,ncol=sum(ncellsper))
    test_data_corr_test80[1,]<-dat_orig_alt_all[j,]
    test_resid_corr_test80<-matrix(nrow=test_size,ncol=sum(ncellsper))
    test_resid_corr_test80[1,]<-residuals_orig_alt_all[j,]

    sum_stats_test_corr_test80[[i]]<-matrix(NA,nrow=test_size,ncol=2)
    sum_stats_test_corr_test80[[i]][1,]<-sum_stats_alt[[j]][1,]
    # gives the null genes that will be used
    #Here
    index4<-sample(setdiff(1:ngenes_null_indep,c(j,index_corr_ref)),test_size,replace=test_size>30)
    temp<-numeric(length=test_size)

    # Gives the correlated genes that will be used
    # sample all possibilities (1:nruns) but will use less than all of them
    # depending on desired composition of test set of null and alternative genes
    #Here
    index4a<-sample(1:nruns,test_size-1,replace=test_size>30)
    temp[1]<-LR_stat_orig_alt_all[j]
    for(d in 2:(test_size)){
      if(d<=.8*test_size){
        g<-j
        h<-index4a[d]
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_alt_all[[g]][h]
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_alt,"/residuals_large_noise_sim",sim_number_alt,"_",g,"_",h,".RData"))
        test_resid_corr_test80[d,]<-residuals_large_noise
        sum_stats_test_corr_test80[[i]][d,]<-sum_stats_alt[[g]][h,]
      }
      if(d>.8*test_size){ # use null data
        g<-index4[d]
        h<-sample(1:nruns,1)
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_null_all[[g]][h]
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/residuals_large_noise_sim",sim_number_null,"_",g,"_",h,".RData"))
        test_resid_corr_test80[d,]<-residuals_large_noise
        sum_stats_test_corr_test80[[i]][d,]<-sum_stats_null[[g]][h,]
      }
      test_data_corr_test80[d,]<-test_dat[h,]
    }
    cor_temp<-cor(t(test_data_corr_test80))
    rho_corr_test80_dat[i]<-mean(cor_temp[upper.tri(cor_temp)])

    test_set_corr_test80[[i]]<-temp


    #Test set that is not perfectly DE
    test_data_corr_test20<-matrix(nrow=test_size,ncol=sum(ncellsper))
    test_data_corr_test20[1,]<-dat_orig_alt_all[j,]
    test_resid_corr_test20<-matrix(nrow=test_size,ncol=sum(ncellsper))
    test_resid_corr_test20[1,]<-residuals_orig_alt_all[j,]

    sum_stats_test_corr_test20[[i]]<-matrix(NA,nrow=test_size,ncol=2)
    sum_stats_test_corr_test20[[i]][1,]<-sum_stats_alt[[j]][1,]
    # gives the null genes that will be used
    #Here
    index4<-sample(setdiff(1:ngenes_null_indep,c(j,index_corr_ref)),test_size,replace = test_size>30)
    temp<-numeric(length=test_size)

    # Gives the correlated genes that will be used
    # sample all possibilities (1:nruns) but will use less than all of them
    # depending on desired composition of test set of null and alternative genes
    index4a<-sample(1:nruns,test_size-1,replace = test_size>30)
    temp[1]<-LR_stat_orig_alt_all[j]
    for(d in 2:(test_size)){
      if(d<=.2*test_size){
        g<-j
        h<-index4a[d]
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_alt_all[[g]][h]
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_alt,"/residuals_large_noise_sim",sim_number_alt,"_",g,"_",h,".RData"))
        test_resid_corr_test20[d,]<-residuals_large_noise
        sum_stats_test_corr_test20[[i]][d,]<-sum_stats_alt[[g]][h,]
      }
      if(d>.2*test_size){# use null data
        g<-index4[d]
        h<-sample(1:nruns,1)
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_null_all[[g]][h]
        load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/residuals_large_noise_sim",sim_number_null,"_",g,"_",h,".RData"))
        test_resid_corr_test20[d,]<-residuals_large_noise
        sum_stats_test_corr_test20[[i]][d,]<-sum_stats_null[[g]][h,]
      }
      test_data_corr_test20[d,]<-test_dat[h,]
    }
    test_set_corr_test20[[i]]<-temp
    cor_temp<-cor(t(test_data_corr_test20))
    rho_corr_test20_dat[i]<-mean(cor_temp[upper.tri(cor_temp)])

    # Test 0 
    
    test_data_corr_test0<-matrix(nrow=test_size,ncol=sum(ncellsper))
    test_data_corr_test0[1,]<-dat_orig_alt_all[j,]
    test_resid_corr_test0<-matrix(nrow=test_size,ncol=sum(ncellsper))
    test_resid_corr_test0[1,]<-residuals_orig_alt_all[j,]
    
    sum_stats_test_corr_test0[[i]]<-matrix(NA,nrow=test_size,ncol=2)
    sum_stats_test_corr_test0[[i]][1,]<-sum_stats_alt[[j]][1,]
    # gives the null genes that will be used
    #Here
    index4<-sample(setdiff(1:ngenes_null_indep,c(j,index_corr_ref)),test_size,replace = test_size>30)
    temp<-numeric(length=test_size)
    
    # Gives the correlated genes that will be used
    # sample all possibilities (1:nruns) but will use less than all of them
    # depending on desired composition of test set of null and alternative genes
    index4a<-sample(1:nruns,test_size-1,replace = test_size>30)
    temp[1]<-LR_stat_orig_null_all[j]
    for(d in 2:(test_size)){
      # use null data
      g<-index4[d]
      h<-sample(1:nruns,1)
      load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
      temp[d]<-LR_stat_large_noise_null_all[[g]][h]
      load(paste0("/pine/scr/e/v/evanbure/temp_sims/raw_dat/sim",sim_number_null,"/residuals_large_noise_sim",sim_number_null,"_",g,"_",h,".RData"))
      test_resid_corr_test0[d,]<-residuals_large_noise
      sum_stats_test_corr_test0[[i]][d,]<-sum_stats_null[[g]][h,]
      
      test_data_corr_test0[d,]<-test_dat[h,]
    }
    test_set_corr_test0[[i]]<-temp
    cor_temp<-cor(t(test_data_corr_test0))
    rho_corr_test0_dat[i]<-mean(cor_temp[upper.tri(cor_temp)])
    
    
    index5a<-sample(setdiff(1:ngenes_alt_indep,c(j,index_alt_indep)),size=test_size/2-1)
    index5<-sample(setdiff(1:ngenes_null_indep,c(j,index_null_indep)),size=test_size/2)

    test_set_indep_test50[[i]]<-c(LR_stat_orig_alt_all[c(j,index5a)],LR_stat_orig_null_all[index5])
    test_data_indep_test50<-rbind(dat_orig_alt_all[c(j,index5a),],dat_orig_null_all[index5,])
    test_resid_indep_test50<-rbind(residuals_orig_alt_all[c(j,index5a),],residuals_orig_null_all[index5,])

    index5a<-sample(setdiff(1:ngenes_alt_indep,c(j,index_alt_indep)),size=.8*test_size-1)
    index5<-sample(setdiff(1:ngenes_null_indep,c(j,index_null_indep)),size=test_size-.8*test_size)

    test_set_indep_test80[[i]]<-c(LR_stat_orig_alt_all[c(j,index5a)],LR_stat_orig_null_all[index5])
    test_data_indep_test80<-rbind(dat_orig_alt_all[c(j,index5a),],dat_orig_null_all[index5,])
    test_resid_indep_test80<-rbind(residuals_orig_alt_all[c(j,index5a),],residuals_orig_null_all[index5,])

    index5a<-sample(setdiff(1:ngenes_alt_indep,c(j,index_alt_indep)),size=.2*test_size-1)
    index5<-sample(setdiff(1:ngenes_null_indep,c(j,index_null_indep)),size=test_size-.2*test_size)

    test_set_indep_test20[[i]]<-c(LR_stat_orig_alt_all[c(j,index5a)],LR_stat_orig_null_all[index5])
    test_data_indep_test20<-rbind(dat_orig_alt_all[c(j,index5a),],dat_orig_null_all[index5,])
    test_resid_indep_test20<-rbind(residuals_orig_alt_all[c(j,index5a),],residuals_orig_null_all[index5,])

    ##################################################################
    # iDEA results
    ##################################################################
    annot<-c(rep(1,test_size),rep(0,ref_size))
    annot2<-1:test_size
    gs<-list(annot2)
    names(gs)<-"Set"

    rm(idea,idea2,dat_page)
    sum_stats_ref50<-rbind(unlist(sum_stats_test_corr[[i]]),unlist(sum_stats_ref_corr_50[[i]]))
    
    suppressMessages({
    tryCatch({
      idea <- CreateiDEAObject(sum_stats_ref50, annot, num_core=1)
      idea <- iDEA.fit(idea,
        fit_noGS=FALSE,
        init_beta=NULL,
        min_degene=0,
        em_iter=150,fit.tol=1e-8,
        mcmc_iter=1000,
        modelVariant = F,
         verbose=FALSE)
      idea2<-iDEA.louis(idea)
      p.val_corr_idea_ref50[i]<-idea2@gsea$pvalue},error=function(e){})})
    
    fgsea_stats<-sum_stats_ref50[,1]/sqrt(sum_stats_ref50[,2])
    names(fgsea_stats)<-1:(test_size+ref_size)
    p.val_corr_fgsea_ref50[i]<-fgsea(pathways = list("path1" = as.character(1:test_size)), 
          stats = fgsea_stats,
          nperm=10000)$pval
    
    dat_page<-data.frame(Coefficient=sum_stats_ref50[,1])
    p.val_corr_page_ref50[i]<-as.numeric(PGSEA(dat_page,cl = gs,range = c(0, 15000), ref = NULL, center = FALSE, p.value = TRUE, weighted = TRUE)$p.results)

    
    rm(idea,idea2,dat_page)
    
    suppressMessages({
    sum_stats_test80_ref50<-rbind(unlist(sum_stats_test_corr_test80[[i]]),unlist(sum_stats_ref_corr_50[[i]]))
    sum_stats_test80_ref50[,1]<-sum_stats_test80_ref50[,1]/1
    
    tryCatch({
      idea <- CreateiDEAObject(sum_stats_test80_ref50, annot, num_core=1)
      idea <- iDEA.fit(idea,
        fit_noGS=FALSE,
        init_beta=NULL,
        min_degene=0,
        em_iter=15,fit.tol=1e-2,
        mcmc_iter=1000,
        modelVariant = F,
         verbose=FALSE)
      idea2<-iDEA.louis(idea)
      p.val_corr_idea_test80_ref50[i]<-idea2@gsea$pvalue},error=function(e){})})
    
    fgsea_stats<-sum_stats_test80_ref50[,1]/sqrt(sum_stats_test80_ref50[,2])
    names(fgsea_stats)<-1:(test_size+ref_size)
    p.val_corr_fgsea_test80_ref50[i]<-fgsea(pathways = list("path1" = as.character(1:test_size)), 
                                     stats = fgsea_stats,
                                     nperm=10000)$pval
    
    dat_page<-data.frame(Coefficient=sum_stats_test80_ref50[,1])
    p.val_corr_page_test80_ref50[i]<-as.numeric(PGSEA(dat_page,cl = gs,range = c(0, 15000), ref = NULL, center = FALSE, p.value = TRUE, weighted = TRUE)$p.results)

    rm(idea,idea2,dat_page)
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_3"),"URL",paste0("Gene",1:test_size)))
    suppressMessages({
    sum_stats_test80_ref20<-rbind(unlist(sum_stats_test_corr_test80[[i]]),unlist(sum_stats_ref_corr_20[[i]]))
    sum_stats_test80_ref20[,1]<-sum_stats_test80_ref20[,1]/1
 
    tryCatch({
      idea <- CreateiDEAObject(sum_stats_test80_ref20, annot, num_core=1)
      idea <- iDEA.fit(idea,
        fit_noGS=FALSE,
        init_beta=NULL,
        min_degene=0,
        em_iter=15,fit.tol=1e-2,
        mcmc_iter=1000,
        modelVariant = F,
         verbose=FALSE)
      idea2<-iDEA.louis(idea)
      p.val_corr_idea_test80_ref20[i]<-idea2@gsea$pvalue},error=function(e){})})
    
    fgsea_stats<-sum_stats_test80_ref20[,1]/sqrt(sum_stats_test80_ref20[,2])
    names(fgsea_stats)<-1:(test_size+ref_size)
    p.val_corr_fgsea_test80_ref20[i]<-fgsea(pathways = list("path1" = as.character(1:test_size)), 
                                            stats = fgsea_stats,
                                            nperm=10000)$pval
    
    dat_page<-data.frame(Coefficient=sum_stats_test80_ref20[,1])
    p.val_corr_page_test80_ref20[i]<-as.numeric(PGSEA(dat_page,cl = gs,range = c(0, 15000), ref = NULL, center = FALSE, p.value = TRUE, weighted = TRUE)$p.results)

    rm(idea,idea2,dat_page)
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_4"),"URL",paste0("Gene",1:test_size)))
    suppressMessages({
    
    sum_stats_test50_ref20<-rbind(unlist(sum_stats_test_corr_test50[[i]]),unlist(sum_stats_ref_corr_20[[i]]))
    sum_stats_test50_ref20[,1]<-sum_stats_test50_ref20[,1]/1
 
    tryCatch({
      idea <- CreateiDEAObject(sum_stats_test50_ref20, annot, num_core=1)
      idea <- iDEA.fit(idea,
        fit_noGS=FALSE,
        init_beta=NULL,
        min_degene=0,
        em_iter=15,fit.tol=1e-2,
        mcmc_iter=1000,
        modelVariant = F,
         verbose=FALSE)
      idea2<-iDEA.louis(idea)
      p.val_corr_idea_test50_ref20[i]<-idea2@gsea$pvalue},error=function(e){})})
    
    fgsea_stats<-sum_stats_test50_ref20[,1]/sqrt(sum_stats_test50_ref20[,2])
    names(fgsea_stats)<-1:(test_size+ref_size)
    p.val_corr_fgsea_test50_ref20[i]<-fgsea(pathways = list("path1" = as.character(1:test_size)), 
                                            stats = fgsea_stats,
                                            nperm=10000)$pval
    
    dat_page<-data.frame(Coefficient=sum_stats_test50_ref20[,1])
    p.val_corr_page_test50_ref20[i]<-as.numeric(PGSEA(dat_page,cl = gs,range = c(0, 15000), ref = NULL, center = FALSE, p.value = TRUE, weighted = TRUE)$p.results)

    rm(idea,idea2,dat_page)

    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_5"),"URL",paste0("Gene",1:test_size)))
    suppressMessages({
    sum_stats_test20_ref20<-rbind(unlist(sum_stats_test_corr_test20[[i]]),unlist(sum_stats_ref_corr_20[[i]]))
    sum_stats_test20_ref20[,1]<-sum_stats_test20_ref20[,1]/1
    tryCatch({
      idea <- CreateiDEAObject(sum_stats_test20_ref20, annot, num_core=1)
      idea <- iDEA.fit(idea,
        fit_noGS=FALSE,
        init_beta=NULL,
        min_degene=0,
        em_iter=15,fit.tol=1e-2,
        mcmc_iter=1000,
        modelVariant = F,
         verbose=FALSE)
      idea2<-iDEA.louis(idea)
      p.val_corr_idea_test20_ref20[i]<-idea2@gsea$pvalue},error=function(e){})})
    
    fgsea_stats<-sum_stats_test20_ref20[,1]/sqrt(sum_stats_test20_ref20[,2])
    names(fgsea_stats)<-1:(test_size+ref_size)
    p.val_corr_fgsea_test20_ref20[i]<-fgsea(pathways = list("path1" = as.character(1:test_size)), 
                                            stats = fgsea_stats,
                                            nperm=10000)$pval
    
    dat_page<-data.frame(Coefficient=sum_stats_test20_ref20[,1])
    p.val_corr_page_test20_ref20[i]<-as.numeric(PGSEA(dat_page,cl = gs,range = c(0, 15000), ref = NULL, center = FALSE, p.value = TRUE, weighted = TRUE)$p.results)

    rm(idea,idea2,dat_page)
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_6"),"URL",paste0("Gene",1:test_size)))
    suppressMessages({
    
    sum_stats_test50_ref50<-rbind(unlist(sum_stats_test_corr_test50[[i]]),unlist(sum_stats_ref_corr_50[[i]]))
    sum_stats_test50_ref50[,1]<-sum_stats_test50_ref50[,1]/1
    tryCatch({
      idea <- CreateiDEAObject(sum_stats_test50_ref50, annot, num_core=1)
      idea <- iDEA.fit(idea,
        fit_noGS=FALSE,
        init_beta=NULL,
        min_degene=0,
        em_iter=15,fit.tol=1e-2,
        mcmc_iter=1000,
        modelVariant = F,
         verbose=FALSE)
      idea2<-iDEA.louis(idea)
      p.val_corr_idea_test50_ref50[i]<-idea2@gsea$pvalue},error=function(e){})})
    
    fgsea_stats<-sum_stats_test50_ref50[,1]/sqrt(sum_stats_test50_ref50[,2])
    names(fgsea_stats)<-1:(test_size+ref_size)
    p.val_corr_fgsea_test50_ref50[i]<-fgsea(pathways = list("path1" = as.character(1:test_size)), 
                                            stats = fgsea_stats,
                                            nperm=10000)$pval
    
    dat_page<-data.frame(Coefficient=sum_stats_test50_ref50[,1])
    p.val_corr_page_test50_ref50[i]<-as.numeric(PGSEA(dat_page,cl = gs,range = c(0, 15000), ref = NULL, center = FALSE, p.value = TRUE, weighted = TRUE)$p.results)
    
    rm(idea,idea2,dat_page)

    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_11"),"URL",paste0("Gene",1:test_size)))
    suppressMessages({
    sum_stats_test0_ref0<-rbind(unlist(sum_stats_test_corr_test0[[i]]),unlist(sum_stats_ref_corr_0[[i]]))
      tryCatch({
        idea <- CreateiDEAObject(sum_stats_test0_ref0, annot, num_core=1)
        idea <- iDEA.fit(idea,
                         fit_noGS=FALSE,
                         init_beta=NULL,
                         min_degene=0,
                         em_iter=15,fit.tol=1e-2,
                         mcmc_iter=1000,
                         modelVariant = F,
                         verbose=FALSE)
        idea2<-iDEA.louis(idea)
        p.val_corr_idea_test0[i]<-idea2@gsea$pvalue},error=function(e){})})
    
    fgsea_stats<-sum_stats_test0_ref0[,1]/sqrt(sum_stats_test0_ref0[,2])
    names(fgsea_stats)<-1:(test_size+ref_size)
    p.val_corr_fgsea_test0[i]<-fgsea(pathways = list("path1" = as.character(1:test_size)), 
                                            stats = fgsea_stats,
                                            nperm=10000)$pval
    
    dat_page<-data.frame(Coefficient=sum_stats_test0_ref0[,1])
    p.val_corr_page_test0[i]<-as.numeric(PGSEA(dat_page,cl = gs,range = c(0, 15000), ref = NULL, center = FALSE, p.value = TRUE, weighted = TRUE)$p.results)

    #
    # # Wilcoxon statistics
    #
    wilcox_stat_corr[i]<-sum(rank(c(test_set_corr[[i]],ref_set_corr_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_corr_test20[i]<-sum(rank(c(test_set_corr_test20[[i]],ref_set_corr_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    wilcox_stat_corr_test0[i]<-sum(rank(c(test_set_corr_test0[[i]],ref_set_corr_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_corr_test20_ref20[i]<-sum(rank(c(test_set_corr_test20[[i]],ref_set_corr_20[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_corr_ref50[i]<-sum(rank(c(test_set_corr[[i]],ref_set_corr_50[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    wilcox_stat_corr_test50[i]<-sum(rank(c(test_set_corr_test50[[i]],ref_set_corr_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_corr_test50_ref50[i]<-sum(rank(c(test_set_corr_test50[[i]],ref_set_corr_50[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_corr_test50_ref20[i]<-sum(rank(c(test_set_corr_test50[[i]],ref_set_corr_20[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_corr_test80_ref20[i]<-sum(rank(c(test_set_corr_test80[[i]],ref_set_corr_20[[i]]))[1:test_size]) - .5*test_size*(test_size+1)


    wilcox_stat_corr_test80[i]<-sum(rank(c(test_set_corr_test80[[i]],ref_set_corr_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_corr_test80_ref50[i]<-sum(rank(c(test_set_corr_test80[[i]],ref_set_corr_50[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    #
    #
    # #####################################################################################################
    # # Reference set here is perfectly competitive alternative
    # ######################################################################################################
    rho<-cor_ind[j]
    rho_corr_ind[i]<-rho
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))

    adj_stat_corr[i]<-(wilcox_stat_corr[i]-(.5*ref_size*test_size))/sqrt(var)
    p.val_corr_wilcox_adjust[i]<-2*(pnorm(-1*abs(adj_stat_corr[i])))

    adj_stat_corr[i]<-(wilcox_stat_corr_ref50[i]-(.5*ref_size*test_size))/sqrt(var)
    p.val_corr_wilcox_ref50_adjust[i]<-2*(pnorm(-1*abs(adj_stat_corr[i])))
    #
    # ################################################################################################################
    #
    # # Combined Test Set
    #
    # ################################################################################################################
    cor_temp<-numeric(length=nind)
    m<-1
    for(y in 1:nind){
      temp<-cor(t(test_resid_corr_test50[,m:(m+ncellsper[y]-1)]))
      cor_temp[y]<-mean(temp[upper.tri(temp)],na.rm = T)
      m<-m+ncellsper[y]
    }
    rho_corr_test50_ind[i]<-mean(cor_temp)
    rho<-mean(cor_temp)
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))
    p.val_corr_wilcox_test50_ref50_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test50_ref50[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_corr_wilcox_test50_ref20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test50_ref20[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_corr_wilcox_test50_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test50[i]-.5*test_size*ref_size)/sqrt(var)))
    #
    cor_temp<-numeric(length=nind)
    m<-1
    for(y in 1:nind){
      temp<-cor(t(test_resid_corr_test80[,m:(m+ncellsper[y]-1)]))
      cor_temp[y]<-mean(temp[upper.tri(temp)],na.rm = T)
      m<-m+ncellsper[y]
    }
    rho_corr_test80_ind[i]<-mean(cor_temp)
    rho<-mean(cor_temp)
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))
    p.val_corr_wilcox_test80_ref50_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test80_ref50[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_corr_wilcox_test80_ref20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test80_ref20[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_corr_wilcox_test80_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test80[i]-.5*test_size*ref_size)/sqrt(var)))
    #
    cor_temp<-numeric(length=nind)
    m<-1
    for(y in 1:nind){
      temp<-cor(t(test_resid_corr_test20[,m:(m+ncellsper[y]-1)]))
      cor_temp[y]<-mean(temp[upper.tri(temp)],na.rm = T)
      m<-m+ncellsper[y]
    }
    rho_corr_test20_ind[i]<-mean(cor_temp)
    rho<-mean(cor_temp)
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))
    p.val_corr_wilcox_test20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test20[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_corr_wilcox_test20_ref20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test20_ref20[i]-.5*test_size*ref_size)/sqrt(var)))

    
    cor_temp<-numeric(length=nind)
    m<-1
    for(y in 1:nind){
      temp<-cor(t(test_resid_corr_test0[,m:(m+ncellsper[y]-1)]))
      cor_temp[y]<-mean(temp[upper.tri(temp)],na.rm = T)
      m<-m+ncellsper[y]
    }
    rho_corr_test20_ind[i]<-mean(cor_temp)
    rho<-mean(cor_temp)
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))
    p.val_corr_wilcox_test0_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test0[i]-.5*test_size*ref_size)/sqrt(var)))
    
    
    if(i%%1==0){
      print(i)
      print(paste("Set of 1 took ",round(proc.time()[3]-time[3],2),"seconds"))
      time<-proc.time()
    }

  }
  assign(paste0("nsims_test",test_size,"_ref_",ref_size,"_sim",sim_number_alt),nsims)
  assign(paste0("gene_selected_test",test_size,"_ref_",ref_size,"_sim",sim_number_alt),gene_selected)

  assign(paste0("rho_indep_ind_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_ind)
  assign(paste0("rho_corr_ind_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_ind)

  assign(paste0("rho_corr_test50_ind_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test50_ind)
  assign(paste0("rho_corr_test20_ind_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test20_ind)
  assign(paste0("rho_corr_test80_ind_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test80_ind)

  assign(paste0("rho_corr_dat_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_dat)

  assign(paste0("rho_corr_test50_dat_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test50_dat)
  assign(paste0("rho_corr_test20_dat_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test20_dat)
  assign(paste0("rho_corr_test80_dat_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test80_dat)

  assign(paste0("rho_indep_test20_ind_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test20_ind)
  assign(paste0("rho_indep_test50_ind_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test50_ind)
  assign(paste0("rho_indep_test80_ind_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test80_ind)

  assign(paste0("rho_indep_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_camera)
  assign(paste0("rho_corr_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_camera)

  assign(paste0("rho_corr_test50_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test50_camera)
  assign(paste0("rho_corr_test20_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test20_camera)
  assign(paste0("rho_corr_test80_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test80_camera)

  assign(paste0("rho_corr_mast_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_mast)

  assign(paste0("rho_corr_test50_mast_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test50_mast)
  assign(paste0("rho_corr_test20_mast_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test20_mast)
  assign(paste0("rho_corr_test80_mast_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_corr_test80_mast)

  assign(paste0("rho_indep_test20_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test20_camera)
  assign(paste0("rho_indep_test50_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test50_camera)
  assign(paste0("rho_indep_test80_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test80_camera)

  for(t in ls(pattern="p.val")){
    obj<-get(t)
    assign(paste0(t,"_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),obj)
  }

  rm(ref_set_corr_0,ref_set_indep_0,ref_set_indep_50,ref_set_corr_50)
  rm(test_set_indep)

  assign(paste0("code_version_geneset_sim",sim_number_alt),"v11")
  rm(list=setdiff(ls(),c("sim_number_alt","test_size","ref_size",ls(pattern="_sim"))))


  rm(list=ls(pattern="residuals"))
  rm(list=ls(pattern="times"))
  setwd("~/sim/array2/pathway/final_pathway")
  #a<-ls(pattern=paste0("p.val"))
  a<-setdiff(ls(pattern="_sim"),ls(pattern="adj_stat"))

  for(i in sim_number_alt){
    #save(c,file=paste0("non_p.val_","_ref_",ref_size,"_sim",i,".RData"))
    #b<-a[grepl(paste0("sim",i),a)]
    save(list=a,file=paste0("p.val_add_idea_add_fgsea_only_testsize_",test_size,"_refsize_",ref_size,"_sim",i,"_v11.RData"))
  }

  print(gc())
}),file=paste0("/pine/scr/e/v/evanbure/output_geneset","_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt,"_allmeth.txt"))


