# set_simulation_code.R

# This script is used to accomplish several set-level tasks
# related to the TWO-SIGMA-G paper simulations comparing
# to MAST, CAMERA, and GSEA

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

# In the simulations we used both reference set size of 30 and 100
test_size<-30
ref_size<-30

# To run GSEA results
library(GSEA)
library(Biobase)
library(ArrayTools)

library(limma) #for camera function
library(MAST)# for MAST function
library(SingleCellExperiment)
library(twosigma)

# Specifies relationship between null and alternative simulations
# as identified by the simulation number
# We need to have access to both to construct our null and alternative scenarios
# For example, to construct T50 R20, take 50% of genes from the alternative
# and 50% from the null simulation in the test set, and similarly for the reference set
# 20% from the alternative and 80% from the null.
capture.output(system.time({
  if(sim_number_alt>=12300 & sim_number_alt<=12399){
    sim_number_null<-sim_number_alt-300
  }
  if(sim_number_alt>=12500 & sim_number_alt<=12599){
    sim_number_null<-sim_number_alt-300
  }
  if(sim_number_alt>=12800 & sim_number_alt<=12899){
    sim_number_null<-sim_number_alt-100
  }
  
  # number of replicates to run in this job. As in the gene-level case
  # we employ distributed computing then combine the saves from many scripts
  # See Supplementary Section S1 for more details
  nsims=10 
  print(paste("Sim ",sim_number_alt, ", nsims = ",nsims))
  print(paste("Sim Number Null",sim_number_null, ", nsims = ",nsims))
  
  # Load precomputed results for null and alternative hypotheses
  load(paste0("results_sim",sim_number_null,"_pathway.RData"))
  load(paste0("residuals_sim",sim_number_null,"_pathway.RData"))
  load(paste0("LRparams_sim",sim_number_null,"_pathway.RData"))
  ngenes_null=length(get(paste0("LR_stat_orig_null_all_sim",sim_number_null)))
  
  load(paste0("results_sim",sim_number_alt,"_pathway.RData"))
  load(paste0("residuals_sim",sim_number_alt,"_pathway.RData"))
  load(paste0("LRparams_sim",sim_number_alt,"_pathway.RData"))
  
  ngenes_alt=length(get(paste0("LR_stat_orig_null_all_sim",sim_number_alt)))
  
  nruns=get(paste0("nruns_sim",sim_number_alt))
  
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
  
  #Resimulate covariates used in TWO-SIGMA for comparison to other methods
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
  
  
  #p-values from independent and correlated gene sets for twosigmag
  p.val_indep_twosigmag<-rep(NA,nsims)
  p.val_indep_twosigmag_adjust<-rep(NA,nsims)

  p.val_indep_twosigmag_ref50<-rep(NA,nsims)
  p.val_indep_twosigmag_ref50_adjust<-rep(NA,nsims)

  p.val_indep_twosigmag_test50<-rep(NA,nsims)
  p.val_indep_twosigmag_test50_adjust<-rep(NA,nsims)
  p.val_indep_twosigmag_test50_ref50<-rep(NA,nsims)
  p.val_indep_twosigmag_test50_ref50_adjust<-rep(NA,nsims)
  p.val_indep_twosigmag_test50_ref20<-rep(NA,nsims)
  p.val_indep_twosigmag_test50_ref20_adjust<-rep(NA,nsims)

  p.val_indep_twosigmag_test80<-rep(NA,nsims)
  p.val_indep_twosigmag_test80_adjust<-rep(NA,nsims)
  p.val_indep_twosigmag_test80_ref50<-rep(NA,nsims)
  p.val_indep_twosigmag_test80_ref50_adjust<-rep(NA,nsims)
  p.val_indep_twosigmag_test80_ref20<-rep(NA,nsims)
  p.val_indep_twosigmag_test80_ref20_adjust<-rep(NA,nsims)
  
  p.val_corr_twosigmag<-rep(NA,nsims)
  p.val_corr_twosigmag_adjust<-rep(NA,nsims)
  
  p.val_corr_twosigmag_ref50<-rep(NA,nsims)
  p.val_corr_twosigmag_ref50_adjust<-rep(NA,nsims)
  
  p.val_indep_twosigmag_test50<-rep(NA,nsims)
  p.val_corr_twosigmag_test50_adjust<-rep(NA,nsims)
  p.val_corr_twosigmag_test50_ref50<-rep(NA,nsims)
  p.val_corr_twosigmag_test50_ref50_adjust<-rep(NA,nsims)
  p.val_corr_twosigmag_test50_ref20<-rep(NA,nsims)
  p.val_corr_twosigmag_test50_ref20_adjust<-rep(NA,nsims)
  
  p.val_corr_twosigmag_test50<-rep(NA,nsims)
  p.val_corr_twosigmag_test80<-rep(NA,nsims)
  p.val_corr_twosigmag_test80_adjust<-rep(NA,nsims)
  p.val_corr_twosigmag_test80_ref50<-rep(NA,nsims)
  p.val_corr_twosigmag_test80_ref50_adjust<-rep(NA,nsims)
  p.val_corr_twosigmag_test80_ref20<-rep(NA,nsims)
  p.val_corr_twosigmag_test80_ref20_adjust<-rep(NA,nsims)
  
  p.val_indep_twosigmag_test20<-rep(NA,nsims)
  p.val_corr_twosigmag_test20<-rep(NA,nsims)
  
  p.val_indep_twosigmag_test20_ref20<-rep(NA,nsims)
  p.val_corr_twosigmag_test20_ref20<-rep(NA,nsims)
  
  p.val_indep_twosigmag_test20_adjust<-rep(NA,nsims)
  p.val_corr_twosigmag_test20_adjust<-rep(NA,nsims)
  
  p.val_indep_twosigmag_test20_ref20_adjust<-rep(NA,nsims)
  p.val_corr_twosigmag_test20_ref20_adjust<-rep(NA,nsims)
  
  #p-values from independent and correlated gene sets for camera
  
  p.val_indep_camera<-rep(NA,nsims)
  p.val_indep_camera_adjust<-rep(NA,nsims)
  p.val_indep_camera_ref50<-rep(NA,nsims)
  p.val_indep_camera_ref50_adjust<-rep(NA,nsims)
  
  p.val_indep_camera_test50<-rep(NA,nsims)
  p.val_indep_camera_test50_adjust<-rep(NA,nsims)
  p.val_indep_camera_test50_ref50<-rep(NA,nsims)
  p.val_indep_camera_test50_ref50_adjust<-rep(NA,nsims)
  p.val_indep_camera_test50_ref20<-rep(NA,nsims)
  p.val_indep_camera_test50_ref20_adjust<-rep(NA,nsims)
  
  p.val_indep_camera_test80<-rep(NA,nsims)
  p.val_indep_camera_test80_adjust<-rep(NA,nsims)
  p.val_indep_camera_test80_ref50<-rep(NA,nsims)
  p.val_indep_camera_test80_ref50_adjust<-rep(NA,nsims)
  p.val_indep_camera_test80_ref20<-rep(NA,nsims)
  p.val_indep_camera_test80_ref20_adjust<-rep(NA,nsims)
  
  p.val_indep_camera_test20<-rep(NA,nsims)
  p.val_indep_camera_test20_adjust<-rep(NA,nsims)
  p.val_indep_camera_test20_ref20<-rep(NA,nsims)
  p.val_indep_camera_test20_ref20_adjust<-rep(NA,nsims)
  
  p.val_corr_camera<-rep(NA,nsims)
  p.val_corr_camera_adjust<-rep(NA,nsims)
  p.val_corr_camera_ref50<-rep(NA,nsims)
  p.val_corr_camera_ref50_adjust<-rep(NA,nsims)
  
  p.val_corr_camera_test50<-rep(NA,nsims)
  p.val_corr_camera_test50_adjust<-rep(NA,nsims)
  p.val_corr_camera_test50_ref50<-rep(NA,nsims)
  p.val_corr_camera_test50_ref50_adjust<-rep(NA,nsims)
  p.val_corr_camera_test50_ref20<-rep(NA,nsims)
  p.val_corr_camera_test50_ref20_adjust<-rep(NA,nsims)
  
  p.val_corr_camera_test80<-rep(NA,nsims)
  p.val_corr_camera_test80_adjust<-rep(NA,nsims)
  p.val_corr_camera_test80_ref50<-rep(NA,nsims)
  p.val_corr_camera_test80_ref50_adjust<-rep(NA,nsims)
  p.val_corr_camera_test80_ref20<-rep(NA,nsims)
  p.val_corr_camera_test80_ref20_adjust<-rep(NA,nsims)
  
  p.val_corr_camera_test20<-rep(NA,nsims)
  p.val_corr_camera_test20_adjust<-rep(NA,nsims)
  
  p.val_corr_camera_test20_ref20<-rep(NA,nsims)
  p.val_corr_camera_test20_ref20_adjust<-rep(NA,nsims)
  
  wilcox_stat_corr<-rep(NA,nsims)
  wilcox_stat_indep<-rep(NA,nsims)
  wilcox_stat_corr_ref50<-rep(NA,nsims)
  wilcox_stat_indep_ref50<-rep(NA,nsims)
  
  wilcox_stat_indep_test20<-rep(NA,nsims)
  wilcox_stat_indep_test20_ref20<-rep(NA,nsims)
  wilcox_stat_indep_test50<-rep(NA,nsims)
  wilcox_stat_indep_test80<-rep(NA,nsims)
  wilcox_stat_indep_test50_ref50<-rep(NA,nsims)
  wilcox_stat_indep_test50_ref20<-rep(NA,nsims)
  wilcox_stat_indep_test80_ref50<-rep(NA,nsims)
  wilcox_stat_indep_test80_ref20<-rep(NA,nsims)
  
  wilcox_stat_corr_test20<-rep(NA,nsims)
  wilcox_stat_corr_test20_ref20<-rep(NA,nsims)
  wilcox_stat_corr_test50<-rep(NA,nsims)
  wilcox_stat_corr_test80<-rep(NA,nsims)
  wilcox_stat_corr_test50_ref50<-rep(NA,nsims)
  wilcox_stat_corr_test50_ref20<-rep(NA,nsims)
  wilcox_stat_corr_test80_ref50<-rep(NA,nsims)
  wilcox_stat_corr_test80_ref20<-rep(NA,nsims)
  
  ref_set_corr_0<-vector('list',length=nsims)
  ref_set_corr_50<-vector('list',length=nsims)
  ref_set_corr_20<-vector('list',length=nsims)
  
  ref_set_indep_0<-vector('list',length=nsims)
  ref_set_indep_50<-vector('list',length=nsims)
  ref_set_indep_20<-vector('list',length=nsims)
  
  test_set_indep<-vector('list',length=nsims)
  test_set_corr<-vector('list',length=nsims)
  
  test_set_indep_test20<-vector('list',length=nsims)
  test_set_indep_test50<-vector('list',length=nsims)
  test_set_indep_test80<-vector('list',length=nsims)
  
  test_set_corr_test20<-vector('list',length=nsims)
  test_set_corr_test50<-vector('list',length=nsims)
  test_set_corr_test80<-vector('list',length=nsims)
  
  rho_indep_ind<-rep(NA,nsims)
  rho_indep_test50_ind<-rep(NA,nsims)
  rho_indep_test20_ind<-rep(NA,nsims)
  rho_indep_test80_ind<-rep(NA,nsims)
  
  rho_corr_ind<-rep(NA,nsims)
  rho_corr_test50_ind<-rep(NA,nsims)
  rho_corr_test20_ind<-rep(NA,nsims)
  rho_corr_test80_ind<-rep(NA,nsims)
  
  rho_corr_dat<-rep(NA,nsims)
  rho_corr_test50_dat<-rep(NA,nsims)
  rho_corr_test20_dat<-rep(NA,nsims)
  rho_corr_test80_dat<-rep(NA,nsims)
  
  rho_indep_camera<-rep(NA,nsims)
  rho_indep_test50_camera<-rep(NA,nsims)
  rho_indep_test20_camera<-rep(NA,nsims)
  rho_indep_test80_camera<-rep(NA,nsims)
  
  rho_corr_camera<-rep(NA,nsims)
  rho_corr_test50_camera<-rep(NA,nsims)
  rho_corr_test20_camera<-rep(NA,nsims)
  rho_corr_test80_camera<-rep(NA,nsims)
  
  # Save information from MAST
  p.val_indep_mast_adjust<-rep(NA,nsims)
  p.val_indep_mast_ref50_adjust<-rep(NA,nsims)
  p.val_indep_mast_test50_adjust<-rep(NA,nsims)
  p.val_indep_mast_test50_ref50_adjust<-rep(NA,nsims)
  p.val_indep_mast_test50_ref20_adjust<-rep(NA,nsims)
  p.val_indep_mast_test80_adjust<-rep(NA,nsims)
  p.val_indep_mast_test80_ref50_adjust<-rep(NA,nsims)
  p.val_indep_mast_test80_ref20_adjust<-rep(NA,nsims)
  p.val_indep_mast_test20_adjust<-rep(NA,nsims)
  p.val_indep_mast_test20_ref20_adjust<-rep(NA,nsims)
  
  p.val_corr_mast_adjust<-rep(NA,nsims)
  p.val_corr_mast_ref50_adjust<-rep(NA,nsims)
  p.val_corr_mast_test50_adjust<-rep(NA,nsims)
  p.val_corr_mast_test50_ref50_adjust<-rep(NA,nsims)
  p.val_corr_mast_test50_ref20_adjust<-rep(NA,nsims)
  p.val_corr_mast_test80_adjust<-rep(NA,nsims)
  p.val_corr_mast_test80_ref50_adjust<-rep(NA,nsims)
  p.val_corr_mast_test80_ref20_adjust<-rep(NA,nsims)
  p.val_corr_mast_test20_adjust<-rep(NA,nsims)
  p.val_corr_mast_test20_ref20_adjust<-rep(NA,nsims)
  
  rho_corr_mast<-matrix(NA,nrow=nsims,ncol=2)
  rho_corr_test50_mast<-matrix(NA,nrow=nsims,ncol=2)
  rho_corr_test20_mast<-matrix(NA,nrow=nsims,ncol=2)
  rho_corr_test80_mast<-matrix(NA,nrow=nsims,ncol=2)
  
  rho_indep_mast<-matrix(NA,nrow=nsims,ncol=2)
  rho_indep_test50_mast<-matrix(NA,nrow=nsims,ncol=2)
  rho_indep_test20_mast<-matrix(NA,nrow=nsims,ncol=2)
  rho_indep_test80_mast<-matrix(NA,nrow=nsims,ncol=2)
  
  if(ncol(design_camera)==2){mast_form<-~t2d_sim}
  if(ncol(design_camera)==4){mast_form<-~t2d_sim+age_sim+cdr_sim}
  
  # Save information from GSEA
  
  p.val_indep_gsea<-rep(NA,nsims)
  p.val_indep_gsea_ref50<-rep(NA,nsims)
  p.val_indep_gsea_test50<-rep(NA,nsims)
  p.val_indep_gsea_test50_ref50<-rep(NA,nsims)
  p.val_indep_gsea_test50_ref20<-rep(NA,nsims)
  p.val_indep_gsea_test80<-rep(NA,nsims)
  p.val_indep_gsea_test80_ref50<-rep(NA,nsims)
  p.val_indep_gsea_test80_ref20<-rep(NA,nsims)
  p.val_indep_gsea_test20<-rep(NA,nsims)
  p.val_indep_gsea_test20_ref20<-rep(NA,nsims)
  
  p.val_corr_gsea<-rep(NA,nsims)
  p.val_corr_gsea_ref50<-rep(NA,nsims)
  p.val_corr_gsea_test50<-rep(NA,nsims)
  p.val_corr_gsea_test50_ref50<-rep(NA,nsims)
  p.val_corr_gsea_test50_ref20<-rep(NA,nsims)
  p.val_corr_gsea_test80<-rep(NA,nsims)
  p.val_corr_gsea_test80_ref50<-rep(NA,nsims)
  p.val_corr_gsea_test80_ref20<-rep(NA,nsims)
  p.val_corr_gsea_test20<-rep(NA,nsims)
  p.val_corr_gsea_test20_ref20<-rep(NA,nsims)
    
  time<-proc.time()
  for(i in 1:nsims){
    # j is the "master gene" and corresponding test set
    # No matter what scenario evaluating we want to use this in the test set
    j<-sample(1:ngenes_alt,1)
    gene_selected[i]<-j
    
    # Randomly pick correlated genes to be in the test set for correlated case
    # will be all if test_size is 30 based on simulation setup
    index<-sample(1:nruns,test_size-1,replace=test_size>30)
    
    index_null_indep<-sample(setdiff(1:ngenes_null,j),ref_size)
    index_alt_indep<-sample(setdiff(1:ngenes_alt,j),ref_size)
    
    #perfectly competitive ref set for independent case
    ref_data_0_indep<-dat_orig_null_all[index_null_indep,]
    ref_set_indep_0[[i]]<-LR_stat_orig_null_all[index_null_indep]
    
    #polluted reference set for independent case
    # Half of genes in reference under alt and half under null
    index_alt_ref50<-sample(index_alt_indep,round(.5*ref_size,0))
    index_null_ref50<-sample(index_null_indep,ref_size-round(.5*ref_size,0))
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
    index_corr_ref<-sample(setdiff(1:ngenes_alt,j),ref_size)
    
    #index_corr_ref0<-sample(index_corr_ref,ref_size)
    ref_data_corr_0<-matrix(nrow=ref_size,ncol=sum(ncellsper))
    temp<-numeric(length=ref_size)
    for(d in 1:ref_size){
      g<-index_corr_ref[d]
      # From the randomly sampled set choose a random gene to use
      # This preserves the general structure from the "correlated" genes but ensures we do not have correlated genes being chosen
      h<-sample(1:nruns,1)
      # Load data
      load(paste0("/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
      temp[d]<-LR_stat_large_noise_null_all[[g]][h]
      ref_data_corr_0[d,]<-test_dat[h,]
    }
    ref_set_corr_0[[i]]<-temp
    
    #Polluted Competitive Reference Set made of genes added with noise
    # would still like to preserve the independence in the reference set though
    # so be sure no correlation exists in reference set
    #index_corr_ref50<-sample(setdiff(1:ngenes_alt,j),ref_size)
    ref_data_50_corr<-matrix(nrow=ref_size,ncol=sum(ncellsper))
    temp<-numeric(length=ref_size)
    for(d in 1:ref_size){
      g<-index_corr_ref[d]
      # For the simulated set, randomly choose a given gene to include
      # Again we want to get data that is simulated like the "correlated" genes
      # But ensure that we do not have correlation in the reference set
      h<-sample(1:nruns,1)
      if(d<=.5*ref_size){
        load(paste0("/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_alt_all[[g]][h]
      }
      if(d>.5*ref_size){
        load(paste0("/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_null_all[[g]][h]
      }
      ref_data_50_corr[d,]<-test_dat[h,]
    }
    ref_set_corr_50[[i]]<-temp
    
    #index_corr_ref20<-sample(setdiff(1:ngenes_alt,j),ref_size)
    ref_data_20_corr<-matrix(nrow=ref_size,ncol=sum(ncellsper))
    temp<-numeric(length=ref_size)
    for(d in 1:ref_size){
      g<-index_corr_ref[d]
      # For the simulated set, randomly choose a given gene to include
      # Again we want to get data that is simulated like the "correlated" genes
      # But ensure that we do not have correlation in the reference set
      h<-sample(1:nruns,1)
      if(d<=.2*ref_size){
        load(paste0("/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_alt_all[[g]][h]
      }
      if(d>.2*ref_size){
        load(paste0("/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_null_all[[g]][h]
      }
      ref_data_20_corr[d,]<-test_dat[h,]
    }
    ref_set_corr_20[[i]]<-temp
    
    #Randomly pick test set for indep case from statstics from independent genes
    # Being sure to exclude original and genes that are in any reference set (index_alt)
    index_alt_test<-sample(setdiff(1:ngenes_alt,c(j,index_alt_indep)),test_size-1)
    
    # First set of statistics combines the original gene's statistic with either independent gene statistics
    
    test_set_indep[[i]]<-LR_stat_orig_alt_all[c(j,index_alt_test)]
    
    test_set_corr[[i]]<-c(LR_stat_orig_alt_all[j],LR_stat_large_noise_alt_all[[j]][index])
    
    max_in_test_set_corr[i]<-max(c(test_set_corr[[i]],ref_set_corr_0[[i]]))%in%test_set_corr[[i]]
    max_in_test_set_indep[i]<-max(c(test_set_indep[[i]],ref_set_corr_0[[i]]))%in%test_set_indep[[i]]
    
    test_data_indep_100<-dat_orig_alt_all[c(j,index_alt_test),]
    test_resid_indep<-residuals_orig_alt_all[c(j,index_alt_test),]
    
    #loading this object will have the name test_dat
    load(paste0("/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",j,".RData"))
    test_data_corr_100<-rbind(matrix(dat_orig_alt_all[j,],nrow=1),test_dat)
    cor_temp<-cor(t(test_data_corr_100))
    rho_corr_dat[i]<-mean(cor_temp[upper.tri(cor_temp)])
    
    #Test set that is not perfectly DE
    test_data_corr_test50<-matrix(nrow=test_size,ncol=sum(ncellsper))
    test_data_corr_test50[1,]<-dat_orig_alt_all[j,]
    
    test_resid_corr_test50<-matrix(nrow=test_size,ncol=sum(ncellsper))
    test_resid_corr_test50[1,]<-residuals_orig_alt_all[j,]
    
    # gives the null genes that will be used 
    # Being sure to remove any gene that could have been used in a reference set
    index4<-sample(setdiff(1:ngenes_null,c(j,index_corr_ref)),test_size)
    temp<-numeric(length=test_size)
    
    # Gives the correlated genes that will be used
    # sample all possibilities (1:nruns) but will use less than all of them
    # depending on desired composition of test set of null and alternative genes
    # always pull 
    index4a<-sample(1:nruns,test_size-1)
    temp[1]<-LR_stat_orig_alt_all[j]
    for(d in 2:(test_size)){
      if(d<=.5*test_size){
        g<-j
        h<-index4a[d]
        load(paste0("/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_alt_all[[g]][h]
        load(paste0("/raw_dat/sim",sim_number_alt,"/residuals_large_noise_sim",sim_number_alt,"_",g,"_",h,".RData"))
        test_resid_corr_test50[d,]<-residuals_large_noise
      }
      if(d>.5*test_size){# use null data
        g<-index4[d]
        h<-sample(1:nruns,1)
        load(paste0("/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_null_all[[g]][h]
        load(paste0("/raw_dat/sim",sim_number_null,"/residuals_large_noise_sim",sim_number_null,"_",g,"_",h,".RData"))
        test_resid_corr_test50[d,]<-residuals_large_noise
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
    
    # gives the null genes that will be used 
    index4<-sample(setdiff(1:ngenes_null,c(j,index_corr_ref)),test_size)
    temp<-numeric(length=test_size)
    
    # Gives the correlated genes that will be used
    # sample all possibilities (1:nruns) but will use less than all of them
    # depending on desired composition of test set of null and alternative genes
    index4a<-sample(1:nruns,test_size-1)
    temp[1]<-LR_stat_orig_alt_all[j]
    for(d in 2:(test_size)){
      if(d<=.8*test_size){
        g<-j
        h<-index4a[d]
        load(paste0("/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_alt_all[[g]][h]
        load(paste0("/raw_dat/sim",sim_number_alt,"/residuals_large_noise_sim",sim_number_alt,"_",g,"_",h,".RData"))
        test_resid_corr_test80[d,]<-residuals_large_noise
      }
      if(d>.8*test_size){ # use null data
        g<-index4[d]
        h<-sample(1:nruns,1)
        load(paste0("/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_null_all[[g]][h]
        load(paste0("/raw_dat/sim",sim_number_null,"/residuals_large_noise_sim",sim_number_null,"_",g,"_",h,".RData"))
        test_resid_corr_test80[d,]<-residuals_large_noise
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
    
    # gives the null genes that will be used 
    index4<-sample(setdiff(1:ngenes_null,c(j,index_corr_ref)),test_size)
    temp<-numeric(length=test_size)
    
    # Gives the correlated genes that will be used
    # sample all possibilities (1:nruns) but will use less than all of them
    # depending on desired composition of test set of null and alternative genes
    index4a<-sample(1:nruns,test_size-1)
    temp[1]<-LR_stat_orig_alt_all[j]
    for(d in 2:(test_size)){
      if(d<=.2*test_size){
        g<-j
        h<-index4a[d]
        load(paste0("/raw_dat/sim",sim_number_alt,"/dat_large_noise_null_all_sim",sim_number_alt,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_alt_all[[g]][h]
        load(paste0("/raw_dat/sim",sim_number_alt,"/residuals_large_noise_sim",sim_number_alt,"_",g,"_",h,".RData"))
        test_resid_corr_test20[d,]<-residuals_large_noise
      }
      if(d>.2*test_size){# use null data
        g<-index4[d]
        h<-sample(1:nruns,1)
        load(paste0("/raw_dat/sim",sim_number_null,"/dat_large_noise_null_all_sim",sim_number_null,"_",g,".RData"))
        temp[d]<-LR_stat_large_noise_null_all[[g]][h]
        load(paste0("/raw_dat/sim",sim_number_null,"/residuals_large_noise_sim",sim_number_null,"_",g,"_",h,".RData"))
        test_resid_corr_test20[d,]<-residuals_large_noise
      }
      test_data_corr_test20[d,]<-test_dat[h,]
    }
    test_set_corr_test20[[i]]<-temp
    cor_temp<-cor(t(test_data_corr_test20))
    rho_corr_test20_dat[i]<-mean(cor_temp[upper.tri(cor_temp)])
    
    index5a<-sample(setdiff(1:ngenes_alt,c(j,index_alt_indep)),size=test_size/2-1)
    index5<-sample(setdiff(1:ngenes_null,c(j,index_null_indep)),size=test_size/2)
    
    test_set_indep_test50[[i]]<-c(LR_stat_orig_alt_all[c(j,index5a)],LR_stat_orig_null_all[index5])
    test_data_indep_test50<-rbind(dat_orig_alt_all[c(j,index5a),],dat_orig_null_all[index5,])
    test_resid_indep_test50<-rbind(residuals_orig_alt_all[c(j,index5a),],residuals_orig_null_all[index5,])
    
    index5a<-sample(setdiff(1:ngenes_alt,c(j,index_alt_indep)),size=.8*test_size-1)
    index5<-sample(setdiff(1:ngenes_null,c(j,index_null_indep)),size=test_size-.8*test_size)
    
    test_set_indep_test80[[i]]<-c(LR_stat_orig_alt_all[c(j,index5a)],LR_stat_orig_null_all[index5])
    test_data_indep_test80<-rbind(dat_orig_alt_all[c(j,index5a),],dat_orig_null_all[index5,])
    test_resid_indep_test80<-rbind(residuals_orig_alt_all[c(j,index5a),],residuals_orig_null_all[index5,])
    
    index5a<-sample(setdiff(1:ngenes_alt,c(j,index_alt_indep)),size=.2*test_size-1)
    index5<-sample(setdiff(1:ngenes_null,c(j,index_null_indep)),size=test_size-.2*test_size)
    
    test_set_indep_test20[[i]]<-c(LR_stat_orig_alt_all[c(j,index5a)],LR_stat_orig_null_all[index5])
    test_data_indep_test20<-rbind(dat_orig_alt_all[c(j,index5a),],dat_orig_null_all[index5,])
    test_resid_indep_test20<-rbind(residuals_orig_alt_all[c(j,index5a),],residuals_orig_null_all[index5,])
###############################################################
#Run GSEA for all scenarios
###############################################################
    tryCatch({
    test_dat<-rbind(test_data_corr_100,ref_data_50_corr)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")

    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_corr_gsea_ref50[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_corr_100,ref_data_corr_0)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")

    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_corr_gsea[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_corr_test80,ref_data_corr_0)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_corr_gsea_test80[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_corr_test80,ref_data_20_corr)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_corr_gsea_test80_ref20[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_corr_test80,ref_data_50_corr)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_corr_gsea_test80_ref50[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_corr_test50,ref_data_corr_0)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_corr_gsea_test50[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_corr_test50,ref_data_20_corr)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_corr_gsea_test50_ref20[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_corr_test20,ref_data_corr_0)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_corr_gsea_test20[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_corr_test20,ref_data_20_corr)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_corr_gsea_test20_ref20[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_corr_test50,ref_data_50_corr)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_corr_gsea_test50_ref50[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_indep_100,ref_data_50_indep)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_indep_gsea_ref50[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_indep_100,ref_data_indep_0)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_indep_gsea[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_indep_test80,ref_data_indep_0)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_indep_gsea_test80[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_indep_test80,ref_data_20_indep)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_indep_gsea_test80_ref20[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_indep_test80,ref_data_50_indep)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_indep_gsea_test80_ref50[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_indep_test50,ref_data_indep_0)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_indep_gsea_test50[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_indep_test50,ref_data_20_indep)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_indep_gsea_test50_ref20[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_indep_test20,ref_data_indep_0)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_indep_gsea_test20[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_indep_test20,ref_data_20_indep)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_indep_gsea_test20_ref20[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    tryCatch({
    test_dat<-rbind(test_data_indep_test50,ref_data_50_indep)
    colnames(test_dat)<-paste0("Cell",1:length(t2d_sim))
    rownames(test_dat)<-paste0("Gene",1:(test_size+ref_size))
    phenDat<-data.frame(t2d_sim)
    rownames(phenDat)<-paste0("Cell",1:length(t2d_sim))
    phenDat2<-AnnotatedDataFrame(phenDat)
    temp<-ExpressionSet(test_dat,phenoData=phenDat2)
    t2d_sim_cat<-as.character(t2d_sim)
    createGSEAFiles(mydir=paste0("/gsea_sim",sim_number_alt),eSet=temp,catVar="t2d_sim")
    
    gs3<-list(c(paste0("Pathway_",sim_number_alt, "_1"),"URL",paste0("Gene",1:test_size)))
    temp<-GSEA(input.ds=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.probe.gct"),input.cls=paste0("/gsea_sim",sim_number_alt,"/Result/Pathway_Analysis/GSEA/t2d_sim.phenotype.cls"),gs.db=gs3,output.directory = "",gsea.type="GSEA",collapse.dataset = FALSE,nperm=200)
    p.val_indep_gsea_test50_ref50[i]<-as.numeric(as.character(temp$report1$`NOM p-val`))},error=function(e){})
    
    ##################################################################
    # Run MAST
    ##################################################################
    
    suppressMessages({
      coldat<-as.data.frame(cbind(t2d_sim,age_sim,cdr_sim))
      rownames(coldat)<-colnames(test_data_corr_100)
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_corr_100,ref_data_corr_0)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      boots_corr <- bootVcov1(mast_model_corr, 25)
      mast_gsea_corr<-gseaAfterBoot(mast_model_corr,boots_corr,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_corr<-summary(mast_gsea_corr)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_corr_mast_adjust[i]<-mast_sum_corr$combined_P[1]
    # Only keep correlation from first case of each test composition
    rho_corr_mast[i,]<-mast_gsea_corr@tests[1,,'avgCor','test']

    suppressMessages({
      coldat<-as.data.frame(cbind(t2d_sim,age_sim,cdr_sim))
      rownames(coldat)<-colnames(test_data_corr_100)
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_corr_100,ref_data_50_corr)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      boots_corr <- bootVcov1(mast_model_corr, 25)
      mast_gsea_corr<-gseaAfterBoot(mast_model_corr,boots_corr,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_corr<-summary(mast_gsea_corr)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_corr_mast_ref50_adjust[i]<-mast_sum_corr$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_corr_test50,ref_data_corr_0)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      boots_corr <- bootVcov1(mast_model_corr, 25)
      mast_gsea_corr<-gseaAfterBoot(mast_model_corr,boots_corr,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_corr<-summary(mast_gsea_corr)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    rho_corr_test50_mast[i,]<-mast_gsea_corr@tests[1,,'avgCor','test']
    p.val_corr_mast_test50_adjust[i]<-mast_sum_corr$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_corr_test50,ref_data_50_corr)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      boots_corr <- bootVcov1(mast_model_corr, 25)
      mast_gsea_corr<-gseaAfterBoot(mast_model_corr,boots_corr,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_corr<-summary(mast_gsea_corr)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_corr_mast_test50_ref50_adjust[i]<-mast_sum_corr$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_corr_test50,ref_data_20_corr)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      boots_corr <- bootVcov1(mast_model_corr, 25)
      mast_gsea_corr<-gseaAfterBoot(mast_model_corr,boots_corr,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_corr<-summary(mast_gsea_corr)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_corr_mast_test50_ref20_adjust[i]<-mast_sum_corr$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_corr_test80,ref_data_corr_0)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      boots_corr <- bootVcov1(mast_model_corr, 25)
      mast_gsea_corr<-gseaAfterBoot(mast_model_corr,boots_corr,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_corr<-summary(mast_gsea_corr)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    rho_corr_test80_mast[i,]<-mast_gsea_corr@tests[1,,'avgCor','test']
    p.val_corr_mast_test80_adjust[i]<-mast_sum_corr$combined_P[1]


    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_corr_test80,ref_data_50_corr)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      boots_corr <- bootVcov1(mast_model_corr, 25)
      mast_gsea_corr<-gseaAfterBoot(mast_model_corr,boots_corr,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_corr<-summary(mast_gsea_corr)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_corr_mast_test80_ref50_adjust[i]<-mast_sum_corr$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_corr_test80,ref_data_20_corr)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      boots_corr <- bootVcov1(mast_model_corr, 25)
      mast_gsea_corr<-gseaAfterBoot(mast_model_corr,boots_corr,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_corr<-summary(mast_gsea_corr)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_corr_mast_test80_ref20_adjust[i]<-mast_sum_corr$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_corr_test20,ref_data_corr_0)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      boots_corr <- bootVcov1(mast_model_corr, 25)
      mast_gsea_corr<-gseaAfterBoot(mast_model_corr,boots_corr,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_corr<-summary(mast_gsea_corr)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    rho_corr_test20_mast[i,]<-mast_gsea_corr@tests[1,,'avgCor','test']
    p.val_corr_mast_test20_adjust[i]<-mast_sum_corr$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_corr_test20,ref_data_20_corr)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_corr<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
        check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_corr<-zlm(mast_form,sca=mast_dat_corr,parallel = F)
      boots_corr <- bootVcov1(mast_model_corr, 25)
      mast_gsea_corr<-gseaAfterBoot(mast_model_corr,boots_corr,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_corr<-summary(mast_gsea_corr)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_corr_mast_test20_ref20_adjust[i]<-mast_sum_corr$combined_P[1]

    suppressMessages({
      coldat<-as.data.frame(cbind(t2d_sim,age_sim,cdr_sim))
      rownames(coldat)<-colnames(test_data_indep_100)
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_indep_100,ref_data_0_indep)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_indep<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
                                           check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_indep<-zlm(mast_form,sca=mast_dat_indep,parallel = F)
      boots_indep <- bootVcov1(mast_model_indep, 25)
      mast_gsea_indep<-gseaAfterBoot(mast_model_indep,boots_indep,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_indep<-summary(mast_gsea_indep)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_indep_mast_adjust[i]<-mast_sum_indep$combined_P[1]
    # Only keep indepelation from first case of each test composition
    rho_indep_mast[i,]<-mast_gsea_indep@tests[1,,'avgCor','test']

    suppressMessages({
      coldat<-as.data.frame(cbind(t2d_sim,age_sim,cdr_sim))
      rownames(coldat)<-colnames(test_data_indep_100)
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_indep_100,ref_data_50_indep)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_indep<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
                                           check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_indep<-zlm(mast_form,sca=mast_dat_indep,parallel = F)
      boots_indep <- bootVcov1(mast_model_indep, 25)
      mast_gsea_indep<-gseaAfterBoot(mast_model_indep,boots_indep,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_indep<-summary(mast_gsea_indep)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_indep_mast_ref50_adjust[i]<-mast_sum_indep$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_indep_test50,ref_data_0_indep)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_indep<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
                                           check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_indep<-zlm(mast_form,sca=mast_dat_indep,parallel = F)
      boots_indep <- bootVcov1(mast_model_indep, 25)
      mast_gsea_indep<-gseaAfterBoot(mast_model_indep,boots_indep,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_indep<-summary(mast_gsea_indep)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    rho_indep_test50_mast[i,]<-mast_gsea_indep@tests[1,,'avgCor','test']
    p.val_indep_mast_test50_adjust[i]<-mast_sum_indep$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_indep_test50,ref_data_50_indep)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_indep<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
                                           check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_indep<-zlm(mast_form,sca=mast_dat_indep,parallel = F)
      boots_indep <- bootVcov1(mast_model_indep, 25)
      mast_gsea_indep<-gseaAfterBoot(mast_model_indep,boots_indep,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_indep<-summary(mast_gsea_indep)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_indep_mast_test50_ref50_adjust[i]<-mast_sum_indep$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_indep_test50,ref_data_20_indep)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_indep<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
                                           check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_indep<-zlm(mast_form,sca=mast_dat_indep,parallel = F)
      boots_indep <- bootVcov1(mast_model_indep, 25)
      mast_gsea_indep<-gseaAfterBoot(mast_model_indep,boots_indep,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_indep<-summary(mast_gsea_indep)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_indep_mast_test50_ref20_adjust[i]<-mast_sum_indep$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_indep_test80,ref_data_0_indep)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_indep<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
                                           check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_indep<-zlm(mast_form,sca=mast_dat_indep,parallel = F)
      boots_indep <- bootVcov1(mast_model_indep, 25)
      mast_gsea_indep<-gseaAfterBoot(mast_model_indep,boots_indep,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_indep<-summary(mast_gsea_indep)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    rho_indep_test80_mast[i,]<-mast_gsea_indep@tests[1,,'avgCor','test']
    p.val_indep_mast_test80_adjust[i]<-mast_sum_indep$combined_P[1]


    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_indep_test80,ref_data_50_indep)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_indep<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
                                           check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_indep<-zlm(mast_form,sca=mast_dat_indep,parallel = F)
      boots_indep <- bootVcov1(mast_model_indep, 25)
      mast_gsea_indep<-gseaAfterBoot(mast_model_indep,boots_indep,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_indep<-summary(mast_gsea_indep)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_indep_mast_test80_ref50_adjust[i]<-mast_sum_indep$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_indep_test80,ref_data_20_indep)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_indep<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
                                           check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_indep<-zlm(mast_form,sca=mast_dat_indep,parallel = F)
      boots_indep <- bootVcov1(mast_model_indep, 25)
      mast_gsea_indep<-gseaAfterBoot(mast_model_indep,boots_indep,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_indep<-summary(mast_gsea_indep)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_indep_mast_test80_ref20_adjust[i]<-mast_sum_indep$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_indep_test20,ref_data_0_indep)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_indep<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
                                           check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_indep<-zlm(mast_form,sca=mast_dat_indep,parallel = F)
      boots_indep <- bootVcov1(mast_model_indep, 25)
      mast_gsea_indep<-gseaAfterBoot(mast_model_indep,boots_indep,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_indep<-summary(mast_gsea_indep)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    rho_indep_test20_mast[i,]<-mast_gsea_indep@tests[1,,'avgCor','test']
    p.val_indep_mast_test20_adjust[i]<-mast_sum_indep$combined_P[1]

    suppressMessages({
      sce<- SingleCellExperiment(assays = list(counts = log2(as.matrix(rbind(test_data_indep_test20,ref_data_20_indep)+1))),colData=coldat)
      rownames(sce)<-paste("Gene",1:(test_size+ref_size))
      mast_dat_indep<-SceToSingleCellAssay(sce, class = "SingleCellAssay",
                                           check_sanity = TRUE)

      mast_time_temp<-Sys.time()
      mast_model_indep<-zlm(mast_form,sca=mast_dat_indep,parallel = F)
      boots_indep <- bootVcov1(mast_model_indep, 25)
      mast_gsea_indep<-gseaAfterBoot(mast_model_indep,boots_indep,list("Set 1" = 1:test_size,"Set 2" = 1:test_size),CoefficientHypothesis("t2d_sim"))
      mast_sum_indep<-summary(mast_gsea_indep)})

    #mast_times[i]<-Sys.time()-mast_time_temp
    p.val_indep_mast_test20_ref20_adjust[i]<-mast_sum_indep$combined_P[1]
    
    
    
###########################################################################################
#Run CAMERA   
###########################################################################################
    
    p.val_indep_camera[i]<-camera(log2(rbind(test_data_indep_100,ref_data_0_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue

    camera_temp<-camera(log2(rbind(test_data_indep_100,ref_data_0_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)

    p.val_indep_camera_adjust[i]<-camera_temp$PValue
    rho_indep_camera[i]<-camera_temp$Corr
    
    
    camera_temp<-camera(log2(rbind(test_data_corr_100,ref_data_corr_0)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)
    p.val_corr_camera_adjust[i]<-camera_temp$PValue
    rho_corr_camera[i]<-camera_temp$Corr
    
    #ref50--polluted competitive alternative
    
    p.val_indep_camera_ref50[i]<-camera(log2(rbind(test_data_indep_100,ref_data_50_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue

    p.val_indep_camera_ref50_adjust[i]<-camera(log2(rbind(test_data_indep_100,ref_data_50_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue
    
    #p.val_corr_camera_ref50[i]<-camera(log2(rbind(test_data_corr_100,ref_data_50_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    
    p.val_corr_camera_ref50_adjust[i]<-camera(log2(rbind(test_data_corr_100,ref_data_50_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue
    
    # Combined Test Set
    
    p.val_indep_camera_test50[i]<-camera(log2(rbind(test_data_indep_test50,ref_data_0_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue

    camera_temp<-camera(log2(rbind(test_data_indep_test50,ref_data_0_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)
    p.val_indep_camera_test50_adjust[i]<-camera_temp$PValue
    rho_indep_test50_camera[i]<-camera_temp$Corr

    p.val_indep_camera_test50_ref50[i]<-camera(log2(rbind(test_data_indep_test50,ref_data_50_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue

    p.val_indep_camera_test50_ref50_adjust[i]<-camera(log2(rbind(test_data_indep_test50,ref_data_50_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue

    p.val_indep_camera_test50_ref20[i]<-camera(log2(rbind(test_data_indep_test50,ref_data_20_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue

    p.val_indep_camera_test50_ref20_adjust[i]<-camera(log2(rbind(test_data_indep_test50,ref_data_20_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue

    p.val_indep_camera_test80_ref20[i]<-camera(log2(rbind(test_data_indep_test80,ref_data_20_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    # 
    # camera_temp<-camera(log2(rbind(test_data_indep_test80,ref_data_20_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)
    p.val_indep_camera_test80_ref20_adjust[i]<-camera_temp$PValue
    rho_indep_test80_camera[i]<-camera_temp$Corr

    p.val_indep_camera_test80[i]<-camera(log2(rbind(test_data_indep_test80,ref_data_0_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue

    p.val_indep_camera_test80_adjust[i]<-camera(log2(rbind(test_data_indep_test80,ref_data_0_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue

    p.val_indep_camera_test80_ref50[i]<-camera(log2(rbind(test_data_indep_test80,ref_data_50_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue

    p.val_indep_camera_test80_ref50_adjust[i]<-camera(log2(rbind(test_data_indep_test80,ref_data_50_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue
    
    #p.val_corr_camera_test50[i]<-camera(log2(rbind(test_data_corr_test50,ref_data_corr_0)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    
    camera_temp<-camera(log2(rbind(test_data_corr_test50,ref_data_corr_0)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)
    p.val_corr_camera_test50_adjust[i]<-camera_temp$PValue
    rho_corr_test50_camera[i]<-camera_temp$Corr
    
    #p.val_corr_camera_test50_ref50[i]<-camera(log2(rbind(test_data_corr_test50,ref_data_50_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    
    p.val_corr_camera_test50_ref50_adjust[i]<-camera(log2(rbind(test_data_corr_test50,ref_data_50_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue
    
    #p.val_corr_camera_test50_ref20[i]<-camera(log2(rbind(test_data_corr_test50,ref_data_20_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    
    p.val_corr_camera_test50_ref20_adjust[i]<-camera(log2(rbind(test_data_corr_test50,ref_data_20_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue
    
    camera_temp<-camera(log2(rbind(test_data_corr_test80,ref_data_20_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)
    p.val_corr_camera_test80_ref20_adjust[i]<-camera_temp$PValue
    rho_corr_test80_camera[i]<-camera_temp$Corr
    
    #p.val_corr_camera_test80_ref20[i]<-camera(log2(rbind(test_data_corr_test80,ref_data_20_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    
    #p.val_corr_camera_test80[i]<-camera(log2(rbind(test_data_corr_test80,ref_data_corr_0)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    
    p.val_corr_camera_test80_adjust[i]<-camera(log2(rbind(test_data_corr_test80,ref_data_corr_0)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue
    
    #p.val_corr_camera_test80_ref50[i]<-camera(log2(rbind(test_data_corr_test80,ref_data_50_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    
    p.val_corr_camera_test80_ref50_adjust[i]<-camera(log2(rbind(test_data_corr_test80,ref_data_50_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue
    
     p.val_indep_camera_test20[i]<-camera(log2(rbind(test_data_indep_test20,ref_data_0_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    # 
     camera_temp<-camera(log2(rbind(test_data_indep_test20,ref_data_0_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)
     p.val_indep_camera_test20_adjust[i]<-camera_temp$PValue
     rho_indep_test20_camera[i]<-camera_temp$Corr
    # 
     p.val_indep_camera_test20_ref20[i]<-camera(log2(rbind(test_data_indep_test20,ref_data_20_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    
    p.val_indep_camera_test20_ref20_adjust[i]<-camera(log2(rbind(test_data_indep_test20,ref_data_20_indep)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue
    
    #p.val_corr_camera_test20[i]<-camera(log2(rbind(test_data_corr_test20,ref_data_corr_0)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    
    camera_temp<-camera(log2(rbind(test_data_corr_test20,ref_data_corr_0)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)
    p.val_corr_camera_test20_adjust[i]<-camera_temp$PValue
    rho_corr_test20_camera[i]<-camera_temp$Corr
    
    #p.val_corr_camera_test20_ref20[i]<-camera(log2(rbind(test_data_corr_test20,ref_data_20_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = 0)$PValue
    
    p.val_corr_camera_test20_ref20_adjust[i]<-camera(log2(rbind(test_data_corr_test20,ref_data_20_corr)+1),index=1:test_size,design=design_camera,contrast=2,inter.gene.cor = NULL)$PValue
    
################################################################################################################
# Run stripped down version of twosigmag
# done in this way because gene-level models were fit previously
# so a waste of computation to duplicate
################################################################################################################
    
    wilcox_stat_corr[i]<-sum(rank(c(test_set_corr[[i]],ref_set_corr_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    wilcox_stat_corr_test20[i]<-sum(rank(c(test_set_corr_test20[[i]],ref_set_corr_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    wilcox_stat_corr_test20_ref20[i]<-sum(rank(c(test_set_corr_test20[[i]],ref_set_corr_20[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    wilcox_stat_corr_ref50[i]<-sum(rank(c(test_set_corr[[i]],ref_set_corr_50[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    
    wilcox_stat_indep[i]<-sum(rank(c(test_set_indep[[i]],ref_set_indep_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_indep_test20[i]<-sum(rank(c(test_set_indep_test20[[i]],ref_set_indep_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_indep_test20_ref20[i]<-sum(rank(c(test_set_indep_test20[[i]],ref_set_indep_20[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_indep_ref50[i]<-sum(rank(c(test_set_indep[[i]],ref_set_indep_50[[i]]))[1:test_size]) - .5*test_size*(test_size+1)


    wilcox_stat_indep_test50[i]<-sum(rank(c(test_set_indep_test50[[i]],ref_set_indep_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_indep_test50_ref50[i]<-sum(rank(c(test_set_indep_test50[[i]],ref_set_indep_50[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_indep_test50_ref20[i]<-sum(rank(c(test_set_indep_test50[[i]],ref_set_indep_20[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_indep_test80_ref20[i]<-sum(rank(c(test_set_indep_test80[[i]],ref_set_indep_20[[i]]))[1:test_size]) - .5*test_size*(test_size+1)


    wilcox_stat_indep_test80[i]<-sum(rank(c(test_set_indep_test80[[i]],ref_set_indep_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)

    wilcox_stat_indep_test80_ref50[i]<-sum(rank(c(test_set_indep_test80[[i]],ref_set_indep_50[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    
    wilcox_stat_corr_test50[i]<-sum(rank(c(test_set_corr_test50[[i]],ref_set_corr_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    wilcox_stat_corr_test50_ref50[i]<-sum(rank(c(test_set_corr_test50[[i]],ref_set_corr_50[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    wilcox_stat_corr_test50_ref20[i]<-sum(rank(c(test_set_corr_test50[[i]],ref_set_corr_20[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    wilcox_stat_corr_test80_ref20[i]<-sum(rank(c(test_set_corr_test80[[i]],ref_set_corr_20[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    
    wilcox_stat_corr_test80[i]<-sum(rank(c(test_set_corr_test80[[i]],ref_set_corr_0[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    wilcox_stat_corr_test80_ref50[i]<-sum(rank(c(test_set_corr_test80[[i]],ref_set_corr_50[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    
    # Wilcoxon p-values for unadjusted
    
    p.val_indep_twosigmag[i]<-2*pnorm(-1*abs((wilcox_stat_indep[i]-.5*test_size*ref_size)/sqrt(var0)))
    
    p.val_indep_twosigmag_test20[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test20[i]-.5*test_size*ref_size)/sqrt(var0)))
    
    p.val_indep_twosigmag_test20_ref20[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test20_ref20[i]-.5*test_size*ref_size)/sqrt(var0)))
    
    p.val_indep_twosigmag_ref50[i]<-2*pnorm(-1*abs((wilcox_stat_indep_ref50[i]-.5*test_size*ref_size)/sqrt(var0)))

    
     p.val_indep_twosigmag_test50_ref50[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test50_ref50[i]-.5*test_size*ref_size)/sqrt(var0)))
    p.val_indep_twosigmag_test50_ref20[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test50_ref20[i]-.5*test_size*ref_size)/sqrt(var0)))
    p.val_indep_twosigmag_test80_ref20[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test80_ref20[i]-.5*test_size*ref_size)/sqrt(var0)))
    p.val_indep_twosigmag_test50[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test50[i]-.5*test_size*ref_size)/sqrt(var0)))
    
    p.val_indep_twosigmag_test80_ref50[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test80_ref50[i]-.5*test_size*ref_size)/sqrt(var0)))
    p.val_indep_twosigmag_test80[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test80[i]-.5*test_size*ref_size)/sqrt(var0)))
    
    
    
    #####################################################################################################
    # Compute correlation estimate using residuals and get twosigmag equivalent p-value
    ######################################################################################################
    
    
    cor_temp<-numeric(length=nind)
    m<-1
    for(y in 1:nind){
      temp<-cor(t(test_resid_indep[,m:(m+ncellsper[y]-1)]))
      cor_temp[y]<-mean(temp[upper.tri(temp)],na.rm = T)
      m<-m+ncellsper[y]
    }
    rho_indep_ind[i]<-mean(cor_temp)
    rho<-mean(cor_temp)
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))
    adj_stat_indep[i]<-(wilcox_stat_indep[i]-(.5*ref_size*test_size))/sqrt(var)
    p.val_indep_twosigmag_adjust[i]<-2*(pnorm(-1*abs(adj_stat_indep[i])))

    adj_stat_indep[i]<-(wilcox_stat_indep_ref50[i]-(.5*ref_size*test_size))/sqrt(var)
    p.val_indep_twosigmag_ref50_adjust[i]<-2*(pnorm(-1*abs(adj_stat_indep[i])))
    
    # Now correlated gene sets
    
    rho<-cor_ind[j]
    rho_corr_ind[i]<-rho
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))
    
    adj_stat_corr[i]<-(wilcox_stat_corr[i]-(.5*ref_size*test_size))/sqrt(var)
    p.val_corr_twosigmag_adjust[i]<-2*(pnorm(-1*abs(adj_stat_corr[i])))
    
    adj_stat_corr[i]<-(wilcox_stat_corr_ref50[i]-(.5*ref_size*test_size))/sqrt(var)
    p.val_corr_twosigmag_ref50_adjust[i]<-2*(pnorm(-1*abs(adj_stat_corr[i])))
    
    cor_temp<-numeric(length=nind)
    m<-1
    for(y in 1:nind){
      temp<-cor(t(test_resid_indep_test50[,m:(m+ncellsper[y]-1)]))
      cor_temp[y]<-mean(temp[upper.tri(temp)],na.rm = T)
      m<-m+ncellsper[y]
    }
    rho_indep_test50_ind[i]<-mean(cor_temp)
    rho<-mean(cor_temp)
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))
    p.val_indep_twosigmag_test50_ref50_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test50_ref50[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_indep_twosigmag_test50_ref20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test50_ref20[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_indep_twosigmag_test50_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test50[i]-.5*test_size*ref_size)/sqrt(var)))
    
    cor_temp<-numeric(length=nind)
    m<-1
    for(y in 1:nind){
      temp<-cor(t(test_resid_indep_test80[,m:(m+ncellsper[y]-1)]))
      cor_temp[y]<-mean(temp[upper.tri(temp)],na.rm = T)
      m<-m+ncellsper[y]
    }
    rho_indep_test80_ind[i]<-mean(cor_temp)
    rho<-mean(cor_temp)
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))
    p.val_indep_twosigmag_test80_ref50_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test80_ref50[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_indep_twosigmag_test80_ref20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test80_ref20[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_indep_twosigmag_test80_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test80[i]-.5*test_size*ref_size)/sqrt(var)))
    # 
    cor_temp<-numeric(length=nind)
    m<-1
    for(y in 1:nind){
      temp<-cor(t(test_resid_indep_test20[,m:(m+ncellsper[y]-1)]))
      cor_temp[y]<-mean(temp[upper.tri(temp)],na.rm = T)
      m<-m+ncellsper[y]
    }
    rho_indep_test20_ind[i]<-mean(cor_temp)
    rho<-mean(cor_temp)
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))
    p.val_indep_twosigmag_test20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test20[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_indep_twosigmag_test20_ref20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_indep_test20_ref20[i]-.5*test_size*ref_size)/sqrt(var)))
    
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
    p.val_corr_twosigmag_test50_ref50_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test50_ref50[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_corr_twosigmag_test50_ref20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test50_ref20[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_corr_twosigmag_test50_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test50[i]-.5*test_size*ref_size)/sqrt(var)))
    
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
    p.val_corr_twosigmag_test80_ref50_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test80_ref50[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_corr_twosigmag_test80_ref20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test80_ref20[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_corr_twosigmag_test80_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test80[i]-.5*test_size*ref_size)/sqrt(var)))
    
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
    p.val_corr_twosigmag_test20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test20[i]-.5*test_size*ref_size)/sqrt(var)))
    p.val_corr_twosigmag_test20_ref20_adjust[i]<-2*pnorm(-1*abs((wilcox_stat_corr_test20_ref20[i]-.5*test_size*ref_size)/sqrt(var)))
    
    if(i%%1==0){
      print(i)
      print(paste("Set of 1 took ",round(proc.time()[3]-time[3],2),"seconds for all methods"))
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
  
  assign(paste0("rho_indep_mast_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_mast)
  assign(paste0("rho_indep_test50_mast_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test50_mast)
  assign(paste0("rho_indep_test20_mast_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test20_mast)
  assign(paste0("rho_indep_test80_mast_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test80_mast)
  
  assign(paste0("rho_indep_test20_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test20_camera)
  assign(paste0("rho_indep_test50_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test50_camera)
  assign(paste0("rho_indep_test80_camera_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),rho_indep_test80_camera)
  
  for(t in ls(pattern="p.val")){
    obj<-get(t)
    assign(paste0(t,"_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt),obj)
  }
  
  rm(ref_set_corr_0,ref_set_indep_0,ref_set_indep_50,ref_set_corr_50)
  rm(test_set_indep)
  
  rm(list=setdiff(ls(),c("sim_number_alt","test_size","ref_size",ls(pattern="_sim"))))
  
  
  rm(list=ls(pattern="residuals"))
  rm(list=ls(pattern="times"))
  a<-setdiff(ls(pattern="_sim"),ls(pattern="adj_stat"))
  
  for(i in sim_number_alt){
    save(list=a,file=paste0("p.vals_all_meth_",test_size,"_refsize_",ref_size,"_sim",i,".RData"))
  }
  
  print(gc())
}),file=paste0("output_geneset","_testsize_",test_size,"_refsize_",ref_size,"_sim",sim_number_alt,".txt"))


