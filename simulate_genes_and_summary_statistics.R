#simulate_genes_and_summary_statistics.R

# This script is used to accomplish several gene-level tasks
# related to the TWO-SIGMA-G paper simulations

# 1. It simulates and saves correlated gene expression counts

# 2. It calculates and saves gene-level summary statistics using 
# twosigma. It is desirable to do so here to avoid saving
# more to disk than is necessary and to avoid repeating computation.

# 3. It saves gene-level likelihood ratio statistics and
# residuals from the gene-level models used to collect summary statistics. 
# These will be used later to conduct set-level inference, and are saved here
# to avoid repeating computation later.


# Arbitrary number used only to distinguish simulation scenarios
sim_num<-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

suppressPackageStartupMessages({
  library(devtools)
  library(R.utils)
  library(pscl)
  library(twosigma)
  library(methods) #Include to prevent conflicts with glmmTMB package
  library(pryr)
  library(MASS)
  library(limma)})

# This file is a list object which contains all objects needed 
#to run each of the simulation scenarios we considered 
load("sim_params_twosigmag.RData")
print(sim_num)
# assign the objects for scenario being computed to the Global Environment
list2env(sim_params_twosigmag[[sim_num]],envir = parent.frame())

assign(paste0("sim.seed_sim",sim_number),sim.seed)
assign(paste0("ncellsper_sim",sim_number),ncellsper)
assign(paste0("alpha_sim",sim_number),alpha)
assign(paste0("beta_sim",sim_number),beta)
assign(paste0("phi_sim",sim_number),phi)
assign(paste0("sigma.a_sim",sim_number),sigma.a)
assign(paste0("sigma.b_sim",sim_number),sigma.b)

# Maximum number of runs allowed
# equal to minimum but framework allowed for flexibility
assign(paste0("nmax_sim",sim_number),nreps)
assign(paste0("nreps_sim",sim_number),nreps)
assign(paste0("nruns_sim",sim_number),nruns)
assign(paste0("nind_sim",sim_number),nind)

assign(paste0("mean_form_sim",sim_number),mean_form)
assign(paste0("zi_form_sim",sim_number),zi_form)

assign(paste0("w_2",sim_number),alpha_lower)

assign(paste0("beta_lower_sim",sim_number),beta_lower)
assign(paste0("beta_upper_sim",sim_number),beta_upper)

assign(paste0("rand_alpha_sim",sim_number),rand_alpha)
assign(paste0("rand_beta_sim",sim_number),rand_beta)

assign(paste0("w_2_sim",sim_number),w_2)

# All related simulations will simulate individual-level covariates
# using the same simulation seed to ensure the cell population 
# is the same. This is needed because, to mimic a real dataset,
# one the pool of cells is fixed, variability comes from
# simulating many genes that are all observed on the fixed 
# pool of cells.
set.seed(sim.seed[1])

id.levels<-1:nind
t2d_ind<-rbinom(nind,1,p=.4)
t2d_sim<-rep(t2d_ind,times=ncellsper)
nind<-length(id.levels)
id<-rep(id.levels,times=ncellsper)
cdr_sim<-rbeta(sum(ncellsper),3,6)
age_sim_ind<-sample(c(20:60),size=nind,replace = TRUE)
age_sim<-rep(age_sim_ind,times=ncellsper)

Z<-cbind(1,scale(t2d_sim),scale(age_sim),scale(cdr_sim))
colnames(Z)<-c("Intercept","t2d_sim","age_sim","cdr_sim")
X<-cbind(1,scale(t2d_sim),scale(age_sim),scale(cdr_sim))
colnames(X)<-c("Intercept","t2d_sim","age_sim","cdr_sim")

#simulate random effects
a <- as.matrix(rnorm(nind,mean=0,sd=sigma.a))
a.rep <- rep(a,times=ncellsper)
b <- as.matrix(rnorm(nind,mean=0,sd=sigma.b))
b.rep <- rep(b,times=ncellsper)

# Simulate the counts according to our model
sim_ZINB_RE_data<-function(ncellsper,alpha,beta,nind
                           ,phi,sigma.a,sigma.b,id.levels
                           ,X,Z,a.rep,b.rep)
{
  # Simulate t2d status randomly, so no relationship between t2d and later simulated counts
  logit.p<-Z%*%alpha +a.rep
  log.mu<-X%*%beta+b.rep
  
  p  <- exp(logit.p)/(1+exp(logit.p)) # inverse logit function
  mu <- exp(log.mu)
  Y<-rep(NA,sum(ncellsper))
  ind.dropout <- rbinom(length(Y), 1, p)
  for (i in 1:length(Y)){
    if(ind.dropout[i] == 1){Y[i]=0}
    if(ind.dropout[i] == 0)
    {
      if(phi==0) # not really going to use
      {
        Y[i]<-rpois(1,lambda=mu[i])
      }
      if(phi>0)
      {
        Y[i]<-rnbinom(1,size=(1/phi),p=(1/(1+phi*mu[i])))
      }
      #Watch out for negative binomial parameterizations
    }
  }
  return(Y)
}

t2d_sim<-scale(t2d_sim)
age_sim<-scale(age_sim)
cdr_sim<-scale(cdr_sim)

# Vectors of gene-level statistics from twosigma
# for the original gene "orig"
# and the genes simulated with nose "large_noise"
LR_stat_orig_null_all<-rep(NA,length=nreps)
chisq_orig_null_all<-rep(NA,length=nreps)
LR_stat_large_noise_null_all<-vector('list',length=1)
chisq_large_noise_null_all<-vector('list',length=1)

#Will save gene-level summary statistics
sum_stats<-vector('list',length=nreps)

# Indices of whether the model converged successfully
# Needed especially for cases with random effects
index_large_noise_null_all<-vector('list',length=1)
index_orig<-rep(FALSE,length=nreps)

#Matrices to save simulated genes, and model residuals from twosigma
dat_orig_null_all<-matrix(nrow=nreps,ncol=sum(ncellsper))
residuals_orig_null_all<-matrix(NA,nrow=nreps,ncol=sum(ncellsper))
dat_large_noise_null_all<-vector('list',length=1)
residuals_large_noise_null_all<-vector('list',length=1)

# Will hold the randomly simulated regression parameters
alpha_matrix<-matrix(nrow=nreps,ncol=length(alpha))
beta_matrix<-matrix(nrow=nreps,ncol=length(beta))

#Will save gene-level summary statistics
sum_stats<-vector('list',length=nreps)

# To facilitate distributed computing, we split settings into 100
# jobs to run simultaneously, which are accessed later using
# numbers 1:100

code_num=1
set.seed(sim.seed[code_num])

# if random effects are present in the model,
# define objects to save thei magnitudes
if(sigma.b>0){
  re_orig_null_all<-rep(NA,length=nreps)
  re_large_noise_null_all<-vector('list',length=1)
  
  re_zi_orig_null_all<-rep(NA,length=nreps)
  re_zi_large_noise_null_all<-vector('list',length=1)
}

count<-NULL
mean_form<-update.formula(mean_form,as.formula(count~.))


capture.output(system.time(
  {print(paste("sim_number_pathway = ",sim_number,sep=""))
    # nreps gives the number of gene sets to simulate
    for(j in 1:nreps){
      starttime<-proc.time()
      sum_stats[[j]]<-matrix(NA,ncol=2,nrow=nruns+1)
      
      # assign coefficient parameter values
      for(z in 1:length(alpha)){
        if(rand_alpha[z]==T){
          alpha[z]<-runif(1,alpha_lower,alpha_upper)
        }
      }
      for(z in 1:length(beta)){
        if(rand_beta[z]==T){
          beta[z]<-runif(1,beta_lower,beta_upper)
        }
      }
      
      alpha_matrix[j,]<-alpha
      beta_matrix[j,]<-beta
      #master "starting" gene
      original<-sim_ZINB_RE_data(ncellsper,alpha,beta,nind,phi,sigma.a,sigma.b,id.levels,X,Z,a.rep,b.rep)
      dat_orig_null_all[j,]<-original

      sim_data_permute<-function(Y_input=NULL)
      {
    
        rand<-runif(n=1,min=w_2,max=w_2)
        Y_out<-round(.5*Y_input+rand*rnbinom(n=length(Y_input),size=(1/2.5),p=(1/(1+2.5*3))),0)
        if(alpha[2]>0 | beta[2]>0){
          if(beta[2]==.03){ # Alternative Hypothesis
            Y_out[t2d_sim>0]<-round(Y_out[t2d_sim>0]+.1*rnbinom(n=length(Y_out[t2d_sim>0]),size=(1/2.5),p=(1/(1+2.5*3))),0)
          }
          if(beta[2]==.06){ # Mixed DE scenarios involving more signal
            Y_out[t2d_sim>0]<-round(Y_out[t2d_sim>0]+.15*rnbinom(n=length(Y_out[t2d_sim>0]),size=(1/2.5),p=(1/(1+2.5*3))),0)
          }
          if(beta[2]>=.1){ # Summary statistic scenarios
            mult<-runif(1,.15,.75)
            Y_out[t2d_sim>0]<-round(Y_out[t2d_sim>0]+rand*mult*rnbinom(n=length(Y_out[t2d_sim>0]),size=(1/2.5),p=(1/(1+2.5*3))),0)
          }
        }
        Y_out[Y_out<0]<-0
        return(list(Y_input=Y_input,Y_out=Y_out))
      }
      
      Y_large_null<-matrix(ncol=sum(ncellsper),nrow=nruns)
      dat_large_noise_null_all[[j]]<-Y_large_null
      
      for(i in 1:nruns){
        # Simulate correlated genes
        Y_large_null[i,]<-sim_data_permute(original,3=3,a_perm=a1_perm,2.5=2.5)$Y_out 
        
        # Below ensures the proportion of zeros stays the same
        prop<-mean(original==0)-mean(Y_large_null[i,]==0)
        index<-sample(which(Y_large_null[i,]!=0),size=floor(prop*length(Y_large_null[i,])))
        Y_large_null[i,index]<-0
      }
      
      dat_large_noise_null_all[[j]]<-Y_large_null
      cor(Y_large_null[1,],original)
      
      LR_stat_orig_null<-NA
      chisq_orig_null<-NA
      
      # Fit twosigma model to master gene to (1) get summary statistics for relevant simulations
      # , (2) collect residuals, and (3) collect
      # gene-level LR statistics
      tryCatch({
        fit_orig_no_t2d_sim_null<-twosigma_custom(count=matrix(original,nrow=1),mean_form=mean_form,zi_form=zi_form,id=id,return_summary_fits = FALSE)
        
        fit_orig_t2d_sim_null<-twosigma_custom(count=matrix(original,nrow=1),mean_form=update(mean_form,~.+t2d_sim)
                                               ,zi_form=update(zi_form,~.+t2d_sim),id=id,return_summary_fits = FALSE)
        
        sum_stats[[j]][1,1]<-summary(fit_orig_t2d_sim_null[[1]])$coefficients$cond['t2d_sim',1]
        sum_stats[[j]][1,2]<-(summary(fit_orig_t2d_sim_null[[1]])$coefficients$cond['t2d_sim',2])^2
        
        residuals_orig_null_all[j,]<-residuals(fit_orig_t2d_sim_null[[1]])
        
        fit<-summary(fit_orig_t2d_sim_null[[1]])
        
        tryCatch({ # if error in model fit or accessing
          index_orig[j]<-(fit_orig_t2d_sim_null[[1]]$sdr$pdHess && fit_orig_t2d_sim_null[[1]]$sdr$pdHess
                             && fit_orig_no_t2d_sim_null[[1]]$fit$objective<=fit_orig_no_t2d_sim_null[[1]]$fit$objective)}
          ,error=function(e){print(paste("Error Indexing orig at",j))})
        
        
        chisq_orig_null<-(fit$coefficients$cond['t2d_sim','z value'])^2 + (fit$coefficients$zi['t2d_sim','z value'])^2
        
        # If random effects save relevant information
        if(sigma.b>0){
          temp<-strsplit(capture.output(fit_orig_t2d_sim_null[[1]]$sdr),"\t")
          temp<-temp[[length(temp)-2]]
          temp<-strsplit(temp," ")
          temp<-lapply(temp,function(x){x[!x ==""]})
          re_orig<-exp(as.numeric(temp[[1]][2]))
          
          temp<-strsplit(capture.output(fit_orig_t2d_sim_null[[1]]$sdr),"\t")
          temp<-temp[[length(temp)-1]]
          temp<-strsplit(temp," ")
          temp<-lapply(temp,function(x){x[!x ==""]})
          re_zi_orig<-exp(as.numeric(temp[[1]][2]))
        }
        
        LR_stat_orig_null<--2*(summary(fit_orig_no_t2d_sim_null[[1]])$logLik-summary(fit_orig_t2d_sim_null[[1]])$logLik)
      }, error=function(e){print(paste("Orig error at rep",j))})
      
      
      # Fit twosigma model to correlated genes to (1) get summary statistics for relevant simulations
      # , (2) collect residuals, and (3) collect
      # gene-level LR statistics
      
      fit_no_t2dsim_large_noise_null<-vector('list',length=nruns)
      fit_t2dsim_large_noise_null<-vector('list',length=nruns)
      
      LR_stat_large_noise_null<-rep(NA,length=nruns)
      if(sigma.b>0){
        re_large_noise<-rep(NA,length=nruns)
        re_zi_large_noise<-rep(NA,length=nruns)
      }
      chisq_large_noise_null<-rep(NA,length=nruns)
      
      index_large_noise_null<-rep(FALSE,length=nruns)
      residuals_large_noise_null_all[[j]]<-matrix(NA,nrow=nruns,ncol=sum(ncellsper))
      
      for(i in 1:nruns)
      {
        tryCatch({
          
          fit_no_t2dsim_large_noise_null[[i]]<-twosigma_custom(count=matrix(Y_large_null[i,],nrow=1),mean_form=update(mean_form,count~.),zi_form=zi_form,id=id,return_summary_fits = FALSE)
          
          fit_t2dsim_large_noise_null[[i]]<-twosigma_custom(count=matrix(Y_large_null[i,],nrow=1),mean_form=update(update(mean_form,~.+t2d_sim),count~.),zi_form=update(zi_form,~.+t2d_sim),id=id,return_summary_fits = FALSE)
          
          sum_stats[[j]][i+1,1]<-summary(fit_t2dsim_large_noise_null[[i]][[1]])$coefficients$cond['t2d_sim',1]
          sum_stats[[j]][i+1,2]<-(summary(fit_t2dsim_large_noise_null[[i]][[1]])$coefficients$cond['t2d_sim',2])^2
          
          residuals_large_noise_null_all[[j]][i,]<-residuals(fit_t2dsim_large_noise_null[[i]][[1]])
          
          LR_stat_large_noise_null[i]<--2*(summary(fit_no_t2dsim_large_noise_null[[i]][[1]])$logLik-summary(fit_t2dsim_large_noise_null[[i]][[1]])$logLik)
          
          if(sigma.b>0){
            temp<-strsplit(capture.output(fit_t2dsim_large_noise_null[[i]][[1]]$sdr),"\t")
            temp<-temp[[length(temp)-2]]
            temp<-strsplit(temp," ")
            temp<-lapply(temp,function(x){x[!x ==""]})
            re_large_noise[i]<-exp(as.numeric(temp[[1]][2]))
            
            temp<-strsplit(capture.output(fit_t2dsim_large_noise_null[[i]][[1]]$sdr),"\t")
            temp<-temp[[length(temp)-1]]
            temp<-strsplit(temp," ")
            temp<-lapply(temp,function(x){x[!x ==""]})
            re_zi_large_noise[i]<-exp(as.numeric(temp[[1]][2]))
          }
          
          fit<-summary(fit_t2dsim_large_noise_null[[i]][[1]])
          chisq_large_noise_null[i]<-(fit$coefficients$cond['t2d_sim','z value'])^2 + (fit$coefficients$zi['t2d_sim','z value'])^2
          print(i)
          
          tryCatch({
            index_large_noise_null[i]<-(fit_t2dsim_large_noise_null[[i]][[1]]$sdr$pdHess && fit_no_t2dsim_large_noise_null[[i]][[1]]$sdr$pdHess
                                           && fit_t2dsim_large_noise_null[[i]][[1]]$fit$objective<=fit_no_t2dsim_large_noise_null[[i]][[1]]$fit$objective)}
            ,error=function(e){print(paste("Error Indexing Corr at",j))})
        }, error=function(e){print(paste("Corr error at rep",j, " run ",i))})
      }
      
      #Save output in a more convenient way for later accessing
      
      LR_stat_orig_null_all[j]<-LR_stat_orig_null
      chisq_orig_null_all[j]<-chisq_orig_null
      
      LR_stat_large_noise_null[!index_large_noise_null]<-NA
      chisq_large_noise_null[!index_large_noise_null]<-NA
      
      LR_stat_large_noise_null_all[[j]]<-LR_stat_large_noise_null
      
      index_large_noise_null_all[[j]]<-index_large_noise_null
      chisq_large_noise_null_all[[j]]<-chisq_large_noise_null
      
      if(sigma.b>0){
        re_orig_null_all[j]<-re_orig
        re_large_noise_null_all[[j]]<-re_large_noise
        
        re_zi_orig_null_all[j]<-re_zi_orig
        re_zi_large_noise_null_all[[j]]<-re_zi_large_noise
      }
      
      print(paste0("Replication" ,j))
      print(gc())
      times[j]<-proc.time()[3]-starttime[3]
      print(paste("Rep took ",times[j],"seconds"))
    }
    
    LR_stat_orig_null_all[!index_orig]<-NA
    chisq_orig_null_all[!index_orig]<-NA
    
    
    assign(paste0("chisq_orig_null_all_sim",sim_number,"_",code_num),c(chisq_orig_null_all))
    assign(paste0("LR_stat_orig_null_all_sim",sim_number,"_",code_num),c(LR_stat_orig_null_all))
    
    if(sigma.b>0){
      assign(paste0("re_orig_null_all_sim",sim_number,"_",code_num),c(re_orig_null_all))
      assign(paste0("re_zi_orig_null_all_sim",sim_number,"_",code_num),c(re_zi_orig_null_all))
    }
    assign(paste0("index_orig_sim",sim_number,"_",code_num),index_orig)
    
    assign(paste0("LR_stat_large_noise_null_all_sim",sim_number,"_",code_num),c(LR_stat_large_noise_null_all))
    
    if(sigma.b>0){
      assign(paste0("re_large_noise_null_all_sim",sim_number,"_",code_num),c(re_large_noise_null_all))
      assign(paste0("re_zi_large_noise_null_all_sim",sim_number,"_",code_num),c(re_zi_large_noise_null_all))
    }
    
    # Save objects with better names for joint accessing and save results
    
    assign(paste0("chisq_large_noise_null_all_sim",sim_number,"_",code_num),c(chisq_large_noise_null_all))
    assign(paste0("index_large_noise_null_all_sim",sim_number,"_",code_num),c(index_large_noise_null_all))
    
    
    assign(paste0("residuals_orig_null_all_sim",sim_number,"_",code_num),residuals_orig_null_all)
    assign(paste0("residuals_large_noise_null_all_sim",sim_number,"_",code_num),residuals_large_noise_null_all)
    assign(paste0("dat_large_noise_null_all_sim",sim_number,"_",code_num),dat_large_noise_null_all)
    assign(paste0("sum_stats_sim",sim_number,"_",code_num),sum_stats)
    assign(paste0("dat_orig_null_all_sim",sim_number,"_",code_num),dat_orig_null_all)
    
    assign(paste0("alpha_matrix_sim",sim_number,"_",code_num),alpha_matrix)
    assign(paste0("beta_matrix_sim",sim_number,"_",code_num),beta_matrix)
    
    assign(paste0("times_sim",sim_number,"_",code_num),times)
    if(sigma.b>0){
      rm(list=setdiff(ls()
                      ,c("sim_number",paste0("LR_stat_orig_null_all_sim",sim_number,"_",code_num)
                         ,paste0("LR_stat_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("chisq_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("chisq_orig_null_all_sim",sim_number,"_",code_num)
                         ,paste0("index_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("dat_orig_null_all_sim",sim_number,"_",code_num)
                         ,paste0("sum_stats_sim",sim_number,"_",code_num)
                         ,paste0("residuals_orig_null_all_sim",sim_number,"_",code_num)
                         ,paste0("residuals_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("index_orig_sim",sim_number,"_",code_num)
                         ,paste0("alpha_matrix_sim",sim_number,"_",code_num)
                         ,paste0("beta_matrix_sim",sim_number,"_",code_num)
                         ,paste0("re_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("re_zi_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("re_orig_null_all_sim",sim_number,"_",code_num)
                         ,paste0("re_zi_orig_null_all_sim",sim_number,"_",code_num)
                         ,paste0("dat_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("times_sim",sim_number,"_",code_num))))
    }else{
      rm(list=setdiff(ls()
                      ,c("sim_number",paste0("LR_stat_orig_null_all_sim",sim_number,"_",code_num)
                         ,paste0("LR_stat_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("chisq_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("chisq_orig_null_all_sim",sim_number,"_",code_num)
                         ,paste0("index_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("dat_orig_null_all_sim",sim_number,"_",code_num)
                         ,paste0("sum_stats_sim",sim_number,"_",code_num)
                         ,paste0("residuals_orig_null_all_sim",sim_number,"_",code_num)
                         ,paste0("residuals_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("index_orig_sim",sim_number,"_",code_num)
                         ,paste0("alpha_matrix_sim",sim_number,"_",code_num)
                         ,paste0("beta_matrix_sim",sim_number,"_",code_num)
                         ,paste0("dat_large_noise_null_all_sim",sim_number,"_",code_num)
                         ,paste0("times_sim",sim_number,"_",code_num))))
    }
    
    save.image(file=paste0("results_sim",sim_number,"_pathway.sim.RData"))
    print(gc())
  }),file=paste0("output_sim",sim_number,"_path.txt"))









