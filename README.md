# TWO-SIGMA-G Paper Code
 Provides code and instructions to reproduce the analyses in the TWO-SIGMA-G Paper:
 
 
**TWO-SIGMA-G bioRxiv preprint**:

TWO-SIGMA-G: A New Competitive Gene Set Testing Framework for scRNA-seq Data Accounting for Inter-Gene and Cell-Cell Correlation: doi: https://doi.org/10.1101/2021.01.24.427979

## Contact

Please contact Eric Van Buren (evb@hsph.harvard.edu) if you have questions or comments about reproducing these simulations or modifying them for your own purposes.

## Intoduction

All simulations for the TWO-SIGMA-G paper were conducted using a slurm-based computing cluster.  This type of environment would be ideal, although not necessary, to reproduce or modify our simulations. This is because we calibrated the scope of our simulations, in terms of number of cells, number of genes, number of gene sets, diversity of scenarios, and the tradeoff between saving on disk and recomputing information based on what is possible using the distributed computing and mass storage available with a computing cluster. Any user without such a resource that wishes to use our simulation structure for their own purposes can likely do so on a smaller scale by prioritizing the settings (as discussed in Supplementary Section S1) that are most relevant and reducing the scope of the simulations in terms of the number of genes and gene sets being simulated. 

## sim_params_twosigmag.RData

This is a list object which contains all needed objects and parameters values to run simulation code and reproduce results. Included are:

- **alpha**: Initial values for zero-inflation component parameters. Non-intercept terms are changed from zero if **rand_alpha** is TRUE (i.e. if other covariates are assumed in the true model).

- **beta**: Initial values for zero-inflation component parameters. Non-intercept terms are changed from zero if **rand_beta** is TRUE (i.e. if other covariates are assumed in the true model).

- **nind**: Number of individuals (samples or donors) used to simulate data. Set to 100 in all simulations.

- **ncellsper**: Vector of length nind which gives the number of cells per individual.

- **sim.seed** Gives the simulation seeds used to (1) simulate the initial cell population of a given simulation scenario and (2) to simulate the random genes used to construct gene sets. The code was implemented such that there were 100 separate calls to the function simulate_genes_and_summary_statistics.R as a way to distribute computation.

- **id.levels**: Vector assigning numbers (e.g. 1 to nind) to the samples. Used to keep track of simulated cells only.

- **nreps**: How many times the simulation procedure is repeated.

-**phi**: Negative binomial overdispersion parameter used in simulation. Set to 0.1 in all simulations.

-**sigma.a**: Variance component for random intercept terms in the zero-inflation component. If 0, no random effects are simulated.

-**sigma.b**: Variance component for random intercept terms in the mean component. If 0, no random effects are simulated.

-**sim_number**: Arbitrary number used to save all simulation outputs simultaneously and organize results.  Numbers only have meaning in the sense that they relate to the scenarios in Supplementary Section S1 in the following way:

  + Comparisons to MAST, CAMERA, and GSEA:
    - Null Hypothesis, No Random Effect Terms: 12000 <= sim_number <= 12099
    - Null Hypothesis, Random Effect Terms Present: 12200 <= sim_number <= 12299
    - Null Hypothesis, Random Effect Terms Incorrectly Absent: 12700 <= sim_number <= 12799
    - Alternative Hypothesis, No Random Effect Terms: 12300 <= sim_number <= 12399
    - Alternative Hypothesis, Random Effect Terms Present: 12500 <= sim_number <= 12599
    - Alternative Hypothesis, Random Effect Terms Incorrectly Absent: 12800 <= sim_number <= 12899
  + Comparisons to fGSEA, iDEA, and PAGE:
    - Null Hypothesis, No Random Effect Terms: 12000 <= sim_number <= 12099
    - Null Hypothesis, Random Effect Terms Present: 12200 <= sim_number <= 12299
    - Null Hypothesis, Random Effect Terms Incorrectly Absent: 12700 <= sim_number <= 12799
    - Alternative Hypothesis, No Random Effect Terms: 14300 <= sim_number <= 14399
    - Alternative Hypothesis, No Random Effect Terms: 14500 <= sim_number <= 14599
    - Alternative Hypothesis, No Random Effect Terms: 14800 <= sim_number <= 14899

-**nruns**: Gives the number of correlated genes to simulate.  Set to 29 in all simulations to correspond to a gene set size of 30.

-**mean_form**: Gives the form of the mean component model, including whether or not additional covariates or random effect terms were included.

-**zi_form**: Gives the form of the zero-inflation component model, including whether or not additional covariates were included.

-**alpha_lower**: Lower bound for simulation of Uniform distributed covariates for the zero-inflation component. Equal to -1.5 for all simulations.

-**alpha_upper**: Upper bound for simulation of Uniform distributed covariates for the zero-inflation component. Equal to 1.5 for all simulations.

-**beta_lower**: Lower bound for simulation of Uniform distributed covariates for the mean component. Equal to -1.5 for all simulations.

-**beta_upper**: Upper bound for simulation of Uniform distributed covariates for the mean component. Equal to 1.5 for all simulations.

-**rand_alpha**: Logical vector indicating whether covariates should be randomly sampled for the zero-inflation component. If FALSE, corresponds to scenarios without additional covariates. If TRUE, corresponds to situations with additional covariates.

-**rand_beta**: Logical vector indicating whether covariates should be randomly sampled for the mean component. If FALSE, corresponds to scenarios without additional covariates. If TRUE, corresponds to situations with additional covariates.

## (1) simulate_genes_and_summary_statisics.R

This code accomplishes the following:

1. It simulates and saves correlated gene expression counts

2. It calculates and saves gene-level summary statistics using 
twosigma. It is desirable to do so here to avoid saving
more to disk than is necessary and to avoid repeating computation.

3. It saves gene-level likelihood ratio statistics and
residuals from the gene-level models used to collect summary statistics. 
These will be used later to conduct set-level inference, and are saved here
to avoid repeating computation later.

## (2) set_level_simulation_code.R and set_level_summary_statistic_simulation_code.R

These two scripts accomplish the following, for the simulations
using the raw data and summary statistics, respectively:

 1. Constructs test and reference sets for each scenario
 by randomly sampling genes from matched null and alternative
 hypotheses (using the arbitrary "sim_number" to ensure a match)
 2. Computes and saves the IGC estimate using the individual-level
 residual correlation matrix for twosigmag
 3. Runs set-level inference for (CAMERA, MAST, GSEA, and twosigmag) or 
 (iDEA, fGSEA,PAGE and twosigmag).  Including all methods in the same script is 
 not ideal for distributed computing but serves as a way of being 
 especially sure that the same test and reference sets are being used for all methods
 4. Saves output
 
## (3) boxplot_functions.R

This script provides the functions used to plot the simulation results
as boxplots. Included are functions needed for boxplots of both type-I error
and power. 

 
