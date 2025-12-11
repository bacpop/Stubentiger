# model definition derived from Corander et al. (2017)

## Definition of the time-step and output as "time"
initial(time) <- 0
update(time) <- (step + 1) * dt

# calculate probabilities based on relative sizes of groups and their fitness (does not necessarily sum to 1):
Pop_size <- sum(Pop[1:GPSC_no,1:sero_no])

# frequency dependent selection (i.e. the fitness of the genotypes are not constant but based on how frequent the genes/loci they contain are in the population)
# based on Corander et al. (2017)
# frequency of each gene at current time:
gene_freq[,] <-  Genotypes[i,j] * sum(Pop[j,1:sero_no])
freq[] <- sum(gene_freq[i,1:GPSC_no]) / Pop_size

# identify genes that are under NFDS
pi_f_freq[,] <- if ((delta[i] <= prop_f * gene_no)) Genotypes[i,j] * (eq[i] - freq[i]) else 0
pi_f_genotypes[] <- sum(pi_f_freq[1:gene_no,i])

# Genotype specific probability to produce offspring
# those are the individuals' probabilities multiplied by the number of individual that have this genotype
probs[,] <- (1 + exp(sigma_f))^pi_f_genotypes[i] * Pop[i,j] * (1- (as.integer(time >= vacc_time) * vaccTypes[j] * v))

# generate the next generation based on the current one
y[,] <- rpois(capacity * (probs[i,j] / sum(probs[1:GPSC_no,1:sero_no])) * (1-exp(m)) )

# m is the migration rate
# fitness of individuals in the community is reduced by this rate
# determining immigration number:
mig_num <- rbinom(capacity, exp(m))
Pop_mig[,] <- rbinom(mig_num, Pop_mig_dist[i,j])

## Core equation for population (assume constant size here):
update(Pop[,]) <- y[i,j] + Pop_mig[i,j]
update(Pop_tot[]) <- sum(y[i,]) + sum(Pop_mig[i,])

initial(Pop[,]) <- Pop_start[i,j] # deterministic, user-based start value
initial(Pop_tot[]) <- sum(Pop_start[i,])

#calculate equilibrium frequencies of genes
gene_eq[,] <- Genotypes[i,j] * Pop_eq[j]
eq[] <- sum(gene_eq[i,1:GPSC_no]) / sum(Pop_eq[1:GPSC_no])

## User defined parameters - default in parentheses:
dt <- user(1)
GPSC_no <- user() #number of species / strains / genotypes in the population
sero_no <- user() # number of serotypes
gene_no <- user() # number of genes in the data set
Pop_start[,] <- user() # initial frequency of Genotypes
Pop_eq[] <- user()
capacity <- user()
sigma_f <- user()
#sigma_w <- user()
prop_f <- user()
delta[] <- user()
m <- user() # migration rate
#GeneFitness[] <- user() # fitness vector for different genes
Genotypes[,] <- user() # each column is a genotype, giving the information which genes are present in that genotype and which are not
vaccTypes[] <- user() # Boolean vector (0/1) whether genotype is affected by vaccine (1) or not (0)
v <- user() # effect size of vaccine on vaccine genotypes
vacc_time <- user() # time when the vaccination happens / starts to have an effect
Pop_mig_dist[,] <- user()

dim(gene_freq) <- c(gene_no, GPSC_no)
dim(freq) <- gene_no
dim(Pop_start) <- c(GPSC_no,sero_no)
dim(Pop_eq) <- GPSC_no
dim(Pop_mig) <- c(GPSC_no,sero_no)
dim(Pop_tot) <- GPSC_no
dim(gene_eq) <- c(gene_no, GPSC_no) #frequency of genes at equilibrium
dim(eq) <- gene_no
dim(pi_f_freq) <- c(gene_no, GPSC_no)
dim(pi_f_genotypes) <- GPSC_no
dim(delta) <- gene_no
dim(Genotypes) <- c(gene_no, GPSC_no) # we have in each column the genes (present/not present, i.e. 1/0) of one genotype
dim(Pop) <- c(GPSC_no,sero_no)
dim(y) <- c(GPSC_no,sero_no)
dim(probs) <- c(GPSC_no,sero_no)
dim(vaccTypes) <- sero_no
dim(Pop_mig_dist) <- c(GPSC_no,sero_no)
