seed = 1

seqfile =  simulate_IM.txt
Imapfile = simple.Imap.txt
jobname = out

# fixed species tree
species&tree = 3 A B C
                 5 5 5 
		  ((A, B), C);


# phased data for population
phase =   0 0 0

# use sequence likelihood
usedata = 1

nloci = 5 
clock = 1
model = HKY

# invgamma(a, b) for root tau & Dirichlet(a) for other tau's
tauprior = gamma 50 100000
thetaprior = gamma 50 100000

# MCMC samples, locusrate, heredityscalars, Genetrees
print = 1 0 0 1 1
burnin = 10000
sampfreq = 8
nsample = 50000
printlocus = 2 1 3

wprior = 15 .01
migration = 2
A B 
B A 

