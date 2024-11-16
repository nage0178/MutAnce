seed = 1

seqfile = simulate_IM.txt
treefile = simulate_trees.txt
Imapfile = simple.Imap.txt

modelparafile = parameters.txt

# fixed number of species/populations 
*speciesdelimitation = 0

# fixed species tree

species&tree = 3  A B C
                  5 5 5  
		  ((A #0.00049, B #0.00044):.00025 #.0005, C #0.0005):.0005 #.0006;


# phased data for population
phase =   0 0 0 

loci&length = 5 1000
clock = 1
locusrate =0
printlocus = 1 1
model = 7
basefreqs = 1 .27 .23 .28 .22
qrates = 1 8 1 1 1 1 8


migration = 2
A B  0.15
B A  0.2

