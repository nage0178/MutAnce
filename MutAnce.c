#include "MutAnce.h"
gsl_rng * r;  /* global generator */

long opt_locus_count; 
char bpp_errmsg[200] = {0};
int bpp_errno;

void freeMig (gnode_t * node) {

        if(node->mi->me) {
		for (int i = 0; i < node->mi->count; i++) {
			free(node->mi->me[i].source);
			free(node->mi->me[i].target);
		}
		free(node->mi->me);
	}
        free(node->mi);
}

void freeTrees(gtree_t * trees) {
  
  for (int i = 0; i < trees->inner_count + trees->tip_count; i++) {
	gnode_t * node = trees->nodes[i];
  	freeMig(node);

  	for (int i = 0; i < 4; i++)
  	        free(node->pmat[i]);

  	free(node->pmat);

  	free(node->pop);
  	free(node->label);
  	free(node);
  }

  free(trees->nodes);
  
  if (trees->instRate) {
        for (int j = 0; j < 4; j++) 
        	  free(trees->instRate[j]);
        free(trees->instRate);
  }
  
  free(trees);
}

/* Writes the starting seed to a file*/
void writeSeed(unsigned int RGSeed) {
  char filename[120] = "seed.txt\0";

  char seed[15];
  sprintf(seed, "%u", RGSeed);

  FILE *file;
  file = fopen(filename, "w");
  if (! file) {
          fprintf(stderr, "Failed to open %s. Exiting the program\n", filename);
          exit(1);
  }
  fputs(seed, file);
  fclose(file);

  return;
}


void setSeed(long opt_seed, gsl_rng *r) {

  long RGSeed = opt_seed;
  if (RGSeed < 1) {
    RGSeed = time(0);
    writeSeed(RGSeed);
  }
  gsl_rng_set(r, RGSeed);

  return;
}


int main (int argc, char * argv[])
{

  if (argc != 2) {
	  fatal("The control file must given as an argument. No other arguments can be included.");
  }

  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  /* Read in arguments/control file */
  char * opt_gfile = NULL; 
  char * opt_mfile = NULL; 
  char * opt_msafile = NULL; 
  char * opt_subfile = NULL; 
  char * opt_outfile = NULL;
  char * opt_cfile = argv[1];

  long opt_seed = 0;
  long locusNum = -1;

  long opt_site = -1;

  load_cfile(opt_cfile, &opt_msafile, &opt_outfile, &opt_gfile, &opt_mfile, &opt_subfile, &opt_seed, &locusNum, &opt_site);

  /* Change to zero indexing */
  locusNum--;
  opt_site--;

  setSeed(opt_seed, r);
  
  long g_linecount, m_linecount, sub_linecount; 
 
  double * params = NULL;
  double  kappa = 1;
  int jc, site; 
  long samples;
  long msa_count;

  double freqs[4] = {.25, .25, .25, .25};
  double ** instRate = NULL;
  gtree_t * tree;
  char * popOrigin;
  int varsite = 0;

  if (opt_subfile) 
	  jc = 0;
  else 
	  jc = 1;

  g_linecount = getlinecount(opt_gfile);
  m_linecount = getlinecount(opt_mfile);
  sub_linecount = getlinecount(opt_subfile);

  if (g_linecount != m_linecount)
	  fatal("Different number of lines in migration history and gene tree files");
  if (g_linecount != sub_linecount - 1)
	  fatal("Different number of samples in substitution model parameter and gene tree files");

  samples = g_linecount;

  /* Prepare to read in gene trees */
  FILE * fp_gtree = xopen(opt_gfile, "r");

  FILE * fp_mig  = xopen(opt_mfile, "r");

  FILE * fp_sub = xopen(opt_subfile, "r");
  loadSubHeader(fp_sub);

  if (!jc)
	  params = (double *)xmalloc(5 * sizeof(double * ));


  /* Read in sequence data */
  phylip_t * fd = phylip_open(opt_msafile, pll_map_fasta);
  assert(fd);

  msa_t ** msa_list;
  msa_t * msa;
  msa_list = phylip_parse_multisequential(fd, &msa_count);
  assert(msa_list);

  //ANNA???
  //opt_locus_count = 2;
  phylip_close(fd);
  //if (opt_locus_count > msa_count)
  //  fatal("Expected %ld loci but found only %ld", opt_locus_count, msa_count);

  if (locusNum > msa_count)
    fatal("Locus to simulate (%d) is ");
  msa_ambiguous_sites(msa_list[locusNum], pll_map_amb);
  //Allow gaps in alignment?

  // This is the number of sequences at a locus
  // printf("%d\n", msa_list[i]->count);
  // This is the number of loci 
  // printf("%ld \n", msa_count);

  msa = msa_list[locusNum];
  findVariant(msa);
  
  if (opt_site > -1)
    msa->variableCount = 1;
  
  
  if (!msa->variableCount || (opt_site > -1 && msa->variable_sites[opt_site] != 1)) {
  
    free(msa->variable_sites);
    
    for (long i = 0; i < msa_count; i++) {
          msa_destroy(msa_list[i]);
    }

    free(msa_list);
    fclose(fp_gtree);
    fclose(fp_mig);
    fclose(fp_sub);

    if (params)
      free(params);

    if (!msa->variableCount) 
    	printf("No variable sites\n");
    else 
    	printf("Site %ld not variable\n", opt_site + 1);
    
    exit(1);
  }
   
//Why isn't this a sizeof ponter mutHist
  mutHist_t ** mutationHist = xcalloc(samples , sizeof(mutHist_t * ));
  //mutHist_t ** mutationHist = xcalloc(samples , sizeof(mutHist_t));
  
  if (jc) {

	  instRate = xmalloc(4 * sizeof (double * ));

	  for (int i = 0; i < 4; i++) {
		instRate[i] = xmalloc(4 * sizeof(double));

	  	for (int j = 0; j < 4; j++) {

			if (i == j)
				instRate[i][j] = -1;
			else
				instRate[i][j] = 1.0 / 3;
		}
	  }
  }

  int startSite = 0;
  int endSite = msa->length;
  if (opt_site > -1) {
  	startSite = opt_site;
  	endSite = opt_site + 1; 
  } 

  for (int i = 0; i < samples; i++) {

	  fflush(stdout);
	  if (samples < 200 || (i + 1) % (samples / 200) == 0) {
	        printf("\r%4.0f%% ", (i + 1.499) / samples * 100.);	
	        fflush(stdout);
	  } 

	tree = loadTree(fp_gtree, samples);
	loadMigration(fp_mig, tree);
	
        /* Allocate memory for transition probability matrices */
	alloc_pmat(tree->root);

  	/* Read in substitution model parameters */
  	if (!jc) {
  		loadSubstitution(fp_sub, params);
  		kappa = params[0];

  		/* Account for the order in substitution model file */
  		/* Based now in alphabetical order in statFreqs */
		freqs[0] = params[3];
		freqs[1] = params[2];
		freqs[2] = params[4];
		freqs[3] = params[1];

  		update_pmat_hky(tree->root, freqs, kappa);

  	} else
  		update_pmat_hky(tree->root, freqs, 1);
  

  	mutationHist[i] = xcalloc(msa->variableCount, sizeof(mutHist_t));

  	matchGtreeSeqs(tree, msa);
	varsite = 0;
	int numMut = 0;
	
	for (site = startSite; site < endSite; site ++ ) {

		if (msa->variable_sites[site]) {

			likelihoodTips(tree, site, msa);
  			cal_condP(tree->root);


			condP_root(tree->root, freqs, r);

			if (!jc) 
				makeInstantaneousRate(kappa, freqs, tree);
			else
				tree->instRate = instRate;
			
			tree->root->mutTime = -1; 
			tree->root->mut = -1; 

			simulate_mut(tree->root->left, tree->root, tree->instRate, r);
			simulate_mut(tree->root->right, tree->root, tree->instRate, r);

			mutationHist[i][varsite].site = site;
			mutationHist[i][varsite].root = tree->root->state;

			/* Mark which branches have mutations that are have descendants */
			find_mut_descendant(tree->root);
			for (int j = 0; j < tree->tip_count + tree->inner_count; j++) {
				if (tree->nodes[j]->mark) {
					mutationHist[i][varsite].numMut++;	
				}
			}

			assert(mutationHist[i][varsite].numMut != 0);

			mutationHist[i][varsite].mutations = xcalloc(mutationHist[i][varsite].numMut, sizeof(mutEvent_t));
			numMut = 0;
			for (int j = 0; j < tree->tip_count + tree->inner_count; j++) {
				if (tree->nodes[j]->mark) {
					
					mutationHist[i][varsite].mutations[numMut].time = tree->nodes[j]->mutTime;
					mutationHist[i][varsite].mutations[numMut].base = tree->nodes[j]->mut;
					popOrigin = find_pop_origin(tree->nodes[j]);
					mutationHist[i][varsite].mutations[numMut].pop = xstrdup(popOrigin);
					numMut++;
				}
			}

			varsite++;
		}


	}
  	freeTrees(tree);
  }


  summary(mutationHist, varsite, samples, msa, opt_site);
  

  /* Clean up memory */

  /* Free things inside of mutHist */
  for (int i = 0; i < samples; i++) {
  	varsite = 0;
	for (site = startSite; site < endSite; site ++ ) {
		if (msa->variable_sites[site]) {

	  		for (int j = 0; j < mutationHist[i][varsite].numMut; j++) {
				free(mutationHist[i][varsite].mutations[j].pop);
			}
  	        	free(mutationHist[i][varsite].mutations);
			varsite++;

		}
  	}

	free(mutationHist[i]);
  }

  free(msa->variable_sites);
  free(mutationHist);

  for (long i = 0; i < msa_count; i++) {
	msa_destroy(msa_list[i]);
  }

  free(msa_list);
  free(params);

  if (jc) {
	  for (int i = 0; i < 4; i++) 
	  	free(instRate[i]);
	  
	  free(instRate);

  } else {
	free(opt_subfile);
	fclose(fp_sub);
  }

  gsl_rng_free (r);

  free(opt_msafile);
  free(opt_gfile);
  free(opt_mfile);
  free(opt_outfile);

  fclose(fp_gtree);
  fclose(fp_mig);
  freeline();
}
