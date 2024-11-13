#include "MutAnce.h"
gsl_rng * r;  /* global generator */


int simulateDNABranch(int base, double instRate[4][4], double branchLength, const gsl_rng * r) {
  double time = 0;        /* The time that has passed in the simulation */
  double totRate, unif;   /* Total rate of substitions, depends on starting state.
                          unif is used with the random number generator to pick the next base */
  double probBase[4];     /* Given an substition occured, the probability the new base is 0-3 (A-T)*/

  /* Until the end of the branch is reached */
  while (time <= branchLength) {
    /* Rate depends on starting state */
    totRate = -instRate[base][base];

    /* Time of next substition*/
    time = time + gsl_ran_exponential (r, 1 / totRate);

    if (time > branchLength) {
      return base;
    }

    /* Finds the probability of each substition */
    for (int i = 0; i < 4; i++) {
      if (i == base) {
        probBase[i] = 0.0;
      } else {
        probBase[i] = instRate[base][i] / totRate;
      }
    }
    //assert(FloatEquals(probBase[0] + probBase[1] + probBase[2] + probBase[3], 1, 1e-4));

    /* Updates the current base */
    unif = gsl_rng_uniform (r);
    if (unif <= probBase[0]) {
      base = 0;
    } else if (unif <= probBase[0] + probBase[1]) {
      base = 1;
    } else if (unif <= probBase[0] + probBase[1] + probBase[2]) {
      base = 2;
    } else {
      base = 3;
    }
  }
  return base;
}

/* Makes the instantaneous rate matrix given the stationary frequencies, substition rate, and rate parameters (in instRate)*/
void makeInstantaneousRate(double instRate[4][4], double statFreq[]) {

  double diag,  weightSum  = 0;
  /* Finds the diagonal elements and the total substition rate*/
  for (int i = 0; i < 4; i++) {
    instRate[i][i] = 0;
    diag = 0;
    for (int j = 0; j < 4; j ++) {
      diag = diag + instRate[i][j];
    }
    weightSum = diag * statFreq[i] + weightSum;
    instRate[i][i] = -diag;

  }
  /*if (weightSum != 1) {
    printf("Rescaling instantaneous rate matrix so the average substitution rate is one.\nThen multiplying by the substituion rate. \n");
  } */
  if (weightSum == 0) {
    fprintf(stderr, "The average substition rate is zero. Exiting. \n");
    exit(1);
  }

  /* Rescales the instantaneous rate matrix so that the average substitution rate is equal to mu */
  /*for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j ++) {
      instRate[i][j] = instRate[i][j] / weightSum *  mu;
    }
  }*/

  return;
}

void simulate_site(gnode_t * node, int p_state, double instRate[4][4], const gsl_rng * r) {
	node->state = simulateDNABranch(p_state, instRate, node->length, r);
	if (node->left) {
		simulate_site(node->left, node->state, instRate, r);
		simulate_site(node->right, node->state, instRate, r);
	}
}

void root_site (gnode_t * node, double *freqs, const gsl_rng * r) {
	double draw = gsl_ran_flat(r, 0, 1);

	int i; 
	double sum = 0;
	for (i = 0; i < 4; i++ ) {
		if (draw < sum + freqs[i])
			break;
		sum += freqs[i];
	}
	assert(i <4);
	node->state = i;
}

void main () {

	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	gnode_t ** nodes;
	gnode_t * root = calloc(1, sizeof(gnode_t));
	root->length = 0 ;

	root->left = calloc(1, sizeof(gnode_t));
	root->left->length = .08;

	root->right = calloc(1, sizeof(gnode_t));
	root->right->length = .15;

	root->left->left = calloc(1, sizeof(gnode_t));
	root->left->left->length = .1;

	root->left->right = calloc(1, sizeof(gnode_t));
	root->left->right->length = .07;

	/* change here */
	root->left->left->left = calloc(1, sizeof(gnode_t));
	root->left->left->left->length = .09;

	root->left->left->right = calloc(1, sizeof(gnode_t));
	root->left->left->right->length = .12;

	nodes = calloc(7, sizeof(gnode_t *));
	nodes[0] = root->left->left->left;
	nodes[1] = root->left->left->right;
	nodes[2] = root->left->right;
	nodes[3] = root->right;
	nodes[4] = root->left->left;
	nodes[5] = root->left;
	nodes[6] = root;

	double statFreq[4] = {.2, .3, .23, .27};
	double kappa = 2;
	/* Set the values of the instantaneous rate matrix. Sets diagonal and rescales to makes the mutation rate in the function. */
	double rateParams[6] = {1, 1, 1, 1, 1, 1};
	rateParams[1] = rateParams[4] = kappa;
	double instRate[4][4] = { {                          0, rateParams[0] * statFreq[1],    rateParams[1] * statFreq[2],  rateParams[2] * statFreq[3]},
				  {rateParams[0] * statFreq[0],                           0,    rateParams[3] * statFreq[2],  rateParams[4] * statFreq[3]},
				  {rateParams[1] * statFreq[0], rateParams[3] * statFreq[1],                              0,  rateParams[5] * statFreq[3]},
				  {rateParams[2] * statFreq[0], rateParams[4] * statFreq[1],    rateParams[5] * statFreq[2],                           0}};
	
	makeInstantaneousRate(instRate, statFreq);
	//int ** sitePatterns[7][10000];
	int pattern[4] = {0,2,0,0};

	// For each tip
	for (int i = 0; i < 4; i++) {
		nodes[i]->state =  pattern[i];
		//set likelihoodTips
		for (int j = 0; j < 4; j++) {
			if (nodes[i]->state == j)  {
				nodes[i]->likelihood[j] = 1;
				}
			else
				nodes[i]->likelihood[j] = 0;

		}

	}

	alloc_pmat(root);
	update_pmat_hky(root, statFreq, kappa);
	cal_condP(root);
	condP_root(root, statFreq, r);

	/* conditional probability */
	printf("cond prob ");
	for (int i = 0; i < 4; i++) {
		printf("%f ", root->condP[i]);
	}
	printf("\n");
	
	printf("cond prob L, root A ");
	prob_recursive(root->left, 0, r);
	for (int i = 0; i < 4; i++) {
		printf("%.8f ", root->left->condP[i]);
	}
	printf("\n");

	printf("cond prob L, root C ");
	prob_recursive(root->left, 1, r);
	for (int i = 0; i < 4; i++) {
		printf("%.8f ", root->left->condP[i]);
	}
	printf("\n");
	 
	printf("cond prob L, root G ");
	prob_recursive(root->left, 2, r);
	for (int i = 0; i < 4; i++) {
		printf("%.8f ", root->left->condP[i]);
	}
	printf("\n");

	printf("cond prob L, root T ");
	prob_recursive(root->left, 3, r);
	for (int i = 0; i < 4; i++) {
		printf("%.8f ", root->left->condP[i]);
	}
	printf("\n");

	printf("cond prob LL, L A ");
	prob_recursive(root->left->left, 0, r);
	for (int i = 0; i < 4; i++) {
		printf("%.8f ", root->left->left->condP[i]);
	}
	printf("\n");

	printf("cond prob LL, L C ");
	prob_recursive(root->left->left, 1, r);
	for (int i = 0; i < 4; i++) {
		printf("%.8f ", root->left->left->condP[i]);
	}
	printf("\n");

	printf("cond prob LL, L G ");
	prob_recursive(root->left->left, 2, r);
	for (int i = 0; i < 4; i++) {
		printf("%.8f ", root->left->left->condP[i]);
	}
	printf("\n");

	printf("cond prob LL, L T ");
	prob_recursive(root->left->left, 3, r);
	for (int i = 0; i < 4; i++) {
		printf("%.8f ", root->left->left->condP[i]);
	}
	printf("\n");

	long int count = 0;
	long *** siteFreq = malloc(4 * sizeof(long **));
	for (int i = 0; i < 4; i++) {
		siteFreq[i] = malloc(4 * sizeof(long *));
		for (int j = 0; j < 4; j++) {
			siteFreq[i][j] = calloc(4, sizeof(long));
		}
	}

	while (count < 100000000) {
	//while (count < 100000) {
		root_site(root, statFreq, r);
		simulate_site (root->left, root->state, instRate, r);
		simulate_site (root->right, root->state, instRate, r);
		
		if (nodes[0]->state != pattern[0] ||
		    nodes[1]->state != pattern[1] ||
		    nodes[2]->state != pattern[2] ||
		    nodes[3]->state != pattern[3]) {
			//printf("%d %d %d %d \n", nodes[0]->state,
			//			 nodes[1]->state,
			//			 nodes[2]->state,
			//			 nodes[3]->state);
			continue;
		} 

		siteFreq[nodes[4]->state][nodes[5]->state][nodes[6]->state]++;
/*		for(int i = 0; i < 7; i++)
			printf("%d ", nodes[i]->state); */

		//printf("\n");
		count++;
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				printf("%d %d %d %ld\n", i, j, k, siteFreq[i][j][k]);
			}
		}
	}


}


