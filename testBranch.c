#include "MutAnce.h"
gsl_rng * r;  /* global generator */

int simulateDNABranchUncond(int base, double ** instRate, double branchLength, const gsl_rng * r) {
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
void makeInstantaneousRateTest(double ** instRate, double statFreq[]) {

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


void main () {

	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	gnode_t ** nodes;
	gnode_t * root = calloc(1, sizeof(gnode_t));
	root->length = 0 ;
	root->time = .27;

	root->left = calloc(1, sizeof(gnode_t));
	root->left->length = .08;
	root->left->time = .19;


	double statFreq[4] = {.2, .3, .23, .27};
	double kappa = 2;
	/* Set the values of the instantaneous rate matrix. Sets diagonal and rescales to makes the mutation rate in the function. */
	double rateParams[6] = {1, 1, 1, 1, 1, 1};
	rateParams[1] = rateParams[4] = kappa;
	double ** instRate = malloc(4 * sizeof(double * ));
	for(int i =0; i < 4; i++) {
		instRate[i] = malloc(4 * sizeof(double));
		instRate[i][i] = 0;
	}

	instRate[0][1] = rateParams[0] * statFreq[1];
	instRate[0][2] = rateParams[1] * statFreq[2];
	instRate[0][3] = rateParams[2] * statFreq[3];
                                                    ;
	instRate[1][0] = rateParams[0] * statFreq[0];
	instRate[1][2] = rateParams[3] * statFreq[2];
	instRate[1][3] = rateParams[4] * statFreq[3];
                                                    ;
	instRate[2][0] = rateParams[1] * statFreq[0];
	instRate[2][1] = rateParams[3] * statFreq[1];
	instRate[2][3] = rateParams[5] * statFreq[3];
                                                    ;
	instRate[3][0] = rateParams[2] * statFreq[0];
	instRate[3][1] = rateParams[4] * statFreq[1];
	instRate[3][2] = rateParams[5] * statFreq[2];
	makeInstantaneousRateTest(instRate, statFreq);

	
	root->state = 0;
	root->left->state = 1;
	long count = 0;
	//printf("%f \n", instRate[0][0] - instRate[1][1]);
	//exit(1);
 	long samples = 10000000;
	double * times = malloc(sizeof(double) * samples);
	while (count < samples) {
		simulateDNABranch(root->left,  root, instRate, r);
		times[count] = root->left->mutTime - root->left->time; 
		printf("%f\n", times[count]);
		//Reset anything
		count++;
	}

	for(int i =0; i < 4; i++) 
		free(instRate[i]);
	free(instRate);
	free(times);
	free(root->left);
	free(root);
	
}


