#include "MutAnce.h"
gsl_rng * r;  /* global generator */

double simulateDNABranch2Reject(int base, int baseFinal, double ** instRate, double branchLength, const gsl_rng * r) {
  double time = 0, draw;  /* The time that has passed in the simulation */
  double totRate, unif;   /* Total rate of substitions, depends on starting state.
                          unif is used with the random number generator to pick the next base */
  double probBase[4];     /* Given an substition occured, the probability the new base is 0-3 (A-T)*/
  int numMut = 0;
  int startBase = base;

  while (numMut != 2 || base != baseFinal) {
    numMut = 0;
    base = startBase;
    time = 0;

    /* Until the end of the branch is reached */
    while (time <= branchLength ) {

      /* Rate depends on starting state */
      totRate = -instRate[base][base];
    
      /* Time of next substition*/
      draw = gsl_ran_exponential (r, 1 / totRate);
    
      if (time + draw > branchLength) {
        break;
      }
    
      time = time + draw;

      numMut++;

      /* Finds the probability of each substition */
      for (int i = 0; i < 4; i++) {
        if (i == base) {
          probBase[i] = 0.0;
        } else {
          probBase[i] = instRate[base][i] / totRate;
        }
      }
    
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
  }
  return time;
}


double simulateDNABranch1Reject(int base, int baseFinal, double ** instRate, double branchLength, const gsl_rng * r) {
  double time = 0, draw;  /* The time that has passed in the simulation */
  double totRate, unif;   /* Total rate of substitions, depends on starting state.
                          unif is used with the random number generator to pick the next base */
  double probBase[4];     /* Given an substition occured, the probability the new base is 0-3 (A-T)*/
  int numMut = 0;
  int startBase = base;

  while (numMut != 1 || base != baseFinal) {
    numMut = 0;
    base = startBase;
    time = 0;

    /* Until the end of the branch is reached */
    while (time <= branchLength ) {

      /* Rate depends on starting state */
      totRate = -instRate[base][base];
    
      /* Time of next substition*/
      draw = gsl_ran_exponential (r, 1 / totRate);
    
      if (time + draw > branchLength) {
        break;
      }
    
      time = time + draw;

      numMut++;

      /* Finds the probability of each substition */
      for (int i = 0; i < 4; i++) {
        if (i == base) {
          probBase[i] = 0.0;
        } else {
          probBase[i] = instRate[base][i] / totRate;
        }
      }
    
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
  }
  return time;
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

  if (weightSum == 0) {
    fprintf(stderr, "The average substition rate is zero. Exiting. \n");
    exit(1);
  }

  /* Rescales the instantaneous rate matrix so that the average substitution rate is equal to mu */
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j ++) {
      instRate[i][j] = instRate[i][j] / weightSum;
    }
  } 

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
  root->time = 1.37;
  root->state=0;
  
  root->left = calloc(1, sizeof(gnode_t));
  root->left->length = .1;
  root->left->time = .37;
  root->left->state = 1;
  
  /*Generic case */
  //double statFreq[4] = {.19, .3, .24, .27};
  //double statFreq[4] = {.2, .3, .23, .27};
  double statFreq[4] = {.25, .25, .25, .25};
 
  double time, time2;
  
  double kappa = 1;
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

  /* for(int i = 0 ; i < 4; i++) {
  for(int j = 0 ; j < 4; j++) {
  	printf("%f\t", instRate[i][j]);
  }
  printf("\n");
  
  }   */

  /* To calculate the transition probabilites and find the probability of more than two mutations */
  alloc_pmat(root);
  update_pmat_hky(root, statFreq, kappa);
  probMultMut(root->left, root, instRate, r );
  exit(1);

  long samples = 1000000;
  for (int i = 0; i < samples; i++) {
	  
    /* Distribution of times with theory */
    time = solve2Mut(root->left, root, instRate,  r);
    
    /* one mutation */
   // time = simulateBranchOneMut(root->left, root, instRate,  r);
   // time = root->time - root->left->mutTime;

    /* Distribution of times by rejection */
    time2 = simulateDNABranch2Reject(root->state, root->left->state, instRate, root->left->length, r);
    /* one mutation */
//    time2 = simulateDNABranch1Reject(root->state, root->left->state, instRate, root->left->length, r);
    printf("%f\t%f\n", time, time2);
  }

}
