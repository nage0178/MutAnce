#include "MutAnce.h"

int simulateDNABranch(gnode_t * node, gnode_t * parent, double  ** instRate, const gsl_rng * r) {
  double time = 0;        /* The time that has passed in the simulation */
  double totRate, unif;   /* Total rate of substitions, depends on starting state.
                          unif is used with the random number generator to pick the next base */
  double probBase[4];     /* Given an substition occured, the probability the new base is 0-3 (A-T)*/

  int base = parent->state;
  double branchLength = node->length; 
  /* Until the end of the branch is reached */
  node->mut = -1; 

  while (node->mut != node->state) {
	/* Reset stuff */

	node->mutTime = -1;
	node->mut = -1; 
	time = 0;

  	while (time <= branchLength) {
  	  /* Rate depends on starting state */
  	  totRate = -instRate[base][base];

  	  /* Time of next substition*/
  	  time = time + gsl_ran_exponential (r, 1 / totRate);

  	  if (time > branchLength) {
  	   	break; 
  	  }


  	  node->mutTime = parent->time - time;
	  
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

	node->mut = base;
	
  }

  return base;
}


double solve2Mut(gnode_t * node, gnode_t* parent, double **instRate, const gsl_rng * r) {

  double p2ijT = 0; 
  double Qii, Qjj, Qkk, Qik, Qkj;
  double T;
  //double denom, num;
  int i, j;


  j = node->state; 
  i = parent->state; 
  T = parent->time - node->time; 
  if (fabs(parent->time -node->time - node->length) > pow(10, -8)) 
	  fatal("length incorrect %.10f %.10f\n", parent->time- node->time, node->length);
  Qii = instRate[i][i]; 
  Qjj = instRate[j][j];
  
  p2ijT = prob2Mut(node, parent, instRate);

  //p2ijT = 0; 
  /* Solve for pij2T */
  /*for (int k = 0; k < 4; k++) {

    if (k == i || k == j ) 
      continue;

    Qkk = -instRate[k][k];
    Qik = instRate[i][k];
    Qkj = instRate[k][j];

    if ( Qkk == Qjj || Qkk == Qii || Qii == Qjj) {
	    printf("%f %f\n", Qii, Qjj);
      fatal("Qkk is equal to Qjj or Qii. Need to solve this case again.");
    }

    denom = (Qii-Qjj) * (Qii - Qkk) * (Qjj - Qkk);
    num = Qik * Qkj * ( exp(T * Qkk )* (Qii -Qjj) + exp(T * Qii) *(Qjj - Qkk) + exp(T * Qjj) * (Qkk - Qii) ); 
    p2ijT = p2ijT +  num/denom;
  } 

   */

  double lbound = 0, ubound = T; 
  double term1 = 1, term2 = 1, sum = 1;
  double u = gsl_ran_flat(r, 0, 1);
  double t2 =(lbound + ubound) /  2.0; 

  while (fabs(sum) > pow(10, -10) && (ubound - lbound) > pow(10, -10)) {
    t2 =(lbound + ubound) /  2.0; 

    sum = -u;
    for (int k = 0; k < 4; k++) {
      
      if (k == i || k == j ) 
        continue;

      Qkk = instRate[k][k];
      Qik = instRate[i][k];
      Qkj = instRate[k][j];


      if (fabs(Qii- Qjj) > pow (10, -10)  && fabs(Qii - Qkk) > pow( 10, -10) && fabs(Qkk - Qjj) > pow(10, -10)) {
        term1 = Qik * Qkj * exp(T * Qjj) / (Qii - Qkk);
      	term2 = (exp(t2 * (Qii - Qjj)) - 1.0) / (Qii -Qjj) - (exp(t2 * (Qkk - Qjj)) - 1.0) / (Qkk -Qjj);

      } else if (fabs(Qii- Qjj) > pow (10, -10)  && fabs(Qii - Qkk) <= pow( 10, -10) && fabs(Qkk - Qjj) > pow(10, -10)) {
      	term2 = (t2 / (Qkk -Qjj) - 1.0 / pow(Qkk-Qjj,2)) * exp(t2 * (Qkk-Qjj)) + 1.0 / pow(Qkk-Qjj, 2);
	term1 = Qik * Qkj * exp(T * Qjj) ;

      } else if (fabs(Qii- Qjj) <= pow (10, -10)  && fabs(Qii - Qkk) > pow( 10, -10) && fabs(Qkk - Qjj) > pow(10, -10)) {
        term1 = Qik * Qkj * exp(T * Qjj) / (Qii - Qkk);
	term2 = t2 - (exp(t2 * (Qkk-Qjj)) - 1.0) / (Qkk -Qjj);

      } else if (fabs(Qii- Qjj) > pow (10, -10)  && fabs(Qii - Qkk) > pow( 10, -10) && fabs(Qkk - Qjj) <= pow(10, -10)) {
        term1 = Qik * Qkj * exp(T * Qjj) / (Qii - Qkk);
	term2 = (exp(t2 * (Qii-Qkk)) - 1.0)/(Qii -Qkk)  - t2;

      } else if (fabs(Qii- Qjj) <= pow (10, -10)  && fabs(Qii - Qkk) <= pow( 10, -10) && fabs(Qkk - Qjj) <= pow(10, -10)) {
        term1 = Qik * Qkj * exp(T * Qjj) ;
	term2 =  t2 * t2 / 2.0;

      }

      sum = sum + term1 * term2 / p2ijT; 
    }
    

    if ( sum  < 0 )
      lbound = t2; 
    else
      ubound = t2; 


  }

  node->mutTime = parent->time - t2;
  node->mut = node->state;
  return t2; 

}


double prob2Mut (gnode_t * node, gnode_t* parent, double ** instRate ) {
  int i = parent->state;
  int j = node->state;
  double T = node->length;

  double Qii, Qjj, Qkk, Qik, Qkj; 
  Qii = instRate[i][i];
  Qjj = instRate[j][j];
  double sum = 0, num = 1, denom =1;
  
  for (int k = 0; k < 4; k++) {
    if (k == i || k ==j)
      continue;
    
    Qkk = instRate[k][k];
    Qik = instRate[i][k];
    Qkj = instRate[k][j];
    
    if (fabs(Qii- Qjj) > pow (10, -10)  && fabs(Qii - Qkk) > pow( 10, -10) && fabs(Qkk - Qjj) > pow(10, -10)) {
      num = Qik * Qkj * (exp(T * Qkk) * (Qii - Qjj) + exp(T * Qii) * (Qjj -Qkk) + exp(T * Qjj) * (Qkk - Qii) ) ;
      denom = (Qii - Qjj) * (Qii - Qkk) * (Qjj - Qkk);
    
    } else if (fabs(Qii- Qjj) <= pow (10, -10)  && fabs(Qii - Qkk) > pow( 10, -10) && fabs(Qkk - Qjj) > pow(10, -10)) {
      denom = Qii - Qkk; 
      num = Qik * Qkj * exp (T * Qjj) * (T - (exp(T * (Qkk - Qjj)) - 1) / (Qkk - Qjj));
    
    
    } else if (fabs(Qii- Qjj) > pow (10, -10)  && fabs(Qii - Qkk) <= pow( 10, -10) && fabs(Qkk - Qjj) > pow(10, -10)) {
      denom = Qkk - Qjj;
      num = Qik * Qkj * exp(T * Qjj) * (T * exp (T * denom) - (exp(T * denom) - 1)/ denom);
    
    
    } else if (fabs(Qii- Qjj) > pow (10, -10)  && fabs(Qii - Qkk) > pow( 10, -10) && fabs(Qkk - Qjj) <= pow(10, -10)) {
      denom = Qii - Qkk; 
      num = Qik * Qkj * exp(T * Qjj) * ((exp(T * denom) - 1)/denom - T);
    
    } else if (fabs(Qii- Qjj) <= pow (10, -10)  && fabs(Qii - Qkk) <= pow( 10, -10) && fabs(Qkk - Qjj) <= pow(10, -10)) {
      denom = 2;
      num = Qik * Qkj  * exp(T * Qjj) *  pow( T, 2);
    
    } else {
        fatal("Check prob2Mut?");
    }

    sum = sum + (num / denom) ;
    
  }

  return sum; 
}

int probMultMut(gnode_t * node, gnode_t* parent, double ** instRate, const gsl_rng * r ) {
  
  int numMut = 0;
  int i = parent->state;
  int j = node->state;
  double a = instRate[i][i] - instRate[j][j];
  double p1, p2;

  if (fabs(a) < pow(10, -10) )
	p1 = node->length * instRate[i][j] * exp(node->length * instRate[j][j]);
  else 
  	p1  = instRate[i][j] * exp(node->length * instRate[j][j]) * (exp(node->length * a) -1) / a;

  p2 = prob2Mut(node, parent, instRate);

  if ((p1 + p2) / node->pmat[i][j] < 0.99) {
	  printf("Probability of more than two  mutations on branch is greater than .01 : %f\n", 1-(p1+ p2)/node->pmat[i][j]);

  }

  /* Draw number of mutations*/
  double draw = gsl_ran_flat(r, 0, 1);
  if (draw < p1 /(p1 + p2)){
  	numMut = 1;
  }
	
  else {
  	numMut = 2;
  }

  return numMut;

}


int simulateBranchOneMut(gnode_t * node, gnode_t * parent, double  ** instRate, const gsl_rng * r) {

  int i = parent->state;
  int j = node->state;
  double a = instRate[i][i] - instRate[j][j];

  node->mut = node->state; 
  if (fabs(a) < pow(10, -10) ){

  	double time = gsl_ran_flat(r, 0, node->length);
  	node->mutTime = parent->time - time;
	return node->state;
  }


  double u = gsl_ran_flat(r, 0, 1);
  double  time = 1/a * log(u * (exp(node->length * a) - 1) + 1);
  node->mutTime = parent->time - time;

  return node->state;
}

void simulate_mut(gnode_t * node, gnode_t * parent, double ** instRate, gsl_rng *r ) {

  //int finalState; 
 
  /* At fewest, no mutations */
  if (node->state == parent->state || node->state == 4) {
	  simulateDNABranch(node, parent, instRate,  r);
  } else {
  
	/* Find the number of mutations */
	int numMut = probMultMut(node, parent, instRate, r);

	/* Two  mutations */
	if (numMut == 2) {
	      solve2Mut(node, parent, instRate, r);

	/* One mutation */
        } else {
          simulateBranchOneMut(node, parent, instRate,  r);

        }

  }
	if (node->left) {
		simulate_mut(node->left, node, instRate, r);
		simulate_mut(node->right, node, instRate, r);
	}
}

/* Makes the instantaneous rate matrix given the stationary frequencies, substition rate, and rate parameters (in instRate)*/
void makeInstantaneousRate(double kappa, double statFreq[], gtree_t * tree) {

  double ** instRate = tree->instRate;

  if (!instRate) {
  	tree->instRate = xmalloc(4 * sizeof (double * ));
	instRate = tree->instRate;

  	for (int i = 0; i < 4; i++)
		instRate[i] = xmalloc(4 * sizeof (double));

  }

  instRate[0][0] = 0;
  instRate[0][1] = statFreq[1];
  instRate[0][2] = kappa * statFreq[2];
  instRate[0][3] = statFreq[3];

  instRate[1][0] = statFreq[0];
  instRate[1][1] = 0;
  instRate[1][2] = statFreq[2];
  instRate[1][3] = kappa * statFreq[3];
       
  instRate[2][0] = kappa * statFreq[0];
  instRate[2][1] = statFreq[1];
  instRate[2][2] = 0;
  instRate[2][3] = statFreq[3];

  instRate[3][0] = statFreq[0];
  instRate[3][1] = kappa * statFreq[1];
  instRate[3][2] = statFreq[2];
  instRate[3][3] = 0; 

  double diag,  weightSum  = 0;
  /* Finds the diagonal elements and the total substition rate*/
  for (int i = 0; i < 4; i++) {
    diag = 0;
    for (int j = 0; j < 4; j ++) {
      diag = diag + instRate[i][j];
    }
    weightSum = diag * statFreq[i] + weightSum;
    instRate[i][i] = -diag;

  } 

  for (int i = 0; i < 4; i++) {
    diag = 0;
    for (int j = 0; j < 4; j ++) {
      instRate[i][j] = instRate[i][j] / weightSum;
    }

  } 
  
}


//ANNA you need to check this 
int find_mut_descendant(gnode_t * node) {
	int left, right;
	if (node->left) {
		left = find_mut_descendant(node->left);
		right = find_mut_descendant(node->right);
	} else {

		node->mark = (node->mutTime == -1) ? 0 : 1;
		return node->mark;
	}

	if (node->mutTime != -1  && (left == 0 || right == 0 )) {
		node->mark = 1;
		return 1;
	} else {

		node->mark = 0;
		return 0;
	}
}

char * find_pop_origin(gnode_t * node) {
	char * pop;
	int count;
	if (node->mi && node->mi->me) {

		count = node->mi->count - 1;
		pop = node->mi->me[count].target;

		while (count >= 0) {
			if (node->mi->me[count].time > node->mutTime) {
				pop = node->mi->me[count].source;
			} else {
				break;
			}
			count--;
		}
		if (count < 0)
			pop = node->pop;

	} else {
		pop = node->pop;
	}

	return pop;
}
