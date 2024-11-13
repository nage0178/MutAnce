#include "MutAnce.h"

static int cb_cmp_double(const void * a, const void * b)
{
  double * x = (double *)a;
  double * y = (double *)b;

  if ( *x > *y) return 1;
  if ( *x < *y) return -1;
  return 0;
}

static void hpd_interval(double * x,
                         long n,
                         double * ltail,
                         double * rtail,
                         double alpha)
{
  long lrow = (long)(n*alpha/2);
  long urow = (long)(n*(1-alpha/2));
  long diffrow = urow - lrow;
  long l,r;

  long left = lrow;


  double w = x[urow] - x[lrow];

  *ltail = x[lrow];
  *rtail = x[urow];

  if (n <= 2) return;

  for (l=0,r=l+diffrow; r < n; l++,r++)
  {
    if (x[r] - x[l] < w)
    {
      left = l;
      w = x[r] - x[l];
    }
  }

  *ltail = x[left];
  *rtail = x[left + diffrow];
}


void summary (mutHist_t ** mutationHist, int varsite, long samples, msa_t * msa, long opt_site) {

  int site; 
  int * variableSites = xmalloc(varsite * sizeof(int));
  varsite = 0;

  if (opt_site > -1) {
    variableSites[0] = opt_site;
    varsite = 1;
  } else {
    for (site = 0; site < msa->length; site ++ ) {
      if (msa->variable_sites[site]) {
        variableSites[varsite] = site;
        varsite++;
      }
    }
  } 

  double ** sum_time = xmalloc(varsite * sizeof(double  *));
  long  ** sum_base = xmalloc(varsite * sizeof(long  *));
  long  ** sum_root = xmalloc(varsite * sizeof(long  *));
  for (int i = 0; i < varsite; i++) {
	  sum_time[i] = xmalloc(samples * sizeof(double));
	  sum_base[i] = xcalloc(4,  sizeof(long));
	  sum_root[i] = xcalloc(4, sizeof(long));
  }

  long * sum_mult_mut = xcalloc(varsite, sizeof(long));
  long * sum_two_mut = xcalloc(varsite, sizeof(long));

  list_t * mutPop = xcalloc(1, sizeof(list_t));
  list_item_t * item;
  for (int i = 0; i < samples; i++) { 
	  for (int j = 0; j < varsite; j++) {
		  if (mutationHist[i][j].numMut == 1) {
			  sum_time[j][i] = mutationHist[i][j].mutations[0].time;
			  sum_root[j][mutationHist[i][j].root]++;
			  sum_base[j][mutationHist[i][j].mutations[0].base]++;

			   item = mutPop->head;
			  long k;
			  for (k = 0; k < mutPop->count; k++) {
			        if (!strcmp(mutationHist[i][j].mutations[0].pop, (char *) item->data))
			        	break;
				item = item->next;
			  }
			  if (k == mutPop->count){
			          char * data = strdup(mutationHist[i][j].mutations[0].pop);
			          list_append(mutPop, data);
			  }

		  } else {
			  if (mutationHist[i][j].numMut == 2) {
			  	sum_two_mut[j]++;
			  }
			  sum_mult_mut[j]++;
			  sum_time[j][i] = 0;

		  }

		  for (int m = 0; m < mutationHist[i][j].numMut; m++) {
			  item = mutPop->head;
			  long k;
			  for (k = 0; k < mutPop->count; k++) {
			        if (!strcmp(mutationHist[i][j].mutations[m].pop, (char *) item->data))
			        	break;
				item = item->next;
			  }
			  if (k == mutPop->count){
			          char * data = strdup(mutationHist[i][j].mutations[m].pop);
			          list_append(mutPop, data);
			  }
		  }
	  }

  }


  //printf("\n");
  int siteStart = 0;
  int siteEnd = msa->length;
  if (opt_site > -1) {
    siteStart = opt_site;
    siteEnd = opt_site +1; 
  }

  for (site = siteStart; site < siteEnd; site ++ ) {
    if (msa->variable_sites[site]) {
      printf("%d ", site + 1);
            for(int i = 0; i < msa->count; i++)
        	    printf("%c", msa->sequence[i][site]);
    printf("\n");
    }
  }

  printf("site ");
  for (site = 0; site < varsite; site++) {
	  printf("%d\t", variableSites[site] + 1 );
  }
  printf("\n");
  printf("\n");
  
  printf("Posterior probability multiple mutations ");
  for (int i = 0; i < varsite; ++i) {
	  printf("%.4f ", (double) sum_mult_mut[i]/ (double) samples);
  }
  printf("\n\nPosterior probability of three or more mutations ");
  for (int i = 0; i < varsite; ++i) {
	  printf("%.4f ", (double) (sum_mult_mut[i]-sum_two_mut[i])/ (double) samples);
  }

  printf("\n");
  printf("\n");

  item = mutPop->head;

  char ** pops = xcalloc(mutPop->count, sizeof(char *));
  for (int i = 0; i < mutPop->count; i++) {
    pops[i] = (char *) item->data;
    item = item->next;
  }

  long ** sum_pop = xcalloc(varsite, sizeof(long *));
  for (int i = 0 ; i < varsite; i++) {
	  sum_pop[i] = xcalloc(mutPop->count, sizeof(long));
  }

  double ** mut2_time1 = xcalloc(varsite, sizeof(double *));
  double ** mut2_time2 = xcalloc(varsite, sizeof(double *));
  long ** baseTwo1 = xcalloc(varsite, sizeof(long *));
  long ** baseTwo2 = xcalloc(varsite, sizeof(long *));
  long *** jointMut2 = xmalloc(sizeof(long *) * varsite);

  long ** rootTwo = xcalloc(varsite, sizeof(long *));
  long ** pop1 = xcalloc(varsite, sizeof(long *));
  long ** pop2 = xcalloc(varsite, sizeof(long *));
  long *** jointPop2 = xmalloc(sizeof(long *) * varsite);

  for (int i = 0; i < varsite; i++) {
    mut2_time1[i] = xcalloc(sum_two_mut[i],  sizeof(double*));
    mut2_time2[i] = xcalloc(sum_two_mut[i],  sizeof(double*));
    baseTwo1[i] = xcalloc(4, sizeof(long));
    baseTwo2[i] = xcalloc(4, sizeof(long));
    rootTwo[i]  = xcalloc(4, sizeof(long));
    pop1[i] = xcalloc(mutPop->count, sizeof(long));
    pop2[i] = xcalloc(mutPop->count, sizeof(long));

    jointMut2[i] = xcalloc(4, sizeof(long *));
    for (int j = 0; j < 4; j++) {
     jointMut2[i][j] = xcalloc(4, sizeof (long));
    }

    jointPop2[i] = xcalloc(mutPop->count, sizeof(long *));
    for (int j = 0; j < mutPop->count; j++) {
      jointPop2[i][j] = xcalloc(mutPop->count, sizeof (long));
    }
  }

  long * mut2_count = xcalloc(varsite, sizeof(long));
  int first, second;
  int pop1Index, pop2Index;

  for (int i = 0; i < samples; i++) { 
    for (int j = 0; j < varsite; j++) {
      long k;
      if (mutationHist[i][j].numMut == 1) {
         item = mutPop->head;

        for (k = 0; k < mutPop->count; k++) {
          if (!strcmp(mutationHist[i][j].mutations[0].pop, pops[k])) {
           	sum_pop[j][k]++; 
		break;
	  }
        }

      } else if (mutationHist[i][j].numMut == 2) {
        rootTwo[j][mutationHist[i][j].root]++;

        if (mutationHist[i][j].mutations[0].time < mutationHist[i][j].mutations[1].time) {
                first = 0;
                second = 1;
        } else {
                first = 1;
                second = 0;
        }

        mut2_time1[j][mut2_count[j]] = mutationHist[i][j].mutations[first].time;
        mut2_time2[j][mut2_count[j]] = mutationHist[i][j].mutations[second].time;

        int base1 = mutationHist[i][j].mutations[first].base;
        int base2 = mutationHist[i][j].mutations[second].base;

        /* Marginal distribution of bases */ 
        baseTwo1[j][base1]++;
        baseTwo2[j][base2]++;

        /* Joint distribution of bases */
        jointMut2[j][base2][base1]++;
        

        /* Marginal distribution of populations */
        for (k = 0; k < mutPop->count; k++) {
          if (!strcmp(mutationHist[i][j].mutations[first].pop, pops[k])){
             pop1[j][k]++; 
             pop1Index = k;
             break;
          }
        }
        assert(k < mutPop->count);

        for (k = 0; k < mutPop->count; k++) {
          if (!strcmp(mutationHist[i][j].mutations[second].pop, pops[k])) {
            pop2[j][k]++; 
            pop2Index = k;
            break;
          }
        }  
        assert(k < mutPop->count);

        /* Joint distribution of populations */
        jointPop2[j][pop2Index][pop1Index]++;

        mut2_count[j]++;

      }
    }
  }

  for (int i = 0; i < varsite; i++) 
  	qsort(sum_time[i], samples, sizeof(double), cb_cmp_double);

  for (int i = 0; i < varsite; i++) {
  	qsort(mut2_time1[i], mut2_count[i], sizeof(double), cb_cmp_double);
  	qsort(mut2_time2[i], mut2_count[i], sizeof(double), cb_cmp_double);
  }

  printf("population of mutation\n");
  for (int k = 0; k < mutPop->count; k++) {
	  printf("1_mut_pop_%s\t", pops[k]);
    for (int j = 0; j < varsite; j++) {
      printf("%f\t", (double) sum_pop[j][k] / (samples - sum_mult_mut[j]));
    }
  printf("\n");
  }
  printf("\n");

  printf("mean_time\t");
  double * mean = xcalloc(varsite , sizeof(double));
  for (int i = 0; i < varsite; ++i)
  {
    double sum = 0;
    for (int j = 0; j < samples; ++j)
      sum += sum_time[i][j];

    mean[i] = sum/(samples - sum_mult_mut[i]);
  if ( (double) sum_mult_mut[i]/ (double) samples != 1) 
    fprintf(stdout, "%.10f\t", mean[i]);
  else 
    fprintf(stdout, "NA\t");
  }
  printf("\n1_mut_02.5_HPD\t"); 

  double * hpd025 = (double *)xcalloc(varsite * 3 , sizeof(double));
  double * hpd975 = (double *)xcalloc(varsite * 3 , sizeof(double));

  for (int i = 0; i < varsite; ++i) {
	  
	  if ( (double) sum_mult_mut[i]/ (double) samples != 1) {
    		hpd_interval(sum_time[i] + sum_mult_mut[i], samples - sum_mult_mut[i], hpd025 + i, hpd975 + i, 0.05);
	  }

	  if ( (double) sum_mult_mut[i]/ (double) samples != 0) {
    hpd_interval(mut2_time1[i] , sum_two_mut[i], hpd025 + varsite +  i, hpd975 + (varsite)  +  i, 0.05);
    hpd_interval(mut2_time2[i] , sum_two_mut[i], hpd025 + (varsite * 2) + i, hpd975 + (varsite * 2) +  i, 0.05);
	  }
  }
  
  for (int i = 0; i < varsite; ++i) {

  if ( (double) sum_mult_mut[i]/ (double) samples != 1) 
    fprintf(stdout, "%.10f\t", hpd025[i]);
  else 
    fprintf(stdout, "NA\t");
  }
 
  printf("\n1_mut_97.5_HPD\t"); 
  for (int i = 0; i < varsite; ++i) {
  if ( (double) sum_mult_mut[i]/ (double) samples != 1) 
    fprintf(stdout, "%.10f\t", hpd975[i]);
  else
    fprintf(stdout, "NA\t");
  }
  printf("\n");
  printf("\n");


  char dna[5] = "ACGT";
  printf("root state\n"); 
  for (int j = 0; j < 4; j++) {
      printf("1_mut_root_%c ", dna[j]);
    for (int i = 0; i < varsite; ++i) {
	    printf("%f\t", (double) sum_root[i][j] / (samples - sum_mult_mut[i]));
    }
    printf("\n");
  }
  
  printf("\n");

  printf("mutation\n"); 
  for (int j = 0; j < 4; j++) {
      printf("1_mut_mutation_%c ", dna[j]);
    for (int i = 0; i < varsite; ++i) {
	    printf("%f\t", (double) sum_base[i][j] / (samples - sum_mult_mut[i]));
    }
    printf("\n");
  }

  printf("\n");
  /* Two mutations */
  printf("2_mut_time_1\t");

  double * mean1 = xcalloc(varsite , sizeof(double));
  double * mean2 = xcalloc(varsite , sizeof(double));

  for (int i = 0; i < varsite; ++i)
  {
    double sum1 = 0;
    double sum2 = 0;
    for (int j = 0; j < sum_two_mut[i]; ++j) {
      sum1 += mut2_time1[i][j];
      sum2 += mut2_time2[i][j];
    }

    mean1[i] = sum1/(sum_two_mut[i]);
    mean2[i] = sum2/(sum_two_mut[i]);
    fprintf(stdout, "%.10f\t", mean2[i]);
  }

  printf("\n2_mut_02.5_HPD_1\t"); 
  for (int i = 0; i < varsite; ++i)
    fprintf(stdout, "  %.10f", hpd025[varsite * 2 + i]);
 
  printf("\n2_mut_97.5_HPD_1\t"); 
  for (int i = 0; i < varsite; ++i)
    fprintf(stdout, "%.10f\t", hpd975[varsite * 2 + i]);
  printf("\n\n");


  printf("\n\n");
  printf("2_mut_time_2\t");
  for (int i = 0; i < varsite; ++i)
  {
    fprintf(stdout, "%.10f\t", mean1[i]);
  }

  printf("\n2_mut_02.5_HPD_2\t"); 
  for (int i = 0; i < varsite; ++i)
    fprintf(stdout, "%.10f\t", hpd025[varsite + i]);
 
  printf("\n2_mut_97.5_HPD_2\t"); 
  for (int i = 0; i < varsite; ++i)
    fprintf(stdout, "%.10f\t", hpd975[varsite + i]);


  printf("\n\nroot state\n");
  for (int j = 0; j < 4; j++) {
    printf("2_mut_root_%c\t", dna[j]);
    for (int i = 0; i < varsite; i++) {
        printf("%f\t", (double) rootTwo[i][j] / sum_two_mut[i]);
    }
    printf("\n");
  }
  printf("\n");

  printf("First mutation\n");
  for (int j = 0; j < 4; j++) {
    printf("2_mut_mutation_1_%c\t", dna[j]);
    for (int i = 0; i < varsite; i++) {
        printf("%f\t", (double) baseTwo2[i][j] / sum_two_mut[i]);
    }
    printf("\n");
  }
  printf("\n");
  
  printf("Second mutation\n");
  for (int j = 0; j < 4; j++) {
    printf("2_mut_mutation_2_%c\t", dna[j]);
    for (int i = 0; i < varsite; i++) {
        printf("%f\t", (double) baseTwo1[i][j] / sum_two_mut[i]);
    }
    printf("\n");
  }
  printf("\n");


  printf("population of mutation 1\n");
  for (int k = 0; k < mutPop->count; k++) {
	  printf("2_mut_pop_1_%s\t", pops[k]);
    for (int j = 0; j < varsite; j++) {
      printf("%f\t", (double) pop2[j][k] / (sum_two_mut[j]));
    }
  printf("\n");
  }
  printf("\n");

  printf("population of mutation 2\n");
  for (int k = 0; k < mutPop->count; k++) {
	  printf("2_mut_pop_2%s\t", pops[k]);
    for (int j = 0; j < varsite; j++) {
      printf("%f\t", (double) pop1[j][k] / (sum_two_mut[j]));
    }
  printf("\n");
  }
  printf("\n");

  printf("Joint probability distribution of mutations to bases conditional on two mutations\n");
  printf("Row: first mutation, Column: second mutation\n");
  for (int k = 0; k < varsite; k++) {
    printf("site %d\n", variableSites[k] + 1 );

    printf("\t%c\t%c\t%c\t%c\n", dna[0], dna[1], dna[2], dna[3]);
    /* First mutation */
    for (int i = 0; i < 4; i++) {
      printf("\t%c\t", dna[i]);

      /* Second mutation (forward time) */
      for (int j = 0; j < 4; j++) {
        printf("%f\t", (double) jointMut2[k][i][j] / (sum_two_mut[k]));
      }
      printf("\n");
    }
    printf("\n\n");
  }

  printf("Joint probability distribution of populatons conditional on two mutations\n");
  printf("Row: first populaton, Column: second population\n");
  for (int k = 0; k < varsite; k++) {
    printf("site %d\n", variableSites[k] + 1 );

    for (int i = 0; i < mutPop->count; i++)
      printf("\t%s", pops[i]);
    printf("\n");

    /* First mutation */
    for (int i = 0; i < mutPop->count; i++) {
      printf("%s\t", pops[i]);

      /* Second mutation (forward time) */
      for (int j = 0; j < mutPop->count; j++) {
        printf("%f\t", (double) jointPop2[k][i][j] / (sum_two_mut[k]));
      }
      printf("\n");
    }
    printf("\n\n");
  }
  /* Clean up memory */

  for (int i = 0; i < varsite; i++) {
	  free(mut2_time1[i]);
	  free(mut2_time2[i]);
	  free(baseTwo1[i]);
	  free(baseTwo2[i]);
	  free(rootTwo[i]);
	  free(sum_time[i]);
	  free(sum_base[i]);
	  free(sum_root[i]);
	  free(sum_pop[i]);
	  free(pop1[i]);
	  free(pop2[i]);

    for (int j = 0; j < 4; j++) 
      free(jointMut2[i][j]);

    for (int j = 0; j < mutPop->count; j++)
      free(jointPop2[i][j]);

    free(jointMut2[i]);
    free(jointPop2[i]);

  }

  free(mut2_time1);
  free(mut2_time2);
  free(mut2_count);
  free(hpd025);
  free(hpd975);
  free(baseTwo1);
  free(baseTwo2);
  free(rootTwo);
  free(sum_base);
  free(sum_root);
  free(sum_time);
  free(sum_pop);
  free(sum_mult_mut);
  free(sum_two_mut);
  free(mean);
  free(mean1);
  free(mean2);
  free(variableSites);
  free(pop1);
  free(pop2);
  free(jointMut2);
  free(jointPop2);

  free(pops);
  list_clear(mutPop, free);
  free(mutPop);
}
