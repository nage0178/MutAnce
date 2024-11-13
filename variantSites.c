#include "MutAnce.h"

void findVariant (msa_t * msa) {

	msa->variable_sites = xcalloc(msa->length, sizeof(int));

	msa->variableCount = 0;
	/* For each site */
	for (int k = 0; k < msa->length; k++) {

		/* For each sequence */
		for (int j = 1; j < msa->count; j++) {

			if (msa->sequence[j][k] != msa->sequence[0][k] ) {
				//Make array/list that defines if sites are invariant 
				msa->variable_sites[k]= 1 ;
				msa->variableCount++;
				break;
			}
		}
	}

}
