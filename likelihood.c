#include "MutAnce.h"

void cal_condP (gnode_t * node) {
	if (node->left) {
		cal_condP(node->left);
		cal_condP(node->right);
	}

	//Calculate likelhood based on left and right
	// States at parent node
	if (node->left) {
		for (int j = 0; j < 4; j++ ) {
			/* From i to j */
			node->likelihood[j] = 0;
			for (int i = 0; i < 4; i ++ ) {
				for (int k = 0; k < 4; k ++ ) {
					node->likelihood[j] += (node->left->likelihood[i] * node->left->pmat[j][i] )* (node->right->likelihood[k] * node->right->pmat[j][k]);
				}
				
			}
		}

	} 
 }

void prob_recursive (gnode_t * node, int p_state, gsl_rng *r) {
	int i;
	double sum = 0;
	
	for (int j = 0; j < 4; j++ ) 
		sum += node->likelihood[j] * node->pmat[p_state][j];
	for (int j = 0; j < 4; j++ ) {
		node->condP[j] = node->likelihood[j] * node->pmat[p_state][j]/sum;

	}
	

	/* Draw  state */
	double state = gsl_ran_flat(r, 0, 1);	
	sum = 0;
	for (i = 0; i < 4; i++) {
		if (state < sum + node->condP[i])
			break;
		 sum += node->condP[i];
	}

	assert(i < 4);
	node->state = i;


	if (node->left && node->left->left) {
		prob_recursive(node->left, node->state, r);
	}
	if (node->right && node->right->right) {
		prob_recursive(node->right, node->state, r);
	}
}

void condP_root (gnode_t * node, double * freq, gsl_rng *r) {

	int i;

	double sum = 0;
	for (int j = 0; j < 4; j++ ) {
		sum += node->likelihood[j] * freq[j];
	}
	for (int j = 0; j < 4; j++ ) {
		node->condP[j] = node->likelihood[j] * freq[j] /sum;
	}

	/* Draw root state */
	double state = gsl_ran_flat(r, 0, 1);	
	sum = 0;
	for (i = 0; i < 4; i++) {
		if (state < sum + node->condP[i])
			break;

		 sum += node->condP[i];
	}
	assert(i < 4);
	node->state = i;

	if (node->left && node->left->left) {
		prob_recursive(node->left, node->state, r);
	}
	if (node->right && node->right->right) {
		prob_recursive(node->right, node->state, r);
	}

}

void alloc_pmat(gnode_t * node) {
	if (node->left) {
		alloc_pmat(node->left);
	}
	if (node->right){
		alloc_pmat(node->right);
	}

        node->pmat = xmalloc(4 * sizeof (double *));
	for (int i = 0; i < 4; i++)
		node->pmat[i] = xcalloc(4, sizeof(double));

}

void update_pmat_hky(gnode_t * root,
			    const double * freqs, 
			    double kappa)
{
  double bt;
  double a1t,a2t;
  double e1,e2,e3;
  double A,C,G,T,Y,R;

  double ** pmat;

  pmat = root->pmat;
  double bl = root->length;

  A = freqs[0];
  C = freqs[1];
  G = freqs[2];
  T = freqs[3];
  Y = T + C;
  R = A + G;

  double mr = 1 / (2*T*C*kappa + 2*A*G*kappa + 2*Y*R);
  bt = bl*mr;
  a1t = a2t = kappa*bt;

  e1 = expm1(-bt);
  e2 = expm1(-(R*a2t + Y*bt));
  e3 = expm1(-(Y*a1t + R*bt));

  /* row column */
  pmat[0][0]  = 1 + Y*A / R*e1 + G / R*e2;
  pmat[0][1]  = -C*e1;
  pmat[0][2]  = Y*G / R*e1 - G / R*e2;
  pmat[0][3]  = -T*e1;

  pmat[1][0]  = -A*e1;
  pmat[1][1]  = 1 + (R*C*e1 + T*e3) / Y;
  pmat[1][2]  = -G*e1;
  pmat[1][3]  = (R*e1 - e3)*T / Y;

  pmat[2][0]  = Y*A / R*e1 - A / R*e2;
  pmat[2][1]  = -C*e1;
  pmat[2][2] = 1 + Y*G / R*e1 + A / R*e2;
  pmat[2][3] = -T*e1;

  pmat[3][0] = -A*e1;
  pmat[3][1] = (R*e1 - e3)*C / Y;
  pmat[3][2] = -G*e1;
  pmat[3][3] = 1 + (R*T*e1 + C*e3) / Y;

  if (!(root->left)) return;

  update_pmat_hky(root->left, freqs, kappa);

  if (!(root->right)) return;
  update_pmat_hky(root->right, freqs, kappa);

}

void matchGtreeSeqs(gtree_t * tree, msa_t  * msa) {

	for (int i = 0; i < msa->count; i ++ ) {

		for (int j = 0; j < tree->tip_count; j ++) {

			if(!strcmp(msa->label[i], tree->nodes[j]->label)) {
				tree->nodes[j]->msa_index = i;
				break;
			}
		}

		assert(i < msa->count);
	}
}

void likelihoodTips (gtree_t * gtree, int site, msa_t * msa) {

	for (int i = 0; i < gtree->tip_count; i++) {
		gnode_t * node = gtree->nodes[i];
		long msa_index = node->msa_index;

		if (msa->sequence[msa_index][site] == 'a' || msa->sequence[msa_index][site] == 'A') {
			node->likelihood[0] = 1;
			node->likelihood[1] = 0;
			node->likelihood[2] = 0;
			node->likelihood[3] = 0;
			node->state = 0;

		} else if (msa->sequence[msa_index][site] == 'c' || msa->sequence[msa_index][site] == 'C') {
			node->likelihood[0] = 0;
			node->likelihood[1] = 1;
			node->likelihood[2] = 0;
			node->likelihood[3] = 0;
			node->state = 1;

		} else if (msa->sequence[msa_index][site] == 'g' || msa->sequence[msa_index][site] == 'G') {
			node->likelihood[0] = 0;
			node->likelihood[1] = 0;
			node->likelihood[2] = 1;
			node->likelihood[3] = 0;
			node->state = 2;

		} else if (msa->sequence[msa_index][site] == 't' || msa->sequence[msa_index][site] == 'T') {
			node->likelihood[0] = 0;
			node->likelihood[1] = 0;
			node->likelihood[2] = 0;
			node->likelihood[3] = 1;
			node->state = 3;
		} else if (msa->sequence[msa_index][site] == '-' || msa->sequence[msa_index][site] == 'N'|| msa->sequence[msa_index][site] == 'n' ){
			node->likelihood[0] = 1.0;
			node->likelihood[1] = 1.0;
			node->likelihood[2] = 1.0;
			node->likelihood[3] = 1.0;
			node->state = 4;
		}
		if (node->state != 4)
			assert(node->likelihood[0] + node->likelihood[1] + node->likelihood[2] + node->likelihood[3] == 1);
	}

}

