#include "MutAnce.h"

static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;
static char buffer[LINEALLOC];

static void reallocline(size_t newmaxsize)
{
  char * temp = (char *)xmalloc((size_t)newmaxsize*sizeof(char));

  if (line)
  {
    memcpy(temp,line,line_size*sizeof(char));
    free(line);
  }
  line = temp;
  line_maxsize = newmaxsize;
}


static char * getnextline(FILE * fp)
{
  size_t len = 0;

  line_size = 0;

  /* read from file until newline or eof */
  while (fgets(buffer, LINEALLOC, fp))
  {
    len = strlen(buffer);

    if (line_size + len > line_maxsize)
      reallocline(line_maxsize + LINEALLOC);

    memcpy(line+line_size,buffer,len*sizeof(char));
    line_size += len;

    if (buffer[len-1] == '\n')
    {
      #if 0
      if (line_size+1 > line_maxsize)
        reallocline(line_maxsize+1);

      line[line_size] = 0;
      #else
        line[line_size-1] = 0;
      #endif

      return line;
    }
  }

  if (!line_size)
  {
    free(line);
    line_maxsize = 0;
    line = NULL;
    return NULL;
  }

  if (line_size == line_maxsize)
    reallocline(line_maxsize+1);

  line[line_size] = 0;
  return line;
}


/* The function parses a string into a tree structure. Every time there's an open
parenthesis, two new nodes are created and the function is called twice. The
function keeps track of which character in the treeString (currChar) it is looking
at. */
int makeTree(gnode_t * node, char* treeString, int currChar) {

  /* Note: You should never actually reach the end of string character. The last
  time you return is after the last closed parenthesis */
  while( treeString[currChar] != '\0') {
    /* If the current character is an open parenthesis, make new node, call
    makeTree recursively. currChar is updated along the way */
    if (treeString[currChar] == '(') {
      node->left = (gnode_t *)xcalloc(1, sizeof(gnode_t));
      currChar = makeTree(node->left, treeString, currChar + 1);
      node->right = (gnode_t *)xcalloc(1, sizeof(gnode_t));
      currChar = makeTree(node->right, treeString, currChar);


    } else {
      /* First, find the node name. Then, find the branch length. Return the
      current character */
      char* ptr; /* Needed for the strtod function */
      char* residualTreeString = malloc(strlen(treeString) + 1); /* Used so the original string isn't modified */
      strcpy(residualTreeString, treeString + currChar);

      regex_t regDecimal, regInt;
      int regCompiled, regMatch;
      /* Regex of decimal number */
      regCompiled = regcomp(&regDecimal, "^([0-9]+)((\\.)([0-9]+))$" , REG_EXTENDED);
      if (regCompiled == 1) {
        fprintf(stderr, "Regular expression did not compile.\n");
        exit(1);
      }
      /* Regex of integer number */
      regCompiled = regcomp(&regInt, "^([0-9]+)$" , REG_EXTENDED);
      if (regCompiled == 1) {
        fprintf(stderr, "Regular expression did not compile.\n");
        exit(1);
      }

      /* Finds the nodeName by looking for the next ":". Convert it into an
      integer, save it in the node structure. Update currChar. */
      if (treeString[currChar] != ':') {
      	char* nodeName = strtok(residualTreeString, ":"); /*Note this function modified residualTreeString */

      	node->label = xstrdup(nodeName);
      	currChar = currChar + strlen(nodeName) + 1;

      	residualTreeString = strcpy(residualTreeString, treeString + currChar);
      } else {
	
      	currChar = currChar + 1;

      	residualTreeString = strcpy(residualTreeString, treeString + currChar);
      }

      /* Finds the branch length, converts it to a double, saves it in the node
      structure */
      char* branchLength = strtok(residualTreeString, ",);");

      regMatch = regexec(&regDecimal, branchLength, 0, NULL, 0);
      if (regMatch != 0) {
        fprintf(stderr, "Problem reading in tree file. Regular expression does not match a decimal number.\n");
        exit(1);
      }
      node->length = strtod(branchLength, &ptr);
      currChar = currChar + strlen(branchLength) + 1 ;

      free(residualTreeString);
      regfree(&regDecimal);
      regfree(&regInt);

      /* Returns the updated current character */
      return(currChar);
    }
  }
  return(currChar);
}


void parseTH(char * THstring, gnode_t * root) {

      char* ptr; /* Needed for the strtod function */
      regex_t regDecimal;
      int regCompiled, regMatch;

      regCompiled = regcomp(&regDecimal, "^ (\\[)TH=([0-9]+)((\\.)([0-9]+)), TL=([0-9]+)((\\.)([0-9]+))(\\])$" , REG_EXTENDED);
      if (regCompiled == 1) {
        fprintf(stderr, "Regular expression did not compile.\n");
        exit(1);
      }

      regMatch = regexec(&regDecimal, THstring, 0, NULL, 0);
      if (regMatch != 0) {
        fprintf(stderr, "Problem reading in tree file. Problem with TH and TL.\n");
        exit(1);
      }
	
      char* TH = strtok(THstring+5, ",");
      root->time = strtod(TH, &ptr);

      regfree(&regDecimal);

      return;
}

void countNodes (gnode_t * node, unsigned *inner, unsigned *tip) {
	if (node->left && node->right){
		countNodes(node->left, inner, tip);

		countNodes(node->right, inner, tip);
		(*inner)++;

	} else if (!node->left && !node->right) {
		(*tip)++;


	} else {
		fatal("Gene tree is not binary?");
	}
}

void fill_gtree(gtree_t * gtree, gnode_t * node, unsigned tip) {

	if (node->left && node->right){

		node->left->time = node->time - node->left->length;
		node->right->time = node->time - node->right->length;

		fill_gtree(gtree, node->left, tip);
		fill_gtree(gtree, node->right, tip);

		gtree->nodes[tip + gtree->inner_count++] = node;

	}
	else if (!node->left && !node->right) {
		gtree->nodes[gtree->tip_count++] = node;

	} else {
		fatal("Gene tree is not binary?");
	}
}

long getlinecount(const char * filename)
{
  long linecount = 0;
  FILE * fp;

  fp = xopen(filename,"r");

  /* read number of lines */
  while (getnextline(fp)) linecount++;

  fclose(fp);

  return linecount;
}

gtree_t * loadTree(FILE *fp, long row) {
  
  unsigned tip, inner;
  gtree_t * gtree = NULL;
      
  if (getnextline(fp)) {
    tip = 0;
    inner = 0;
    gtree = xcalloc(1, sizeof(gtree_t));
    
    gnode_t * root = (gnode_t *)xcalloc(1, sizeof(gnode_t));
    gtree->root = root;
    
    char * tree = xstrdup(line);
    int endChar = makeTree(root, tree, 0);
    
    parseTH(tree + endChar, root);
    countNodes(root, &inner, &tip);
    gtree->nodes = xmalloc((inner+tip) * sizeof(snode_t));
    
    fill_gtree(gtree, root, tip);
    
    assert(tip == gtree->tip_count);
    assert(inner == gtree->inner_count);
    
    free(tree);
    //free(line);

  } else {
    fatal("Wrong number of lines in gene tree file?");
  }
   return gtree;
}

void printNewick(gnode_t* node) {

  if (node->left == NULL && node->right == NULL) {
    printf("%s:%f", node->label, node->length);
    return;

  } else {
    printf("(") ;
    printNewick(node->left);
    printf(",");
    printNewick(node->right);
    printf("):%f", node->length);

  }
  return;
}

int parseMigration(gnode_t * node, char * migString, int currChar) {
	char * time, * source, * dest; 
	char * ptr;

	if (node->left){
		currChar = parseMigration(node->left, migString, currChar);
		currChar = parseMigration(node->right, migString, currChar);
	} 

      	char* residualString = malloc(strlen(migString+currChar) + 1); /* Used so the original string isn't modified */
      	strcpy(residualString, migString + currChar);

      	regex_t regDecimal;
      	int regCompiled, regMatch;
      	/* Regex of decimal number */
      	regCompiled = regcomp(&regDecimal, "^([0-9]+)((\\.)([0-9]+))$" , REG_EXTENDED);
      	if (regCompiled == 1) {
      	  fprintf(stderr, "Regular expression did not compile.\n");
      	  exit(1);
      	}

	if (strchr(residualString, ':') == NULL) {
        	fprintf(stderr, "Problem reading in migration file. Check the file is correct.\n");
        	exit(1);
	}
      
      	char* pop = strtok(residualString, ":"); /*Note this function modified residualTreeString */


      	node->pop = xstrdup(pop);
	// Could you move this inside the else? 
	node->mi = xcalloc(1, sizeof(miginfo_t));

	int currChar2 = strlen(pop)+1;	

	if ((residualString + currChar2)[0] == ';') {
		currChar2++;
		free(residualString);
      		regfree(&regDecimal);
		return currChar + currChar2;

	} else {
		int index = 0;
		int numME = 1;
		// Count the number of migration events
		while ((residualString + currChar2 + index)[0] != ';')  {
			if ((residualString + currChar2 + index)[0] == '$')
				numME++;
			index++;
		}
		if (numME == 0 ) {
      			fprintf(stderr, "Problem reading in migation file. Format incorrect. Check file.\n");
			exit(1);
		}

		node->mi->me = xcalloc(numME, sizeof(migevent_t));
		node->mi->count = numME;

		int i = 0;
		while (currChar2 != '\0' &&  i < numME) {
			// Make migration event
			time = strtok(residualString + currChar2, "*");

      			regMatch = regexec(&regDecimal, time, 0, NULL, 0);
      			if (regMatch != 0) {
      			  fprintf(stderr, "Problem reading in migation file. Time does not match the regular expression of a decimal number. Check file.\n");
      			  exit(1);
      			}

			//Make migration event

      			node->mi->me[i].time = strtod(time, &ptr);
      	                currChar2 = currChar2 + strlen(time) + 1 ;

			if (strchr(residualString + currChar2, '&') == NULL) {
		        	fprintf(stderr, "Problem reading in migration file. Check the file is correct.\n");
		        	exit(1);
			}
			
			dest = strtok(residualString + currChar2, "&");
			assert(! node->mi->me[i].target);

			node->mi->me[i].target = xmalloc((strlen(dest) + 1) * sizeof(char));
			strcpy(node->mi->me[i].target, dest);
      	                currChar2 = currChar2 + strlen(dest) + 1 ;

			if (strchr(residualString + currChar2, '$') == NULL && strchr(residualString + currChar2, ';') == NULL ) {
		        	fprintf(stderr, "Problem reading in migration file. Check the file is correct.\n");
		        	exit(1);
			}

			source = strtok(residualString + currChar2, "$;");
			node->mi->me[i].source = xmalloc((strlen(source) + 1) * sizeof(char));
			strcpy(node->mi->me[i].source, source);
      	                currChar2 = currChar2 + strlen(source) + 1 ;

			i++;


			
		}
		if (i != numME) {
			fprintf(stderr, "Problem reading in migration file. Check the file is correct.\n");
			exit(1);
		}
	}

	free(residualString);
      	regfree(&regDecimal);
	return currChar + currChar2; 
}

void loadMigration(FILE * fp, gtree_t * tree) {
    
    if (getnextline(fp)) {

	char * mig = xstrdup(line);
	parseMigration(tree->root, line, 0);
	
	free(mig);
    //	free(line);

    } else {
	fatal("Wrong number of lines in migration  file?");
    }

    return;
}

void loadSubHeader(FILE * fp) {
        char  header[] = "kappa\tpi_T\tpi_C\tpi_A\tpi_G";

        getnextline(fp);
        int cmp = strcmp(header, line);
	if (cmp != 0 ) {
		fatal("Header of substitution model parameter file is formatted incorrectly.");
	}
}

void loadSubstitution(FILE * fp, double * params) {
        char * ret, * ptr;
        int numParams = 5;

        /* Read in MCMC */

        if (getnextline(fp)) {
		ret = strtok(line, "\t");
                params[0] = strtod(ret, &ptr);
                for (int i = 1; i < numParams; i++) {
                        ret = strtok(NULL, "\t");
			if (ret)
                        	params[i] = strtod(ret, &ptr);
			else
				fatal("Incorrect formatting of substitution model file");
                }

        } else 
    		fatal("Wrong number of lines in substitution model  file?");
		

        return;
}


double ** loadTau(char * opt_mcmcfile, long samples, int * numTau, char *** tauLabels) {
	printf("Loading taus... ");
	FILE * fp = xopen(opt_mcmcfile, "r");
	char * ret, * token, * ptr;
	int taus = 0, i = 0, row = 0; 


	getnextline(fp);
	ret = strstr(line, "\ttau");
	while (ret  && ret+1) {
		taus++;
		ret = strstr(ret + 1, "\ttau");
	}

	/*Save tau labels */
	char ** tau_labels = (char **)xmalloc((size_t) taus * sizeof(char *));

	ret = strstr(line, "\ttau") + 1;
	while (i < taus) {
		token = strtok(ret, "\t");
		tau_labels[i] = xstrdup(ret);
		ret += strlen(token) + 1;
		i++;
	}

	/* Allocate memory for MCMC */
	double ** tau = (double **) xmalloc((size_t) taus  * sizeof(double *));
	for (int i =0; i < taus; i++) 
		tau[i] = (double *) xmalloc((size_t) samples * sizeof(double));	

	/* Read in MCMC */
	int preceedCol = 2 * taus + 2;

	while (getnextline(fp) && (row < samples)) {
		ret = strtok(line, "\t");
		for (int i = 1; i < preceedCol; i++) {
			ret = strtok(NULL, "\t");
		}
		for (int i = 0; i < taus; i++) {
			ret = strtok(NULL, "\t");
			tau[i][row] = strtod(ret, &ptr);
		}
		
		row++;
	}

	assert(row == samples);
	fclose(fp);

	*tauLabels = tau_labels;
	*numTau = taus;

	printf("Done\n");
	return tau; 
}

void freeline () {
	free(line);
}
