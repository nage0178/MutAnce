#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#define LINEALLOC 2048

#define BPP_FAILURE  0
#define BPP_SUCCESS  1

/* error codes */

#define ERROR_PHYLIP_SYNTAX            106
#define ERROR_PHYLIP_LONGSEQ           107
#define ERROR_PHYLIP_NONALIGNED        108
#define ERROR_PHYLIP_ILLEGALCHAR       109
#define ERROR_PHYLIP_UNPRINTABLECHAR   110
#define ERROR_PARSE_MORETHANEXPECTED   111
#define ERROR_PARSE_LESSTHANEXPECTED   112
#define ERROR_PARSE_INCORRECTFORMAT    113

#define BPP_DNA_MODEL_JC69              0
#define BPP_DNA_MODEL_K80               1
#define BPP_DNA_MODEL_F81               2
#define BPP_DNA_MODEL_HKY               3
#define BPP_DNA_MODEL_T92               4
#define BPP_DNA_MODEL_TN93              5
#define BPP_DNA_MODEL_F84               6
#define BPP_DNA_MODEL_GTR               7


extern long opt_locus_count;
extern int bpp_errno;
extern  char bpp_errmsg[200];
extern const unsigned int pll_map_fasta[256];
extern const unsigned int pll_map_amb[256];
extern const unsigned int pll_map_nt[256];


typedef struct mutEvent_s
{
  double time; 
  char * pop;
  int base;
} mutEvent_t;

typedef struct mutHist_s
{
  int site;
  int numMut;
  int root;
  mutEvent_t * mutations;
} mutHist_t;

typedef struct list_item_s
{
  void * data;
  struct list_item_s * next;
} list_item_t;

typedef struct list_s
{
  list_item_t * head;
  list_item_t * tail;
  long count;
} list_t; 

typedef struct snode_s
{
  char * label;
  double tau;
  struct snode_s * left;
  struct snode_s * right;
  struct snode_s * parent;

} snode_t;

typedef struct migevent_s
{
  double time;
  char * source;
  char * target;

} migevent_t;

typedef struct miginfo_s
{
  long alloc_size;
  long count;

  migevent_t * me;
} miginfo_t;

typedef struct gtree_s 
{
  unsigned int tip_count;
  unsigned int inner_count;

  struct gnode_s ** nodes;
  struct gnode_s * root;

  double ** instRate;
} gtree_t;

typedef struct gnode_s 
{
  char * label;
  double length;
  double time;
  struct gnode_s * left;
  struct gnode_s * right;
  //struct gnode_s * parent;
  double likelihood[4];
  double condP[4];
  double ** pmat;
  int msa_index;
  int state; 

  int mut;
  double mutTime;

//  snode_t * pop;
  char * pop;
  int mark;

  miginfo_t * mi;
} gnode_t;


typedef struct msa_s
{
  int count;
  int length;

  char ** sequence;
  char ** label;

  int amb_sites_count;
  int original_length;

  double * freqs;

  int dtype;
  int model;
  int original_index;

  int * variable_sites;
  int variableCount;

} msa_t;

typedef struct phylip_s
{
  FILE * fp;
  char * line;
  size_t line_size;
  size_t line_maxsize;
  char buffer[LINEALLOC];
  const unsigned int * chrstatus;
  long no;
  long filesize;
  long lineno;
  long stripped_count;
  long stripped[256];
} phylip_t;

void freeTree(gnode_t * node);

/* util */
FILE * xopen(const char * filename, const char * mode);
void * xmalloc(size_t size);
void * xcalloc(size_t nmemb, size_t size);
char * xstrdup(const char * s);
char * xstrchrnul(char *s, int c);
void fatal(const char * format, ...);


/*parse */
gtree_t * loadTree(FILE *fp, long row);
//static char * getnextline(FILE * fp);
void printNewick(gnode_t* node) ;
void fill_gtree(gtree_t * gtree, gnode_t * node, unsigned tip);
long getlinecount(const char * filename);
void loadMigration(FILE * fp, gtree_t * tree);
double ** loadTau(char * opt_mcmcfile, long samples, int * numTau, char *** tauLabels);
void loadSubstitution(FILE * fp, double * params);
void loadSubHeader(FILE *fp);
void freeline();


/* phylip */

phylip_t * phylip_open(const char * filename,
                       const unsigned int * map);
msa_t ** phylip_parse_multisequential(phylip_t * fd, long * count);
void phylip_close(phylip_t * fd);
void msa_destroy(msa_t * msa);
void msa_ambiguous_sites(msa_t * msa, const unsigned int * map);

void findVariant (msa_t *msa);

/* likelihood */
void alloc_pmat(gnode_t * node);
void update_pmat_hky(gnode_t * root,
                            const double * freqs,
                            double kappa);
void cal_condP (gnode_t * node);
void matchGtreeSeqs(gtree_t * tree, msa_t  * msa);
void likelihoodTips (gtree_t * gtree, int site, msa_t * msa);
void condP_root (gnode_t * node, double * freq, gsl_rng *r);
void prob_recursive (gnode_t * node, int p_state, gsl_rng *r);

/* Simulate */
void makeInstantaneousRate(double kappa, double statFreq[], gtree_t * tree);
void simulate_mut(gnode_t * node, gnode_t * parent, double ** instRate, gsl_rng *r );
char * find_pop_origin(gnode_t * node);
int find_mut_descendant(gnode_t * node);
int simulateDNABranch(gnode_t * node, gnode_t * parent, double  ** instRate, const gsl_rng * r);
double solve2Mut(gnode_t * node, gnode_t* parent, double **instRate, const gsl_rng * r);
double prob2Mut (gnode_t * node, gnode_t* parent, double ** instRate );
int probMultMut(gnode_t * node, gnode_t* parent, double ** instRate, const gsl_rng * r );
int simulateBranchOneMut(gnode_t * node, gnode_t * parent, double  ** instRate, const gsl_rng * r) ;

void load_cfile(char * opt_cfile, char ** opt_msafile,
                char ** opt_outfile,
                char ** opt_gfile, char ** opt_mfile,
                char ** opt_subfile,
                long * opt_seed, long * locusSim, long * opt_site);

/* list */
void list_append(list_t * list, void * data);
void list_clear(list_t * list, void (*cb_dealloc)(void *));

/*summary */
void summary (mutHist_t ** mutationHist, int varsite, long samples, msa_t * msa, long opt_site);

