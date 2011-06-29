/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses
 *  with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#include <assert.h>
#include <stdint.h>



#define nmlngth        256         /* number of characters in species name */



#define NUM_BRANCHES 10


#define TRUE             1
#define FALSE            0


#define MIN(x,y)  (((x)<(y)) ?    (x)  : (y))

#define programName        "Parsimonator"
#define programVersion     "1.0.2"
#define programDate        "June 2011"




#define M_GTRCAT         1




#define TIP_TIP     0
#define TIP_INNER   1
#define INNER_INNER 2


#define DNA_DATA         1










typedef  int boolean;




typedef unsigned int parsimonyNumber;



#define PCF 32






typedef struct ratec
{
  double accumulatedSiteLikelihood;
  double rate;
}
  rateCategorize;




struct noderec;













typedef  struct noderec
{
  
  struct noderec  *next;
  struct noderec  *back;
  int              number;
  char             x;
}
  node, *nodeptr;





typedef  struct
{
  int              numsp;
  int              sites;
  unsigned char             **y;
  unsigned char             *y0;
  unsigned char             *yBUF;
  int              *wgt;
} rawdata;

typedef  struct {
  int             *alias;       /* site representing a pattern */
  int             *aliaswgt;    /* weight by pattern */
  int             *rateCategory;
  int              endsite;     /* # of sequence patterns */
  double          *patrat;      /* rates per pattern */
  double          *patratStored; 
} cruncheddata;








typedef struct 
{
  int left;
  int right;
  double likelihood;
} lhEntry;


typedef struct 
{
  int count;
  int size;
  lhEntry *entries;
} lhList;


typedef struct List_{
  void *value; 			
  struct List_ *next; 
} List;

struct stringEnt
{
  int nodeNumber;
  char *word;
  struct stringEnt *next;
};

typedef struct stringEnt stringEntry;

typedef unsigned int hashNumberType;

typedef struct
{
  hashNumberType tableSize;
  stringEntry **table;
}
  stringHashtable;

typedef  struct  {
 
 

  parsimonyNumber **parsimonyState_A;
  parsimonyNumber **parsimonyState_C;
  parsimonyNumber **parsimonyState_G;
  parsimonyNumber **parsimonyState_T;
  unsigned int *parsimonyScore; 
  int *ti;
  unsigned int compressedWidth;

  unsigned char             **yVector;
  

 
 
  stringHashtable  *nameHash;

 
  node           **nodep;
  node            *start;
  int              mxtips;
   
  int              ntips;
  int              nextnode;
  
 
  rawdata         *rdta;
  cruncheddata    *cdta;

  char **nameList;
  char *tree_string;
  int *nodesInTree;

  int treeStringLength;
  unsigned int bestParsimony;
  
  nodeptr removeNode;
  nodeptr insertNode;

} tree;


/***************************************************************/



/**************************************************************/






typedef  struct {


  long           parsimonySeed;
  boolean        restart;
  int            numberOfTrees;
} analdef;


/****************************** FUNCTIONS ****************************************************/


extern void makeParsimonyTreeFastDNA(tree *tr, analdef *adef);
extern void printBothOpen(const char* format, ... );
extern double gettime(void);
extern unsigned int precomputed16_bitcount (unsigned int n);
extern FILE *myfopen(const char *path, const char *mode);
extern void treeReadLen (FILE *fp, tree *tr);
