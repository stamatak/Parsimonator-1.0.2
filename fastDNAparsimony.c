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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with 
 *  thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */


#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>  
#endif

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>

#ifdef __SIM_SSE3

#include <xmmintrin.h>
#include <pmmintrin.h>
  
#endif

#ifdef __AVX

#include <xmmintrin.h>
#include <immintrin.h>

#endif


#include "axml.h"

extern char run_id[128];
extern char  seq_file[1024];
extern char  tree_file[1024];
extern char  resultFileName[1024];
extern char  infoFileName[1024]; 
extern char  randomFileName[1024];

const unsigned int mask32[32] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 
					262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 
					268435456, 536870912, 1073741824, 2147483648U};

/* vector-specific stuff */



#ifdef __SIM_SSE3

#define INTS_PER_VECTOR 4
#define INT_TYPE __m128i
#define CAST __m128i*
#define SET_ALL_BITS_ONE _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define VECTOR_LOAD _mm_load_si128
#define VECTOR_BIT_AND _mm_and_si128
#define VECTOR_BIT_OR  _mm_or_si128
#define VECTOR_STORE  _mm_store_si128
#define VECTOR_AND_NOT _mm_andnot_si128
#define BYTE_ALIGNMENT 16

#endif

#ifdef __AVX

#define INTS_PER_VECTOR 8
#define INT_TYPE __m256d
#define CAST double*
#define SET_ALL_BITS_ONE (__m256d)_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define VECTOR_LOAD _mm256_load_pd
#define VECTOR_BIT_AND _mm256_and_pd
#define VECTOR_BIT_OR  _mm256_or_pd
#define VECTOR_STORE  _mm256_store_pd
#define VECTOR_AND_NOT _mm256_andnot_pd
#define BYTE_ALIGNMENT 32

#endif

#if !(defined(__SIM_SSE3) || defined(__AVX))
#define BYTE_ALIGNMENT 16
#endif

static void *malloc_aligned(size_t size, size_t align) 
{
  void *ptr = (void *)NULL;  
  int res;
  

#if defined (__APPLE__)
  /* 
     presumably malloc on MACs always returns 
     a 16-byte aligned pointer
  */

  ptr = malloc(size);
  
  if(ptr == (void*)NULL) 
   assert(0);

#else
  res = posix_memalign( &ptr, align, size );

  if(res != 0) 
    assert(0);
#endif 
   
  return ptr;
}


/********************************DNA FUNCTIONS *****************************************************************/


static void hookupParsimony(nodeptr p, nodeptr q)
{
  p->back = q;
  q->back = p;
}

static void getxnodeLocal (nodeptr p)
{
  nodeptr  s;

  if((s = p->next)->x || (s = s->next)->x)
    {
      p->x = s->x;
      s->x = 0;
    }
}

static void computeTraversalInfoParsimony(nodeptr p, int *ti, int *counter, int maxTips, boolean full)
{        
  nodeptr 
    q = p->next->back,
    r = p->next->next->back;
  
  if(! p->x)
    getxnodeLocal(p);  
  
  if(full)
    {
       if(q->number > maxTips) 
	 computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  else
    {
      if(q->number > maxTips && !q->x) 
	computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips && !r->x) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  
  
  ti[*counter]     = p->number;
  ti[*counter + 1] = q->number;
  ti[*counter + 2] = r->number;
  *counter = *counter + 4;
}




#define BIT_COUNT(x)  precomputed16_bitcount(x)


#if (defined(__SIM_SSE3) || defined(__AVX))

static inline unsigned int populationCount(INT_TYPE v_N)
{
  unsigned int
    res[INTS_PER_VECTOR] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  
  unsigned int 
    i,
    a = 0;
    
  VECTOR_STORE((CAST)res, v_N);
    
  for(i = 0; i < INTS_PER_VECTOR; i++)
    a += BIT_COUNT(res[i]);
    
  return a;	   
}

#else

static inline unsigned int populationCount(unsigned int n)
{
  return BIT_COUNT(n);
}

#endif



#if (defined(__SIM_SSE3) || defined(__AVX))


static void newviewParsimonyIterativeFast(tree *tr)
{    
  INT_TYPE
    allOne = SET_ALL_BITS_ONE;
  
  int     
    *ti = tr->ti,
    count = ti[0],
    index;

  unsigned int  
    width = tr->compressedWidth;

 

  for(index = 4; index < count; index += 4)
    {        
      int 
	pNumber = ti[index],
	qNumber = ti[index + 1],
	rNumber = ti[index + 2];
      
      parsimonyNumber      
	*leftState_A  = tr->parsimonyState_A[qNumber],
	*rightState_A = tr->parsimonyState_A[rNumber],
	*thisState_A  = tr->parsimonyState_A[pNumber],
	
	*leftState_C  = tr->parsimonyState_C[qNumber],
	*rightState_C = tr->parsimonyState_C[rNumber],
	*thisState_C  = tr->parsimonyState_C[pNumber],

	*leftState_G  = tr->parsimonyState_G[qNumber],
	*rightState_G = tr->parsimonyState_G[rNumber],
	*thisState_G  = tr->parsimonyState_G[pNumber],

	*leftState_T  = tr->parsimonyState_T[qNumber],
	*rightState_T = tr->parsimonyState_T[rNumber],
	*thisState_T  = tr->parsimonyState_T[pNumber]; 	           
      
      
      
      unsigned int	
	i,            
	ts = 0;      
                 

#pragma omp parallel for reduction(+:ts)
       for(i = 0; i < width; i += INTS_PER_VECTOR)
	{	 	  
	  INT_TYPE
	    s_r, s_l, v_N,
	    l_A, l_C, l_G, l_T,
	    v_A, v_C, v_G, v_T;	    	 

	  s_l = VECTOR_LOAD((CAST)(&leftState_A[i]));
	  s_r = VECTOR_LOAD((CAST)(&rightState_A[i]));
	  l_A = VECTOR_BIT_AND(s_l, s_r);
	  v_A = VECTOR_BIT_OR(s_l, s_r);

	  s_l = VECTOR_LOAD((CAST)(&leftState_C[i]));
	  s_r = VECTOR_LOAD((CAST)(&rightState_C[i]));
	  l_C = VECTOR_BIT_AND(s_l, s_r);
	  v_C = VECTOR_BIT_OR(s_l, s_r);
	  
	  s_l = VECTOR_LOAD((CAST)(&leftState_G[i]));
	  s_r = VECTOR_LOAD((CAST)(&rightState_G[i]));
	  l_G = VECTOR_BIT_AND(s_l, s_r);
	  v_G = VECTOR_BIT_OR(s_l, s_r);

	  s_l = VECTOR_LOAD((CAST)(&leftState_T[i]));
	  s_r = VECTOR_LOAD((CAST)(&rightState_T[i]));
	  l_T = VECTOR_BIT_AND(s_l, s_r);
	  v_T = VECTOR_BIT_OR(s_l, s_r);
	  	   	    
	  v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));	  	 	    	  
	 
	  VECTOR_STORE((CAST)(&thisState_A[i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
	  VECTOR_STORE((CAST)(&thisState_C[i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));
	  VECTOR_STORE((CAST)(&thisState_G[i]), VECTOR_BIT_OR(l_G, VECTOR_AND_NOT(v_N, v_G)));
	  VECTOR_STORE((CAST)(&thisState_T[i]), VECTOR_BIT_OR(l_T, VECTOR_AND_NOT(v_N, v_T)));	  	 	 	  	  	  	
	  	 
	  v_N = VECTOR_AND_NOT(v_N, allOne);
	    
	  ts += populationCount(v_N); 
	}	     
	  
      tr->parsimonyScore[pNumber] = ts + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];          	 	    	      	  	
    }
}


static unsigned int evaluateParsimonyIterativeFast(tree *tr)
{
  INT_TYPE 
    allOne = SET_ALL_BITS_ONE;

  int   
    pNumber = tr->ti[1],
    qNumber = tr->ti[2];
  
  unsigned int 
    bestScore = tr->bestParsimony,
    width = tr->compressedWidth,
    i,   
    sum;

  parsimonyNumber
    *rightState_A = tr->parsimonyState_A[pNumber], 
    *leftState_A  = tr->parsimonyState_A[qNumber],
    
    *rightState_C = tr->parsimonyState_C[pNumber], 
    *leftState_C  = tr->parsimonyState_C[qNumber],

    *rightState_G = tr->parsimonyState_G[pNumber], 
    *leftState_G  = tr->parsimonyState_G[qNumber],

    *rightState_T = tr->parsimonyState_T[pNumber], 
    *leftState_T  = tr->parsimonyState_T[qNumber];
  
            
  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr);       
  
  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

#pragma omp parallel for reduction(+:sum)  
  for(i = 0; i < width; i += INTS_PER_VECTOR)
    {                                       
      INT_TYPE      
	l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&leftState_A[i])), VECTOR_LOAD((CAST)(&rightState_A[i]))),
	l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&leftState_C[i])), VECTOR_LOAD((CAST)(&rightState_C[i]))),
	l_G = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&leftState_G[i])), VECTOR_LOAD((CAST)(&rightState_G[i]))),
	l_T = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&leftState_T[i])), VECTOR_LOAD((CAST)(&rightState_T[i]))),
	v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));     

     
      v_N = VECTOR_AND_NOT(v_N, allOne);

      sum += populationCount(v_N);

#ifndef _OPENMP
      if(sum >= bestScore)
	break; 
#endif
    }   
  
  return sum;      
}


#else



static void newviewParsimonyIterativeFast(tree *tr)
{  
  int
    count = tr->ti[0],
    *ti   = tr->ti,
    index;

  unsigned int
    width = tr->compressedWidth;

  for(index = 4; index < count; index += 4)
    {        
      int 
	pNumber = ti[index],
	qNumber = ti[index + 1],
	rNumber = ti[index + 2];
      
      parsimonyNumber      
	*leftState_A  = tr->parsimonyState_A[qNumber],
	*rightState_A = tr->parsimonyState_A[rNumber],
	*thisState_A  = tr->parsimonyState_A[pNumber],
	
	*leftState_C  = tr->parsimonyState_C[qNumber],
	*rightState_C = tr->parsimonyState_C[rNumber],
	*thisState_C  = tr->parsimonyState_C[pNumber],

	*leftState_G  = tr->parsimonyState_G[qNumber],
	*rightState_G = tr->parsimonyState_G[rNumber],
	*thisState_G  = tr->parsimonyState_G[pNumber],

	*leftState_T  = tr->parsimonyState_T[qNumber],
	*rightState_T = tr->parsimonyState_T[rNumber],
	*thisState_T  = tr->parsimonyState_T[pNumber];
      
      register parsimonyNumber
	o_A,
	o_C,
	o_G,
	o_T,
	t_A,
	t_C,
	t_G,
	t_T,	
	t_N;
      
      unsigned int 
	i,	   
	ts = 0;             	         	  

      for(i = 0; i < width; i++)
	{		 	    
	  t_A = leftState_A[i] & rightState_A[i];
	  t_C = leftState_C[i] & rightState_C[i];
	  t_G = leftState_G[i] & rightState_G[i];	  
	  t_T = leftState_T[i] & rightState_T[i];

	  o_A = leftState_A[i] | rightState_A[i];
	  o_C = leftState_C[i] | rightState_C[i];
	  o_G = leftState_G[i] | rightState_G[i];	  
	  o_T = leftState_T[i] | rightState_T[i];

	  t_N = ~(t_A | t_C | t_G | t_T);	  

	  thisState_A[i] = t_A | (t_N & o_A);
	  thisState_C[i] = t_C | (t_N & o_C);
	  thisState_G[i] = t_G | (t_N & o_G);
	  thisState_T[i] = t_T | (t_N & o_T);	 	 	  
	  
	  ts += populationCount(t_N);   	   		      
	}	
	  
      tr->parsimonyScore[pNumber] = ts + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];          	 	    	      	  	
    }
}


static unsigned int evaluateParsimonyIterativeFast(tree *tr)
{
  int   
    pNumber = tr->ti[1],
    qNumber = tr->ti[2];  
  
  unsigned int 
    bestScore = tr->bestParsimony,
    width = tr->compressedWidth,
    i,   
    sum, 
    t_A,
    t_C,
    t_G,
    t_T,
    t_N;

  parsimonyNumber
    *rightState_A = tr->parsimonyState_A[pNumber], 
    *leftState_A  = tr->parsimonyState_A[qNumber],
    
    *rightState_C = tr->parsimonyState_C[pNumber], 
    *leftState_C  = tr->parsimonyState_C[qNumber],

    *rightState_G = tr->parsimonyState_G[pNumber], 
    *leftState_G  = tr->parsimonyState_G[qNumber],

    *rightState_T = tr->parsimonyState_T[pNumber], 
    *leftState_T  = tr->parsimonyState_T[qNumber];
            
  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr);       
  
  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];
  
  for(i = 0; i < width; i++)
    {                
       t_A = leftState_A[i] & rightState_A[i];
       t_C = leftState_C[i] & rightState_C[i];
       t_G = leftState_G[i] & rightState_G[i];	  
       t_T = leftState_T[i] & rightState_T[i];

       t_N = ~(t_A | t_C | t_G | t_T);

       sum += populationCount(t_N);     
       if(sum >= bestScore)
	 break;
    }         

  return sum;      
}

#endif






static unsigned int evaluateParsimony(tree *tr, nodeptr p, boolean full)
{
  volatile unsigned int result;
  nodeptr q = p->back;
  int
    *ti = tr->ti,
    counter = 4;
  
  ti[1] = p->number;
  ti[2] = q->number;

  if(full)
    {
      if(p->number > tr->mxtips)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }
  else
    {
      if(p->number > tr->mxtips && !p->x)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips && !q->x)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }
  
 

  ti[0] = counter;

  
  
  result = evaluateParsimonyIterativeFast(tr);

 

  return result;
}


static void newviewParsimony(tree *tr, nodeptr  p)
{     
  if(p->number <= tr->mxtips)
    return;

  {
    int 
      counter = 4;     
           
    computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
    tr->ti[0] = counter;            
    
    newviewParsimonyIterativeFast(tr);      
  }
}





/****************************************************************************************************************************************/

static void insertParsimony (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupParsimony(p->next,       q);
  hookupParsimony(p->next->next, r); 
   
  newviewParsimony(tr, p);     
} 

static nodeptr buildNewTip (tree *tr, nodeptr p)
{ 
  nodeptr  q;

  q = tr->nodep[(tr->nextnode)++];
  hookupParsimony(p, q);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
} 

static void buildSimpleTree (tree *tr, int ip, int iq, int ir)
{    
  nodeptr  p, s;
  int  i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  hookupParsimony(p, tr->nodep[iq]);
  s = buildNewTip(tr, tr->nodep[ir]);
  insertParsimony(tr, s, p);
}

static void testInsertParsimony (tree *tr, nodeptr p, nodeptr q)
{ 
  unsigned int mp;
 
  nodeptr  r = q->back;   
    
  insertParsimony(tr, p, q);   
  
  mp = evaluateParsimony(tr, p->next->next, FALSE);          
      
  if(mp < tr->bestParsimony)
    {
      tr->bestParsimony = mp;
      tr->insertNode = q;
      tr->removeNode = p;
    }
  
  hookupParsimony(q, r);
  p->next->next->back = p->next->back = (nodeptr) NULL;
       
  return;
} 


static void restoreTreeParsimony(tree *tr, nodeptr p, nodeptr q)
{ 
  nodeptr
    r = q->back;
  
  int counter = 4;
  
  hookupParsimony(p->next,       q);
  hookupParsimony(p->next->next, r);
  
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
  tr->ti[0] = counter;
    
  newviewParsimonyIterativeFast(tr); 
}


static void addTraverseParsimony (tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav, boolean doAll)
{        
  if (doAll || (--mintrav <= 0))               
    testInsertParsimony(tr, p, q);	                 

  if (((q->number > tr->mxtips)) && ((--maxtrav > 0) || doAll))
    {	      
      addTraverseParsimony(tr, p, q->next->back, mintrav, maxtrav, doAll);	      
      addTraverseParsimony(tr, p, q->next->next->back, mintrav, maxtrav, doAll);              	     
    }
}


static nodeptr findAnyTipFast(nodeptr p, int numsp)
{ 
  return  (p->number <= numsp)? p : findAnyTipFast(p->next->back, numsp);
} 

static double randum(long  *seed)
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}

static void makePermutationFast(int *perm, int n, analdef *adef)
{    
  int  i, j, k;


  assert(adef->parsimonySeed != 0);
       
  

  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {      
      double d =  randum(&adef->parsimonySeed);

      k =  (int)((double)(n + 1 - i) * d);
      
      j        = perm[i];

      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }
}

static nodeptr  removeNodeParsimony (nodeptr p)
{ 
  nodeptr  q, r;         

  q = p->next->back;
  r = p->next->next->back;   
    
  hookupParsimony(q, r);

  p->next->next->back = p->next->back = (node *) NULL;
  
  return  q;
}

static int rearrangeParsimony(tree *tr, nodeptr p, int mintrav, int maxtrav, boolean doAll)  
{   
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2; 
           
  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3; 

  
  if(maxtrav < mintrav)
    return 0;

  q = p->back;

  if(p->number > tr->mxtips) 
    {     
      p1 = p->next->back;
      p2 = p->next->next->back;
      
      if ((p1->number > tr->mxtips) || (p2->number > tr->mxtips)) 
	{	  	  
	  removeNodeParsimony(p);	  	 

	  if ((p1->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, p, p1->next->back, mintrav, maxtrav, doAll);         
	      addTraverseParsimony(tr, p, p1->next->next->back, mintrav, maxtrav, doAll);          
	    }
	 
	  if ((p2->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, p, p2->next->back, mintrav, maxtrav, doAll);
	      addTraverseParsimony(tr, p, p2->next->next->back, mintrav, maxtrav, doAll);          
	    }
	    
	   
	  hookupParsimony(p->next,       p1); 
	  hookupParsimony(p->next->next, p2);	   	    	    

	  newviewParsimony(tr, p);
	}
    }  
       
  if ((q->number > tr->mxtips) && maxtrav > 0) 
    {
      q1 = q->next->back;
      q2 = q->next->next->back;

      if (
	  (
	   (q1->number > tr->mxtips) && 
	   ((q1->next->back->number > tr->mxtips) || (q1->next->next->back->number > tr->mxtips))
	   )
	  ||
	  (
	   (q2->number > tr->mxtips) && 
	   ((q2->next->back->number > tr->mxtips) || (q2->next->next->back->number > tr->mxtips))
	   )
	  )
	{	   

	  removeNodeParsimony(q);
	  
	  mintrav2 = mintrav > 2 ? mintrav : 2;
	  
	  if ((q1->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, q, q1->next->back, mintrav2 , maxtrav, doAll);
	      addTraverseParsimony(tr, q, q1->next->next->back, mintrav2 , maxtrav, doAll);         
	    }
	 
	  if ((q2->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, q, q2->next->back, mintrav2 , maxtrav, doAll);
	      addTraverseParsimony(tr, q, q2->next->next->back, mintrav2 , maxtrav, doAll);          
	    }	   
	   
	  hookupParsimony(q->next,       q1); 
	  hookupParsimony(q->next->next, q2);
	   
	  newviewParsimony(tr, q);
	}
    }

  return 1;
} 


static void restoreTreeRearrangeParsimony(tree *tr)
{    
  removeNodeParsimony(tr->removeNode);  
  restoreTreeParsimony(tr, tr->removeNode, tr->insertNode);  
}

static boolean isInformative(tree *tr, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = 15;

  unsigned char
    nucleotide,
    target = 0;
  	
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {	   
      nucleotide = tr->yVector[j][site];	    
      check[nucleotide] =  check[nucleotide] + 1;      	           
    }
  
  
  if(check[1] > 1)
    {
      informativeCounter++;    
      target = target | 1;
    }
  if(check[2] > 1)
    {
      informativeCounter++; 
      target = target | 2;
    }
  if(check[4] > 1)
    {
      informativeCounter++; 
      target = target | 4;
    }
  if(check[8] > 1)
    {
      informativeCounter++; 
      target = target | 8;
    }
	  
  if(informativeCounter >= 2)
    return TRUE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
	{
	  if(j == 3 || j == 5 || j == 6 || j == 7 || j == 9 || j == 10 || j == 11 || 
	     j == 12 || j == 13 || j == 14)
	    {
	      if(check[j] > 1)
		{
		  if(!(target & j))
		    return TRUE;
		}
	    }
	} 
    }
     
  return FALSE;	     
}


static void determineUninformativeSites(tree *tr, int *informative)
{
  int 
    i,
    number = 0;

  /* 
     Not all characters are useful in constructing a parsimony tree. 
     Invariant characters, those that have the same state in all taxa, 
     are obviously useless and are ignored by the method. Characters in 
     which a state occurs in only one taxon are also ignored. 
     All these characters are called parsimony uninformative.

     Alternative definition: informative columns contain at least two types
     of nucleotides, and each nucleotide must appear at least twice in each 
     column. Kind of a pain if we intend to check for this when using, e.g.,
     amibiguous DNA encoding.
  */

  for(i = 0; i < tr->cdta->endsite; i++)
    {
      if(isInformative(tr, i)) 
	informative[i] = 1;
      else
	{
	  informative[i] = 0;
	  number++;
	}            
    }
  
 
  /* printf("Uninformative Patterns: %d\n", number); */
}


static void reorderNodes(tree *tr, nodeptr *np, nodeptr p, int *count)
{
  int i, found = 0;

  if((p->number <= tr->mxtips))    
    return;
  else
    {              
      for(i = tr->mxtips + 1; (i <= (tr->mxtips + tr->mxtips - 1)) && (found == 0); i++)
	{
	  if (p == np[i] || p == np[i]->next || p == np[i]->next->next)
	    {
	      if(p == np[i])			       
		tr->nodep[*count + tr->mxtips + 1] = np[i];		 		
	      else
		{
		  if(p == np[i]->next)		  
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next;		     	   
		  else		   
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next->next;		    		    
		}

	      found = 1;	      	     
	      *count = *count + 1;
	    }
	}            
     
      reorderNodes(tr, np, p->next->back, count);     
      reorderNodes(tr, np, p->next->next->back, count);                
    }
}

static void nodeRectifierFast(tree *tr)
{
  nodeptr *np = (nodeptr *)malloc(2 * (size_t)tr->mxtips * sizeof(nodeptr));
  int i;
  int count = 0;
  
  tr->start       = tr->nodep[1]; 
  
  for(i = tr->mxtips + 1; i <= (tr->mxtips + tr->mxtips - 1); i++)
    np[i] = tr->nodep[i];           
  
  reorderNodes(tr, np, tr->start->back, &count); 
 
  free(np);
}
  
static parsimonyNumber *compressDNA(tree *tr, int *informative)
{
  size_t
    i,
    entries = 0,   
    compressedEntries,
    compressedEntriesPadded;

  parsimonyNumber
    *compressedScratch;
  
  for(i = 0; i < (size_t)tr->cdta->endsite; i++)    
    if(informative[i])
      entries += 1;     
  
  compressedEntries = entries / PCF;

  if(entries % PCF != 0)
    compressedEntries++;
  
  
  /* printf("compression %d -> %d\n", entries, compressedEntries); */

#if (defined(__SIM_SSE3) || defined(__AVX))
  if(compressedEntries % INTS_PER_VECTOR != 0)
    compressedEntriesPadded = compressedEntries + (INTS_PER_VECTOR - (compressedEntries % INTS_PER_VECTOR));
  else
    compressedEntriesPadded = compressedEntries;
#else
  compressedEntriesPadded = compressedEntries;
#endif

  /* printf("padded %d\n", compressedEntriesPadded); */

  compressedScratch = (parsimonyNumber *)malloc_aligned((size_t)compressedEntriesPadded * 8 * (size_t)tr->mxtips * sizeof(parsimonyNumber), BYTE_ALIGNMENT);
     
  for(i = 0; i < compressedEntriesPadded * 8 * tr->mxtips; i++)      
    compressedScratch[i] = 0;    
      

  for(i = 0; i < (size_t)tr->mxtips; i++)
    {
      parsimonyNumber
	*compressedTip_A = &compressedScratch[(compressedEntriesPadded * 4 * (i + 1))],
	*compressedTip_C = &compressedScratch[(compressedEntriesPadded * 4 * (i + 1)) + compressedEntriesPadded * 1],
	*compressedTip_G = &compressedScratch[(compressedEntriesPadded * 4 * (i + 1)) + compressedEntriesPadded * 2],
	*compressedTip_T = &compressedScratch[(compressedEntriesPadded * 4 * (i + 1)) + compressedEntriesPadded * 3],
	compressedValue_A = 0,
	compressedValue_C = 0,
	compressedValue_G = 0,
	compressedValue_T = 0;
      
      size_t
	w = 0,
	compressedIndex = 0,
	compressedCounter = 0,
	index = 0;
	      
      for(index = 0; index < (size_t)tr->cdta->endsite; index++)
	{
	  if(informative[index])
	    {
	      parsimonyNumber value = (parsimonyNumber)(tr->yVector[i + 1][index]);	  
	      
	      for(w = 0; w < 1; w++)
		{	    	     	    
		  if(value & 1)
		    compressedValue_A |= mask32[compressedCounter];
		  if(value & 2)
		    compressedValue_C |= mask32[compressedCounter];
		  if(value & 4)
		    compressedValue_G |= mask32[compressedCounter];
		  if(value & 8)
		    compressedValue_T |= mask32[compressedCounter];
	          
		  compressedCounter++;
		  
		  if(compressedCounter == PCF)
		    {
		      compressedTip_A[compressedIndex] = compressedValue_A;
		      compressedTip_C[compressedIndex] = compressedValue_C;
		      compressedTip_G[compressedIndex] = compressedValue_G;
		      compressedTip_T[compressedIndex] = compressedValue_T;
		      
		      compressedValue_A = 0;
		      compressedValue_C = 0;
		      compressedValue_G = 0;
		      compressedValue_T = 0;
		      
		      compressedCounter = 0;
		      compressedIndex++;
		    }
		}
	    }
	}
                           
      for(;compressedIndex < compressedEntriesPadded; compressedIndex++)
	{	
	  for(;compressedCounter < PCF; compressedCounter++)
	    {	        	      	      
	      compressedValue_A |= mask32[compressedCounter];		 
	      compressedValue_C |= mask32[compressedCounter];		 
	      compressedValue_G |= mask32[compressedCounter];		  
	      compressedValue_T |= mask32[compressedCounter];	      	     
	    }
	  
	   compressedTip_A[compressedIndex] = compressedValue_A;
	   compressedTip_C[compressedIndex] = compressedValue_C;
	   compressedTip_G[compressedIndex] = compressedValue_G;
	   compressedTip_T[compressedIndex] = compressedValue_T;
	   
	   compressedValue_A = 0;
	   compressedValue_C = 0;
	   compressedValue_G = 0;
	   compressedValue_T = 0;
	   
	   compressedCounter = 0;
	}	 	
    }
    
  tr->parsimonyState_A = (parsimonyNumber**)malloc_aligned(sizeof(parsimonyNumber*) * 2 * (size_t)tr->mxtips, BYTE_ALIGNMENT);
  tr->parsimonyState_C = (parsimonyNumber**)malloc_aligned(sizeof(parsimonyNumber*) * 2 * (size_t)tr->mxtips, BYTE_ALIGNMENT);
  tr->parsimonyState_G = (parsimonyNumber**)malloc_aligned(sizeof(parsimonyNumber*) * 2 * (size_t)tr->mxtips, BYTE_ALIGNMENT);
  tr->parsimonyState_T = (parsimonyNumber**)malloc_aligned(sizeof(parsimonyNumber*) * 2 * (size_t)tr->mxtips, BYTE_ALIGNMENT);
  
  tr->parsimonyScore = (unsigned int*)malloc_aligned(sizeof(unsigned int) * 2 * (size_t)tr->mxtips, BYTE_ALIGNMENT);  
    
  for(i = 0; i < 2 * (size_t)tr->mxtips; i++) 
    {
      tr->parsimonyState_A[i] = &compressedScratch[i * 4 * compressedEntriesPadded];       
      tr->parsimonyState_C[i] = &compressedScratch[i * 4 * compressedEntriesPadded + compressedEntriesPadded * 1]; 
      tr->parsimonyState_G[i] = &compressedScratch[i * 4 * compressedEntriesPadded + compressedEntriesPadded * 2]; 
      tr->parsimonyState_T[i] = &compressedScratch[i * 4 * compressedEntriesPadded + compressedEntriesPadded * 3];          
 
      tr->parsimonyScore[i] = 0;
    }
  
  tr->compressedWidth = (unsigned int)compressedEntriesPadded;
  
  return compressedScratch; 
}



static void stepwiseAddition(tree *tr, nodeptr p, nodeptr q)
{            
  nodeptr 
    r = q->back;

  unsigned int 
    mp,
    bestParsimony = tr->bestParsimony;
  
  int 
    counter = 4;
  
  p->next->back = q;
  q->back = p->next;

  p->next->next->back = r;
  r->back = p->next->next;
   
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
  tr->ti[0] = counter;
    
  newviewParsimonyIterativeFast(tr);      
  
#if (defined(__SIM_SSE3) || defined(__AVX))
  {
    INT_TYPE 
      allOne = SET_ALL_BITS_ONE;
    
    int   
      pNumber = p->number,
      qNumber = p->back->number;
  
    unsigned int       
      width = tr->compressedWidth,
      i;
    
    parsimonyNumber
      *rightState_A = tr->parsimonyState_A[pNumber], 
      *leftState_A  = tr->parsimonyState_A[qNumber],
      
      *rightState_C = tr->parsimonyState_C[pNumber], 
      *leftState_C  = tr->parsimonyState_C[qNumber],
      
      *rightState_G = tr->parsimonyState_G[pNumber], 
      *leftState_G  = tr->parsimonyState_G[qNumber],
      
      *rightState_T = tr->parsimonyState_T[pNumber], 
      *leftState_T  = tr->parsimonyState_T[qNumber];
        
    mp = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

#pragma omp parallel for reduction(+:mp)
    for(i = 0; i < width; i += INTS_PER_VECTOR)
      {     	
	INT_TYPE      
	  l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&leftState_A[i])), VECTOR_LOAD((CAST)(&rightState_A[i]))),
	  l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&leftState_C[i])), VECTOR_LOAD((CAST)(&rightState_C[i]))),
	  l_G = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&leftState_G[i])), VECTOR_LOAD((CAST)(&rightState_G[i]))),
	  l_T = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&leftState_T[i])), VECTOR_LOAD((CAST)(&rightState_T[i]))),
	  v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));     		

	v_N = VECTOR_AND_NOT(v_N, allOne);

	mp += populationCount(v_N);

#ifndef _OPENMP
	if(mp >= bestParsimony)
	  goto SKIP;	
#endif

      }           
  }
#else
  {
    int   
      pNumber = p->number,
      qNumber = p->back->number;  
    
    unsigned int 
      width = tr->compressedWidth,
      i,         
      t_A,
      t_C,
      t_G,
      t_T,
      t_N;
    
    parsimonyNumber
      *rightState_A = tr->parsimonyState_A[pNumber], 
      *leftState_A  = tr->parsimonyState_A[qNumber],
      
      *rightState_C = tr->parsimonyState_C[pNumber], 
      *leftState_C  = tr->parsimonyState_C[qNumber],
      
      *rightState_G = tr->parsimonyState_G[pNumber], 
      *leftState_G  = tr->parsimonyState_G[qNumber],
      
      *rightState_T = tr->parsimonyState_T[pNumber], 
      *leftState_T  = tr->parsimonyState_T[qNumber];
    
    mp = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];
    
    for(i = 0; i < width; i++)
      {                
	t_A = leftState_A[i] & rightState_A[i];
	t_C = leftState_C[i] & rightState_C[i];
	t_G = leftState_G[i] & rightState_G[i];	  
	t_T = leftState_T[i] & rightState_T[i];
	
	t_N = ~(t_A | t_C | t_G | t_T);
	
	mp += populationCount(t_N);
	
	if(mp >= bestParsimony)
	  goto SKIP;
      }            
  }
#endif


#ifdef _OPENMP
  if(mp < tr->bestParsimony)
    {
#endif      
      tr->bestParsimony = mp;
      tr->insertNode = q;     
#ifdef _OPENMP
    }
#endif
  
 SKIP: 
  q->back = r;
  r->back = q;
   
  if(q->number > tr->mxtips && tr->parsimonyScore[q->number] > 0)
    {	      
      stepwiseAddition(tr, p, q->next->back);	      
      stepwiseAddition(tr, p, q->next->next->back);              	     
    }
}


static char *Tree2StringREC(char *treestr, tree *tr, nodeptr p)
{     
  if(p->number > 0 && p->number <= tr->mxtips) 
    {
      sprintf(treestr, "%s", tr->nameList[p->number]);   	
      while (*treestr) treestr++;
    }
  else 
    {                 	 
      *treestr++ = '(';
      treestr = Tree2StringREC(treestr, tr, p->next->back);
      *treestr++ = ',';
      treestr = Tree2StringREC(treestr, tr, p->next->next->back);
      if(p == tr->start->back) 
	{
	  *treestr++ = ',';
	  treestr = Tree2StringREC(treestr, tr, p->back);
	}
      *treestr++ = ')';                    
    }

  if(p == tr->start->back)     
    sprintf(treestr, ";\n");	 	  	   
  else     
    sprintf(treestr, "%s", "\0");	         
  
  while (*treestr) treestr++;
  return  treestr;
}





void makeParsimonyTreeFastDNA(tree *tr, analdef *adef)
{   
  nodeptr  
    p, 
    f;    
  
  int 
    treeCounter = 0,
    i, 
    nextsp,
    *perm        = (int *)malloc((size_t)(tr->mxtips + 1) * sizeof(int)),
    *informative = (int *)malloc(sizeof(int) * (size_t)tr->cdta->endsite);  

  unsigned int 
    randomMP, 
    startMP;        

  parsimonyNumber
    *parsimonySpace;

  determineUninformativeSites(tr, informative);     

  parsimonySpace = compressDNA(tr, informative);

  free(informative); 

  tr->ti = (int*)malloc(sizeof(int) * 4 * (size_t)tr->mxtips);  

  if(adef->restart)
    tr->nodesInTree = (int*)calloc((size_t)(tr->mxtips + 1), sizeof(int));

  while(treeCounter < adef->numberOfTrees)
    {
      double t = gettime();

      if(adef->restart)
	{	  
	  int j = 0;

	  FILE *treeFile = fopen(tree_file, "rb");	 

	  treeReadLen(treeFile, tr);

	  printBothOpen("Read Starting tree with %d tips, total %d\n", tr->ntips, tr->mxtips);

	  fclose(treeFile);
	  
	  tr->start = findAnyTipFast(tr->start, tr->rdta->numsp);

	  tr->bestParsimony = INT_MAX;

	  evaluateParsimony(tr, tr->start->back, TRUE);


	  assert(tr->start);
	  assert(adef->parsimonySeed != 0);	 
	  
	  j = tr->ntips + 1;
	  
	  for(i = 1; i <= tr->mxtips; i++)      
	    if(tr->nodesInTree[i] == 0) 
	      perm[j++] = i;	  

	  for(i = tr->ntips + 1; i <= tr->mxtips; i++) 
	    {	     
	      int k, j;
	    
	      k =  (int)((double)(tr->mxtips + 1 - i) * randum(&adef->parsimonySeed));

	      assert(i + k <= tr->mxtips);
	      j        = perm[i];
	      perm[i]     = perm[i + k];
	      perm[i + k] = j;
	    }    
	  
	  f = tr->start;
	}
      else
	{
	  makePermutationFast(perm, tr->mxtips, adef);
  
	  tr->ntips = 0;    
      
	  tr->nextnode = tr->mxtips + 1;       
 
	  buildSimpleTree(tr, perm[1], perm[2], perm[3]);      

	  f = tr->start;
	}     

      while(tr->ntips < tr->mxtips) 
	{	
	  nodeptr q;
	  
	  tr->bestParsimony = INT_MAX;
	  nextsp = ++(tr->ntips);             
	  p = tr->nodep[perm[nextsp]];                 
	  q = tr->nodep[(tr->nextnode)++];
	  p->back = q;
	  q->back = p;
	  
	  
	  stepwiseAddition(tr, q, f->back);      	  	 
	  
	  {
	    nodeptr	  
	      r = tr->insertNode->back;
	    
	    int counter = 4;
	    
	    hookupParsimony(q->next,       tr->insertNode);
	    hookupParsimony(q->next->next, r);
	    
	    computeTraversalInfoParsimony(q, tr->ti, &counter, tr->mxtips, FALSE);              
	    tr->ti[0] = counter;
	    
	    newviewParsimonyIterativeFast(tr);	
	  }
	}    
      
      
      /*printf("ADD: %d\n", tr->bestParsimony); */
    
      
      nodeRectifierFast(tr);
 
      randomMP = tr->bestParsimony;        
  
      do
	{
	  startMP = randomMP;
	  nodeRectifierFast(tr);
	  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	    {
	      rearrangeParsimony(tr, tr->nodep[i], 1, 20, FALSE);
	      if(tr->bestParsimony < randomMP)
		{		
		  restoreTreeRearrangeParsimony(tr);
		  randomMP = tr->bestParsimony;
		}
	    }      		  	   
	}
      while(randomMP < startMP);
      
      
      /* printf("OPT: %d\n", tr->bestParsimony); */
      
      {
	FILE 
	  *f;
	
	char 
	  buf[64],
	  fileName[2048];

	sprintf(buf, "%d", treeCounter);

	strcpy(fileName, resultFileName);
	strcat(fileName, ".");
	strcat(fileName, buf);

	f = myfopen(fileName, "wb");
	
	Tree2StringREC(tr->tree_string, tr, tr->start->back);
	
	fprintf(f, "%s", tr->tree_string);
	
	fclose(f);
	
	printBothOpen("Parsimony tree [%d] with length %u computed in %f seconds written to file: %s\n", treeCounter, tr->bestParsimony, (gettime() - t), fileName);
      }

      treeCounter++;
    }

  free(parsimonySpace);
} 



 
