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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands 
 *  of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef WIN32
#include <unistd.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "axml.h"

static const double MNBRAK_GOLD =    1.618034;
static const double MNBRAK_TINY =      1.e-20;
static const double MNBRAK_GLIMIT =     100.0;
static const double BRENT_ZEPS  =      1.e-5;
static const double BRENT_CGOLD =   0.3819660;

double assignmentEvalTicks = 0;

extern int optimizeRatesInvocations;
extern int optimizeRateCategoryInvocations;
extern int optimizeAlphaInvocations;
extern int optimizeInvarInvocations;
extern double masterTime;
extern char ratesFileName[1024];
extern char workdir[1024];
extern char run_id[128];
extern char lengthFileName[1024];
extern char lengthFileNameModel[1024];
extern char *protModels[20];
extern int protEmpiricalFreqs;
extern int assertionError;


#ifdef _USE_PTHREADS
extern volatile int             NumberOfThreads;
extern volatile double          *reductionBuffer;
// [JH]
#else
volatile int NumberOfThreads=1;
#endif



/*********************FUNCTIONS FOOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/


static void setRateModel(tree *tr, int model, double rate, int position)
{
  int
    states   = tr->partitionData[model].states,
    numRates = (states * states - states) / 2;

  if(tr->partitionData[model].dataType == DNA_DATA)
    assert(position >= 0 && position < (numRates - 1));
  else
    assert(position >= 0 && position < numRates);

  assert(tr->partitionData[model].dataType != BINARY_DATA); 

  assert(rate >= RATE_MIN && rate <= RATE_MAX);

  if(tr->partitionData[model].nonGTR)
    {    
      int 
	i, 
	k = tr->partitionData[model].symmetryVector[position];

      assert(tr->partitionData[model].dataType == SECONDARY_DATA ||
	     tr->partitionData[model].dataType == SECONDARY_DATA_6 ||
	     tr->partitionData[model].dataType == SECONDARY_DATA_7);

      if(k == -1)
	tr->partitionData[model].substRates[position] = 0.0;
      else
	{
	  if(k == tr->partitionData[model].symmetryVector[numRates - 1])
	    {
	      for(i = 0; i < numRates - 1; i++)
		if(tr->partitionData[model].symmetryVector[i] == k)
		  tr->partitionData[model].substRates[position] = 1.0;
	    }
	  else
	    {
	      for(i = 0; i < numRates - 1; i++)
		{
		  if(tr->partitionData[model].symmetryVector[i] == k)
		    tr->partitionData[model].substRates[i] = rate; 
		}	      	     
	    }
	}
    }
  else
    tr->partitionData[model].substRates[position] = rate;
}





static linkageList* initLinkageList(int *linkList, tree *tr)
{
  int 
    k,
    partitions,
    numberOfModels = 0,
    i,
    pos;
  linkageList* ll = (linkageList*)malloc(sizeof(linkageList));
      
  for(i = 0; i < tr->NumberOfModels; i++)
    {
      if(linkList[i] > numberOfModels)
	numberOfModels = linkList[i];
    }

  numberOfModels++;
  
  ll->entries = numberOfModels;
  ll->ld      = (linkageData*)malloc(sizeof(linkageData) * numberOfModels);


  for(i = 0; i < numberOfModels; i++)
    {
      ll->ld[i].valid = TRUE;
      partitions = 0;

      for(k = 0; k < tr->NumberOfModels; k++)	
	if(linkList[k] == i)
	  partitions++;	    

      ll->ld[i].partitions = partitions;
      ll->ld[i].partitionList = (int*)malloc(sizeof(int) * partitions);
      
      for(k = 0, pos = 0; k < tr->NumberOfModels; k++)	
	if(linkList[k] == i)
	  ll->ld[i].partitionList[pos++] = k;
    }

  return ll;
}


static linkageList* initLinkageListGTR(tree *tr)
{
  int
    i,
    *links = (int*)malloc(sizeof(int) * tr->NumberOfModels),
    firstAA = tr->NumberOfModels + 2,
    countGTR = 0,
    countOtherModel = 0;
  linkageList* ll;

  for(i = 0; i < tr->NumberOfModels; i++)
    {     
      if(tr->partitionData[i].dataType == AA_DATA)
	{
	  if(tr->partitionData[i].protModels == GTR)
	    {
	      if(i < firstAA)
		firstAA = i;
	      countGTR++;
	    }
	  else
	    countOtherModel++;
	}
    }
  
  assert((countGTR > 0 && countOtherModel == 0) || (countGTR == 0 && countOtherModel > 0) ||  (countGTR == 0 && countOtherModel == 0));

  if(countGTR == 0)
    {
      for(i = 0; i < tr->NumberOfModels; i++)
	links[i] = i;
    }
  else
    {
      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  switch(tr->partitionData[i].dataType)
	    {	   
	    case DNA_DATA:
	    case BINARY_DATA:
	    case GENERIC_32:
	    case GENERIC_64:
	    case SECONDARY_DATA:
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7: 
	      links[i] = i;
	      break;
	    case AA_DATA:	  
	      links[i] = firstAA;
	      break;
	    default:
	      assert(0);
	    }
	}
    }
  

  ll = initLinkageList(links, tr);

  free(links);
  
  return ll;
}



static void freeLinkageList( linkageList* ll)
{
  int i;    

  for(i = 0; i < ll->entries; i++)    
    free(ll->ld[i].partitionList);         

  free(ll->ld);
  free(ll);   
}

#define ALPHA_F 0
#define INVAR_F 1
#define RATE_F  2




static void evaluateChange(tree *tr, int rateNumber, double *value, double *result, boolean* converged, analdef *adef, int whichFunction, int numberOfModels, linkageList *ll)
{ 
  int i, k, pos;

  switch(whichFunction)
    {    
    case RATE_F:
      for(i = 0, pos = 0; i < ll->entries; i++)
	{
	  if(ll->ld[i].valid)
	    {
	      if(converged[pos])
		{
		  for(k = 0; k < ll->ld[i].partitions; k++)
		    tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;
		}
	      else
		{
		  for(k = 0; k < ll->ld[i].partitions; k++)
		    {
		      int index = ll->ld[i].partitionList[k];		  	      
		      setRateModel(tr, index, value[pos], rateNumber);  
		      initReversibleGTR(tr, adef, index);		 
		    }
		}
	      pos++;
	    }
	  else
	    {
	      for(k = 0; k < ll->ld[i].partitions; k++)
		tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;	     
	    }
	 
	}

      assert(pos == numberOfModels);

#ifdef _USE_PTHREADS
      {
	volatile double result;
	
	masterBarrier(THREAD_OPT_RATE, tr);
	if(tr->NumberOfModels == 1)
	  {
	    for(i = 0, result = 0.0; i < NumberOfThreads; i++)    	  
	      result += reductionBuffer[i];  	        
	    tr->perPartitionLH[0] = result;
	  }
	else
	  {
	    int j;
	    volatile double partitionResult;
	
	    result = 0.0;

	    for(j = 0; j < tr->NumberOfModels; j++)
	      {
		for(i = 0, partitionResult = 0.0; i < NumberOfThreads; i++)          	      
		  partitionResult += reductionBuffer[i * tr->NumberOfModels + j];
		result +=  partitionResult;
		tr->perPartitionLH[j] = partitionResult;
	      }
	  }
      }
#else
#ifdef _FINE_GRAIN_MPI
      masterBarrierMPI(THREAD_OPT_RATE, tr);     
#else
      evaluateGenericInitrav(tr, tr->start);      
#endif
#endif     
      
      for(i = 0, pos = 0; i < ll->entries; i++)	
	{
	  if(ll->ld[i].valid)
	    {
	      result[pos] = 0.0;
	      for(k = 0; k < ll->ld[i].partitions; k++)
		{
		  int index = ll->ld[i].partitionList[k];

		  assert(tr->perPartitionLH[index] <= 0.0);

		  result[pos] -= tr->perPartitionLH[index];
		  
		}
	      pos++;
	    }
	   for(k = 0; k < ll->ld[i].partitions; k++)
	     {
	       int index = ll->ld[i].partitionList[k];
	       tr->executeModel[index] = TRUE;
	     }	  
	}

      assert(pos == numberOfModels);
      break;
    case ALPHA_F:
      for(i = 0; i < ll->entries; i++)
	{
	  if(converged[i])
	    {
	      for(k = 0; k < ll->ld[i].partitions; k++)
		tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;
	    }
	  else
	    {
	      for(k = 0; k < ll->ld[i].partitions; k++)
		{
		  int index = ll->ld[i].partitionList[k];
		  tr->executeModel[index] = TRUE;
		  tr->partitionData[index].alpha = value[i];
		  makeGammaCats(tr->partitionData[index].alpha, tr->partitionData[index].gammaRates, 4);
		}
	    }
	}
#ifdef _USE_PTHREADS   
      {
	volatile double result;
	
	masterBarrier(THREAD_OPT_ALPHA, tr);
	if(tr->NumberOfModels == 1)
	  {
	    for(i = 0, result = 0.0; i < NumberOfThreads; i++)    	  
	      result += reductionBuffer[i];  	        
	    tr->perPartitionLH[0] = result;
	  }
	else
	  {
	    int j;
	    volatile double partitionResult;
	
	    result = 0.0;

	    for(j = 0; j < tr->NumberOfModels; j++)
	      {
		for(i = 0, partitionResult = 0.0; i < NumberOfThreads; i++)          	      
		  partitionResult += reductionBuffer[i * tr->NumberOfModels + j];
		result +=  partitionResult;
		tr->perPartitionLH[j] = partitionResult;
	      }
	  }
      }
#else
#ifdef _FINE_GRAIN_MPI
      masterBarrierMPI(THREAD_OPT_ALPHA, tr);     
#else
      evaluateGenericInitrav(tr, tr->start);
#endif
#endif
            
      for(i = 0; i < ll->entries; i++)	
	{	  
	  result[i] = 0.0;
	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    {
	      int index = ll->ld[i].partitionList[k];
	      	      
	      assert(tr->perPartitionLH[index] <= 0.0);		
	      
	      result[i] -= tr->perPartitionLH[index];	            
	      tr->executeModel[index] = TRUE;
	    }
	}
      break;
    default:
      assert(0);	
    }

}



static void brentGeneric(double *ax, double *bx, double *cx, double *fb, double tol, double *xmin, double *result, int numberOfModels, 
			 int whichFunction, int rateNumber, analdef *adef, tree *tr, linkageList *ll, double lim_inf, double lim_sup)
{
  int iter, i;
  double 
    *a     = (double *)malloc(sizeof(double) * numberOfModels),
    *b     = (double *)malloc(sizeof(double) * numberOfModels),
    *d     = (double *)malloc(sizeof(double) * numberOfModels),
    *etemp = (double *)malloc(sizeof(double) * numberOfModels),
    *fu    = (double *)malloc(sizeof(double) * numberOfModels),
    *fv    = (double *)malloc(sizeof(double) * numberOfModels),
    *fw    = (double *)malloc(sizeof(double) * numberOfModels),
    *fx    = (double *)malloc(sizeof(double) * numberOfModels),
    *p     = (double *)malloc(sizeof(double) * numberOfModels),
    *q     = (double *)malloc(sizeof(double) * numberOfModels),
    *r     = (double *)malloc(sizeof(double) * numberOfModels),
    *tol1  = (double *)malloc(sizeof(double) * numberOfModels),
    *tol2  = (double *)malloc(sizeof(double) * numberOfModels),
    *u     = (double *)malloc(sizeof(double) * numberOfModels),
    *v     = (double *)malloc(sizeof(double) * numberOfModels),
    *w     = (double *)malloc(sizeof(double) * numberOfModels),
    *x     = (double *)malloc(sizeof(double) * numberOfModels),
    *xm    = (double *)malloc(sizeof(double) * numberOfModels),
    *e     = (double *)malloc(sizeof(double) * numberOfModels);
  boolean *converged = (boolean *)malloc(sizeof(boolean) * numberOfModels);
  boolean allConverged;
  
  for(i = 0; i < numberOfModels; i++)    
    converged[i] = FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      e[i] = 0.0;
      d[i] = 0.0;
    }

  for(i = 0; i < numberOfModels; i++)
    {
      a[i]=((ax[i] < cx[i]) ? ax[i] : cx[i]);
      b[i]=((ax[i] > cx[i]) ? ax[i] : cx[i]);
      x[i] = w[i] = v[i] = bx[i];
      fw[i] = fv[i] = fx[i] = fb[i];
    }

  for(i = 0; i < numberOfModels; i++)
    {      
      assert(a[i] >= lim_inf && a[i] <= lim_sup);
      assert(b[i] >= lim_inf && b[i] <= lim_sup);
      assert(x[i] >= lim_inf && x[i] <= lim_sup);
      assert(v[i] >= lim_inf && v[i] <= lim_sup);
      assert(w[i] >= lim_inf && w[i] <= lim_sup);
    }
  
  

  for(iter = 1; iter <= ITMAX; iter++)
    {
      allConverged = TRUE;

      for(i = 0; i < numberOfModels && allConverged; i++)
	allConverged = allConverged && converged[i];

      if(allConverged)
	{
	  free(converged);
	  free(a);
	  free(b);
	  free(d);
	  free(etemp);
	  free(fu);
	  free(fv);
	  free(fw);
	  free(fx);
	  free(p);
	  free(q);
	  free(r);
	  free(tol1);
	  free(tol2);
	  free(u);
	  free(v);
	  free(w);
	  free(x);
	  free(xm);
	  free(e);
	  return;
	}     

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {	     	      
	      assert(a[i] >= lim_inf && a[i] <= lim_sup);
	      assert(b[i] >= lim_inf && b[i] <= lim_sup);
	      assert(x[i] >= lim_inf && x[i] <= lim_sup);
	      assert(v[i] >= lim_inf && v[i] <= lim_sup);
	      assert(w[i] >= lim_inf && w[i] <= lim_sup);
  
	      xm[i] = 0.5 * (a[i] + b[i]);
	      tol2[i] = 2.0 * (tol1[i] = tol * fabs(x[i]) + BRENT_ZEPS);
	  
	      if(fabs(x[i] - xm[i]) <= (tol2[i] - 0.5 * (b[i] - a[i])))
		{		 
		  result[i] =  -fx[i];
		  xmin[i]   = x[i];
		  converged[i] = TRUE;		  
		}
	      else
		{
		  if(fabs(e[i]) > tol1[i])
		    {		     
		      r[i] = (x[i] - w[i]) * (fx[i] - fv[i]);
		      q[i] = (x[i] - v[i]) * (fx[i] - fw[i]);
		      p[i] = (x[i] - v[i]) * q[i] - (x[i] - w[i]) * r[i];
		      q[i] = 2.0 * (q[i] - r[i]);
		      if(q[i] > 0.0)
			p[i] = -p[i];
		      q[i] = fabs(q[i]);
		      etemp[i] = e[i];
		      e[i] = d[i];
		      if((fabs(p[i]) >= fabs(0.5 * q[i] * etemp[i])) || (p[i] <= q[i] * (a[i]-x[i])) || (p[i] >= q[i] * (b[i] - x[i])))
			d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i] : b[i] - x[i]));
		      else
			{
			  d[i] = p[i] / q[i];
			  u[i] = x[i] + d[i];
			  if( u[i] - a[i] < tol2[i] || b[i] - u[i] < tol2[i])
			    d[i] = SIGN(tol1[i], xm[i] - x[i]);
			}
		    }
		  else
		    {		     
		      d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i]: b[i] - x[i]));
		    }
		  u[i] = ((fabs(d[i]) >= tol1[i]) ? (x[i] + d[i]): (x[i] +SIGN(tol1[i], d[i])));
		}

	      if(!converged[i])
		assert(u[i] >= lim_inf && u[i] <= lim_sup);
	    }
	}
                 
      evaluateChange(tr, rateNumber, u, fu, converged, adef, whichFunction, numberOfModels, ll);

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {
	      if(fu[i] <= fx[i])
		{
		  if(u[i] >= x[i])
		    a[i] = x[i];
		  else
		    b[i] = x[i];
		  
		  SHFT(v[i],w[i],x[i],u[i]);
		  SHFT(fv[i],fw[i],fx[i],fu[i]);
		}
	      else
		{
		  if(u[i] < x[i])
		    a[i] = u[i];
		  else
		    b[i] = u[i];
		  
		  if(fu[i] <= fw[i] || w[i] == x[i])
		    {
		      v[i] = w[i];
		      w[i] = u[i];
		      fv[i] = fw[i];
		      fw[i] = fu[i];
		    }
		  else
		    {
		      if(fu[i] <= fv[i] || v[i] == x[i] || v[i] == w[i])
			{
			  v[i] = u[i];
			  fv[i] = fu[i];
			}
		    }	    
		}
	      
	      assert(a[i] >= lim_inf && a[i] <= lim_sup);
	      assert(b[i] >= lim_inf && b[i] <= lim_sup);
	      assert(x[i] >= lim_inf && x[i] <= lim_sup);
	      assert(v[i] >= lim_inf && v[i] <= lim_sup);
	      assert(w[i] >= lim_inf && w[i] <= lim_sup);
	      assert(u[i] >= lim_inf && u[i] <= lim_sup);
	    }
	}
    }

  free(converged);
  free(a);
  free(b);
  free(d);
  free(etemp);
  free(fu);
  free(fv);
  free(fw);
  free(fx);
  free(p);
  free(q);
  free(r);
  free(tol1);
  free(tol2);
  free(u);
  free(v);
  free(w);
  free(x);
  free(xm);
  free(e);

  printf("\n. Too many iterations in BRENT !");
  assert(0);
}



static int brakGeneric(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, 
		       double *fc, double lim_inf, double lim_sup, 
		       int numberOfModels, int rateNumber, analdef *adef, int whichFunction, tree *tr, linkageList *ll)
{
  double 
    *ulim = (double *)malloc(sizeof(double) * numberOfModels),
    *u    = (double *)malloc(sizeof(double) * numberOfModels),
    *r    = (double *)malloc(sizeof(double) * numberOfModels),
    *q    = (double *)malloc(sizeof(double) * numberOfModels),
    *fu   = (double *)malloc(sizeof(double) * numberOfModels),
    *dum  = (double *)malloc(sizeof(double) * numberOfModels), 
    *temp = (double *)malloc(sizeof(double) * numberOfModels);
  
  int 
    i,
    *state    = (int *)malloc(sizeof(int) * numberOfModels),
    *endState = (int *)malloc(sizeof(int) * numberOfModels);

  boolean *converged = (boolean *)malloc(sizeof(boolean) * numberOfModels);
  boolean allConverged;

  for(i = 0; i < numberOfModels; i++)
    converged[i] = FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      state[i] = 0;
      endState[i] = 0;

      u[i] = 0.0;

      param[i] = ax[i];

      if(param[i] > lim_sup) 	
	param[i] = ax[i] = lim_sup;
      
      if(param[i] < lim_inf) 
	param[i] = ax[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
   
  
  evaluateChange(tr, rateNumber, param, fa, converged, adef, whichFunction, numberOfModels, ll);


  for(i = 0; i < numberOfModels; i++)
    {
      param[i] = bx[i];
      if(param[i] > lim_sup) 
	param[i] = bx[i] = lim_sup;
      if(param[i] < lim_inf) 
	param[i] = bx[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
  
  evaluateChange(tr, rateNumber, param, fb, converged, adef, whichFunction, numberOfModels, ll);

  for(i = 0; i < numberOfModels; i++)  
    {
      if (fb[i] > fa[i]) 
	{	  
	  SHFT(dum[i],ax[i],bx[i],dum[i]);
	  SHFT(dum[i],fa[i],fb[i],dum[i]);
	}
      
      cx[i] = bx[i] + MNBRAK_GOLD * (bx[i] - ax[i]);
      
      param[i] = cx[i];
      
      if(param[i] > lim_sup) 
	param[i] = cx[i] = lim_sup;
      if(param[i] < lim_inf) 
	param[i] = cx[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
  
 
  evaluateChange(tr, rateNumber, param, fc, converged, adef, whichFunction, numberOfModels,  ll);

   while(1) 
     {       
       allConverged = TRUE;

       for(i = 0; i < numberOfModels && allConverged; i++)
	 allConverged = allConverged && converged[i];

       if(allConverged)
	 {
	   for(i = 0; i < numberOfModels; i++)
	     {	       
	       if(ax[i] > lim_sup) 
		 ax[i] = lim_sup;
	       if(ax[i] < lim_inf) 
		 ax[i] = lim_inf;

	       if(bx[i] > lim_sup) 
		 bx[i] = lim_sup;
	       if(bx[i] < lim_inf) 
		 bx[i] = lim_inf;
	       
	       if(cx[i] > lim_sup) 
		 cx[i] = lim_sup;
	       if(cx[i] < lim_inf) 
		 cx[i] = lim_inf;
	     }

	   free(converged);
	   free(ulim);
	   free(u);
	   free(r);
	   free(q);
	   free(fu);
	   free(dum); 
	   free(temp);
	   free(state);   
	   free(endState);
	   return 0;
	   
	 }

       for(i = 0; i < numberOfModels; i++)
	 {
	   if(!converged[i])
	     {
	       switch(state[i])
		 {
		 case 0:
		   endState[i] = 0;
		   if(!(fb[i] > fc[i]))		         
		     converged[i] = TRUE;		       		     
		   else
		     {
		   
		       if(ax[i] > lim_sup) 
			 ax[i] = lim_sup;
		       if(ax[i] < lim_inf) 
			 ax[i] = lim_inf;
		       if(bx[i] > lim_sup) 
			 bx[i] = lim_sup;
		       if(bx[i] < lim_inf) 
			 bx[i] = lim_inf;
		       if(cx[i] > lim_sup) 
			 cx[i] = lim_sup;
		       if(cx[i] < lim_inf) 
			 cx[i] = lim_inf;
		       
		       r[i]=(bx[i]-ax[i])*(fb[i]-fc[i]);
		       q[i]=(bx[i]-cx[i])*(fb[i]-fa[i]);
		       u[i]=(bx[i])-((bx[i]-cx[i])*q[i]-(bx[i]-ax[i])*r[i])/
			 (2.0*SIGN(MAX(fabs(q[i]-r[i]),MNBRAK_TINY),q[i]-r[i]));
		       
		       ulim[i]=(bx[i])+MNBRAK_GLIMIT*(cx[i]-bx[i]);
		       
		       if(u[i] > lim_sup) 
			 u[i] = lim_sup;
		       if(u[i] < lim_inf) 
			 u[i] = lim_inf;
		       if(ulim[i] > lim_sup) 
			 ulim[i] = lim_sup;
		       if(ulim[i] < lim_inf) 
			 ulim[i] = lim_inf;
		       
		       if ((bx[i]-u[i])*(u[i]-cx[i]) > 0.0)
			 {
			   param[i] = u[i];
			   if(param[i] > lim_sup) 			     
			     param[i] = u[i] = lim_sup;
			   if(param[i] < lim_inf)
			     param[i] = u[i] = lim_inf;
			   endState[i] = 1;
			 }
		       else 
			 {
			   if ((cx[i]-u[i])*(u[i]-ulim[i]) > 0.0) 
			     {
			       param[i] = u[i];
			       if(param[i] > lim_sup) 
				 param[i] = u[i] = lim_sup;
			       if(param[i] < lim_inf) 
				 param[i] = u[i] = lim_inf;
			       endState[i] = 2;
			     }		  	       
			   else
			     {
			       if ((u[i]-ulim[i])*(ulim[i]-cx[i]) >= 0.0) 
				 {
				   u[i] = ulim[i];
				   param[i] = u[i];	
				   if(param[i] > lim_sup) 
				     param[i] = u[i] = ulim[i] = lim_sup;
				   if(param[i] < lim_inf) 
				     param[i] = u[i] = ulim[i] = lim_inf;
				   endState[i] = 0;
				 }		  		
			       else 
				 {		  
				   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
				   param[i] = u[i];
				   endState[i] = 0;
				   if(param[i] > lim_sup) 
				     param[i] = u[i] = lim_sup;
				   if(param[i] < lim_inf) 
				     param[i] = u[i] = lim_inf;
				 }
			     }	  
			 }
		     }
		   break;
		 case 1:
		   endState[i] = 0;
		   break;
		 case 2:
		   endState[i] = 3;
		   break;
		 default:
		   assert(0);
		 }
	       assert(param[i] >= lim_inf && param[i] <= lim_sup);
	     }
	 }
             
       evaluateChange(tr, rateNumber, param, temp, converged, adef, whichFunction, numberOfModels, ll);

       for(i = 0; i < numberOfModels; i++)
	 {
	   if(!converged[i])
	     {	       
	       switch(endState[i])
		 {
		 case 0:
		   fu[i] = temp[i];
		   SHFT(ax[i],bx[i],cx[i],u[i]);
		   SHFT(fa[i],fb[i],fc[i],fu[i]);
		   state[i] = 0;
		   break;
		 case 1:
		   fu[i] = temp[i];
		   if (fu[i] < fc[i]) 
		     {
		       ax[i]=(bx[i]);
		       bx[i]=u[i];
		       fa[i]=(fb[i]);
		       fb[i]=fu[i]; 
		       converged[i] = TRUE;		      
		     } 
		   else 
		     {
		       if (fu[i] > fb[i]) 
			 {
			   assert(u[i] >= lim_inf && u[i] <= lim_sup);
			   cx[i]=u[i];
			   fc[i]=fu[i];
			   converged[i] = TRUE;			  
			 }
		       else
			 {		   
			   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
			   param[i] = u[i];
			   if(param[i] > lim_sup) {param[i] = u[i] = lim_sup;}
			   if(param[i] < lim_inf) {param[i] = u[i] = lim_inf;}	  
			   state[i] = 1;		 
			 }		  
		     }
		   break;
		 case 2: 
		   fu[i] = temp[i];
		   if (fu[i] < fc[i]) 
		     {		     
		       SHFT(bx[i],cx[i],u[i], cx[i]+MNBRAK_GOLD*(cx[i]-bx[i]));
		       state[i] = 2;
		     }	   
		   else
		     {
		       state[i] = 0;
		       SHFT(ax[i],bx[i],cx[i],u[i]);
		       SHFT(fa[i],fb[i],fc[i],fu[i]);
		     }
		   break;	   
		 case 3:		  
		   SHFT(fb[i],fc[i],fu[i], temp[i]);
		   SHFT(ax[i],bx[i],cx[i],u[i]);
		   SHFT(fa[i],fb[i],fc[i],fu[i]);
		   state[i] = 0;
		   break;
		 default:
		   assert(0);
		 }
	     }
	 }
    }
   

   assert(0);
   free(converged);
   free(ulim);
   free(u);
   free(r);
   free(q);
   free(fu);
   free(dum); 
   free(temp);
   free(state);   
   free(endState);

  

   return(0);
}


/**********************************************************************************************************/
/* ALPHA PARAM ********************************************************************************************/







static void optAlpha(tree *tr, double modelEpsilon, linkageList *ll)
{
  int 
    i, 
    k,
    numberOfModels = ll->entries;
  
  double 
    lim_inf     = ALPHA_MIN,
    lim_sup     = ALPHA_MAX;
  double
    *startLH    = (double *)malloc(sizeof(double) * numberOfModels),
    *startAlpha = (double *)malloc(sizeof(double) * numberOfModels),
    *endAlpha   = (double *)malloc(sizeof(double) * numberOfModels),
    *_a         = (double *)malloc(sizeof(double) * numberOfModels),
    *_b         = (double *)malloc(sizeof(double) * numberOfModels),
    *_c         = (double *)malloc(sizeof(double) * numberOfModels),
    *_fa        = (double *)malloc(sizeof(double) * numberOfModels),
    *_fb        = (double *)malloc(sizeof(double) * numberOfModels),
    *_fc        = (double *)malloc(sizeof(double) * numberOfModels),
    *_param     = (double *)malloc(sizeof(double) * numberOfModels),
    *result     = (double *)malloc(sizeof(double) * numberOfModels),
    *_x         = (double *)malloc(sizeof(double) * numberOfModels);   

#if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
   int revertModel = 0;
#endif   

  evaluateGenericInitrav(tr, tr->start);
   /* 
     at this point here every worker has the traversal data it needs for the 
     search, so we won't re-distribute it he he :-)
  */

  for(i = 0; i < numberOfModels; i++)
    {
      assert(ll->ld[i].valid);

      startAlpha[i] = tr->partitionData[ll->ld[i].partitionList[0]].alpha;
      _a[i] = startAlpha[i] + 0.1;
      _b[i] = startAlpha[i] - 0.1;      
      if(_b[i] < lim_inf) 
	_b[i] = lim_inf;

      startLH[i] = 0.0;
      
      for(k = 0; k < ll->ld[i].partitions; k++)	
	{
	  startLH[i] += tr->perPartitionLH[ll->ld[i].partitionList[k]];
	  assert(tr->partitionData[ll->ld[i].partitionList[0]].alpha ==  tr->partitionData[ll->ld[i].partitionList[k]].alpha);
	}
    }					  

  brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, -1, (analdef*)NULL, ALPHA_F, tr, ll);       
  brentGeneric(_a, _b, _c, _fb, modelEpsilon, _x, result, numberOfModels, ALPHA_F, -1, (analdef*)NULL, tr, ll, lim_inf, lim_sup);

  for(i = 0; i < numberOfModels; i++)
    endAlpha[i] = result[i];
  
  for(i = 0; i < numberOfModels; i++)
    {
      if(startLH[i] > endAlpha[i])
	{    	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    {	      
	      tr->partitionData[ll->ld[i].partitionList[k]].alpha = startAlpha[i];
	      makeGammaCats(tr->partitionData[ll->ld[i].partitionList[k]].alpha, tr->partitionData[ll->ld[i].partitionList[k]].gammaRates, 4); 		
	    }
#if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
	  revertModel++;
#endif
	}  
    }

#ifdef _USE_PTHREADS
  if(revertModel > 0)
    masterBarrier(THREAD_COPY_ALPHA, tr);
#endif

#ifdef _FINE_GRAIN_MPI
  if(revertModel > 0)
    masterBarrierMPI(THREAD_COPY_ALPHA, tr);
#endif

  
  free(startLH);
  free(startAlpha);
  free(endAlpha);
  free(result);
  free(_a);
  free(_b);
  free(_c);
  free(_fa);
  free(_fb);
  free(_fc);
  free(_param);
  free(_x);  

}



static void optRates(tree *tr, analdef *adef, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{
  int 
    i, 
    k, 
    j, 
    pos,
    numberOfRates = ((states * states - states) / 2) - 1;
    
  double lim_inf = RATE_MIN;
  double lim_sup = RATE_MAX;
  double 
    *startRates;
  double 
    *startLH= (double *)malloc(sizeof(double) * numberOfModels),
    *endLH  = (double *)malloc(sizeof(double) * numberOfModels),
    *_a     = (double *)malloc(sizeof(double) * numberOfModels),
    *_b     = (double *)malloc(sizeof(double) * numberOfModels),
    *_c     = (double *)malloc(sizeof(double) * numberOfModels),
    *_fa    = (double *)malloc(sizeof(double) * numberOfModels),
    *_fb    = (double *)malloc(sizeof(double) * numberOfModels),
    *_fc    = (double *)malloc(sizeof(double) * numberOfModels),
    *_param = (double *)malloc(sizeof(double) * numberOfModels),
    *result = (double *)malloc(sizeof(double) * numberOfModels),
    *_x     = (double *)malloc(sizeof(double) * numberOfModels); 

#if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
   int revertModel = 0;
#endif

   assert(states != -1);

  startRates = (double *)malloc(sizeof(double) * numberOfRates * numberOfModels);

  evaluateGenericInitrav(tr, tr->start);
  /* 
     at this point here every worker has the traversal data it needs for the 
     search 
  */

  for(i = 0, pos = 0; i < ll->entries; i++)
    {
      if(ll->ld[i].valid)
	{
	  endLH[pos] = unlikely;
	  startLH[pos] = 0.0;

	  for(j = 0; j < ll->ld[i].partitions; j++)
	    {
	      int index = ll->ld[i].partitionList[j];
	      
	      startLH[pos] += tr->perPartitionLH[index];
	      for(k = 0; k < numberOfRates; k++)
		startRates[pos * numberOfRates + k] = tr->partitionData[index].substRates[k];      
	    }
	  pos++;
	}
    }  

  assert(pos == numberOfModels);
  
  for(i = 0; i < numberOfRates; i++)
    {     
      for(k = 0, pos = 0; k < ll->entries; k++)
	{
	  if(ll->ld[k].valid)
	    {
	      int index = ll->ld[k].partitionList[0];
	      _a[pos] = tr->partitionData[index].substRates[i] + 0.1;
	      _b[pos] = tr->partitionData[index].substRates[i] - 0.1;
	      
	      if(_a[pos] < lim_inf) _a[pos] = lim_inf;
	      if(_a[pos] > lim_sup) _a[pos] = lim_sup;
	      
	      if(_b[pos] < lim_inf) _b[pos] = lim_inf;
	      if(_b[pos] > lim_sup) _b[pos] = lim_sup;    
	      pos++;
	    }
	}                       	     

      assert(pos == numberOfModels);

      brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, i, adef, RATE_F, tr, ll);
      
      for(k = 0; k < numberOfModels; k++)
	{
	  assert(_a[k] >= lim_inf && _a[k] <= lim_sup);
	  assert(_b[k] >= lim_inf && _b[k] <= lim_sup);	  
	  assert(_c[k] >= lim_inf && _c[k] <= lim_sup);	    
	}      

      brentGeneric(_a, _b, _c, _fb, modelEpsilon, _x, result, numberOfModels, RATE_F, i, adef, tr,  ll, lim_inf, lim_sup);
	
      for(k = 0; k < numberOfModels; k++)
	endLH[k] = result[k];	           
      
      for(k = 0, pos = 0; k < ll->entries; k++)
	{
#if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
	  revertModel = 0;
#endif
	  if(ll->ld[k].valid)
	    { 
	      if(startLH[pos] > endLH[pos])
		{		  
		  for(j = 0; j < ll->ld[k].partitions; j++)
		    {
		      int index = ll->ld[k].partitionList[j];
		      tr->partitionData[index].substRates[i] = startRates[pos * numberOfRates + i];	             	  
		      initReversibleGTR(tr, adef, index);
		    }
#if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
		  revertModel++;
#endif
		}
	      pos++;
	    }
	}

#ifdef _USE_PTHREADS      
      if(revertModel > 0)
	masterBarrier(THREAD_COPY_RATES, tr);
#endif
#ifdef _FINE_GRAIN_MPI
      if(revertModel > 0)
	masterBarrierMPI(THREAD_COPY_RATES, tr);
#endif
      
      assert(pos == numberOfModels);
    }

 
  free(startLH);
  free(endLH);
  free(result);
  free(_a);
  free(_b);
  free(_c);
  free(_fa);
  free(_fb);
  free(_fc);
  free(_param);
  free(_x);  
  free(startRates);
}

static boolean AAisGTR(tree *tr)
{
  int i, count = 0;

  for(i = 0; i < tr->NumberOfModels; i++)   
    {
      if(tr->partitionData[i].dataType == AA_DATA)
	{
	  count++;
	  if(tr->partitionData[i].protModels != GTR)
	    return FALSE;
	}
    }

  if(count == 0)
    return FALSE;

  return TRUE;
}

static void optRatesGeneric(tree *tr, analdef *adef, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    dnaPartitions = 0,
    aaPartitions  = 0,
    secondaryPartitions = 0,
    secondaryModel = -1,
    states = -1;

  /* assumes homogeneous super-partitions, that either contain DNA or AA partitions !*/
  /* does not check whether AA are all linked */

  /* first do DNA */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case DNA_DATA:	
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  ll->ld[i].valid = TRUE;
	  dnaPartitions++;  
	  break;
	case BINARY_DATA:
	case AA_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }   

  if(dnaPartitions > 0)
    optRates(tr, adef, modelEpsilon, ll, dnaPartitions, states);
  

  /* then SECONDARY */

   for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	  /* we only have one type of secondary data models in one analysis */
	case SECONDARY_DATA_6:
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  secondaryModel = SECONDARY_DATA_6;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case SECONDARY_DATA_7: 
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  secondaryModel = SECONDARY_DATA_7;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case SECONDARY_DATA:
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  secondaryModel = SECONDARY_DATA;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case BINARY_DATA:
	case AA_DATA:	
	case DNA_DATA:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }

  
   
   if(secondaryPartitions > 0)
     {
       assert(secondaryPartitions == 1);

       switch(secondaryModel)
	 {
	 case SECONDARY_DATA:
	   optRates(tr, adef, modelEpsilon, ll, secondaryPartitions, states);
	   break;
	 case SECONDARY_DATA_6:
	   optRates(tr, adef, modelEpsilon, ll, secondaryPartitions, states);
	   break;
	 case SECONDARY_DATA_7:
	   optRates(tr, adef, modelEpsilon, ll, secondaryPartitions, states);
	   break; 
	 default:
	   assert(0);
	 }
     }

  /* then AA for GTR */

  if(AAisGTR(tr))
    {
      for(i = 0; i < ll->entries; i++)
	{
	  switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	    {
	    case AA_DATA:
	      states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	      ll->ld[i].valid = TRUE;
	      aaPartitions++;
	      break;
	    case DNA_DATA:	    
	    case BINARY_DATA:
	    case SECONDARY_DATA:	
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7:
	      ll->ld[i].valid = FALSE;
	      break;
	    default:
	      assert(0);
	    }	 
	}

      assert(aaPartitions == 1);     
      
      optRates(tr, adef, modelEpsilon, ll, aaPartitions, states);
    }
  
  /* then multi-state */

  /* 
     now here we have to be careful, because every multi-state partition can actually 
     have a distinct number of states, so we will go to every multi-state partition separately,
     parallel efficiency for this will suck, but what the hell .....
  */

  if(tr->multiStateModel == GTR_MULTI_STATE)
    {     
      for(i = 0; i < ll->entries; i++)
	{
	  switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	    {
	    case GENERIC_32:
	      {
		int k;
		
		states = tr->partitionData[ll->ld[i].partitionList[0]].states;			      

		ll->ld[i].valid = TRUE;
		
		for(k = 0; k < ll->entries; k++)
		  if(k != i)
		    ll->ld[k].valid = FALSE;
		
		optRates(tr, adef, modelEpsilon, ll, 1, states);
	      }
	      break;
	    case AA_DATA:	    
	    case DNA_DATA:	    
	    case BINARY_DATA:
	    case SECONDARY_DATA:	
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7:
	    case GENERIC_64:
	      break;
	    default:
	      assert(0);
	    }	 
	}           
    }

  for(i = 0; i < ll->entries; i++)
    ll->ld[i].valid = TRUE;
}





/*********************FUNCTIONS FOR APPROXIMATE MODEL OPTIMIZATION ***************************************/






static int catCompare(const void *p1, const void *p2)
{
 rateCategorize *rc1 = (rateCategorize *)p1;
 rateCategorize *rc2 = (rateCategorize *)p2;

  double i = rc1->accumulatedSiteLikelihood;
  double j = rc2->accumulatedSiteLikelihood;
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}


static void categorizePartition(tree *tr, rateCategorize *rc, int model, int lower, int upper)
{
  int
    zeroCounter,
    i, 
    k;
  
  double 
    diff, 
    min;

  for (i = lower, zeroCounter = 0; i < upper; i++, zeroCounter++) 
      {
	double
	  temp = tr->cdta->patrat[i];

	int
	  found = 0;
	
	for(k = 0; k < tr->partitionData[model].numberOfCategories; k++)
	  {
	    if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
	      {
		found = 1;
		tr->cdta->rateCategory[i] = k;				
		break;
	      }
	  }
	
	if(!found)
	  {
	    min = fabs(temp - rc[0].rate);
	    tr->cdta->rateCategory[i] = 0;

	    for(k = 1; k < tr->partitionData[model].numberOfCategories; k++)
	      {
		diff = fabs(temp - rc[k].rate);

		if(diff < min)
		  {
		    min = diff;
		    tr->cdta->rateCategory[i] = k;
		  }
	      }
	  }
      }

  for(k = 0; k < tr->partitionData[model].numberOfCategories; k++)
    tr->partitionData[model].perSiteRates[k] = rc[k].rate; 
}


#if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))

void optRateCatPthreads(tree *tr, double lower_spacing, double upper_spacing, double *lhs, int n, int tid)
{
  int 
    model, 
    localIndex, 
    i;

  for(model = 0; model < tr->NumberOfModels; model++)
    {      
      int 
	localIndex = 0;

      boolean 
	execute = ((tr->manyPartitions && isThisMyPartition(tr, tid, model, n)) || (!tr->manyPartitions));

      if(execute)
	for(i = tr->partitionData[model].lower;  i < tr->partitionData[model].upper; i++)
	  {
	    if(tr->manyPartitions || (i % n == tid))
	      {
	      
		double initialRate, initialLikelihood, 
		  leftLH, rightLH, leftRate, rightRate, v;
		const double epsilon = 0.00001;
		int k;	      
		
		tr->cdta->patrat[i] = tr->cdta->patratStored[i];     
		initialRate = tr->cdta->patrat[i];
		
		initialLikelihood = evaluatePartialGeneric(tr, localIndex, initialRate, model); /* i is real i ??? */
		
		
		leftLH = rightLH = initialLikelihood;
		leftRate = rightRate = initialRate;
		
		k = 1;
		
		while((initialRate - k * lower_spacing > 0.0001) && 
		      ((v = evaluatePartialGeneric(tr, localIndex, initialRate - k * lower_spacing, model)) 
		       > leftLH) && 
		      (fabs(leftLH - v) > epsilon))  
		  {	  
#ifndef WIN32
		    if(isnan(v))
		      assert(0);
#endif
		    
		    leftLH = v;
		    leftRate = initialRate - k * lower_spacing;
		    k++;	  
		  }      
		
		k = 1;
		
		while(((v = evaluatePartialGeneric(tr, localIndex, initialRate + k * upper_spacing, model)) > rightLH) &&
		      (fabs(rightLH - v) > epsilon))    	
		  {
#ifndef WIN32
		    if(isnan(v))
		      assert(0);
#endif     
		    rightLH = v;
		    rightRate = initialRate + k * upper_spacing;	 
		    k++;
		  }           
		
		if(rightLH > initialLikelihood || leftLH > initialLikelihood)
		  {
		    if(rightLH > leftLH)	    
		      {	     
			tr->cdta->patrat[i] = rightRate;
			lhs[i] = rightLH;
		      }
		    else
		      {	      
			tr->cdta->patrat[i] = leftRate;
			lhs[i] = leftLH;
		      }
		  }
		else
		  lhs[i] = initialLikelihood;
		
		tr->cdta->patratStored[i] = tr->cdta->patrat[i];
		localIndex++;
	      }
	  }
      assert(localIndex == tr->partitionData[model].width);    
    }
}



#else


static void optRateCatModel(tree *tr, int model, double lower_spacing, double upper_spacing, double *lhs)
{
  int lower = tr->partitionData[model].lower;
  int upper = tr->partitionData[model].upper;
  int i;
  for(i = lower; i < upper; i++)
    {
      double initialRate, initialLikelihood, 
	leftLH, rightLH, leftRate, rightRate, v;
      const double epsilon = 0.00001;
      int k;
      
      tr->cdta->patrat[i] = tr->cdta->patratStored[i];     
      initialRate = tr->cdta->patrat[i];
      
      initialLikelihood = evaluatePartialGeneric(tr, i, initialRate, model); 
      
      
      leftLH = rightLH = initialLikelihood;
      leftRate = rightRate = initialRate;
      
      k = 1;
      
      while((initialRate - k * lower_spacing > 0.0001) && 
	    ((v = evaluatePartialGeneric(tr, i, initialRate - k * lower_spacing, model)) 
	     > leftLH) && 
	    (fabs(leftLH - v) > epsilon))  
	{	  
#ifndef WIN32
	  if(isnan(v))
	    assert(0);
#endif
	  
	  leftLH = v;
	  leftRate = initialRate - k * lower_spacing;
	  k++;	  
	}      
      
      k = 1;
      
      while(((v = evaluatePartialGeneric(tr, i, initialRate + k * upper_spacing, model)) > rightLH) &&
	    (fabs(rightLH - v) > epsilon))    	
	{
#ifndef WIN32
	  if(isnan(v))
	    assert(0);
#endif     
	  rightLH = v;
	  rightRate = initialRate + k * upper_spacing;	 
	  k++;
	}           
  
      if(rightLH > initialLikelihood || leftLH > initialLikelihood)
	{
	  if(rightLH > leftLH)	    
	    {	     
	      tr->cdta->patrat[i] = rightRate;
	      lhs[i] = rightLH;
	    }
	  else
	    {	      
	      tr->cdta->patrat[i] = leftRate;
	      lhs[i] = leftLH;
	    }
	}
      else
	lhs[i] = initialLikelihood;
      
      tr->cdta->patratStored[i] = tr->cdta->patrat[i];
    }

}


#endif



/* 
   set scaleRates to FALSE everywhere such that 
   per-site rates are not scaled to obtain an overall mean rate 
   of 1.0
*/

void updatePerSiteRates(tree *tr, boolean scaleRates)
{
  int 
    i,
    model;

  if(tr->multiBranch)
    {            
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int 	       
	    lower = tr->partitionData[model].lower,
	    upper = tr->partitionData[model].upper;
	  
	  if(scaleRates)
	    {
	      double 
		scaler = 0.0,       
		accRat = 0.0; 

	      int 
		accWgt     = 0;
	      
	      for(i = lower; i < upper; i++)
		{
		  int 
		    w = tr->cdta->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->cdta->rateCategory[i]];
		  
		  assert(0 <= tr->cdta->rateCategory[i] && tr->cdta->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}	   
	  
	      accRat /= ((double)accWgt);
	  
	      scaler = 1.0 / ((double)accRat);
	  	  
	      for(i = 0; i < tr->partitionData[model].numberOfCategories; i++)
		tr->partitionData[model].perSiteRates[i] *= scaler;	    

	      accRat = 0.0;	 
	      
	      for(i = lower; i < upper; i++)
		{
		  int 
		    w = tr->cdta->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->cdta->rateCategory[i]];
		  
		  assert(0 <= tr->cdta->rateCategory[i] && tr->cdta->rateCategory[i] < tr->maxCategories);	      
		  
		  accRat += (w * rate);
		}	         

	      accRat /= ((double)accWgt);	  

	      assert(ABS(1.0 - accRat) < 1.0E-5);
	    }
	  else
	    {
	      double 		   
		accRat = 0.0; 

	      int 
		accWgt     = 0;
	      
	      for(i = lower; i < upper; i++)
		{
		  int 
		    w = tr->cdta->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->cdta->rateCategory[i]];
		  
		  assert(0 <= tr->cdta->rateCategory[i] && tr->cdta->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}	   
	  
	      accRat /= ((double)accWgt);
	      
	      assert(ABS(1.0 - accRat) < 1.0E-5);
	    }

	  for(i = lower; i < upper; i++)
	    {
	      double
		w = ((double)tr->cdta->aliaswgt[i]);	      

	       double
		 wtemp,
		 temp = tr->partitionData[model].perSiteRates[tr->cdta->rateCategory[i]];

	       wtemp = temp * w;

	       tr->wr[i]  = wtemp;
	       tr->wr2[i] = temp * wtemp;
	    }	            	  
	  
#if !(defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
	  {
	    int 
	      localCount = 0;
	    
	    for(i = lower, localCount = 0; i < upper; i++, localCount++)
	      {	    	      
		tr->partitionData[model].wr[localCount]  = tr->wr[i];
		tr->partitionData[model].wr2[localCount] = tr->wr2[i];
	      }
	  }
#endif
	}
    }
  else
    {
      int
	accWgt = 0;

      double 
	scaler = 0.0,       
	accRat = 0.0; 

      if(scaleRates)
	{
	  for(model = 0, accRat = 0.0, accWgt = 0; model < tr->NumberOfModels; model++)
	    {
	      int 
		localCount = 0,
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper;
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{
		  int 
		    w = tr->cdta->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->cdta->rateCategory[i]];
		  
		  assert(0 <= tr->cdta->rateCategory[i] && tr->cdta->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}
	    }
	  
	  accRat /= ((double)accWgt);
	  
	  scaler = 1.0 / ((double)accRat);
	  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      for(i = 0; i < tr->partitionData[model].numberOfCategories; i++)
		tr->partitionData[model].perSiteRates[i] *= scaler;
	    }

	  for(model = 0, accRat = 0.0; model < tr->NumberOfModels; model++)
	    {
	      int 
		localCount = 0,
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper;
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{
		  int 
		    w = tr->cdta->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->cdta->rateCategory[i]];
		  
		  assert(0 <= tr->cdta->rateCategory[i] && tr->cdta->rateCategory[i] < tr->maxCategories);	      
		  
		  accRat += (w * rate);
		}
	    }           

	  accRat /= ((double)accWgt);	  

	  assert(ABS(1.0 - accRat) < 1.0E-5);
	}
      else
	{
	  for(model = 0, accRat = 0.0, accWgt = 0; model < tr->NumberOfModels; model++)
	    {
	      int 
		localCount = 0,
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper;
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{
		  int 
		    w = tr->cdta->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->cdta->rateCategory[i]];
		  
		  assert(0 <= tr->cdta->rateCategory[i] && tr->cdta->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}
	    }
	  
	  accRat /=  (double)accWgt;

	  assert(ABS(1.0 - accRat) < 1.0E-5);
	}
         
       for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int 
	    localCount = 0,
	    lower = tr->partitionData[model].lower,
	    upper = tr->partitionData[model].upper;

	  for(i = lower, localCount = 0; i < upper; i++, localCount++)
	    {
	      double
		w = ((double)tr->cdta->aliaswgt[i]);	      

	       double
		 wtemp,
		 temp = tr->partitionData[model].perSiteRates[tr->cdta->rateCategory[i]];

	       wtemp = temp * w;

	       tr->wr[i]  = wtemp;
	       tr->wr2[i] = temp * wtemp;
	    }
	}         
#if !(defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
      for(model = 0; model < tr->NumberOfModels; model++)
	{   	  	  	 
	  int 
	    localCount,
	    lower = tr->partitionData[model].lower,
	    upper = tr->partitionData[model].upper;	  
	  
	  for(i = lower, localCount = 0; i < upper; i++, localCount++)
	    {	    	      
	      tr->partitionData[model].wr[localCount]  = tr->wr[i];
	      tr->partitionData[model].wr2[localCount] = tr->wr2[i];
	    }
	}
#endif
    }
  
#ifdef _FINE_GRAIN_MPI
  masterBarrierMPI(THREAD_COPY_RATE_CATS, tr);
#endif  
   
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_COPY_RATE_CATS, tr);
#endif               
}

static void optimizeRateCategories(tree *tr, int _maxCategories)
{
  assert(_maxCategories > 0);

  if(_maxCategories > 1)
    {
      double  
	temp,  
	lower_spacing, 
	upper_spacing,
	initialLH = tr->likelihood,	
	*ratStored = (double *)malloc(sizeof(double) * tr->cdta->endsite),
	*lhs =       (double *)malloc(sizeof(double) * tr->cdta->endsite),
	**oldCategorizedRates = (double **)malloc(sizeof(double *) * tr->NumberOfModels);

      int  
	i,
	k,
	maxCategories = _maxCategories,
	*oldCategory =  (int *)malloc(sizeof(int) * tr->cdta->endsite),
	model,
	*oldNumbers = (int *)malloc(sizeof(int) * tr->NumberOfModels);
  
      assert(isTip(tr->start->number, tr->rdta->numsp));         

      if(tr->multiGene)
	determineFullTraversalMulti(tr->start, tr);
      else
	determineFullTraversal(tr->start, tr);

      if(optimizeRateCategoryInvocations == 1)
	{
	  lower_spacing = 0.5 / ((double)optimizeRateCategoryInvocations);
	  upper_spacing = 1.0 / ((double)optimizeRateCategoryInvocations);
	}
      else
	{
	  lower_spacing = 0.05 / ((double)optimizeRateCategoryInvocations);
	  upper_spacing = 0.1 / ((double)optimizeRateCategoryInvocations);
	}
      
      if(lower_spacing < 0.001)
	lower_spacing = 0.001;
      
      if(upper_spacing < 0.001)
	upper_spacing = 0.001;
      
      optimizeRateCategoryInvocations++;

      memcpy(oldCategory, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);	     
      memcpy(ratStored,   tr->cdta->patratStored, sizeof(double) * tr->cdta->endsite);

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  oldNumbers[model]          = tr->partitionData[model].numberOfCategories;

	  oldCategorizedRates[model] = (double *)malloc(sizeof(double) * tr->maxCategories);
	  
	  memcpy(oldCategorizedRates[model], tr->partitionData[model].perSiteRates, tr->maxCategories * sizeof(double));	  	 	  
	}      
      
#ifdef _USE_PTHREADS
      tr->lhs = lhs;
      tr->lower_spacing = lower_spacing;
      tr->upper_spacing = upper_spacing;
      masterBarrier(THREAD_RATE_CATS, tr);
#else
#ifdef _FINE_GRAIN_MPI
      tr->lhs = lhs;
      tr->lower_spacing = lower_spacing;
      tr->upper_spacing = upper_spacing;
      masterBarrierMPI(THREAD_RATE_CATS, tr);
#else
      for(model = 0; model < tr->NumberOfModels; model++)      
	optRateCatModel(tr, model, lower_spacing, upper_spacing, lhs);
#endif     
#endif

      for(model = 0; model < tr->NumberOfModels; model++)
	{     
	  int 
	    where = 1,
	    found = 0,
	    width = tr->partitionData[model].upper -  tr->partitionData[model].lower,
	    upper = tr->partitionData[model].upper,
	    lower = tr->partitionData[model].lower;
	    
	  rateCategorize 
	    *rc = (rateCategorize *)malloc(sizeof(rateCategorize) * width);		 
	
	  for (i = 0; i < width; i++)
	    {
	      rc[i].accumulatedSiteLikelihood = 0.0;
	      rc[i].rate = 0.0;
	    }  
	
	  rc[0].accumulatedSiteLikelihood = lhs[lower];
	  rc[0].rate = tr->cdta->patrat[lower];
	
	  tr->cdta->rateCategory[lower] = 0;
	
	  for (i = lower + 1; i < upper; i++) 
	    {
	      temp = tr->cdta->patrat[i];
	      found = 0;
	    
	      for(k = 0; k < where; k++)
		{
		  if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
		    {
		      found = 1;						
		      rc[k].accumulatedSiteLikelihood += lhs[i];	
		      break;
		    }
		}
	    
	      if(!found)
		{	    
		  rc[where].rate = temp;	    
		  rc[where].accumulatedSiteLikelihood += lhs[i];	    
		  where++;
		}
	    }
	
	  qsort(rc, where, sizeof(rateCategorize), catCompare);
	
	  if(where < maxCategories)
	    {
	      tr->partitionData[model].numberOfCategories = where;
	      categorizePartition(tr, rc, model, lower, upper);
	    }
	  else
	    {
	      tr->partitionData[model].numberOfCategories = maxCategories;	
	      categorizePartition(tr, rc, model, lower, upper);
	    }
	
	  free(rc);
	}
        	
      updatePerSiteRates(tr, TRUE);	

      evaluateGenericInitrav(tr, tr->start);
      
      if(tr->likelihood < initialLH)
	{	 		  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      tr->partitionData[model].numberOfCategories = oldNumbers[model];
	      memcpy(tr->partitionData[model].perSiteRates, oldCategorizedRates[model], tr->maxCategories * sizeof(double));
	    }	      
	  
	  memcpy(tr->cdta->patratStored, ratStored, sizeof(double) * tr->cdta->endsite);
	  memcpy(tr->cdta->rateCategory, oldCategory, sizeof(int) * tr->cdta->endsite);	     
	  
	  updatePerSiteRates(tr, FALSE);
	  
	  evaluateGenericInitrav(tr, tr->start);

	  /* printf("REVERT: %1.40f %1.40f\n", initialLH, tr->likelihood); */

	  assert(initialLH == tr->likelihood);
	}
          
      for(model = 0; model < tr->NumberOfModels; model++)
	free(oldCategorizedRates[model]);
                   
      free(oldCategorizedRates);
      free(oldCategory);
      free(ratStored);       
      free(lhs); 
      free(oldNumbers);
    }
}
  






/*****************************************************************************************************/

void resetBranches(tree *tr)
{
  nodeptr  p, q;
  int  nodes, i;
  
  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {   
      for(i = 0; i < tr->numBranches; i++)
	p->z[i] = defaultz;
	
      q = p->next;
      while(q != p)
	{	
	  for(i = 0; i < tr->numBranches; i++)
	    q->z[i] = defaultz;	    
	  q = q->next;
	}
      p++;
    }
}


static void printAAmatrix(tree *tr, double epsilon)
{
  if(AAisGTR(tr))
    {
      int model;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->partitionData[model].dataType == AA_DATA) 
	    {
	      char gtrFileName[1024];
	      char epsilonStr[1024];
	      FILE *gtrFile;
	      double *rates = tr->partitionData[model].substRates;
	      double *f     = tr->partitionData[model].frequencies;
	      double q[20][20];
	      int    r = 0;
	      int i, j;

	      assert(tr->partitionData[model].protModels == GTR);

	      sprintf(epsilonStr, "%f", epsilon);

	      strcpy(gtrFileName, workdir);
	      strcat(gtrFileName, "RAxML_proteinGTRmodel.");
	      strcat(gtrFileName, run_id);
	      strcat(gtrFileName, "_");
	      strcat(gtrFileName, epsilonStr);

	      gtrFile = myfopen(gtrFileName, "wb");

	      for(i = 0; i < 20; i++)
		for(j = 0; j < 20; j++)
		  q[i][j] = 0.0;

	      for(i = 0; i < 19; i++)
		for(j = i + 1; j < 20; j++)
		  q[i][j] = rates[r++];

	      for(i = 0; i < 20; i++)
		for(j = 0; j <= i; j++)
		  {
		    if(i == j)
		      q[i][j] = 0.0;
		    else
		      q[i][j] = q[j][i];
		  }
	   
	      for(i = 0; i < 20; i++)
		{
		  for(j = 0; j < 20; j++)		
		    fprintf(gtrFile, "%1.80f ", q[i][j]);
		
		  fprintf(gtrFile, "\n");
		}
	      for(i = 0; i < 20; i++)
		fprintf(gtrFile, "%1.80f ", f[i]);
	      fprintf(gtrFile, "\n");

	      fclose(gtrFile);

	      printBothOpen("\nPrinted intermediate AA substitution matrix to file %s\n\n", gtrFileName);
	      
	      break;
	    }

	}	  
    }
}


static void autoProtein(tree *tr, analdef *adef)
{
  int 
    countAutos = 0,
    i,
    model;
    
   for(model = 0; model < tr->NumberOfModels; model++)	      
     if(tr->partitionData[model].protModels == AUTO)
       countAutos++;
  
  if(countAutos > 0)
    {
      int 
	numProteinModels = AUTO,
	*bestIndex = (int*)malloc(sizeof(int) * tr->NumberOfModels);

      double
	*bestScores = (double*)malloc(sizeof(double) * tr->NumberOfModels);

      /*printf("Entry: %f\n", tr->likelihood);*/
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  bestIndex[model] = -1;
	  bestScores[model] = unlikely;
	}

      for(i = 0; i < numProteinModels; i++)
	{
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {	   
	      if(tr->partitionData[model].protModels == AUTO)
		{
		  tr->partitionData[model].autoProtModels = i;
		  initReversibleGTR(tr, adef, model);  
		}
	    }
	  
#ifdef _USE_PTHREADS	
	  masterBarrier(THREAD_BROADCAST_RATE, tr);	   
#endif
#ifdef _FINE_GRAIN_MPI
	  masterBarrierMPI(THREAD_BROADCAST_RATE, tr);     
#endif
	  
	  resetBranches(tr);
	  evaluateGenericInitrav(tr, tr->start);  
	  treeEvaluate(tr, 0.5);     

	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->partitionData[model].protModels == AUTO)
		{
		  if(tr->perPartitionLH[model] > bestScores[model])
		    {
		      bestScores[model] = tr->perPartitionLH[model];
		      bestIndex[model] = i;		      
		    }
		}

	    }       
	}

      printBothOpen("\n\n");
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{	   
	  if(tr->partitionData[model].protModels == AUTO)
	    {
	      tr->partitionData[model].autoProtModels = bestIndex[model];
	      initReversibleGTR(tr, adef, model);  
	      printBothOpen("Partition: %d best-scoring AA model: %s likelihood %f\n", model, protModels[tr->partitionData[model].autoProtModels], bestScores[model]);
	    }	 
	}
      
      printBothOpen("\n\n");
            
#ifdef _USE_PTHREADS	
	  masterBarrier(THREAD_BROADCAST_RATE, tr);	   
#endif
#ifdef _FINE_GRAIN_MPI
	  masterBarrierMPI(THREAD_BROADCAST_RATE, tr);     
#endif

      resetBranches(tr);
      evaluateGenericInitrav(tr, tr->start); 
      treeEvaluate(tr, 0.5);
      
      /*printf("Exit: %f\n", tr->likelihood);*/
      
      free(bestIndex);
      free(bestScores);
    }
}


static void calcPerSiteProteinModels(tree *tr)
{
    
  determineFullTraversal(tr->start, tr);
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_OPTIMIZE_PER_SITE_AA, tr);  
#else
#ifdef _FINE_GRAIN_MPI
  masterBarrierMPI(THREAD_OPTIMIZE_PER_SITE_AA, tr); 
#else
  {
    int
      s,
      p;
    
    double  
      *bestScore = (double *)malloc(tr->cdta->endsite * sizeof(double));	  
    
    for(s = 0; s < tr->cdta->endsite; s++)	    
      bestScore[s] = unlikely;
    for(p = 0; p < NUM_PROT_MODELS - 2; p++)
      {
	int 
	  model;
	
	for(model = 0; model < tr->NumberOfModels; model++)
	  {
	    double
	      lh;
	    
	    int
	      counter,
	      i,
	      lower = tr->partitionData[model].lower,
	      upper = tr->partitionData[model].upper;
	    	    
	    memcpy(tr->partitionData[model].EIGN,        tr->siteProtModel[p].EIGN,        sizeof(double) * 19);
	    memcpy(tr->partitionData[model].EV,          tr->siteProtModel[p].EV,          sizeof(double) * 400);                
	    memcpy(tr->partitionData[model].EI,          tr->siteProtModel[p].EI,          sizeof(double) * 380);
	    memcpy(tr->partitionData[model].substRates,  tr->siteProtModel[p].substRates,  sizeof(double) * 190);        
	    memcpy(tr->partitionData[model].frequencies, tr->siteProtModel[p].frequencies, sizeof(double) * 20);
	    memcpy(tr->partitionData[model].tipVector,   tr->siteProtModel[p].tipVector,   sizeof(double) * 460);
	    
	    for(i = lower, counter = 0; i < upper; i++, counter++)
	      {
		lh = evaluatePartialGeneric(tr, i, 0.0, model);
		
		if(lh > bestScore[i])
		  {
		    bestScore[i] = lh;		    
		    tr->partitionData[model].perSiteAAModel[counter] = p;
		  }
	      }
	  }	     	           
      }
    
    free(bestScore);
  }
#endif
#endif
}


void modOpt(tree *tr, analdef *adef, boolean resetModel, double likelihoodEpsilon, boolean testGappedImplementation)
{ 
  int i, model, catOpt = 0; 
  double 
    currentLikelihood,
    modelEpsilon = 0.0001;
  linkageList *alphaList;
  linkageList *invarList;
  linkageList *rateList; 
  
  int *unlinked = (int *)malloc(sizeof(int) * tr->NumberOfModels);
      
  modelEpsilon = 0.0001;
    
  for(i = 0; i < tr->NumberOfModels; i++)
    unlinked[i] = i;
  
  alphaList = initLinkageList(unlinked, tr);
  invarList = initLinkageList(unlinked, tr);
  rateList  = initLinkageListGTR(tr);
   
  tr->start = tr->nodep[1];
                 
  do
    {           
      currentLikelihood = tr->likelihood;     
     
      optRatesGeneric(tr, adef, modelEpsilon, rateList);
           
      onlyInitrav(tr, tr->start);                                       

      autoProtein(tr, adef);

      treeEvaluate(tr, 0.0625);      
      
     

      switch(tr->rateHetModel)
	{
	case GAMMA:      
	  optAlpha(tr, modelEpsilon, alphaList); 
	  onlyInitrav(tr, tr->start); 	 	 
	  treeEvaluate(tr, 0.1);	  	 
	  break;
	case CAT:
	  if(catOpt < 3)
	    {	      	     	     
	      optimizeRateCategories(tr, adef->categories);	      	     	      	      
	      catOpt++;
	    }
	  break;	  
	default:
	  assert(0);
	}                   
      
       if(tr->estimatePerSiteAA)
	{	 	     
	  calcPerSiteProteinModels(tr);	  

	  /*for(s = 0; s < tr->cdta->endsite; s++)
	    printf("Site %d model %s\n", s, protModels[tr->modelAssignment[s]]);
	    printf("\n\n");*/
	  
	  onlyInitrav(tr, tr->start); 	 	 
	 	  
	  treeEvaluate(tr, 0.1);	  	 
	}


      printAAmatrix(tr, fabs(currentLikelihood - tr->likelihood));    
    }
  while(fabs(currentLikelihood - tr->likelihood) > likelihoodEpsilon);  
  
  free(unlinked);
  freeLinkageList(alphaList);
  freeLinkageList(rateList);
  freeLinkageList(invarList);  
}

#ifdef _JOERG



/* loop and apply numerical optimization routines for branch lengths and alpha until the difference
 in log likelihood improvement gets smaller than likelihoodEpsilon log likelihood units.
 Note that, the actual pergormance run time depends heavily on likelihoodEpsilon since this influences
 the number of inner do while lopp iterations */
void optimize(tree *tr, linkageList *alphaList) {
	double likelihoodEpsilon = 5.0, currentLikelihood, modelEpsilon = 0.0001;
//    clock_t begin, end;
//
//	begin = clock();

#ifdef _FINE_GRAIN_MPI
		masterBarrierMPI(THREAD_COPY_INIT_MODEL, tr);
#endif

#ifdef _USE_PTHREADS
		masterBarrier(THREAD_COPY_INIT_MODEL, tr);
#endif

		evaluateGenericInitrav(tr, tr->start);

	do {
		/* remember the current likelihood */
		currentLikelihood = tr->likelihood;

		/* entirely re-traverse the tree using Felsensteins pruning algorithm */
		onlyInitrav(tr, tr->start);

		/* optimize the branch lengths a bit */
		treeEvaluate(tr, 0.0625);

		switch (tr->rateHetModel) {
		case GAMMA:
			/* optimize the alpha shape parameter for each partition individually and independently
			 using brent's numerical optimization algorithm. For details on the general algorithm
			 see the book: Numerical recipees in C */
			optAlpha(tr, modelEpsilon, alphaList);

			/* re-traverse the entire tree using Felsenstein's algorithm */
			onlyInitrav(tr, tr->start);

			/* optimize branch lengths a bit more */
			treeEvaluate(tr, 0.1);
			break;
		default:
			assert(0);
			break;
		}
	} while (fabs(currentLikelihood - tr->likelihood) > likelihoodEpsilon);

//	end = clock();
//
//	assignmentEvalTicks += ((double)(end - begin)) / (CLOCKS_PER_SEC * NumberOfThreads);
}

void init(tree *tr, analdef *adef) {
	int model;

	resetBranches(tr);

	for (model = 0; model < tr->NumberOfModels; model++) {
		/* set the alpha shape parameter to its default value */
		tr->partitionData[model].alpha = 1.0;

		/* initialize the discretization of the GAMMA function --- we use 4 discrete rates
		 to approximate the integral over GAMMA we actually want to compute */
		makeGammaCats(tr->partitionData[model].alpha, tr->partitionData[model].gammaRates, 4);

//		use empirical protein frequencies if specified
		if(!adef->protEmpiricalFreqs)
			tr->partitionData[model].protFreqs = 0;
		else tr->partitionData[model].protFreqs = 1;
	}
}

assignment* mallocAssignment(int models) {
	assignment *a;

	a = (assignment*) malloc(sizeof(assignment));
	a->partitionLH = (double *) malloc(models * sizeof(double));
	a->partitionModel = (int *) malloc(models * sizeof(int));
	a->chosen = 0;

	return a;
}

void writeTest(mtest *test) {
	int i, model;

	FILE *f = fopen("testDump", "w");

	fwrite(test, sizeof(mtest), 1, f);
	for(i = 0; i < test->nrRuns; i++) {
		fwrite(test->run[i]->partitionModel, sizeof(int), test->nrModels, f);
		fwrite(test->run[i]->partitionLH, sizeof(double), test->nrModels, f);
		fwrite(&test->run[i]->overallLH, sizeof(double), 1, f);
	}
	fclose(f);
}

mtest* readTest() {
	FILE *f;
	mtest *test = (mtest*) malloc(sizeof(mtest));
	size_t s;
	int i;

	f = fopen("testDump", "r");
	s = fread(test, sizeof(mtest), 1, f);

	test->run = (assignment**) malloc(test->nrRuns * sizeof(assignment*));

	for(i = 0; i < test->nrRuns; i++) {
		test->run[i] = mallocAssignment(test->nrModels);
		s = fread(test->run[i]->partitionModel, sizeof(int), test->nrModels, f);
		s = fread(test->run[i]->partitionLH, sizeof(double), test->nrModels, f);
		s = fread(&test->run[i]->overallLH, sizeof(double), 1, f);
	}
	fclose(f);

	return test;
}

// computes the likelihood values for the given assignment
void evaluateAssignment(tree *tr, analdef *adef, assignment *ass,  linkageList *alphaList) {
	int i, model, catOpt = 0;
//    clock_t begin, end;
//
//	begin = clock();

	init(tr, adef);

	for (model = 0; model < tr->NumberOfModels; model++) {
		tr->partitionData[model].protModels = ass->partitionModel[model];

		initReversibleGTR(tr, adef, model);
//		if(assertionError) {
//			printf("The following assignment of models caused evaluation problems\n");
//			printAssignment(ass, tr->NumberOfModels);
//			assertionError = 0;
//		}
	}

	optimize(tr, alphaList);

	for (model = 0; model < tr->NumberOfModels; model++) {
		ass->partitionLH[model] = tr->perPartitionLH[model];
		ass->partitionModel[model] = tr->partitionData[model].protModels;
	}

	ass->overallLH = tr->likelihood;
//	end = clock();
//	assignmentEvalTicks += ((double)(end - begin)) / (CLOCKS_PER_SEC * NumberOfThreads);
}

void printTime(double p) {
	printf("\r%d%% processed. Overall time estimate is %f sec. or %3.1f min....", (int) p, (assignmentEvalTicks * 100 / p), ((assignmentEvalTicks * 100 / p) /60));
	fflush(stdout);
}



assignment* partitionWiseOptimum(mtest *test) {
	int i, model;

	assignment *result = mallocAssignment(test->nrModels);

	for(model = 0; model < test->nrModels; model++)
		result->partitionLH[model] = unlikely;
	result->overallLH = 0;

	for(i = 0; i < test->nrRuns; i++) {
		for (model = 0; model < test->nrModels; model++) {
			if (test->run[i]->partitionLH[model] > result->partitionLH[model]) {
				result->partitionLH[model] = test->run[i]->partitionLH[model];
				result->partitionModel[model] = test->run[i]->partitionModel[model];
			}
		}
	}

	for (model = 0; model < test->nrModels; model++)
		result->overallLH += result->partitionLH[model];

	return result;
}

assignment* overallOptimum(mtest *test) {
	int i, model;

	assignment *result = mallocAssignment(test->nrModels);
	result->overallLH = unlikely;

	for(i = 0; i < test->nrRuns; i++) {
		if(test->run[i]->overallLH > result->overallLH) {
			result->overallLH = test->run[i]->overallLH;
			for (model = 0; model < test->nrModels; model++) {
				result->partitionLH[model] = test->run[i]->partitionLH[model];
				result->partitionModel[model] = test->run[i]->partitionModel[model];
			}
		}
	}

	return result;
}

void sortUnlinked(mtest *test) {
	int n = test->nrRuns, model, i, swap, tmpModel;
	assignment *tmp;
	double tmpLH, sum;

	do{
		swap = 0;

		for(i = 0; i < test->nrRuns - 1; i++) {
			for(model = 0; model < test->nrModels; model++) {
				if(test->run[i]->partitionLH[model] < test->run[i + 1]->partitionLH[model]) {
					swap = 1;
					tmpLH = test->run[i]->partitionLH[model];
					tmpModel = test->run[i]->partitionModel[model];

					test->run[i]->partitionLH[model] = test->run[i + 1]->partitionLH[model];
					test->run[i]->partitionModel[model] = test->run[i + 1]->partitionModel[model];

					test->run[i + 1]->partitionLH[model] = tmpLH;
					test->run[i + 1]->partitionModel[model] = tmpModel;
				}
			}
		}
		n--;
	} while (swap && n > 0);

	for(i = 0; i < test->nrRuns; i++) {
		sum = 0;
		for(model = 0; model < test->nrModels; model++)
			sum += test->run[i]->partitionLH[model];
		test->run[i]->overallLH = sum;
	}

}

//TODO this is simple bubble sort, replace by quicksort or something else
void sort(mtest *test) {
	int n = test->nrRuns, i, swap;
	assignment *tmp;

	do{
		swap = 0;

		for(i = 0; i < test->nrRuns - 1; i++) {
			if(test->run[i]->overallLH < test->run[i + 1]->overallLH) {
				tmp = test->run[i];
				test->run[i] = test->run[i + 1];
				test->run[i + 1] = tmp;
				swap = 1;
			}
		}
		n--;
	} while (swap && n > 0);
}

// mutates n random partitions of an assignment with m partitions
//TODO check whether every partition is assigned a new model if n = m = nrModels or not
assignment* mutate(assignment *start, int m, int n) {
	int i, model, partition, *parts = malloc(sizeof(int) * n), nrAAModels = NUM_PROT_MODELS - 2, in, im;
	long randomSeed = rand();

	assignment *res = mallocAssignment(m);

	for(i = 0; i < m; i++) {
		res->partitionModel[i] = start->partitionModel[i];
		res->partitionLH[i] = unlikely;
	}
	res->overallLH = unlikely;

	im = 0;
	for (in = 0; in < m && im < n; ++in) {
	  int rn = m - in;
	  int rm = n - im;
	  if (rand() % rn < rm)
		parts[im++] = in;
	}

//	printf("partitions: ");
//	for(i = 0; i < n; i++)
//		printf("%d ", parts[i]);

	for(i = 0; i < n; i++) {
		do {
			res->partitionModel[parts[i]] = (int)(randum(&randomSeed) * nrAAModels);
		} while(res->partitionModel[parts[i]] == start->partitionModel[parts[i]]);
	}

	return res;
}

// enable disable joint branch lengths
void jbl(tree *tr, int enable) {
	if(enable) {
		printf("enabling JBL.\n");
		tr->multiBranch = 0;
		tr->numBranches = 1;

#ifdef _USE_PTHREADS
		masterBarrier(THREAD_CHANGE_NUM_BRANCHES, tr);
#endif
	}
	else {
		printf("disabling JBL(-M)\n");
		tr->multiBranch = 1;
		tr->numBranches = tr->NumberOfModels;

#ifdef _USE_PTHREADS
		masterBarrier(THREAD_CHANGE_NUM_BRANCHES, tr);
#endif
	}
}

// simply consider partitions as unlinked
mtest* simple(tree *tr, analdef *adef, linkageList *alphaList, int loops) {

	int i, model, nrAAModels = NUM_PROT_MODELS - 2;

	mtest *test = malloc(sizeof(mtest));

	// diable joint branchlengths
	jbl(tr, 0);

	test->run = (assignment**) malloc(sizeof(assignment*) * nrAAModels);
	if(loops > 0) {
		printf("loops will be reduced to %d\n", loops);
		test->nrRuns = loops;
	}
	else
		test->nrRuns = nrAAModels;
	test->nrModels = tr->NumberOfModels;

	printf("Simple Test, number of partitions: %d, available AA models: %d\n\n", tr->NumberOfModels, (int) nrAAModels);


	for (i = 0; i < test->nrRuns; i++) {
		test->run[i] = mallocAssignment(tr->NumberOfModels);

		// update models
		for (model = 0; model < tr->NumberOfModels; model++)
			test->run[i]->partitionModel[model] = i;

		evaluateAssignment(tr, adef, test->run[i], alphaList);
	}

	sortUnlinked(test);
//	printf("evaluating one assignment took %f sec. in average\n", assignmentEvalTicks/nrAAModels);
	test->result = partitionWiseOptimum(test);

	// enable joint branchlengths again
	jbl(tr, 1);

	return test;
}


// inputs: base, digits, value
// output: gray
// Convert a value to a graycode with the given base and digits.
// Iterating through a sequence of values would result in a sequence
// of Gray codes in which only one digit changes at a time.
unsigned* to_gray(unsigned base, unsigned digits, unsigned value) //, unsigned gray[digits])
{
        unsigned baseN[digits]; // Stores the ordinary base-N number, one digit per entry
        unsigned i;             // The loop variable
        unsigned *gray = malloc(sizeof(unsigned) * digits);

        // Put the normal baseN number into the baseN array. For base 10, 109
        // would be stored as [9,0,1]
        for (i = 0; i < digits; i++) {
                baseN[i] = value % base;
                value /= base;
        }

        // Convert the normal baseN number into the graycode equivalent. Note that
        // the loop starts at the most significant digit and goes down.
        unsigned shift = 0;
        while (i--) {
                // The gray digit gets shifted down by the sum of the higher
                // digits.
                gray[i] = (baseN[i] + shift) % base;
                shift += base - gray[i];        // Subtract from base so shift is positive
        }
        return gray;
}

unsigned* to_binary(unsigned base, unsigned digits, unsigned value) {
    unsigned *baseN = malloc(sizeof(unsigned) * digits);
    unsigned i;

    // Put the normal baseN number into the baseN array. For base 10, 109
    // would be stored as [9,0,1]
    for (i = 0; i < digits; i++) {
            baseN[i] = value % base;
            value /= base;
    }

    return baseN;
}


// exhaustively tests all possible model assignments
mtest* linkedExhaustive(tree *tr, analdef *adef, linkageList *alphaList, int loops) {
	int i, j, model, nrAAModels = NUM_PROT_MODELS - 2;

	unsigned *val;

	mtest *test = malloc(sizeof(mtest));

	test->method = "Exhaustive";
	test->run = (assignment**) malloc(sizeof(assignment*) * pow(nrAAModels, tr->NumberOfModels));
	// if loops is given use number of loops
	if(loops > 0)
		test->nrRuns = loops;
	else
		test->nrRuns = (int) pow(nrAAModels, tr->NumberOfModels);
	test->nrModels = tr->NumberOfModels;

	printf("Exhaustive search, number of partitions: %d, available AA models: %d, resulting combinations: %d\n\n", tr->NumberOfModels, nrAAModels, test->nrRuns);

	for (i = 0; i < test->nrRuns; i++) {
		test->run[i] = mallocAssignment(tr->NumberOfModels);
		val = to_gray(nrAAModels, test->nrModels, i);

		for(model = 0; model < test->nrModels; model++)
			test->run[i]->partitionModel[model] = val[model];

		evaluateAssignment(tr, adef, test->run[i], alphaList);
//		if(i % (2 * NumberOfThreads) == 0) printTime(i * 100 / test->nrRuns);
	}
	printf("time needed for optimization %f\n", assignmentEvalTicks);

//	printf("\nEvaluating one assignment took %f sec", assignmentEvalTicks);
	test->result = overallOptimum(test);
	return test;
}


// test some random assignments
mtest* randomTest(tree *tr, analdef *adef, linkageList *alphaList, int loops) {
	int i, model, nrAAModels = NUM_PROT_MODELS - 2;

	long randomSeed = 54321;

	mtest *test = malloc(sizeof(mtest));

	assignment *s = mallocAssignment(tr->NumberOfModels);

	for (model = 0; model < tr->NumberOfModels; model++)
		s->partitionModel[model] = -1;

	test->method = "RANDOM";
	test->run = (assignment**) malloc(sizeof(assignment*) * loops);
	test->nrRuns = loops; test->nrModels = tr->NumberOfModels;

	printf("Random Test, number of partitions: %d, number of iterations: %d\n\n", tr->NumberOfModels, loops);

	for (i = 0; i < loops; i++) {
		test->run[i] = mallocAssignment(tr->NumberOfModels);

		test->run[i] = mutate(s, test->nrModels, test->nrModels);

		// update models
		evaluateAssignment(tr, adef, test->run[i], alphaList);
//		printf("One Assignment took %f sec.\n", (assignmentEvalTicks / i));
//		printf("%d%% processed. Overall time estimate is %f sec. One Assignment took %f sec...\r", (int) i * 100 / loops, (assignmentEvalTicks * 100 / (i * 100 / loops)), (assignmentEvalTicks / i));
//		fflush(stdout);
//		if(i % (2 * NumberOfThreads) == 0)
//		printTime(i * 100 / loops);
	}
	printf("\nEvaluating one assignment took %f sec", assignmentEvalTicks/loops);
	test->result = overallOptimum(test);
	return test;
}

mtest* crossover(mtest *init, int npar) {
	int i, j, model, p, par;
	long randomSeed = rand();
	mtest *res = malloc(sizeof(mtest));
	res->run = (assignment**) malloc(sizeof(assignment*) * init->nrRuns);	
	res->nrModels = init->nrModels; res->nrRuns = init->nrRuns;


	for(i = 0; i < init->nrRuns; i++) {
		p = (int)(randum(&randomSeed) * init->nrModels);
		par = (int)(randum(&randomSeed) * npar);
		res->run[i] = mallocAssignment(init->nrModels);

		for(j = 0; j < p; j++) {
			res->run[i]->partitionModel[j] = init->run[par]->partitionModel[j];
		}

		par = (int)(randum(&randomSeed) * npar);
		for(;j < init->nrModels; j++) {
			res->run[i]->partitionModel[j] = init->run[par]->partitionModel[j];
		}
	}

	free(init);
	return res;
}

// starts with an arbitrary state and optimizes until no optimization is possible anymore
mtest* hillClimbing(tree *tr, analdef *adef, linkageList *alphaList, assignment *start) {
	int i, j, model, nrAAModels = NUM_PROT_MODELS - 2, changed = 1;

	mtest *test = malloc(sizeof(mtest));
	assignment *s = mallocAssignment(tr->NumberOfModels), *opt = mallocAssignment(tr->NumberOfModels);

	test->method = "Hill Climbing";
	test->run = (assignment**) malloc(sizeof(assignment*) * nrAAModels * tr->NumberOfModels);
	test->nrRuns = nrAAModels * tr->NumberOfModels;
	test->nrModels = tr->NumberOfModels;

	printf("Hill Climbing, number of partitions: %d, available AA models: %d, resulting combinations: %d\n\n", tr->NumberOfModels, nrAAModels, test->nrRuns);

	for (model = 0; model < tr->NumberOfModels; model++)
		s->partitionModel[model] = -1;

	opt = mutate(s, test->nrModels, test->nrModels);

	evaluateAssignment(tr, adef, opt, alphaList);

	for (model = 0; model < test->nrModels; model++) {
		for (j = 0; j < nrAAModels; j++) {
			test->run[model * nrAAModels + j] = mallocAssignment(test->nrModels);
			for(i=0; i < test->nrModels; i++) {
				if(model != i) test->run[model * nrAAModels + j]->partitionModel[i] = opt->partitionModel[i];
				else test->run[model * nrAAModels + j]->partitionModel[i] = j;
			}
			evaluateAssignment(tr, adef, test->run[model * nrAAModels + j], alphaList);
			if(test->run[model * nrAAModels + j]->overallLH > opt->overallLH) {
				opt = test->run[model * nrAAModels + j];
				changed = 1;
			}
		}
	}

	test->result = opt;

	return test;
}

// evaluates only the last partition of the partition file under the first given models
mtest* greedy(tree *tr, analdef *adef, linkageList *alphaList, assignment* start) {
	int i, j, model, nrAAModels = NUM_PROT_MODELS - 2;

	mtest *test = malloc(sizeof(mtest));

	test->method = "Greedy";
	test->run = (assignment**) malloc(sizeof(assignment*) * nrAAModels);
	test->nrRuns = (int) nrAAModels;
	test->nrModels = tr->NumberOfModels;

	printf("Part of greedy assignment building...\n");

	for (j = 0; j < nrAAModels; j++) {
		test->run[j] = mallocAssignment(tr->NumberOfModels);

		for (model = 0; model < test->nrModels; model++) {
			if(model < test->nrModels - 1) test->run[j]->partitionModel[model] = tr->partitionData[model].protModels;
			else  test->run[j]->partitionModel[model] = j;
		}
		evaluateAssignment(tr, adef, test->run[j], alphaList);
	}

	test->result = overallOptimum(test);

	return test;
}


mtest* simmulatedAnnealing(tree *tr, analdef *adef, linkageList *alphaList, assignment * start, int loops) {
	int model, i, j, numAccepted = 0;

	mtest *test = malloc(sizeof(mtest)), tmp;
	assignment *y;

	//init
	test->method = "SA";
	test->run = (assignment**) malloc(loops * sizeof(assignment*));
	test->run[0] = start; //mutate(start, tr->NumberOfModels, tr->NumberOfModels);
	test->nrModels = tr->NumberOfModels; test->nrRuns = loops;

	printf("Simmulated annealing: loops is %d....\n", loops);

	evaluateAssignment(tr, adef, test->run[0], alphaList);
	printf("started with: \n");
	printAssignment(test->run[0], test->nrModels);
	printf("\n");

	for(i = 0; i < loops - 1; i++) {
		y = mallocAssignment(test->nrModels);
		y = mutate(test->run[i], test->nrModels, 1);
		evaluateAssignment(tr,adef, y, alphaList);

//		printf("diff: %f  x: %f  y: %f  P: %f  RAND: %f\n", y->overallLH - test->run[i]->overallLH, test->run[i]->overallLH, y->overallLH, exp((y->overallLH - test->run[i]->overallLH) / (loops - i)), rand() / (double)RAND_MAX);
//		printAssignment(y, test->nrModels);
//		printf("diff: %f P: %f\n", y->overallLH - test->run[i]->overallLH, exp((y->overallLH - test->run[i]->overallLH) / (loops - i)));

		if(1){ //exp((y->overallLH - test->run[i]->overallLH) / (loops - i)) > rand() / (double)RAND_MAX) {
			numAccepted++;
			test->run[i + 1] = y;
		} else {
			test->run[i + 1] = test->run[i];
		}
	}

	printf("\naccepted branch was chosen %d times\n", numAccepted);
	test->result = overallOptimum(test);
	return test;
}


mtest* geneticAlgo(tree *tr, analdef *adef, linkageList *alphaList, int loops, int popcount) {
	int model, i, j, fraq = 20, num; // loops = 5;

	mtest *test = malloc(sizeof(mtest)), *population;
	assignment *tmp = mallocAssignment(test->nrModels);

	num = (int) (popcount * fraq / 100);

	printf("Genetic algorithm: %d of each population with %d assignments will breed next population %d times...\n"
			"in summary %d assignments will be evaluated\n", num, popcount, loops, loops * popcount);

	test->method = "GA";
	test->nrModels = tr->NumberOfModels; test->nrRuns = popcount * loops; //num * loops;
	// reservers amount of memory needed, to store loops times best breed plus some to seperate them
	test->run = (assignment**) malloc(sizeof(assignment*) * (popcount * loops)); //(num * loops));

	population = malloc(sizeof(mtest));
	population->run = (assignment**) malloc(sizeof(assignment*) * popcount);
	population->nrRuns = popcount;

	for(j = 0; j < test->nrModels; j++) {
		tmp->partitionModel[j] = -1;
		tmp->partitionLH[j] = unlikely;
	}

	for(j = 0; j < popcount; j++)
		population->run[j] = mutate(tmp, test->nrModels, test->nrModels);

	population->nrRuns = popcount;
	population->nrModels = test->nrModels;

	//loops times redo population generation
	for(i = 0; i < loops; i++) {

		population = crossover(population, num);
		//mutate and evaluate
		for(j = 0; j < popcount; j++) {
			population->run[j] = mutate(population->run[j], test->nrModels, 1);
			evaluateAssignment(tr, adef, population->run[j], alphaList);
		}

		sort(population);
		for(j = 0; j < popcount; j++) {
			if(j < num) population->run[j]->chosen = 1;
			test->run[i * popcount + j] = population->run[j];
		}

//		for(j = 0; j < num; j++)
//			test->run[i * num + j] = population->run[j];
	}

	test->result = overallOptimum(test);
	return test;
}

mtest* finishUnlinked(mtest *t, mtest *y) {
	int i, model;
	double sum;

	mtest *tUnlinked;

	tUnlinked  = (mtest*) malloc(sizeof(mtest));
	tUnlinked->method = "unlinked";
	tUnlinked->run = (assignment**) malloc(sizeof(assignment*) * t->nrRuns);
	tUnlinked->nrRuns = t->nrRuns; tUnlinked->nrModels = t->nrModels;

	for(i = 0; i < t->nrRuns; i++) {
		tUnlinked->run[i] = mallocAssignment(t->nrModels);
		for(model = 0; model < t->nrModels; model++) {
			tUnlinked->run[i]->partitionModel[model] = t->run[i]->partitionModel[model];
			// get the LH(model, partition) of unlinked run
			tUnlinked->run[i]->partitionLH[model] = y->run[t->run[i]->partitionModel[model]]->partitionLH[model];
		}
	}

//	finish unlinked LH values
	for(i = 0; i < tUnlinked->nrRuns; i++) {
		sum = 0;
		for(model = 0; model < tUnlinked->nrModels; model++)
			sum += tUnlinked->run[i]->partitionLH[model];
		tUnlinked->run[i]->overallLH = sum;
	}

	return tUnlinked;
}


void modOptJoerg(tree *tr, analdef *adef, rawdata* rdta) {
	int i, model, *unlinked = (int *) malloc(sizeof(int) * tr->NumberOfModels), loops;
	double sum;

	unsigned *gray;

	printf("starting to test models \n");
	if(adef->protEmpiricalFreqs) protEmpiricalFreqs = 1;


	assignment *tmp, *start;
	mtest *t, *y, *tUnlinked;
	linkageList *alphaList;

	/* this unlinked thing here just tells RAxML that the other relevant model
	 parameter, the alpha shape parameter will be estimated separately for each
	 partition */
	for (model = 0; model < tr->NumberOfModels; model++)
		unlinked[model] = model;

	/* initialize the linkage list for the alpha parameter of rate heterogeneity */
	alphaList = initLinkageList(unlinked, tr);

	/* just an arbitrarily selected node of the tree where we will start
	 the recursion of the Felsenstein pruning algorithm */
	tr->start = tr->nodep[1];

	// example to switch -M on the fly
//	if(1) {
//		if(tr->multiBranch == 0)
//			printf("started for JBL\n");
//		else
//			printf("started with -M\n");
//
//		printf("run JBL:\n");
//		jbl(tr, 1);
//		t = simple(tr, adef, alphaList, 10);
//		printModelTest(t);
//		printf("\n");
//
//		printf("run -M:\n");
//		jbl(tr, 0);
//		t = simple(tr, adef, alphaList, 10);
//		printModelTest(t);
//		printf("\n");
//
//		exit(0);
//	}

	/* start testing protein model assignments */
//	if (tr->allCombinations) {
		loops = 100;
		printf("starting model search... \n");

//		create assignment of the models given in partition file
	start = mallocAssignment(tr->NumberOfModels);
	for(model = 0; model < tr->NumberOfModels; model++) {
		start->partitionModel[model] = tr->partitionData[model].protModels;
	}

	switch(tr->modelAssignment)
	{
		case EXHAUSTIVE:
			t = linkedExhaustive(tr, adef, alphaList, loops);
			break;
		case RANDOM:
			t = randomTest(tr, adef, alphaList, loops);
			break;
		case GA:
		{
			int popcount = 10;
			t = geneticAlgo(tr, adef, alphaList, (int) (loops / popcount), popcount);
		} break;
		case SA:
			t = simmulatedAnnealing(tr, adef, alphaList, y->result, loops);
			break;
		case HILL:
			t = hillClimbing(tr, adef, alphaList, start);
			break;
		case GREEDY:
			t = greedy(tr, adef, alphaList, start);
			break;
		case NAIVE:
			t = simple(tr, adef, alphaList, 0);
			break;
		default:
			errorExit(-1);
			break;
	}

	// if specified compare search path with unlinked lh values
//	if(adef->perGeneBranchLengths) {
//		printf("running initial -M\n");
//		y = simple(tr, adef, alphaList, 0);
//		printModelTest(y);
//		printf("\n");
//		tUnlinked = finishUnlinked(t, y);
//		printSearch(t, tUnlinked);
//		printModelTestFile(tUnlinked, "_unlinked");
//	}
//	else {
//		printf("simple heuristic called\n", tr->NumberOfModels);
//		t = simple(tr, adef, alphaList, 0);
//		printf("creating binary file testDump\n");
//		writeTest(t);
//		sortUnlinked(t);
//	}

	printModelTestFile(t, "");
	printAssignmentFile(t->result, t->nrModels);

	free(unlinked); free(t);
	freeLinkageList(alphaList);
}

#endif


