#ifndef __MYMVTNORM_H__
#define __MYMVTNORM_H__

#include <cmath>
#include <iostream>
#include <stdio.h>
#include "MathMatrix.h"
#include "StringHash.h"
#define MATHLIB_STANDALONE 
#include <Rmath.h>

void GetConditionalDist(Vector & score,Matrix & cov,IntArray & comb, Vector & result)
{
   result.Dimension(2);
   result.Zero();

   //chop score vector and sigma matrix
   Vector mu2;
   Matrix sub_cov;
   int dim = comb.Length()-1;
   mu2.Dimension(dim);
   sub_cov.Dimension(dim,dim);
   for(int i=0;i<dim;i++)
   {
      int idx1 = comb[i+1];
      mu2[i] = score[idx1];
      for(int j=0;j<dim;j++)
      {
         int idx2 = comb[j+1];
         sub_cov[i][j] = cov[idx1][idx2];
      }
   }

   Matrix inv;
   inv.Copy(sub_cov);
   SVD svd;
   svd.InvertInPlace(inv);
   Vector sigma12;
   sigma12.Dimension(dim);
   for(int i=0;i<dim;i++)
   {
      int idx1 = comb[0];
      int idx2 = comb[i+1];
      sigma12[i] = cov[idx1][idx2];
   }
   Vector tmp;
   tmp.Dimension(dim);
   tmp.Zero();
   for(int i=0;i<dim;i++)
   {
      tmp[i] += sigma12.InnerProduct(inv[i]);
   }

   //Calculate Mu and Sigma here
   result[0] = tmp.InnerProduct(mu2);
   result[1] = 1.0 - tmp.InnerProduct(sigma12);

   if(result[1]<0)
      result[1] = fabs(result[1]);
}

double CalculateMVTPvalue(Vector & score,Matrix & cov,double t_max)
{
   double pvalue=0.0;
   int dim = score.Length();
   double chisq = t_max * t_max;

   StringDoubleHash jointProbHash;

   if(dim==1)
   {
      pvalue = pchisq(chisq,1,0,0);
      //printf("pvalue is: %g\n",pvalue);
      return pvalue;
   }
   double uni = pchisq(chisq,1,0,0);
   pvalue += dim*uni;
   //printf("pvalue so far is:%g\n",pvalue);

   //Initializations for combination
   vector<int> indx;
   vector<int> alpha (dim,0);
   for(int i=0;i<dim;i++)
      alpha[i] = i;
   int n=dim;
   for(int r=2;r<=dim;r++)
   {
      int j=r;
      int k=r;

IntArray comb;
      Vector par;
      par.Dimension(2);

      for(int twk=j;twk<=k;twk++)
      {
         int r=twk;
         bool done=true;
         for(int iwk=0;iwk<r;iwk++)indx.push_back(iwk);
         while(done)
         {
            done=false;
            for(int owk=0;owk<r;owk++)
            {
               //A new combination has been generated!
               comb.Push(alpha[indx[owk]]);
            }
/*
            printf("New combination:\n");
            for(int i=0;i<comb.Length();i++)
               printf("%g ",comb[i]);
            printf("\n");
*/
            //Calculate P(A|BC)
            GetConditionalDist(score,cov,comb,par);
//            printf("pars are:%g %g\n",par[0],par[1]);
            double chisq,condProb,prob;
            if(par[1]==0.0)
            {
               condProb = 0.0;
            }
            else
            {
               chisq = (t_max-par[0])*(t_max-par[0])/par[1];
               if(chisq<0.0)
                  chisq = -chisq;
               condProb = pchisq(chisq,1,0,0);
            }
//          printf("conditional prob is:%g\n",condProb);
            String hashKey;

	    if(r==2)
	    {
	       int tmp = comb[0];
	       hashKey += tmp;
	       tmp = comb[1];
	       hashKey += tmp;
	       //              printf("hashKey is: %s\n",hashKey.c_str());
	       prob = condProb*uni;
	       //prob = exp(log(condProb) + log(uni));
	       //printf("joint prob is:%g\n",prob);
	       jointProbHash.SetDouble(hashKey,prob);
	       hashKey.Clear();
	    }
	    else
	    {
	       for(int i=1;i<r;i++)
	       {
		  int tmp = comb[i];
		  hashKey += tmp;
	       }
	       //              printf("hashKey is: %s\n",hashKey.c_str());
	       prob = jointProbHash.Double(hashKey);
//	       printf("joint prob is:%g\n",prob);
	       //if(prob == _NAN_)
	       //{
	       //         printf("Error! Joint probability %s is not hashed in properly.\n",hashKey);
	       //     }
	       prob *= condProb; //This is the new joint prob.
	       String newKey;
	       int tmp = comb[0];
	       newKey += tmp;
	       newKey += hashKey;
	       //              printf("newKey is: %s\n",newKey.c_str());
	       jointProbHash.SetDouble(newKey,prob);
	       hashKey.Clear();
	    }

	    pvalue -= prob;
//	    printf("pvalue so far is:%g\n",pvalue);
	    comb.Clear();
	    for(int iwk=r-1;iwk>=0;iwk--)
	    {
	       if(indx[iwk]<=(n-1)-(r-iwk))
	       {
		  indx[iwk]++;
		  for(int swk=iwk+1;swk<r;swk++)
		  {
		     indx[swk]=indx[swk-1]+1;
		  }
		  iwk=-1;
		  done=true;
	       }
	    }
	 }
	 indx.clear();
      }
   }
   // printf("pvalue is:%g\n",pvalue);
   return pvalue;
}

#endif
