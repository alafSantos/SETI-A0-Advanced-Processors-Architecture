#define MAIN
#include <stdio.h>
#include "cache.h"

#ifndef N
#  define N 64	/* Taille  des matrices et vecteurs*/
#endif

#if CUSTOM_MEMORY_LAYOUT
double tt[N*N+2*N];
double *x,*y,*z,s=0;
#else
double x[N*N], y[N], z[N], s=0;
#endif 
FILE *results;



int main()
{
  int i,j,k;
  results = fopen( "res-matvec", "w+" );

  /** Initialisation des matrices**/
#if CUSTOM_MEMORY_LAYOUT
  x=&tt[0]; y=&tt[N*N]; z=&y[N];
#endif
 
  for (i=0;i<N;i++)
    {	
      y[i]=i;
      for (j=0;j<N;j++)
	x[i*N+j]=i+j;
    }
  
 fprintf(results, "PRODUIT MATRICE VECTEUR (N=%d)\n \n",N);
 fprintf(results, "&x[0]=%llx \t &y[0] =%llx \n",(unsigned long long)&x[0],(unsigned long long) &y[0]);

 fprintf(results, "SIZE\tWAY\tLINE\tN\tEchec\tOblig.\tConfl.\tCapac.\n");

   LOOP_CACHE_CONFIG
     {
       initcache();
       
       for (i=0;i<N;i++)
	 {
	   s=0;
	   for (k=0;k<N;k++)
	     {
	       ac(&x[i*N+k]);ac(&y[k]);
	       s = s + x[i*N+k]*y[k];
	     }
	   ac(&z[i]);
	   z[i]=s;
	 }
	 
       print_cache_stats(results, N);

     }
 
}

