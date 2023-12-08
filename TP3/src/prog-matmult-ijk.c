#define MAIN
#include <stdio.h>
#include "cache.h"

#ifndef N
#  define N 64	/* Taille  des matrices et vecteurs*/
#endif

#if CUSTOM_MEMORY_LAYOUT
double tt[3*N*N], *x, *y, *z, s=0;
#else
double x[N*N], y[N*N], z[N*N], s=0;
#endif

FILE *results;



int main()
{
  int i,j,k;
  results = fopen( "res-matmult-ijk", "w+" );

  /** Initialisation des matrices**/

#if CUSTOM_MEMORY_LAYOUT
  x=&tt[0]; y=&x[N*N]; z=&y[N*N];
#endif
  
  fprintf(results, "Multiplication de matrices ijk   (N=%d)\n \n",N);
  fprintf(results, "&x[0]=%p \t &y[0] =%p &z[0]=%p \n",&x[0],&y[0],&z[0]);
  fprintf(results, "SIZE\tWAY\tLINE\tN\tEchec\tOblig.\tConfl.\tCapac.\n");


  for (i=0;i<N;i++)
    {for (j=0;j<N;j++)
	{
	  x[i*N+j]=i+j;
	  y[i*N+j]=i-j;
	}
    }

   LOOP_CACHE_CONFIG
	{
	  initcache();
	  
	  for (i=0;i<N;i++)
	    {
	      for (j=0;j<N;j++)
		{
		  s=0;
		  for (k=0;k<N;k++)
		    {
		
		      ac(&x[i*N+k]);ac(&y[k*N+j]);
		      s = s + x[i*N+k]*y[k*N+j];
		    }
		  ac(&z[i*N+j]);
		  z[i*N+j]=s;
		
		}
	    }
		

          print_cache_stats(results,N);
	}

}
