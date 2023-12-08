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
  results = fopen( "res-matmult-ikj", "w+" );

#if CUSTOM_MEMORY_LAYOUT
  x=&tt[0]; y=&x[N*N]; z=&y[N*N];
#endif

  /** Initialisation des matrices**/

  fprintf(results, "Multiplication de matrices IKJ  (%d) \n \n",N);
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
	      for (k=0;k<N;k++)
		{  ac(&x[i*N+k]);
		  s= x[i*N+k];
		  for (j=0;j<N;j++)
		    {
		      ac(&z[i*N+j]);ac(&y[k*N+j]);
		      z[i*N+j]+=s*y[k*N+j];
		      ac(&z[i*N+j]);	
		    }
		}
	    }
		

          print_cache_stats(results,N);

	}

}
