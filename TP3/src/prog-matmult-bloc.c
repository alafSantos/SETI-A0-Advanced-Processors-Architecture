#define MAIN
#include <stdio.h>
#include "cache.h"
#ifndef N
#  define N 500	/* Taille  des matrices et vecteurs*/
#endif


#if CUSTOM_MEMORY_LAYOUT
double tt[3*N*N];
double *x,*y,*z,s=0;
#else
double x[N*N], y[N*N], z[N*N], s=0;
#endif 


FILE *results;

static int min(int a, int b){
  return a<b ? a : b ;
}

int main()
{
  int i,j,k,jj,kk,B;
  results = fopen( "res-matmult-bloc", "w+" );

#if CUSTOM_MEMORY_LAYOUT
  x=&tt[0]; y=&x[N*N]; z=&y[N*N];
#endif
  
  /** Initialisation des matrices**/
  fprintf(results, "Multiplication de matrices par blocs\n");
  fprintf(results, "&x[0]=%p \t &y[0] =%p &z[0]=%p \n",&x[0],&y[0],&z[0]);
  fprintf(results, "SIZE\tWAY\tLINE\tN\tB\tEchec\tOblig.\tConfl.\tCapac.\n");
  for (i=0;i<N;i++)
    {for (j=0;j<N;j++)
	{
	  x[i*N+j]=i+j;
	  y[i*N+j]=i-j;
	}
    }
  size = 16*1024;// taille cache fixée

  for (line=16;line <100;line+=line)
    for (way=1;way<8;way+=way)
      for (B=4;B<130;B+=B)
	{
	  initcache();
	  for (jj=0;jj<N;jj+=B)
	    for (kk=0;kk<N;kk+=B)
	      for (i=0;i<N;i++)
		{
		  for (j=jj;j<min(jj+B-1,N);j++)
		    {
		      s=0;
		      for (k=kk;k<min(kk+B-1,N);k++)
			{
			  ac(&x[i*N+k]);ac(&y[k*N+j]);
			  s = s + x[i*N+k]*y[k*N+j];
			}
		      ac(&z[i*N+j]);ac(&z[i*N+j]);
		      z[i*N+j]=s+z[i*N+j];
		    }
		}

		
	  fprintf(results,"%d\t%d\t%d\t%d\t%d\t%f\t%f %f %f\n",size,way, line, N,B,dc/(double)am,dcob/(double)am,dcco/(double)am,dcca/(double)am);

	}

}
