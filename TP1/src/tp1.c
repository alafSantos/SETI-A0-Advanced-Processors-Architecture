#include <stdio.h>
#include <stdint.h>


#ifndef N
#define N 1000		// la taille des matrices [N][N]/vecteurs[N^2]


#endif
#ifndef M
#define M 16		// nombre d'itérations de chaque fonction
#endif
#ifndef TYPE
#define TYPE float	// ou double ou int, etc. On peut essayer les
			// différents types de données
#endif

// pour éviter des conversions int-float inutiles
#if TYPE==int || TYPE==short || TYPE==char || TYPE==uint64_t
#  define ZERO 0
#elif  TYPE==float
#  define ZERO 0.0f
#else
#  define ZERO 0.0
#endif


#ifndef BL
#define BL 16		// La taille des blocs pour la mult de matrices par blocs
#endif

#define STR1(x) #x
#define STR(x) STR1(x)

double dt ; // le nombre de cycles écoulés en double

// our functions
double min_res(double nb_iter);
double median_res(double nb_iter);

// https://www.tutorialspoint.com/learn_c_by_examples/median_program_in_c.htm
void swap(double *p,double *q) {
   double t;
   t=*p; 
   *p=*q; 
   *q=t;
}

void sort(double a[],int n) { 
   int i,j,temp;
   for(i = 0;i < n-1;i++) {
      for(j = 0;j < n-i-1;j++) {
         if(a[j] > a[j+1])
            swap(&a[j],&a[j+1]);
      }
   }
}




// Permet de lire le compteur interne du nombre de cycles du
// processeur (TSC timestamp counter).
// Le TSC renvoie le nombre de cycles codé sur 64 bits (long long)
// start_timer() initie la mesure et stop_timer() la termine
// dtime() convertit en double la différence pour simplifier les calculs ultérieurs
//      et retranche le temps (approximatif) bécessaire à la mesure du TSC

static unsigned long long start_timer() ;
static unsigned long long stop_timer() ;
static double dtime(long long debut, long long fin);


// le fichier où sont les fonctions d'accès au TSC
#include "tsc.h"


TYPE AF[N][N], YF[N][N], XF[N][N], YT[N][N],	// matrices [N][N]
  BF[N*N], CF[N*N],				// vecteurs N^2
  SF;						// accumulateur

double resultats[M];                            // les temps de calcul en double

static inline int min (int a,int b)
{
  if (a<b) return a;
  else return b; 
}

long long debut, fin;
double benchtime;

// pour les calculs, on met les temps obtenu pour chaque essai dans un tableau.
// Le contenu du tableau est ensuite affiché et/ou traité.
void add_res(double res, int where)
{
  resultats[where]=res;
}

//Gives min
double min_res(double nb_iter){
  double normalize = 1.0/nb_iter;
  double tmp=resultats[0]*normalize;
  for(int i=1; i<M; i++)
    tmp=(tmp<=resultats[i]*normalize)?tmp:resultats[i]*normalize;
  return tmp;
}

double median_res(double nb_iter){
  sort(resultats,M);
  int n = (M+1) / 2 - 1;
  return (double) resultats[n]/nb_iter;
}

// la fonction d'affichage des résutats. Peut être redéfinie en fonction des besoins
// pour faire des statistiques, enlever les point aberrants, etc
// Le deuxième paramètre est le nombre d'itérations pour normaliser
void print_res(char *funcname, double nb_iter){
  double normalize = 1.0/nb_iter;
  printf("%s\t", funcname);
  for(int i=0; i<M; i++)
    printf("%g%s",resultats[i]*normalize, (i==M-1 ? "\n" : "\t"));
}

// entre deux mesures d'un même programme. Rédéfinir si nécessaire
void separateur(){
  printf("\n");
}

//Mise à zéro d'un vecteur
void zero(){

  int i, j, m ;
  
  for (m=0;m<M;m++)
    {
      debut=start_timer();
      for (i=0;i<(N*N);i++)
	      BF[i]=ZERO;
      benchtime=dtime(debut, stop_timer());
      add_res(benchtime,m);
    }
  print_res ("ZERO",(N*N));
  separateur();
}

//copie de matrices
void copy_ij(){
  int i, j, m ;
  
  for (m=0;m<M;m++)
    {
      debut=start_timer();

      for (i=0;i<N;i++)
	for (j=0;j<N;j++)
	  AF[i][j]=YF[i][j];
      benchtime=dtime(debut, stop_timer());
      add_res(benchtime,m);
    }
  print_res ("COPY_ij", (N*N));
  separateur();
}

void copy_ji(){
  int i, j, m ;
  
  for (m=0;m<M;m++)
    {
      debut=start_timer();

      for (j=0;j<N;j++)
	for (i=0;i<N;i++)
	  AF[i][j]=YF[i][j];
      benchtime=dtime(debut, stop_timer());
      add_res(benchtime,m);
    }
  print_res ("COPY_ji", (double)(N*N));
  separateur();
}

// addition de deux matrices
void add_ij(){
  int i, j, m ;

  for (m=0;m<M;m++)
    {
      debut=start_timer();

      for (i=0;i<N;i++)
	for (j=0;j<N;j++)
	  AF[i][j]+=YF[i][j];
      benchtime=dtime(debut, stop_timer());
      add_res(benchtime,m);
    }
  print_res ("ADD_ij", (double) (N*N));
  separateur();
}

void add_ji(){
  int i, j, m ;

  for (m=0;m<M;m++)
    {
      debut=start_timer();

      for (j=0;j<N;j++)
	for (i=0;i<N;i++)
	  AF[i][j]+=YF[i][j];
      benchtime=dtime(debut, stop_timer());
      add_res(benchtime,m);

    }
  print_res ("ADD_ji", (double) (N*N));
  separateur();
}


// Produit scalaire 
void ps()
{
  int i, j, k, m;
  
  for (m=0;m<M;m++)
    {
      debut=start_timer();

      SF=ZERO;
      for (i=0;i<(N*N);i++){
	SF+=BF[i]*CF[i]; 
      }
      benchtime=dtime(debut, stop_timer());
      add_res(benchtime,m);

    }
  print_res ("PS", (double) (N*N));
  separateur();
}



// Multiplication de matrices ijk
void mm_ijk()
{
  int i, j, k, m;
  
  for (m=0;m<M;m++)
    {
      debut=start_timer();

      for (i=0;i<N;i++)
	for (j=0;j<N;j++)
	  {
	    SF=ZERO;
	    for (k=0; k<N; k++)
	      SF+=AF[i][k]*XF[k][j];
	    YF[i][j]=SF;
	  }
      benchtime=dtime(debut, stop_timer());
      add_res(benchtime,m);

    }
  print_res ("MM_ijk", ((double) N*N*N));
  separateur();
}


// Multiplication de matrices ikj
void mm_ikj()
{
  int i, j, k, m;
  for (m=0;m<M;m++)
    {
      debut=start_timer();
      for (i=0;i<N;i++)
	for (k=0;k<N;k++)
	  {
	    SF=AF[i][k];
	    for (j=0; j<N; j++)
	      YF[i][j]+=SF*XF[k][j];
	  }
      benchtime=dtime(debut, stop_timer());
      add_res(benchtime,m);

    }
  print_res ("MM_ikj", ((double) N*N*N));
  separateur();
}


// Multiplication de matrices par blocs
void mm_b_ijk(){
  int i, j, k, m, ii, jj, kk;
  
  for (m=0;m<M;m++)
    {
      debut=start_timer();
      for (jj=0;jj<N;jj+=BL)
	      for (kk=0;kk<N;kk+=BL)
	        for (i=0;i<N;i++)
	        {
	          for (j=jj;j<min(jj+BL-1,N);j++)
            {
              SF=ZERO;
              for (k=kk;k<min(kk+BL-1,N);k++)
                SF += AF[i][k]*XF[k][j];
              YF[i][j]=SF;
            }
	    }

      benchtime=dtime(debut, stop_timer());
      add_res(benchtime,m);
    }
  print_res ("MM_B_ijk", ((double) N*N*N));

  separateur();
}


int main()
{
  // Commenter et décommenter les appels de fonctions suivant les questions du TP.
  printf("Evaluation : N=%d, type="STR(TYPE)"\n",N);

  copy_ij();
  // dt=median_res(N*N);
  // printf("copy_ij: %lf \n \n",dt);
  
  copy_ji();
  // dt=median_res(N*N);
  // printf("copy_ji: %lf \n \n",dt);

  // add_ij();
  // dt=median_res(N*N);
  // printf("add_ij: %lf \n \n",dt);

  // add_ji();
  // dt=median_res(N*N);
  // printf("add_ji: %lf \n \n",dt);

  // ps();
  // dt=median_res(N*N);
  // printf("ps: %lf \n \n",dt);

  mm_ijk();
  // dt=min_res(N*N*N);
  // printf("mm_ijk: %lf \n \n",dt);

  mm_ikj();
  // dt=min_res(N*N*N);
  // printf("mm_ijk: %lf \n \n",dt);

  mm_b_ijk();
  // dt=min_res(N*N*N);
  // printf("mm_b_ijk: %lf \n \n",dt);

    // Commenter et décommenter les appels de fonctions suivant les questions du TP.
  // printf("Evaluation : N=%d, type="STR(TYPE)"\n",N);
  //zero();
  //copy_ij();
  //copy_ji();
  //add_ij();
  //add_ji();
  //ps();
  //mm_ijk();
  // mm_ikj();
  //mm_b_ijk();

  return 0;
} 

