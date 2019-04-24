#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "timing.h"
#define n 1500
#define N n*n /* matrix size */
#define E_h 300
#define E_d 500
#define NONE 0
#define H 1



int main (int argc, char *argv[])
{
int i, j, k;
timing_t tstart, tend;
double  M[N][N+1],           /* matrix A to be multiplied */
b[N][N],           /* matrix B to be multiplied */
c[N][N];           /* result matrix C */

get_time(&tstart);

C = (N * E_h) + E_d

for (i=0; i<N; i++)
 for (j=0; j<N+1; j++)
   if (j == N) {
     M[i][j]= E_d / C;
   }
   M[i][j]= E_h / C;


for (i=0; i<n; i++)
 for (j=0; j<n; j++)
   b[i][j] = NONE;



for (i=0; i<N; i++)
  {
  for(j=0; j<N; j++)
    for (k=0; k<N; k++)
      c[i][j] += a[i][k] * b[k][j];
  }

get_time(&tend);


 printf("Elapsed time: %g s\n", timespec_diff(tstart, tend));

  printf("*****************************************************\n");

  printf("Result Matrix:\n");
  for (i=0; i<N; i++)
  {
    for (j=0; j<N; j++)
      printf("%6.2f   ", c[i][j]);
    printf("\n");
  }
  printf("******************************************************\n");
  printf ("Done.\n");
printf("Elapsed time: %g s\n", timespec_diff(tstart, tend));
}



int add_particle (double probs[]) // Return index of next particle to add
{
  int i, double cumprobs[];
  for (i=0; i<sizeof(probs); i++) {
    cumprobs[i] = probs[i]; // Place in current probability value at this index
    if (i>0)
      cumprobs[i] += cumprobs[i-1]; // Accumulate probability
  }

  randn = (double) rand() / (double) RAND_MAX; // Random number from 0 to 1
  for (i=1, i<sizeof(probs); i++) {
    if (randn >= cumprobs[i]) & (randn < cumprobs[i+1])
      return i;
  }

}








