#include <stdio.h>
#include <math.h>

int main(void)
{
      int N = 32;
      double d2dx2[N*N], psi[N], x[N], energy[N];
      double a = -3.0, b = 3.0;
      double h = (b-a)/(double)(N);
      double coeff[] = {-205/72.0, 8/5.0, -1/5.0, 8/315.0, -1/560.0};
      int i, j, info;
      FILE *file;
      file = fopen("eign", "wt");
      //construct the 9 diagonal derivative matrix
      for(i=0; i<N; i++)
      {
            for(j=0; j<N; j++)
            {
                  if(abs(j-i)<5 && (i>3) && (i<(N-4)))
                        d2dx2[i+N*j] = - coeff[abs(i-j)]/pow(h, 2);
                  else
                        d2dx2[i+N*j] = 0;
            }
      }
      //construct the (1-x^2) function
      for(i=0; i<N; i++)
      {
            x[i] = 1 - pow(a+i*h,2);
      }

      spotrf('U', N, &x[0], N, &info);

      dsygst(1, 'U', N, &d2dx2[0], N, &x[0], N, &info);
      
      return 0;
}
