#include <stdio.h>
#include <stdlib.h> /* because of the constant EXIT_FAILURE*/
#include <time.h>

#include "OMs.h"

// reads representatives from the Finschi's list and stores it in our format.

int R; // rank
int N; // the number of elements

int B;       // the number of bases
int nr_ints; // the number of integers needed to store the plus (resp. minus) of
             // a chirotope

int main() {

  R = 2;
  N = 8;
  makebases();

  int i, j, k;

  FILE *in, *out;
  char text[300];

  sprintf(text, "Finschi_%d_%d.txt", R,
          N); // read uniform representatives given in the Homepage of oriented
              // matroids by Finschi
  in = fopen(text, "r");
  if (in == NULL) {
    fprintf(stderr, "error fopen():  Could not open the file %s.\n", text);
    exit(EXIT_FAILURE);
  }

  sprintf(text, "uniform_representatives_rank%d_%delements.txt", R,
          N); // the output gets stored here
  out = fopen(text, "w");
  if (out == NULL) {
    fprintf(stderr, "error fopen():  Could not open the file %s.\n", text);
    exit(EXIT_FAILURE);
  }

  fprintf(out,
          "All uniform representatives (reorientations+permutations) of rank "
          "%d on %d elements:\n",
          R, N);

  char b[B][R];

  int ubound, lbound;

  lbound = 12; // the first lbound characters of every line should be ignored
  ubound = lbound + B;

  for (i = 0; i < R;
       i++) // read the bases, as written in the Homepage of oriented matroids
  {
    fgets(text, 300, in);
    for (j = lbound; j < ubound; j++)
      b[j - lbound][i] = text[j] - '0' - 1;
  }

  for (i = 0; i < B; i++) {
    for (j = 0; j < R; j++)
      printf("%d", b[i][j]);
    putchar('\n');
  }

  OM X = makeOM();

  while (fgets(text, 300, in) != NULL) // translate every OM
  {
    for (i = 0; i < nr_ints; i++) {
      X.plus[i] = 0;
      X.minus[i] = 0;
    }

    for (j = 0; j < B; j++) {
      k = ind(b[j]);
      if (text[j + lbound] == '+')
        X.plus[k >> 5] += 1 << (k & 31);
      else if (text[j + lbound] == '-')
        X.minus[k >> 5] += 1 << (k & 31);
      putchar(text[j + lbound]);
    }
    putchar('\n');

    writeOM(X, out);
  }

  fclose(in);
  fclose(out);

  return 0;
}
