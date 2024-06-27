#include </home/mi/palic/OMs.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// This code takes all oriented matroids in the lower cones of uniform oriented
// matroids that are representatives of reorientation and permutation classes,
// as found by Finschi, and constructs all elements of their classes. The output
// are all oriented matroids of the given rank and number of elements.

int R; // rank
int N; // the number of elements

int B;       // the number of bases
int nr_ints; // the number of integers needed to store the plus (resp. minus) of
             // a chirotope

char **bases; // the list of bases

int w;       // number of permutations
char **perm; // the list of permutations

void reorient(struct OM M, struct OM *chirotopes, int *count,
              int *c) // constructs all reorientations of the oriented matroid M
{
  struct OM X, Y;

  int i, j, k, l, s;
  long long int m;

  s = *count; // we compare only OMs that are in the reorientation class of M,
              // since we know that M is not in the permutation-reorientation
              // class of previous OMs

  makeOM(&X);

  char elements[N - 1]; // elements to be reoriented - all nonloops, but the
                        // smallest one
  int nr_elements = -1; // the number of elements to be reoriented

  for (i = 0; i < N; i++) // makes the list of elements to be reoriented
  {
    for (j = 0; j < B; j++) {
      for (k = 0; k < R; k++) {
        if (bases[j][k] == i) {
          l = j >> 5;
          m = 1l << (j & 31);
          if ((M.plus[l] & m) || (M.minus[l] & m)) {
            if (nr_elements >= 0)
              elements[nr_elements] = i;
            nr_elements++;
            k = R;
            j = B;
          }
        }
      }
    }
  }

  int x = (1 << nr_elements) - 1; // x encodes which elements change sign

  while (x >= 0) {
    copyOM(M, &X);

    for (i = 0; i < nr_elements; i++)
      if (x & (1 << i)) // elements[i] has to be reoriented
        for (j = 0; j < B; j++)
          for (k = 0; k < R; k++) {

            if (bases[j][k] == elements[i]) {
              l = j >> 5;
              m = 1l << (j & 31);
              if (X.plus[l] & m) {
                X.plus[l] -= m;
                X.minus[l] += m;
              } else if (X.minus[l] & m) {
                X.minus[l] -= m;
                X.plus[l] += m;
              }
            } else if (bases[j][k] > elements[i])
              k = R;
          }

    for (i = s; i < (*count);
         i++) // checks whether we have already constructed this OM
      if (isequal(X, chirotopes[i]))
        i = *count + 1;

    if (i == *count) // this is the first time to construct X, so store it
    {
      makeOM(&chirotopes[*count]);
      copyOM(X, &(chirotopes[*count]));
      *count = *count + 1;
      if (*count == *c) // overflow - change c
      {
        printf("Overflow: count=%d, c=%d\n", *count, *c);
        removeOM(&X);
        return;
      }
    }

    x--;
  }

  removeOM(&X);
}

int permutereorient(struct OM M, FILE *out) // permutes elements of an OM
{
  int i, j, k, count;
  struct OM X;

  int c = (1 << (N - 1)) * factorial(N) +
          1; // how many elements are there in a class?
  struct OM chirotopes[c];

  count = 0; // counts elements in the class of M

  makeOM(&X);

  for (i = 0; i < w; i++) // first permute the elements
  {
    X = permute(M, perm[i]);
    for (j = 0; j < count; j++)
      if (isequal(X, chirotopes[j])) // this chirotope has already been found
        j = count + 1;

    if (j == count) // the chirotope has not been constructed yet, thus find all
                    // its reorientations
      reorient(X, chirotopes, &count, &c);
  }

  for (i = 0; i < count; i++) {
    standardizeOM(&chirotopes[i]);
    writeOM(chirotopes[i], out);
  }

  for (i = 0; i < count; i++)
    removeOM(&chirotopes[i]);

  removeOM(&X);

  return count;
}

int makeallchirotopes(
    struct OM M) // makes all OMs in the permutation/reorientation class of M
{
  FILE *tfile;
  char text[300];
  struct OM X;

  int c, k;

  makeOM(&X);

  k = countbases(M);

  sprintf(text, "/storage/mi/palic/all_OMs_rank%d_%delements_%dbases.txt", R, N,
          k); // the OMs are stored separately, sorted by the number of bases

  tfile = fopen(
      text,
      "r"); // reads all uniform OMs that are representatives of their classes
  if (tfile != NULL) // if the file already exists, check whether we have
                     // already constructed M
  {
    while (readOM(&X, tfile) != 0)
      if (isequal(
              X,
              M)) // this OM has already been constructed, no need to continue
      {
        fclose(tfile);
        removeOM(&X);
        return 0;
      }
    fclose(tfile);
  }

  sprintf(text, "/storage/mi/palic/all_OMs_rank%d_%delements_%dbases.txt", R, N,
          k);

  tfile = fopen(text, "a"); // this file contains all OMs with k bases
  if (tfile == NULL) {
    fprintf(stderr, "error fopen():  Could not open text file %s.txt.\n", text);
    exit(EXIT_FAILURE);
  }

  c = permutereorient(M, tfile);

  fclose(tfile);
  removeOM(&X);

  return c;
}

int main(int argc, char *argv[]) {
  int step, nr_cones;

  R = 3; // The code works only for B<=64!
  N = 7;
  makebases();
  makepermutations();

  if (R == 2) // nr_cones encodes the number of lower cones
    nr_cones = 1;
  else if (R == 3)
    if (N <= 5)
      nr_cones = 1;
    else if (N == 6)
      nr_cones = 4;
    else if (N == 7)
      nr_cones = 11;
    else
      nr_cones = 0;

  FILE *in;
  char text[300];
  struct OM M;
  makeOM(&M);
  int i = 0;

  for (step = 0; step < nr_cones; step++) // make sure to read every lower cone
  {
    printf("step=%d\n", step);
    sprintf(text, "lower_cones_rank%d_%delements_%d.txt", R, N, step);

    in = fopen(
        text,
        "r"); // reads all uniform OMs that are representatives of their classes
    if (in == NULL) {
      fprintf(stderr, "error fopen():  Could not open text file %s.\n", text);
      exit(EXIT_FAILURE);
    }
    fgets(text, 300, in);

    while (readOM(&M, in) != 0) // reads elements of lower cones and makes their
                                // permutation/reorientaion classes
      i += makeallchirotopes(M);

    fclose(in);
  }

  removeOM(&M);

  printf("%d chirotopes\n", i);

  return 0;
}
