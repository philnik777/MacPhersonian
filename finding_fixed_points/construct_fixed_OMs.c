#include </home/mi/palic/OMs.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// This program makes a list of tests for checking whether an oriented matroid
// is fixed under the action of a group Z_3 for N=9. It makes only first size
// many signs, the rest is made by the program whole_l.c. Finally, that list is
// used by construct_fixed_OMs.c in order to construct all fixed oriented
// matroids under the action of Z_3

int B; // the number of bases
int N; // the number of elements
int R; // dimension
int nr_ints;

char **bases; // the list of bases

int w;       // number of permutations
char **perm; // the list of permutations

int sizeofgroup;
char **action;

int **orbits;      // the list of orbits we make here
char **orbitsigns; // the signs of permutations used for making orbits
int orbitnr;       // the total number of orbits
int *orind;

// struct OM chirotopes[3000];
int count;

int size;

/*
void makegroupaction()						//transitive
action of Z_2+Z_2+Z_2
{
        sizeofgroup=8;
        action=(char**) malloc(sizeofgroup*sizeof(char*));
        int s;
        for (s=0;s<sizeofgroup;s++)
                action[s]=(char*) malloc(N*sizeof(char));

        char i,j,k,l;

        for (i=0;i<2;i++)
        for (j=0;j<2;j++)
        for (k=0;k<2;k++)
                for (l=0;l<8;l++)
                        action[i+(j<<1)+(k<<2)][l]=l+(1-((l&1)<<1))*i+(2-((l&2)<<1))*j+(4-((l&4)<<1))*k;

        for (i=0;i<sizeofgroup;i++)				//the action for
N>8 for (j=sizeofgroup;j<N;j++) action[i][j]=j;
}

/*
void makegroupaction()			//action of a cyclic group
{
        sizeofgroup=7;
        action=(char**) malloc(sizeofgroup*sizeof(char*));
        int s;
        for (s=0;s<sizeofgroup;s++)
                action[s]=(char*) malloc(N*sizeof(char));

        char i,j,k,l;

        for (i=0;i<sizeofgroup;i++)
        {
                for (j=0;j<sizeofgroup;j++)
                        action[i][j]=(i+j)%sizeofgroup;
                for (j=sizeofgroup;j<N;j++)
                        action[i][j]=j;
        }

}

/*
void makegroupaction()						//transitive
action of Z_2+Z_2+Z_2
{
        sizeofgroup=9;
        action=(char**) malloc(sizeofgroup*sizeof(char*));
        int s;
        for (s=0;s<sizeofgroup;s++)
                action[s]=(char*) malloc(N*sizeof(char));

        char i,j,k,l;

        for (i=0;i<3;i++)
        for (j=0;j<3;j++)
                for (l=0;l<sizeofgroup;l++)
                        action[3*i+j][l]=((l+j)%3+(l/3)*3+3*i)%9;

        for (i=0;i<sizeofgroup;i++)				//the action for
N>9 for (j=sizeofgroup;j<N;j++) action[i][j]=j;


}

*/
void makegroupaction() // the group that fixes RS(8)
{
  action = (char **)malloc(4 * sizeof(char *));
  int s;
  for (s = 0; s < 4; s++)
    action[s] = (char *)malloc(N * sizeof(char));

  int i, j;
  struct OM X, M;
  makeOM(&X);
  makeOM(&M);

  // char text[300]="-----+-++---+++++++-------+++++----+++++++--+---++-+++----------++++-+";				//The OM EFM
  // char
  // text[300]="----+++---+++++---+----++------+++++--++--+-++---++--++++++++----------+++---++++++-";
  // //The OM for r=3,N=9 that Richter-Gebert gave me
  char text[300] = "+++++++-++----+-+--+++--++--+-+++-+----+++-+++-+-++++--+---"
                   "-+++---+--+"; // RS(8)

  for (i = 0; i < B; i++)
    if (text[i] == '+')
      X.plus[i >> 5] += 1 << (i % 32);
    else
      X.minus[i >> 5] += 1 << (i % 32);

  sizeofgroup = 0;
  for (i = 0; i < w; i++) {
    M = permute(X, perm[i]);
    if (isequal(M, X)) {
      for (j = 0; j < N; j++)
        action[sizeofgroup][j] = perm[i][j];
      sizeofgroup++;
    }
  }

  removeOM(&X);
  removeOM(&M);
}
/*
void makegroupaction()				//the group Z_3 that fixes EFM
(subgroup of Z_2+Z_3)
{
        action=(char**) malloc(3*sizeof(char*));
        int s;
        for (s=0;s<3;s++)
                action[s]=(char*) malloc(N*sizeof(char));

        int i,j;
        struct OM X,M;
        makeOM(&X);
        makeOM(&M);

        char
text[300]="-----+-++---+++++++-------+++++----+++++++--+---++-+++----------++++-+";
//The OM EFM
        //char
text[300]="----+++---+++++---+----++------+++++--++--+-++---++--++++++++----------+++---++++++-";
//The OM for r=3,N=9 that Richter-Gebert gave me


        for (i=0;i<B;i++)
                if (text[i]=='+')
                        X.plus[i>>5]+=1<<(i & 31);
                else X.minus[i>>5]+=1<<(i & 31);


        sizeofgroup=0;
        for (i=0;i<w;i++)
        {
                M=permute(X,perm[i]);
                if (isequal(M,X))
                {
                        for (j=0;j<N;j++)
                                action[sizeofgroup][j]=perm[i][j];
                        sizeofgroup++;

                }
                if (sizeofgroup==3)
                        break;
        }

        removeOM(&X);
        removeOM(&M);

}
*/
/*
void makegroupaction()						//transitive
action of Z_3+Z_3
{
        sizeofgroup=9;
        action=(char**) malloc(sizeofgroup*sizeof(char*));
        int s;
        for (s=0;s<sizeofgroup;s++)
                action[s]=(char*) malloc(N*sizeof(char));

        char i,j,k,l;

        for (i=0;i<3;i++)
        for (j=0;j<3;j++)
                for (l=0;l<sizeofgroup;l++)
                        action[3*i+j][l]=((l+j)%3+(l/3)*3+3*i)%9;

        for (i=0;i<sizeofgroup;i++)				//the action for
N>9 for (j=sizeofgroup;j<N;j++) action[i][j]=j;


}
*/
/*
void makegroupaction()			//action of Z_N
{
        sizeofgroup=N;
        action=(char**) malloc(sizeofgroup*sizeof(char*));
        int s;
        for (s=0;s<sizeofgroup;s++)
                action[s]=(char*) malloc(N*sizeof(char));

        char i,j,k,l;

        for (i=0;i<sizeofgroup;i++)
                for (j=0;j<sizeofgroup;j++)
                        action[i][j]=(i+j)%sizeofgroup;


}
*/
void initializeprogram() {
  makebases();
  makepermutations();
  makegroupaction();
}

void closeprogram() {
  removebases();
  removepermutations();
  removegroupaction();
}

struct OM permutewithsign(
    struct OM M,
    char s[]) // given an OM, it transforms it into a new one - permutes the
              // labels of the elements, but does not check whether p<m s[] is
              // an array of length N that stores the permutation
{
  int i, j;
  int b[B];     // b[i] is the index of the i-th basis after permutation
  char sign[B]; // sign[i] stores the sign of the permutation on the basis i
  char x[R];
  long long int t;
  char limit;
  char p, q;
  for (i = 0; i < B; i++) {
    for (j = 0; j < R; j++)
      x[j] = s[bases[i][j]];
    sign[i] = sort(x);
    b[i] = ind(x);
  }

  struct OM X;
  makeOM(&X);

  for (i = 0; i < nr_ints; i++) {
    if (i == nr_ints - 1)
      limit = B & 31;
    else
      limit = 32;

    for (j = 0; j < limit; j++) {
      t = 1 << j;

      p = b[(i << 5) + j] >> 5;
      q = b[(i << 5) + j] & 31;

      if (((M.plus[i] & t) && sign[(i << 5) + j] == 1) ||
          ((M.minus[i] & t) &&
           sign[(i << 5) + j] ==
               -1)) // if the b[i]-th basis was positive in the old chirotope,
                    //  the i-th basis is positive in the new chirotope
        X.plus[p] += 1 << q;
      if (((M.plus[i] & t) && sign[(i << 5) + j] == -1) ||
          ((M.minus[i] & t) && sign[(i << 5) + j] == 1))
        X.minus[p] += 1 << q;
    }
  }

  return X;
}

void makeorbits() // this function makes for every basis all its permutations
                  // under the group action, but stores each basis only once
{
  struct OM M, X;
  int i, j, s, p, q, k, limit_j, limit_q, found;

  int t[sizeofgroup];
  int sign[sizeofgroup];

  orbits = (int **)malloc(B * sizeof(int *));
  for (i = 0; i < B; i++)
    orbits[i] = (int *)malloc(sizeofgroup * sizeof(int));

  orbitsigns = (char **)malloc(B * sizeof(char *));
  for (i = 0; i < B; i++)
    orbitsigns[i] = (char *)malloc(sizeofgroup * sizeof(char));

  orind = (int *)malloc(B * sizeof(int));

  orbitnr = 0;
  makeOM(&M);
  makeOM(&X);

  limit_j = 32;
  for (i = 0; i < nr_ints; i++) {
    if (i == nr_ints - 1)
      limit_j = B & 31;
    for (j = 0; j < limit_j; j++) //(i<<5)+j is the index of a basis
    {
      found = 0;
      for (p = 0; p < orbitnr; p++)
        for (q = 0; q < sizeofgroup; q++)
          if (orbits[p][q] == (i << 5) + j) // we have already got this basis as
                                            // a permutation of some other basis
          {
            p = orbitnr;
            q = sizeofgroup;
            found = 1;
          }

      if (found == 0) // this is the first basis in its orbit that we found
      {
        M.plus[i] = 1u << j;

        orbits[orbitnr][0] = (i << 5) + j;
        orbitsigns[orbitnr][0] = 1;
        orind[(i << 5) + j] = orbitnr;
        for (k = 1; k < sizeofgroup; k++) // find the whole orbit
        {
          X = permutewithsign(M, action[k]);

          limit_q = 32;
          for (p = 0; p < nr_ints; p++) {
            if (p == nr_ints - 1)
              limit_q = B & 31;
            for (q = 0; q < limit_q; q++) {

              if (X.plus[p] & 1u << q) {
                orbits[orbitnr][k] = (p << 5) + q;
                orbitsigns[orbitnr][k] = 1;
                orind[(p << 5) + q] = orbitnr;
                p = nr_ints;
                q = 32;
              } else if (X.minus[p] & 1u << q) {
                orbits[orbitnr][k] = (p << 5) + q;
                orbitsigns[orbitnr][k] = -1;
                orind[(p << 5) + q] = orbitnr;
                p = nr_ints;
                q = 32;
              }
            }
          }
        }

        orbitnr++;
      }
    }
    M.plus[i] = 0;
  }

  orbits = (int **)realloc(orbits, orbitnr * sizeof(int *));
  orbitsigns = (char **)realloc(orbitsigns, orbitnr * sizeof(char *));

  removeOM(&M);
  removeOM(&X);
}

void construct(char *l) // given l, it tries all signs for the bases (given by
                        // l) and checks which of them yield oriented matroids
{
  struct OM M;
  int i, k, p, q;
  long long int t;

  makeOM(&M);

  int x;

  x = (1 << sizeofgroup - 1) -
      1; // used to say which bases will change signs. Note that if the j-th
         // bases of one orbit changes the signs, the j-th bases of every orbit
         // has to change the sign.

  while (x >= 0) {
    for (p = 0; p < nr_ints; p++) {
      M.plus[p] = 0;
      M.minus[p] = 0;
    }

    for (i = 0; i < orbitnr;
         i++) // makes an OM for this choice of basis orientations
      if (l[i])
        for (k = 0; k < sizeofgroup; k++) {
          t = 1 << (orbits[i][k] & 31);
          q = orbits[i][k] >> 5;
          if (((M.plus[q] | M.minus[q]) & t) == 0)
            if ((l[i] * orbitsigns[i][k] == 1 && (x & 1 << k) == 0) ||
                (l[i] * orbitsigns[i][k] == -1 && (x & 1 << k)))
              M.plus[q] += t;
            else
              M.minus[q] += t;
          else if ((M.plus[q] & t &&
                    ((l[i] == orbitsigns[i][k] && x & 1 << k) ||
                     (l[i] != orbitsigns[i][k] && (x & 1 << k) == 0))) ||
                   (M.minus[q] & t &&
                    ((l[i] != orbitsigns[i][k] && x & 1 << k) ||
                     (l[i] == orbitsigns[i][k] &&
                      (x & 1 << k) ==
                          0)))) // contaridiction within the same orbit
          {
            k = sizeofgroup;
            i = orbitnr + 1;
          }
        }

    if (i == orbitnr && ischirotope(M)) // no contradiction and M is a chirotope
    {
      if (isfixed(M)) // actually true by construction, gets called rarely, so
                      // it can remain here
      {
        standardizeOM(&M);
        FILE *tfile;
        char text[100];
        sprintf(text, "q-fixed_points%d_%d_Z%d.txt", R, N,
                sizeofgroup); // output that stores all fixed points
        tfile = fopen(text, "a");
        if (tfile == NULL) {
          fprintf(stderr, "error fopen():  Could not open the file %s.\n",
                  text);
          exit(EXIT_FAILURE);
        }

        for (p = nr_ints - 1; p >= 0; p--)
          fprintf(tfile, "%15u", M.plus[p]);
        fprintf(tfile, "\t\t");
        for (p = nr_ints - 1; p >= 0; p--)
          fprintf(tfile, "%15u", M.minus[p]);
        fputc('\n', tfile);

        fclose(tfile);

        count++;
      }
    }

    x--;
  }

  removeOM(&M);
}

void constructall(char *l) // decides on signs of rows of triples (orbits)
{
  int i, j, k;
  char L[orbitnr];
  int used[orbitnr];
  int nrls = 0;

  for (i = 0; i < orbitnr; i++)
    if (l[i]) {
      used[nrls] = i;
      nrls++;
    }

  if (nrls == 0)
    return;

  int y = (1u << nrls - 1) -
          1; // encodes the signs of orbits, we can take all bases of an orbit
             // to be positive, or to be negative

  for (i = 0; i < orbitnr; i++)
    printf("%3d", l[i]);
  putchar('\n');

  while (y >= 0) {
    for (i = 0; i < orbitnr; i++)
      L[i] = l[i];
    for (i = 0; i < nrls - 1; i++)
      if (y & (1 << i))
        L[used[i]] = -l[used[i]];

    construct(L);
    y--;
  }
}

main() {
  R = 4;
  N = 8;
  initializeprogram();

  int i, j, k, p, q;

  for (i = 0; i < sizeofgroup; i++) {
    for (j = 0; j < N; j++)
      printf("%d ", action[i][j]);
    putchar('\n');
  }

  makeorbits();

  printf("%d orbits\n", orbitnr);
  for (i = 0; i < orbitnr; i++) {
    for (p = 0; p < sizeofgroup; p++) {
      if (orbitsigns[i][p] == 1)
        printf("+(");
      else
        printf("-(");
      for (q = 0; q < R - 1; q++)
        printf("%d,", bases[orbits[i][p]][q] + 1);
      printf("%d) ", bases[orbits[i][p]][R - 1] + 1);
    }
    putchar('\n');
  }

  FILE *tfile;
  char text[300];

  sprintf(text, "q-fixed_points%d_%d_Z%d.txt", R, N,
          sizeofgroup); // output that stores all fixed points
  tfile = fopen(text, "w");
  if (tfile == NULL) {
    fprintf(stderr, "error fopen():  Could not open the file %s.\n", text);
    exit(EXIT_FAILURE);
  }

  fprintf(tfile,
          "Fixed points in MacP(%d,%d) under the action of Z_%d (writen for "
          "the code for rank 4):\n",
          R, N, sizeofgroup);

  fclose(tfile);

  sprintf(text, "q-possible_chirotopes_tests%d_%d_Z%d_length%d.txt", R, N,
          sizeofgroup, orbitnr); // output that stores all fixed points
  tfile = fopen(text, "r");
  if (tfile == NULL) {
    fprintf(stderr, "error fopen():  Could not open the file %s.\n", text);
    exit(EXIT_FAILURE);
  }

  char l[orbitnr];
  count = 0;

  for (i = 0; i < orbitnr; i++)
    k = fscanf(tfile, "%d", &l[i]);

  while (k == 1) {
    constructall(l);
    for (i = 0; i < orbitnr; i++)
      k = fscanf(tfile, "%d", &l[i]);
  }

  fclose(tfile);

  printf("%d fixed points for MacP(%d,%d) under the action of Z%d\n", count, R,
         N, sizeofgroup);

  closeprogram();
}
