#include </home/mi/palic/OMs.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// This program makes a list of orbits of the action of the group on bases. It
// makes only the first part, the rest is done by the program q-longer_l.c

int R; // rank
int N; // the number of elements

int B;       // the number of bases
int nr_ints; // the number of integers needed to store the plus (resp. minus) of
             // a chirotope

char **bases; // the list of bases

int w;       // number of permutations
char **perm; // the list of permutations

int sizeofgroup; // the number of elements of the group that acts on [N] (and
                 // thus on MacP)
char **action; // the list of elements of the group, seen as a subgroup of the
               // permutation group

int **orbits;      // the list of orbits we make here
char **orbitsigns; // the signs of permutations used for making orbits
int orbitnr;       // the total number of orbits
int *orind;        // for every basis stores the orbit it is in

int size; // the length of test combinations we make here

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

*//*
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
(subgroup of D_3)
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

int isgood(char *l, int z, int in1,
           int in2) // used in ispossible(), checks whether two bases in1 and
                    // in2 can be as they are given by l
{
  int i, j, k, found, s1, s2, t1, t2, p;
  found = 0;

  char x[R];
  char y[R];

  for (j = 0; j < R;
       j++) // different permutations of the first basis (chirotopes are
            // alternating and we do not store all the bases)
  {
    for (i = 0; i < R; i++) // axiom B2'
    {

      for (k = 0; k < R; k++) // make x_i's and y_i's as in the axiom B2'
      {
        x[k] = bases[in1][k];
        y[k] = bases[in2][k];
      }
      p = x[0];
      x[0] = x[j];
      x[j] = p;

      p = x[0];
      x[0] = y[i];
      y[i] = p;

      s1 = sort(x);
      s2 = sort(y);

      if (s1 && s2) // the x's and the y's are non-zero bases, so we want to
                    // find an i from Axiom B2'
      {

        t1 = orind[ind(x)];
        t2 = orind[ind(y)];

        if (t1 >= z && t2 >= z) // we don't know the signs of these bases - i
                                // could be among them
        {
          found++;
          break;
        }
        if ((t1 >= z && s2) ||
            (t2 >= z &&
             s1)) // one basis is not zero and we do not know anything about the
                  // other one, that one could be i
        {
          found++;
          break;
        }
        if (t1 < z && t2 < z && l[t1] &&
            l[t2]) // both bases are not zero, we found i (possibly, we still
                   // don't know the sign)
        {
          found++;
          break;
        }
      }
    }

    if (found <= j) // for one pair of bases Axiom B2' is not fullfilled
      return 0;
  }

  return 1;
}

int ispossible(
    char *l,
    int z) // checks for the list of signs of orbits of length z whether it can
           // be extended to a full list of signs. Most sign combination do not
           // lead to a chirotope and that is what we want to exclude here.
{
  int i, j, p, q, limit;

  for (i = 0; i < z; i++) // orbits
    if (l[i])
      for (j = 0; j < sizeofgroup; j++) // bases in the orbit
      {
        limit = sizeofgroup;
        for (p = 0; p <= i; p++) // previous orbits
        {
          if (p == i)
            limit = j;
          for (q = 0; q < limit; q++)
            if (l[orind[orbits[p][q]]] != 0)
              if (isgood(l, z, orbits[i][j], orbits[p][q]) ==
                  0) // checking partially chirotope axioms
                return 0;
        }
      }

  return 1;
}

long long int made;

void makepossible(char *l, int i,
                  int z) // z is the length of l that we are making, in final
                         // version z=orbitnr make first z signs of tests that
                         // do not lead to a contradiction
{
  if (i == z) {
    if (ispossible(l, z)) {
      FILE *tfile;
      char text[300];
      int j;
      sprintf(text, "q-possible_chirotopes_tests%d_%d_Z%d_length%d.txt", R, N,
              sizeofgroup, size);

      tfile = fopen(text, "a");
      if (tfile == NULL) {
        fprintf(stderr, "error fopen():  Could not open the file %s.\n", text);
        exit(EXIT_FAILURE);
      }

      for (j = 0; j < z; j++)
        fprintf(tfile, "%2d", l[j]);
      fputc('\n', tfile);
      fclose(tfile);

      made++;
    }

    return;
  }

  l[i] = 0;
  makepossible(l, i + 1, z);
  l[i] = 1;
  makepossible(l, i + 1, z);
}

main() {
  R = 4;
  N = 8;
  initializeprogram();

  printf("B=%d, nr_ints=%d\n", B, nr_ints);

  int i, j, k, s, limit;

  makeorbits();

  printf("%d orbits\n", orbitnr);

  for (i = 0; i < sizeofgroup; i++) {
    for (j = 0; j < N; j++)
      printf("%3d", action[i][j] + 1);
    putchar('\n');
  }

  size = 3;

  made = 0;

  char l[orbitnr];

  FILE *tfile;
  char text[300];
  sprintf(text, "q-possible_chirotopes_tests%d_%d_Z%d_length%d.txt", R, N,
          sizeofgroup, size);
  tfile = fopen(text, "w");
  if (tfile == NULL) {
    fprintf(stderr, "error fopen():  Could not open the file %s.\n", text);
    exit(EXIT_FAILURE);
  }
  fclose(tfile);

  clock_t begin, end;
  begin = clock();

  makepossible(l, 0, size);
  end = clock();
  printf("%.2f min, %llu l's, R=%d, N=%d, sizeofgroup=%d, size=%d\n",
         ((double)(end - begin)) / CLOCKS_PER_SEC / 60, made, R, N, sizeofgroup,
         size);

  closeprogram();
}
