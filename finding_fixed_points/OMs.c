#include </home/mi/palic/OMs.h>
#include <stdio.h>
#include <stdlib.h>

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

void makebases() // makes the list of all posible bases a chirotope could have,
                 // bases are 012, 013,...
{
  int i, k, s;

  B = 1;
  for (i = N; i > N - R; i--)
    B *= i;
  for (i = 2; i <= R; i++)
    B /= i;

  nr_ints = B >> 5;

  if (B & 31)
    nr_ints++;

  bases = (char **)malloc(B * sizeof(char *));
  for (s = 0; s < B; s++)
    bases[s] = (char *)malloc(R * sizeof(char));

  for (i = 0; i < R; i++)
    bases[0][i] = i;
  s = 1;

  while (1) {
    k = R - 1;
    while (k >= 0 && bases[s - 1][k] >= N - R + k)
      k--;
    if (k == -1)
      break;
    for (i = 0; i < k; i++)
      bases[s][i] = bases[s - 1][i];
    for (i = k; i < R; i++)
      bases[s][i] = bases[s - 1][k] + i - k + 1;
    s++;
  }
}

void removebases() // frees the memory used by bases
{
  int s;
  for (s = 0; s < B; s++)
    free(bases[s]);
  free(bases);
}

void makeOM(struct OM *M) // allocates memory for an oriented matroid of rank R
                          // on N elements
{
  M->plus = malloc(nr_ints * sizeof(int));
  M->minus = malloc(nr_ints * sizeof(int));

  int i;
  for (i = 0; i < nr_ints; i++) {
    M->plus[i] = 0;
    M->minus[i] = 0;
  }
}

void removeOM(struct OM *M) // frees the memory used by an oriented matroid
{
  free(M->plus);
  free(M->minus);
}

void showbits(unsigned int *plus) // prints a list if integers in the binary
                                  // representation (smallest bit on the right)
{

  int i, j;

  for (i = (B & 31) - 1; i >= 0; i--)
    (plus[nr_ints - 1] & (1u << i)) ? putchar('1') : putchar('0');

  for (j = nr_ints - 2; j >= 0; j--)
    for (i = 31; i >= 0; i--)
      (plus[j] & (1u << i)) ? putchar('1') : putchar('0');
}

void showchirotope(struct OM M) // prints a chirotope
{
  int i, j;

  for (j = 0; j < nr_ints - 1; j++)
    for (i = 0; i < 32; i++)
      if (M.plus[j] & (1u << i))
        putchar('+');
      else if (M.minus[j] & (1u << i))
        putchar('-');
      else
        putchar('0');

  for (i = 0; i < (B & 31); i++)
    if (M.plus[nr_ints - 1] & (1u << i))
      putchar('+');
    else if (M.minus[nr_ints - 1] & (1u << i))
      putchar('-');
    else
      putchar('0');

  putchar('\n');
}

int countbases(struct OM M) // counts the number of bases of a chirotope
{
  int c = 0;
  int i, j;
  for (i = (B & 31) - 1; i >= 0; i--)
    if ((M.plus[nr_ints - 1] & 1 << i) || (M.minus[nr_ints - 1] & 1 << i))
      c++;
  for (j = nr_ints - 2; j >= 0; j--)
    for (i = 31; i >= 0; i--)
      if ((M.plus[j] & 1l << i) || (M.minus[j] & 1l << i))
        c++;
  return c;
}

int weakmap(
    struct OM M1,
    struct OM M2) // checks whether there is a weak map M_1 \wm M_2
                  // this actually checks weak maps for oriented matroids, since
                  // we make only one chirotope from each pair chi, -chi
{
  int i;
  int x, y;

  int good = 0;

  for (i = 0; i < nr_ints; i++) {
    x = M1.plus[i] ^ M2.plus[i];
    y = M1.minus[i] ^ M2.minus[i];
    if (((x & M1.plus[i]) == x) && ((y & M1.minus[i]) == y)) // chi_1 >= chi_2
      good++;
    else
      break;
  }

  if (good == nr_ints)
    return 1;

  good = 0;
  for (i = 0; i < nr_ints; i++) {
    x = M1.plus[i] ^ M2.minus[i];
    y = M1.minus[i] ^ M2.plus[i];
    if (((x & M1.plus[i]) == x) && ((y & M1.minus[i]) == y)) // chi_1 >= -chi_2
      good++;
    else
      break;
  }

  return good == nr_ints;
}

int isequal(
    struct OM M1,
    struct OM M2) // checks whether the oriented matroids M1 and M2 are the same
{
  int i, good;
  good = 0;
  for (i = 0; i < nr_ints; i++) // chi_1==chi_2
    if (M1.plus[i] == M2.plus[i] && M1.minus[i] == M2.minus[i])
      good++;
    else
      break;
  if (good == nr_ints)
    return 1;

  good = 0;
  for (i = 0; i < nr_ints; i++) // chi_1==-chi_2
    if (M1.plus[i] == M2.minus[i] && M1.minus[i] == M2.plus[i])
      good++;
    else
      break;
  if (good == nr_ints)
    return 1;
  else
    return 0;
}

int ind(char *a) // returns the index of the basis (a[0],a[1],...,a[R-1]) in the
                 // array bases[][], assumes that a[0]<a[1]<...<a[R]
{
  if (R == 3 && N <= 9) // quick code for common used parameters
  {
    if (N == 9) {
      if (a[0] == 0)
        if (a[1] == 1)
          return a[2] - 2;
        else if (a[1] == 2)
          return a[2] + 4;
        else if (a[1] == 3)
          return a[2] + 9;
        else if (a[1] == 4)
          return a[2] + 13;
        else if (a[1] == 5)
          return a[2] + 16;
        else
          return a[1] + a[2] + 12;
      if (a[0] == 1)
        if (a[1] == 2)
          return a[2] + 25;
        else if (a[1] == 3)
          return a[2] + 30;
        else if (a[1] == 4)
          return a[2] + 34;
        else if (a[1] == 5)
          return a[2] + 37;
        else
          return a[1] + a[2] + 33;
      if (a[0] == 2)
        if (a[1] == 3)
          return a[2] + 45;
        else if (a[1] == 4)
          return a[2] + 49;
        else if (a[1] == 5)
          return a[2] + 52;
        else
          return a[1] + a[2] + 48;
      if (a[0] == 3)
        if (a[1] == 4)
          return a[2] + 59;
        else if (a[1] == 5)
          return a[2] + 62;
        else
          return a[1] + a[2] + 58;
      if (a[0] == 4)
        if (a[1] == 5)
          return a[2] + 68;
        else
          return a[1] + a[2] + 64;
      return a[0] + a[1] + a[2] + 62;
    }
    if (N == 8) {
      if (a[0] == 0)
        if (a[1] == 1)
          return a[2] - 2;
        else if (a[1] == 2)
          return a[2] + 3;
        else if (a[1] == 3)
          return a[2] + 7;
        else if (a[1] == 4)
          return a[2] + 10;
        else
          return a[1] + a[2] + 7;
      if (a[0] == 1)
        if (a[1] == 2)
          return a[2] + 18;
        else if (a[1] == 3)
          return a[2] + 22;
        else if (a[1] == 4)
          return a[2] + 25;
        else
          return a[1] + a[2] + 22;
      if (a[0] == 2)
        if (a[1] == 3)
          return a[2] + 32;
        else if (a[1] == 4)
          return a[2] + 35;
        else
          return a[1] + a[2] + 32;
      if (a[0] == 3)
        if (a[1] == 4)
          return a[2] + 41;
        else
          return a[1] + a[2] + 38;
      return a[0] + a[1] + a[2] + 37;
    }
    if (N == 7) {
      if (a[0] == 0)
        if (a[1] == 1)
          return a[2] - 2;
        else if (a[1] == 2)
          return a[2] + 2;
        else if (a[1] == 3)
          return a[2] + 5;
        else
          return a[1] + a[2] + 3;
      if (a[0] == 1)
        if (a[1] == 2)
          return a[2] + 12;
        else if (a[1] == 3)
          return a[2] + 15;
        else
          return a[1] + a[2] + 13;
      if (a[0] == 2)
        if (a[1] == 3)
          return a[2] + 21;
        else
          return a[1] + a[2] + 19;
      return a[0] + a[1] + a[2] + 19;
    }

    if (N == 6) {
      if (a[0] == 0)
        if (a[1] == 1)
          return a[2] - 2;
        else if (a[1] == 2)
          return a[2] + 1;
        else
          return a[1] + a[2];
      if (a[0] == 1)
        if (a[1] == 2)
          return a[2] + 7;
        else
          return a[1] + a[2] + 6;
      return a[0] + a[1] + a[2] + 7;
    }

    if (N == 5) {
      if (a[0] == 0)
        if (a[1] == 1)
          return a[2] - 2;
        else
          return a[1] + a[2] - 2;
      return a[0] + a[1] + a[2];
    }
    if (N == 4) {
      return a[0] + a[1] + a[2] - 3;
    }
  }
  if (R == 4 && N <= 9) {
    if (N == 9) {
      if (a[0] == 0)
        if (a[1] == 1)
          if (a[2] == 2)
            return a[3] - 3;
          else if (a[2] == 3)
            return a[3] + 2;
          else if (a[2] == 4)
            return a[3] + 6;
          else if (a[2] == 5)
            return a[3] + 9;
          else
            return a[2] + a[3] + 5;
        else if (a[1] == 2)
          if (a[2] == 3)
            return a[3] + 17;
          else if (a[2] == 4)
            return a[3] + 21;
          else if (a[2] == 5)
            return a[3] + 24;
          else
            return a[2] + a[3] + 20;
        else if (a[1] == 3)
          if (a[2] == 4)
            return a[3] + 31;
          else if (a[2] == 5)
            return a[3] + 34;
          else
            return a[2] + a[3] + 30;
        else if (a[1] == 4)
          if (a[2] == 5)
            return a[3] + 40;
          else
            return a[2] + a[3] + 36;
        else
          return a[1] + a[2] + a[3] + 34;
      if (a[0] == 1)
        if (a[1] == 2)
          if (a[2] == 3)
            return a[3] + 52;
          else if (a[2] == 4)
            return a[3] + 56;
          else if (a[2] == 5)
            return a[3] + 59;
          else
            return a[2] + a[3] + 55;
        else if (a[1] == 3)
          if (a[2] == 4)
            return a[3] + 66;
          else if (a[2] == 5)
            return a[3] + 69;
          else
            return a[2] + a[3] + 65;
        else if (a[1] == 4)
          if (a[2] == 5)
            return a[3] + 75;
          else
            return a[2] + a[3] + 71;
        else
          return a[1] + a[2] + a[3] + 69;
      if (a[0] == 2)
        if (a[1] == 3)
          if (a[2] == 4)
            return a[3] + 86;
          else if (a[2] == 5)
            return a[3] + 89;
          else
            return a[2] + a[3] + 85;
        else if (a[1] == 4)
          if (a[2] == 5)
            return a[3] + 95;
          else
            return a[2] + a[3] + 91;
        else
          return a[1] + a[2] + a[3] + 89;
      if (a[0] == 3)
        if (a[1] == 4)
          if (a[2] == 5)
            return a[3] + 105;
          else
            return a[2] + a[3] + 101;
        else
          return a[1] + a[2] + a[3] + 99;
      return a[0] + a[1] + a[2] + a[3] + 99;
    }

    if (N == 8) {
      if (a[0] == 0)
        if (a[1] == 1)
          if (a[2] == 2)
            return a[3] - 3;
          else if (a[2] == 3)
            return a[3] + 1;
          else if (a[2] == 4)
            return a[3] + 4;
          else
            return a[2] + a[3] + 1;
        else if (a[1] == 2)
          if (a[2] == 3)
            return a[3] + 11;
          else if (a[2] == 4)
            return a[3] + 14;
          else
            return a[2] + a[3] + 11;
        else if (a[1] == 3)
          if (a[2] == 4)
            return a[3] + 20;
          else
            return a[2] + a[3] + 17;
        else
          return a[1] + a[2] + a[3] + 16;
      if (a[0] == 1)
        if (a[1] == 2)
          if (a[2] == 3)
            return a[3] + 31;
          else if (a[2] == 4)
            return a[3] + 34;
          else
            return a[2] + a[3] + 31;
        else if (a[1] == 3)
          if (a[2] == 4)
            return a[3] + 40;
          else
            return a[2] + a[3] + 37;
        else
          return a[1] + a[2] + a[3] + 36;
      if (a[0] == 2)
        if (a[1] == 3)
          if (a[2] == 4)
            return a[3] + 50;
          else
            return a[2] + a[3] + 47;
        else
          return a[1] + a[2] + a[3] + 46;
      return a[0] + a[1] + a[2] + a[3] + 47;
    }

    if (N == 7) {
      if (a[0] == 0)
        if (a[1] == 1)
          if (a[2] == 2)
            return a[3] - 3;
          else if (a[2] == 3)
            return a[3];
          else
            return a[2] + a[3] - 2;
        else if (a[1] == 2)
          if (a[2] == 3)
            return a[3] + 6;
          else
            return a[2] + a[3] + 4;
        else
          return a[1] + a[2] + a[3] + 4;
      if (a[0] == 1)
        if (a[1] == 2)
          if (a[2] == 3)
            return a[3] + 16;
          else
            return a[2] + a[3] + 14;
        else
          return a[1] + a[2] + a[3] + 14;
      return a[0] + a[1] + a[2] + a[3] + 16;
    }

    if (N == 6) {
      if (a[0] == 0)
        if (a[1] == 1)
          if (a[2] == 2)
            return a[3] - 3;
          else
            return a[2] + a[3] - 4;
        else
          return a[1] + a[2] + a[3] - 3;
      return a[0] + a[1] + a[2] + a[3];
    }
    if (N == 5)
      return a[0] + a[1] + a[2] + a[3] - 6;
  }

  int l, u, m, i; // code that works for every parameter, slow...
  l = 0;
  u = B - 1;
  m = (l + u) >> 1;

  for (i = 0; i < R; i++) {
    while (a[i] != bases[m][i])
      if (a[i] > bases[m][i]) {
        l = m + 1;
        if (l == u)
          return l;
        m = (u + l) >> 1;

      } else if (a[i] < bases[m][i]) {
        u = m - 1;
        if (l == u)
          return l;
        m = (u + l) >> 1;
      }
    l = m;
    u = m;
    while (l >= 0 && a[i] == bases[l][i])
      l--;
    l++;
    while (u < B && a[i] == bases[u][i])
      u++;
    u--;
    if (l == u)
      return l;
    m = (u + l) >> 1;
  }
}

int sort(char a[]) // sorts integers in the array and returns the sign of the
                   // permutation
{
  char p, q;
  if (R == 3) // quick code for three integers
    if (a[1] < a[0])
      if (a[2] < a[1]) {
        p = a[2];
        a[2] = a[0];
        a[0] = p;
        return -1;
      } else if (a[2] == a[1])
        return 0;
      else if (a[0] < a[2]) {
        p = a[1];
        a[1] = a[0];
        a[0] = p;
        return -1;
      } else if (a[0] == a[2])
        return 0;
      else {
        p = a[0];
        a[0] = a[1];
        a[1] = a[2];
        a[2] = p;
        return 1;
      }
    else if (a[1] == a[0])
      return 0;
    else if (a[2] < a[1])
      if (a[0] < a[2]) {
        p = a[1];
        a[1] = a[2];
        a[2] = p;
        return -1;
      } else if (a[0] == a[2])
        return 0;
      else {
        p = a[0];
        a[0] = a[2];
        a[2] = a[1];
        a[1] = p;
        return 1;
      }
    else if (a[2] == a[1])
      return 0;
    else
      return 1;

  if (R == 4) // quick code for four integers
    if (a[1] < a[0])
      if (a[2] < a[1])
        if (a[3] < a[2]) {
          p = a[0];
          a[0] = a[3];
          a[3] = p;
          p = a[1];
          a[1] = a[2];
          a[2] = p;
          return 1;
        } else if (a[3] == a[2])
          return 0;
        else if (a[3] < a[1]) {
          p = a[0];
          a[0] = a[2];
          a[2] = a[1];
          a[1] = a[3];
          a[3] = p;
          return -1;
        } else if (a[1] == a[3])
          return 0;
        else if (a[3] < a[0]) {
          p = a[0];
          a[0] = a[2];
          a[2] = a[3];
          a[3] = p;
          return 1;
        } else if (a[3] == a[0])
          return 0;
        else {
          p = a[0];
          a[0] = a[2];
          a[2] = p;
          return -1;
        }
      else if (a[2] == a[1])
        return 0;
      else if (a[3] < a[1])
        if (a[2] < a[0]) {
          p = a[0];
          a[0] = a[3];
          a[3] = p;
          return -1;
        } else if (a[2] == a[0])
          return 0;
        else {
          p = a[0];
          a[0] = a[3];
          a[3] = a[2];
          a[2] = p;
          return 1;
        }
      else if (a[3] == a[1])
        return 0;
      else if (a[3] < a[2])
        if (a[2] < a[0]) {
          p = a[0];
          a[0] = a[1];
          a[1] = a[3];
          a[3] = p;
          return 1;
        } else if (a[2] == a[0])
          return 0;
        else if (a[3] < a[0]) {
          p = a[0];
          a[0] = a[1];
          a[1] = a[3];
          a[3] = a[2];
          a[2] = p;
          return -1;
        } else if (a[3] == a[0])
          return 0;
        else {
          p = a[0];
          a[0] = a[1];
          a[1] = p;
          p = a[2];
          a[2] = a[3];
          a[3] = p;
          return 1;
        }

      else if (a[3] == a[2])
        return 0;
      else if (a[0] < a[2]) {
        p = a[0];
        a[0] = a[1];
        a[1] = p;
        return -1;
      } else if (a[2] == a[0])
        return 0;
      else if (a[0] < a[3]) {
        p = a[0];
        a[0] = a[1];
        a[1] = a[2];
        a[2] = p;
        return 1;
      } else if (a[0] == a[3])
        return 0;
      else {
        p = a[0];
        a[0] = a[1];
        a[1] = a[2];
        a[2] = a[3];
        a[3] = p;
        return -1;
      }
    else if (a[1] == a[0])
      return 0;
    else if (a[3] < a[1])
      if (a[0] < a[3])
        if (a[2] < a[0]) {
          p = a[0];
          a[0] = a[2];
          a[2] = a[3];
          a[3] = a[1];
          a[1] = p;
          return -1;
        } else if (a[2] == a[0])
          return 0;
        else if (a[2] < a[3]) {
          p = a[1];
          a[1] = a[2];
          a[2] = a[3];
          a[3] = p;
          return 1;
        } else if (a[2] == a[3])
          return 0;
        else if (a[2] < a[1]) {
          p = a[1];
          a[1] = a[3];
          a[3] = p;
          return -1;
        } else if (a[2] == a[1])
          return 0;
        else {
          p = a[1];
          a[1] = a[3];
          a[3] = a[2];
          a[2] = p;
          return 1;
        }
      else if (a[0] == a[3])
        return 0;
      else if (a[2] < a[3]) {
        p = a[0];
        a[0] = a[2];
        a[2] = p;
        p = a[1];
        a[1] = a[3];
        a[3] = p;
        return 1;
      } else if (a[2] == a[3])
        return 0;
      else if (a[2] < a[0]) {
        p = a[0];
        a[0] = a[3];
        a[3] = a[1];
        a[1] = a[2];
        a[2] = p;
        return -1;
      } else if (a[2] == a[0])
        return 0;
      else if (a[2] < a[1]) {
        p = a[0];
        a[0] = a[3];
        a[3] = a[1];
        a[1] = p;
        return 1;
      } else if (a[2] == a[1])
        return 0;
      else {
        p = a[0];
        a[0] = a[3];
        a[3] = a[2];
        a[2] = a[1];
        a[1] = p;
        return -1;
      }
    else if (a[3] == a[1])
      return 0;
    else if (a[2] < a[0]) {
      p = a[0];
      a[0] = a[2];
      a[2] = a[1];
      a[1] = p;
      return 1;
    } else if (a[2] == a[0])
      return 0;
    else if (a[2] < a[1]) {
      p = a[1];
      a[1] = a[2];
      a[2] = p;
      return -1;
    } else if (a[2] == a[1])
      return 0;
    else if (a[2] < a[3])
      return 1;
    else if (a[2] == a[3])
      return 0;
    else {
      p = a[2];
      a[2] = a[3];
      a[3] = p;
      return -1;
    }

  int sorted, min, ind_min, i, j, sign; // slow code, but works for all R

  sorted = 0;
  sign = 1;
  while (sorted < R) {
    min = N;
    for (i = sorted; i < R; i++)
      if (a[i] < min) {
        min = a[i];
        ind_min = i;
      }
    if (ind_min != sorted) {
      j = a[sorted];
      a[sorted] = min;
      a[ind_min] = j;
      sign = -sign;
    }
    sorted++;
  }

  for (i = 1; i < R; i++)
    if (a[i - 1] == a[i])
      return 0;
  return sign;
}

char axB2(struct OM M, char sign, char s1, char s2, int in1, int in2)
// used in b2prime to check Axiom B2' in ischirotope, returns 1 if
// \chi(y1,x2,x3)*\chi(x1,y2,y3) and \chi(x1,x2,x3)*\chi(y1,y2,y3) have the same
// sign
{
  struct OM M1, M2;
  makeOM(&M1);
  makeOM(&M2);
  char i, j, p1, p2;
  long long int i1, i2;
  char res;
  if (s1 == 1) // if sign(y1,x2,x3)=1, then chi(y1,x2,x3) is exactly what we
               // give by plus,minus;
    for (i = 0; i < nr_ints; i++) {
      M1.plus[i] = M.plus[i];
      M1.minus[i] = M.minus[i];
    }
  else // if sign(y1,x2,x3)=-1, then we have to switch the signs of the bases
       // (chirotopes are alternating)
    for (i = 0; i < nr_ints; i++) {
      M1.plus[i] = M.minus[i];
      M1.minus[i] = M.plus[i];
    }

  if (s2 == 1) // the same as for s1
    for (i = 0; i < nr_ints; i++) {
      M2.plus[i] = M.plus[i];
      M2.minus[i] = M.minus[i];
    }
  else // if sign(y1,x2,x3)=-1, then we have to switch the signs of the bases
       // (chirotopes are alternating)
    for (i = 0; i < nr_ints; i++) {
      M2.plus[i] = M.minus[i];
      M2.minus[i] = M.plus[i];
    }

  p1 = in1 >> 5;
  p2 = in2 >> 5;
  in1 = in1 & 31;
  in2 = in2 & 31;

  i1 = 1l << in1;
  i2 = 1l << in2;

  res = -1;
  if (sign == -1)
    res = (((M1.plus[p1] & i1) && (M2.minus[p2] & i2)) ||
           ((M1.minus[p1] & i1) &&
            (M2.plus[p2] & i2))); //\chi(y1,x2,x3)*\chi(x1,y2,y3)
  else if (sign == 1)
    res = (((M1.plus[p1] & i1) && (M2.plus[p2] & i2)) ||
           ((M1.minus[p1] & i1) &&
            (M2.minus[p2] & i2))); //\chi(y1,x2,x3)*\chi(x1,y2,y3)

  removeOM(&M1);
  removeOM(&M2);

  return res;
}

char b2prime(struct OM M, char sign, char *X, char *Y) // checks Axiom B2'
{
  int s1, s2, in1, in2, i, j, q, sx, sy;
  char *x, *y;

  x = (char *)malloc(R * sizeof(char));
  y = (char *)malloc(R * sizeof(char));

  for (i = 0; i < R; i++) {
    x[i] = X[i];
    y[i] = Y[i];
  }

  for (j = 0; j < R; j++) {
    for (i = 0; i < R; i++) {
      x[i] = X[i];
      y[i] = Y[i];
    }
    x[0] = Y[j];
    y[j] = X[0];

    s1 = sort(x); // checks \chi(y1,x2,x3)*\chi(x1,y2,y3)
    s2 = sort(y);

    if (s1 != 0 &&
        s2 != 0) // s1==0 means that two of y1,x2,x3 are the same, its chirotope
                 // value is 0 and we want \chi(y1,x2,x3)*\chi(x1,y2,y3)<0
    {

      if (axB2(M, sign, s1, s2, ind(x), ind(y))) {
        free(x);
        free(y);
        return 1;
      }
    }
  }

  free(x);
  free(y);
  return 0;
}

char ischirotope(struct OM M) // checks chirotope axioms, see "Oriented
                              // matroids" BLSWZ, Definition 3.5.3 returns 0 if
                              // it is not a chirotope, 1 if it is a chirotope
{
  char i, j, k, l, p, q;
  long long int h;

  for (k = 0; k < nr_ints; k++)
    if (M.plus[k] &
        M.minus[k]) // if the same basis is both positive and negative
      return 0;

  for (k = 0; k < nr_ints; k++)
    if (M.plus[k] != 0 || M.minus[k] != 0) //(B0)
      break;

  if (k == nr_ints)
    return 0;

  char *x, *y; //(B2') Lemma 3.5.4
  x = (char *)malloc(R * sizeof(char));
  y = (char *)malloc(R * sizeof(char));

  char sign, limit_i, limit_j;
  int pi, pj, mi, mj;

  for (k = 0; k < nr_ints; k++) {
    if (k == nr_ints)
      limit_i = B & 31;
    else
      limit_i = 32;
    for (i = 0; i < limit_i; i++) {

      h = 1u << i;
      pi = M.plus[k] & h;  // pi!=0 iff \chi(bases[i])=1
      mi = M.minus[k] & h; // mi!=0 iff \chi(bases[i])=-1
      if ((pi == 0) &&
          (mi ==
           0)) //\chi(bases[i])=0, we do not have to worry about this basis
        continue;

      for (l = k; l < nr_ints; l++) {

        if (l == k)
          j = i + 1;
        else
          j = 0;
        if (l == nr_ints - 1)
          limit_j = B & 31;
        else
          limit_j = 32;
        for (j; j < limit_j; j++) {
          h = 1u << j;
          pj = M.plus[l] & h;
          mj = M.minus[l] & h;
          if (pj == 0 && mj == 0)
            continue;

          if ((pi && mj) ||
              (mi && pj)) //\chi(x_1,x_2,x_3)* \chi(y_1,y_2,y_3)=-1
            sign = -1;
          else if ((pi && pj) ||
                   (mi && mj)) //\chi(x_1,x_2,x_3)* \chi(y_1,y_2,y_3)=1
            sign = 1;

          for (p = 0; p < R; p++) // we have to check B2' for all permutations
                                  // of x1,x2,x3, but it suffices to have all
                                  // entries of x in the first position
          {
            for (q = 0; q < R; q++) {
              x[q] = bases[i + (k << 5)][q];
              y[q] = bases[j + (l << 5)][q];
            }
            x[0] = bases[i + (k << 5)][p];
            x[p] = bases[i + (k << 5)][0];

            if (p == 1)
              sign = -sign;

            if (b2prime(M, sign, x, y) == 0) // checks B2'
            {
              free(x);
              free(y);
              return 0;
            }
          }
        }
      }
    }
  }

  free(x);
  free(y);
  return 1;
}

void standardizeOM(struct OM *M) // we store OMs in such a way that the
                                 // lexicographically largest basis is positive
{
  char i;
  int x;

  for (i = nr_ints - 1; i >= 0; i--)
    if (M->plus[i] < M->minus[i])
      i = -2;
    else if (M->minus[i] < M->plus[i])
      i = -1;

  if (i != -2)
    for (i = 0; i < nr_ints; i++) {
      x = M->plus[i];
      M->plus[i] = M->minus[i];
      M->minus[i] = x;
    }
}

void permutations(
    char *p,
    int l) // makes all permutations of N elements and stores them in perm
{
  char i;
  if (l == N) {
    for (i = 0; i < N; i++)
      perm[w][i] = p[i];
    w++;
    return;
  }
  char j, x;
  for (j = l; j < N; j++) {
    x = p[l];
    p[l] = p[j];
    p[j] = x;
    permutations(p, l + 1);
    x = p[l];
    p[l] = p[j];
    p[j] = x;
  }
  return;
}

int factorial(int n) {
  int i, R;
  R = 1;
  for (i = 2; i <= n; i++)
    R *= i;
  return R;
}

void makepermutations() // makes all permutations on N
{
  char i;
  int s;

  perm = (char **)malloc(factorial(N) * sizeof(char *));
  for (s = 0; s < factorial(N); s++)
    perm[s] = (char *)malloc(N * sizeof(char));

  char p[N];
  for (i = 0; i < N; i++)
    p[i] = i;
  w = 0;
  permutations(&p[0], 0);
}

void removepermutations() // frees the memory previously allocated for
                          // permutations
{
  int s;
  for (s = 0; s < w; s++)
    free(perm[s]);
  free(perm);
}

struct OM permute(struct OM M,
                  char s[]) // given an OM, it transforms it into a new one -
                            // permutes the labels of the elements s[] is an
                            // array of length N that stores the permutation
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

  standardizeOM(&X); // for convinience, we always store chirotopes s.t. the
                     // lexicographically largest basis is positive

  return X;
}

int isfixed(
    struct OM M) // returns 1 if the OM M is fixed under the given group action
{
  int i;
  struct OM X;
  makeOM(&X);

  for (i = 1; i < sizeofgroup; i++) {
    X = permute(M, action[i]);
    if (isequal(X, M) == 0) {
      removeOM(&X);
      return 0;
    }
  }
  removeOM(&X);
  return 1;
}

void removegroupaction() {
  int i;
  for (i = 0; i < sizeofgroup; i++)
    free(action[i]);
  free(action);
}
