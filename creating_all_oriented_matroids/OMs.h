#ifndef _OMs_H_
#define _OMs_H_

#include <mdspan>
#include <memory>
#include <stdio.h>

// Oriented matroids of rank R on N elements

// The code works only for B<=64!
inline constexpr int R = 2; // rank
inline constexpr int N = 8; // the number of elements

constexpr int calculate_bases_count() {
  int bases_count = 1;
  for (int i = N; i > N - R; i--)
    bases_count *= i;
  for (int i = 2; i <= R; i++)
    bases_count /= i;

  return bases_count;
}

inline constexpr int B = calculate_bases_count(); // the number of bases

constexpr int calculate_nr_ints() {
  int nr_ints = B >> 5;

  if (B & 31)
    nr_ints++;

  return nr_ints;
}

// the number of integers needed to store the plus (resp.
// minus) of a chirotope
inline constexpr int nr_ints = calculate_nr_ints();

static_assert(nr_ints <= 4, "Requires more ints than available in OM");

// makes the list of all posible bases a chirotope could have, bases are 012,
// 013,...
constexpr std::array<char, B * R> makebases() {
  std::array<char, B * R> ret;

  auto bases = std::mdspan{ret.data(), B, R};

  for (int i = 0; i < R; i++)
    bases[0, i] = i;
  int s = 1;

  while (1) {
    int k = R - 1;
    while (k >= 0 && bases[s - 1, k] >= N - R + k)
      k--;
    if (k == -1)
      break;
    for (int i = 0; i < k; i++)
      bases[s, i] = bases[s - 1, i];
    for (int i = k; i < R; i++)
      bases[s, i] = bases[s - 1, k] + i - k + 1;
    s++;
  }

  return ret;
}

constexpr inline auto bases_backing = makebases();

constexpr inline std::mdspan<const char, std::dextents<size_t, 2>> bases =
    std::mdspan{bases_backing.data(), B, R}; // the list of bases

struct OM {
  OM() {
    for (int i = 0; i < nr_ints; i++) {
      plus[i] = 0;
      minus[i] = 0;
    }
  }

  unsigned int plus[4];
  unsigned int minus[4];
};

// prints a list if integers in the binary representation (smallest bit on the
// right)
void showbits(unsigned int *plus);

// prints a chirotope
void showchirotope(const OM &M, FILE *out = stdout);

// counts the number of bases of a chirotope
int countbases(struct OM M);

// checks whether there is a weak map M_1 \wm M_2
// this actually checks weak maps for oriented matroids, since we make only one
// chirotope from each pair chi, -chi
int weakmap(struct OM M1, struct OM M2);

// returns 1 if M1 and M2 are the same as oriented matroids
int isequal(const OM &M1, const OM &M2);

// returns the index of the basis (a[0],a[1],...,a[R-1]) in the array bases[][],
// assumes that a[0]<a[1]<...<a[R]
int ind(char *a);

// sorts integers in the array and returns the sign of the permutation
int sort(char *a);

// used in b2prime to check Axiom B2' in ischirotope
char axB2(const OM &M, char sign, char s1, char s2, int in1, int in2);

// checks Axiom B2' of BLSWZ, Lemma 3.5.4
char b2prime(const OM &M, char sign, char *X, char *Y);

// checks chirotope axioms, see "Oriented matroids" BLSWZ, Definition 3.5.3
char ischirotope(const OM &M);

// we store OMs in such a way that the largest basis is positive
void standardizeOM(struct OM *M);

// recursively makes all permutations of N elements and stores them in perm
void permutations(char *p, int l);

// makes all permutations on N
void makepermutations();

// computes n!
int factorial(int n);

// given an OM, it transforms it into a new one - permutes the labels of the
// elements, the permutation is given by s
OM permute(const OM &M, char s[]);

// checks whether this OM is fixed under the group action
int isfixed(struct OM M);

// frees the memory that is allocated for the group action
void removegroupaction();

void writeOM(const OM &, FILE *);
int readOM(struct OM *, FILE *);

#endif
