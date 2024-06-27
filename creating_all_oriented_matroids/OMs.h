#ifndef _OMs_H_
#define _OMs_H_

#include <stdio.h>

// Oriented matroids of rank R on N elements

extern int R; // rank
extern int N; // the number of elements

extern int B;       // the number of bases
extern int nr_ints; // the number of integers needed to store the plus (resp. minus) of
             // a chirotope

extern char **bases; // the list of bases

struct OM {
  unsigned int *plus;
  unsigned int *minus;
};

// makes the list of all posible bases a chirotope could have, bases are 012,
// 013,...
void makebases();

// frees the memory used by bases
void removebases();

// allocates memory for an oriented matroid of rank R on N elements
void makeOM(struct OM *M);

// frees the memory used by an oriented matroid
void removeOM(struct OM *M);

// prints a list if integers in the binary representation (smallest bit on the
// right)
void showbits(unsigned int *plus);

// prints a chirotope
void showchirotope(struct OM M, FILE* out = stdout);

// counts the number of bases of a chirotope
int countbases(struct OM M);

// checks whether there is a weak map M_1 \wm M_2
// this actually checks weak maps for oriented matroids, since we make only one
// chirotope from each pair chi, -chi
int weakmap(struct OM M1, struct OM M2);

// returns 1 if M1 and M2 are the same as oriented matroids
int isequal(struct OM M1, struct OM M2);

// returns the index of the basis (a[0],a[1],...,a[R-1]) in the array bases[][],
// assumes that a[0]<a[1]<...<a[R]
int ind(char *a);

// sorts integers in the array and returns the sign of the permutation
int sort(char *a);

// used in b2prime to check Axiom B2' in ischirotope
char axB2(struct OM M, char sign, char s1, char s2, int in1, int in2);

// checks Axiom B2' of BLSWZ, Lemma 3.5.4
char b2prime(struct OM M, char sign, char *X, char *Y);

// checks chirotope axioms, see "Oriented matroids" BLSWZ, Definition 3.5.3
char ischirotope(struct OM M);

// we store OMs in such a way that the largest basis is positive
void standardizeOM(struct OM *M);

// recursively makes all permutations of N elements and stores them in perm
void permutations(char *p, int l);

// makes all permutations on N
void makepermutations();

// frees the memory previously allocated for permutations
void removepermutations();

// computes n!
int factorial(int n);

// given an OM, it transforms it into a new one - permutes the labels of the
// elements, the permutation is given by s
struct OM permute(struct OM M, char s[]);

// checks whether this OM is fixed under the group action
int isfixed(struct OM M);

// frees the memory that is allocated for the group action
void removegroupaction();

void writeOM(struct OM, FILE *);
int readOM(struct OM *, FILE *);

#endif
