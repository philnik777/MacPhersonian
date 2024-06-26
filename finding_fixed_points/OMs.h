#ifndef _OMs_H_
#define _OMs_H_
#endif

// Oriented matroids of rank R on N elements

struct OM {
  unsigned int *plus;
  unsigned int *minus;
};

void makebases(); // makes the list of all posible bases a chirotope could have,
                  // bases are 012, 013,...

void removebases(); // frees the memory used by bases

void makeOM(struct OM *M); // allocates memory for an oriented matroid of rank R
                           // on N elements

void removeOM(struct OM *M); // frees the memory used by an oriented matroid

void showbits(unsigned int *plus); // prints a list if integers in the binary
                                   // representation (smallest bit on the right)

void showchirotope(struct OM M); // prints a chirotope

int countbases(struct OM M); // counts the number of bases of a chirotope

int weakmap(
    struct OM M1,
    struct OM M2); // checks whether there is a weak map M_1 \wm M_2
                   // this actually checks weak maps for oriented matroids,
                   // since we make only one chirotope from each pair chi, -chi

int isequal(
    struct OM M1,
    struct OM M2); // returns 1 if M1 and M2 are the same as oriented matroids

int ind(char *a); // returns the index of the basis (a[0],a[1],...,a[R-1]) in
                  // the array bases[][], assumes that a[0]<a[1]<...<a[R]

int sort(char *a); // sorts integers in the array and returns the sign of the
                   // permutation

char axB2(struct OM M, char sign, char s1, char s2, int in1,
          int in2); // used in b2prime to check Axiom B2' in ischirotope

char b2prime(struct OM M, char sign, char *X,
             char *Y); // checks Axiom B2' of BLSWZ, Lemma 3.5.4

char ischirotope(struct OM M); // checks chirotope axioms, see "Oriented
                               // matroids" BLSWZ, Definition 3.5.3

void standardizeOM(struct OM *M); // we store OMs in such a way that the largest
                                  // basis is positive

void permutations(char *p, int l); // recursively makes all permutations of N
                                   // elements and stores them in perm

void makepermutations(); // makes all permutations on N

void removepermutations(); // frees the memory previously allocated for
                           // permutations

int factorial(int n); // computes n!

struct OM
permute(struct OM M,
        char s[]); // given an OM, it transforms it into a new one - permutes
                   // the labels of the elements, the permutation is given by s

int isfixed(
    struct OM M); // checks whether this OM is fixed under the group action

void removegroupaction(); // frees the memory that is allocated for the group
                          // action
