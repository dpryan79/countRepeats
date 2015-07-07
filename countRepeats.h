#include <gsl/gsl_vector.h>

void pushIncCntTable(cntTable *ct, char *s, int len);
char *us2char(cntTable *ct, hashTable *ht, uniqueSet *us, int len);
int getVecLen(cntTable *ct);
uint32_t uniqueVal(cntTable *ct, int l, int32_t i);

/*
EM methods
\alpha_r must be logit transformed with some min/max value

-LL(alpha_r)= -log(\sum\limits_{e \in E} c_e \sum\limits_{r\in e} \frac{\alpha_r}{l_r})

dLL(alpha_r) = -\frac{log(\sum\limits_{e \in E} c_e \sum\limits_{r\in e} \frac{1}{l_r})}{\sum\limits_{e\in E} c_e \sum\limits_{r\in e} \frac{\alpha_r}{l_r}}

\alpha_r is the logit transformed probability
e is the row in the matrix E
r is the repeat name/class/family
l_r is the combined length of a repeat name/class/family
N is the total number of alignments/groups
\alpha_r N is the expected count

Can stop the optimization when the change is <1% or so
*/
gsl_vector *performEM(cntTable *ct, cntTable *ct2, gsl_vector *l_r, uint32_t N);
typedef struct {
    char *c_Name, *c_Class, *c_Family;
    int l_Name, l_Class, l_Family;
} overlapsRepeats_struct;

