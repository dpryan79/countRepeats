#include "gtf.h"
#include <assert.h>
#include <gsl/gsl_multimin.h>
//The user should be able to change these
#define MINPROB 0.02
#define MINEPS 1e-3
#define MAXITER 100

double logit(double prob) {
    if(prob<MINPROB) prob = MINPROB;
    else if(prob>1-MINPROB) prob = 1-MINPROB;
    return -log((1/prob)-1);
}

double invlogit(double prob) {
    if(prob<logit(MINPROB)) return 0;
    else if(prob>logit(1-MINPROB)) return 1;
    return exp(prob)/(1+exp(prob));
}

//Is the ith bit set?
int bitSet(char *s, int i) {
    uint8_t byte = s[i>>3];
    return byte&(1<<(i&7));
}

//ct is the cntHashTable, ct2 the associated hashTable, aka the coefficients
//This could be sped up a bit, since the vectors are sparse, but we'll leave that for later
gsl_matrix *cntHash2matrix(cntTable *ct, cntTable *ct2) {
    char *s;
    int i,j;
    gsl_matrix *m = gsl_matrix_alloc(ct->ht->l, ct2->ht->l);
    assert(m);
    gsl_matrix_set_all(m, logit(MINPROB));
    for(i=0; i<ct->ht->l; i++) {
        s = val2strHT(ct->ht, i);
        for(j=0; j<ct2->ht->l; j++) {
            if(bitSet(s, j)) {
                gsl_matrix_set(m, i, j, logit(ct->cnts[i]));
            }
        }
    }
    return m;
}

typedef struct {
    gsl_matrix *m;
    gsl_vector *l_r;
} min_params_struct;

//The minus log likelihood
/*
params needs to hold:
  c_e, the count matrix
  l_r, the vector of lengths
*/
double mLL(const gsl_vector *v, void *params) {
    min_params_struct *mps = (min_params_struct *) params;
    int i, j;
    double out = 0.0, val;
    gsl_vector *v2 = gsl_vector_alloc(v->size);
    assert(v2);
    gsl_vector_memcpy(v2, v);

    for(i=0; i<v->size; i++) {
        for(j=0; j<mps->m->size2; j++) {
            val = gsl_matrix_get(mps->m, i, j);
            if(val) {
                out += val*invlogit(gsl_vector_get(v,i))/gsl_vector_get(mps->l_r,i);
            }
        }
    }
    gsl_vector_free(v2);
    if(out <= 0.0) return GSL_NAN;
    return -log(out);
}

//The gradient of above
void dmLL(const gsl_vector *v, void *params, gsl_vector *df) {
    int i, j;
    double cnt, num = 0.0, denom = 0.0;
    min_params_struct *mps = (min_params_struct *) params;

    for(i=0; i<v->size; i++) {
        num = denom = 0.0;
        for(j=0; j<mps->m->size2; j++) {
            cnt = gsl_matrix_get(mps->m, i, j);
            if(cnt) {
                num += cnt/gsl_vector_get(mps->l_r, i);
                denom += cnt*invlogit(gsl_vector_get(v,i))/gsl_vector_get(mps->l_r,i);
            }
        }
        if(num <= 0.0) gsl_vector_set(df, i, GSL_NAN);
        else gsl_vector_set(df, i, -log(num)/denom);
    }
}

//The minimization driver function
void mLL_dmLL(const gsl_vector *x, void *params, double *f, gsl_vector *df) {
    *f = mLL(x, params);
    dmLL(x, params, df);
}

void initializeCoefs(gsl_vector *coefs, cntTable *ct, uint32_t N, int l) {
    int32_t i;
    uint32_t val;
    gsl_vector_set_all(coefs, logit(0));
    for(i=0; i<coefs->size; i++) { //This is kind of slow...
        val = uniqueVal(ct, l, i);
        //get the value
        if(val) gsl_vector_set(coefs, i, logit(val/(double) N));
    }
}

//ct is the cntHashTable, ct2 the associated cntTable, l_r the vector of lengths
//and N the number of counted clusters
gsl_vector *performEM(cntTable *ct, cntTable *ct2, gsl_vector *l_r, uint32_t N) {
    int i = 0, rv;
    gsl_multimin_function_fdf EM_func;
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    gsl_matrix *m = cntHash2matrix(ct, ct2);
    assert(m);
    gsl_vector *coefs = gsl_vector_alloc(ct2->ht->l);
    assert(coefs);
    min_params_struct *mps = malloc(sizeof(min_params_struct));
    assert(mps);
    mps->m = m;
    mps->l_r = l_r;

    //Initialize the coefficients
    initializeCoefs(coefs, ct, N, getVecLen(ct));

    EM_func.n = l_r->size;
    EM_func.f = &mLL;
    EM_func.df = &dmLL;
    EM_func.fdf = &mLL_dmLL;
    EM_func.params = (void *) mps;

    //This should be checked
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, l_r->size);
    gsl_multimin_fdfminimizer_set(s, &EM_func, coefs, 0.01, 1e-4);

      //The actual optimization
    do {
        rv = gsl_multimin_fdfminimizer_iterate(s);
        if(rv) break;
        rv = gsl_multimin_test_gradient(s->gradient, MINEPS);
        i++;
    } while (rv == GSL_CONTINUE && i < MAXITER);

    //Free up space
    gsl_matrix_free(m);
    free(mps);
    gsl_multimin_fdfminimizer_free(s);

    //Convert the coefs to expected counts
    for(i=0; i<coefs->size; i++) {
        gsl_vector_set(coefs, i, \
            invlogit(gsl_vector_get(coefs,i))*N);
    }
    return coefs;
}

//Get l_r and N
