#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "murmur3.h"
#include "gtf.h"

//These are just modified versions of the functions in libGTF

uint64_t hashString2(char *s, int len) {
    uint64_t hash_val[2];
    uint32_t seed = 0xAAAAAAAA;
#if UINTPTR_MAX == 0xffffffff
    MurmurHash3_x86_128((void *) s, len, seed, (void *) &hash_val);
#else
    MurmurHash3_x64_128((void *) s, len, seed, (void *) &hash_val);
#endif
    return hash_val[0];
}

static void rehashElement2(hashTable *ht, hashTableElement *e, int len) {
    hashTableElement *next = e->next;
    if(!e) return;
    uint64_t hash = hashString2(ht->str[e->val], len);
    e->next = NULL;
    insertHTelement(ht, e, hash);
    if(next) rehashElement2(ht, next, len);
}

static void rehashHT2(hashTable *ht, int len) {
    int32_t i;
    hashTableElement *e;
    for(i=0; i<ht->l; i++) {
        if(ht->elements[i]) {
            e = ht->elements[i];
            ht->elements[i] = NULL;
            rehashElement2(ht, e, len);
        }
    }
}

static void growHT2(hashTable *ht, int len) {
    int i;
    ht->m = ht->l+1;
    kroundup32(ht->m);
    ht->str = realloc(ht->str, ht->m*sizeof(char*));
    assert(ht->str);
    ht->elements = realloc(ht->elements, ht->m*sizeof(hashTableElement*));

    for(i=ht->l; i<ht->m; i++) {
        ht->str[i] = NULL;
        ht->elements[i] = NULL;
    }
    rehashHT2(ht, len);
}

int32_t addHTelement2(hashTable *ht, char *s, int len) {
    uint64_t hash = hashString2(s, len);
    int32_t val = ht->l++;
    if(ht->l >= ht->m) growHT2(ht, len);
    ht->str[val] = strndup(s, len);

    hashTableElement *e = calloc(1, sizeof(hashTableElement));
    assert(e);
    e->val = val;
    insertHTelement(ht, e, hash);
    return val;
}

int strExistsHT2(hashTable *ht, char *s, int len) {
    uint64_t h = hashString2(s, len);
    hashTableElement *curr = ht->elements[h%ht->m];
    while(curr) {
        if(strncmp(ht->str[curr->val], s, len) == 0) return 1;
        curr = curr->next;
    }
    return 0;
}

int32_t str2valHT2(hashTable *ht, char *s, int len) {
    uint64_t h = hashString2(s, len);
    hashTableElement *curr = ht->elements[h%ht->m];
    while(curr) {
        if(strncmp(ht->str[curr->val], s, len) == 0) return curr->val;
        curr = curr->next;
    }
    return -1;
}

void growCntTable(cntTable *ct, int len) {
    growHT2(ct->ht, len);
    ct->cnts = realloc(ct->cnts, sizeof(uint32_t)*(ct->ht->m));
    assert(ct->cnts);
    int i;
    for(i=ct->ht->l; i<ct->ht->m; i++) ct->cnts[i] = 0;
}

void pushIncCntTable(cntTable *ct, char *s, int len) {
    int32_t idx;
    if(strExistsHT2(ct->ht, s, len)) {
        idx = str2valHT2(ct->ht, s, len);
        ct->cnts[idx]++;
    } else {
        //Grow if needed
        if(ct->ht->l+1 >= ct->ht->m) growCntTable(ct, len);
        idx = addHTelement2(ct->ht, s, len);
        ct->cnts[idx] = 1;
    }
}

int32_t str2idx(cntTable *ct, char *s) {
    if(!s) return -1;
    uint64_t h = hashString(s);
    hashTableElement *curr = ct->ht->elements[h%ct->ht->m];
    while(curr) {
        if(strcmp(ct->ht->str[curr->val], s) == 0) return curr->val;
        curr = curr->next;
    }
    return -1;
}

//A uniqueSet object holds the index in a hashTable to the string representation
//We can take that string, get its index in a cntTable and then set the
//appropriate bit
char *us2char(cntTable *ct, hashTable *ht, uniqueSet *us, int len) {
    char *o = (char *) calloc(len, 1);
    int byte, bit, i, idx;

    for(i=0; i<us->l; i++) {
        idx = str2idx(ct, val2strHT(ht, us->IDs[i]));
        if(idx == -1) continue;
        byte = idx>>3;
        assert(byte<len);
        bit = idx%8;
        o[byte] |= 1<<bit;
    }
    return o;
}

int getVecLen(cntTable *ct) {
    return ct->ht->l>>3+(ct->ht->l&7)?1:0;
}

//l is the output of getVecLen()
uint32_t uniqueVal(cntTable *ct, int l, int32_t i) {
    char *o = (char *) calloc(l, 1);
    uint32_t rv = 0;
    int byte = i>>3;
    int bit = i%8;
    o[byte] |= 1<<bit;

    if(strExistsHT2(ct->ht, o, l)) {
        rv = ct->cnts[str2valHT2(ct->ht, o, l)];
    }
    free(o);
    return rv;
}
