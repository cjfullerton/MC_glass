#ifndef RANDOM_H
#define RANDOM_H

/* random number stuff up here */

/* MT Period parameters */  
#define def_N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
class JRand_dyn {

  private:
    int N;
    unsigned long *mt;  /* the array for the state vector  */
    int mti;                  
                             /* mti==N+1 means mt[N] is not initialized */
    void init_genrand(unsigned long s);
    unsigned long genrand_int32(void);

  public:
    JRand_dyn() { N=def_N; mti=N; mt=new unsigned long[N]; }
    void srandom(int seed) { init_genrand(seed); }
    int    random(void) { return genrand_int32() & 0x7fffffff; }
    double prob_gen(void) { return random()/(double)0x7fffffff; }
    char   random_bin(void) { return (random() & 8 ? -1 : 0); }

};
#undef def_N

#endif
