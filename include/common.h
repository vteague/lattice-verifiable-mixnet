#include <cstddef>

#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <nfl.hpp>

#include "bench.h"
#include <sys/random.h>

using namespace std;

/* Parameter v in the commitment  scheme (laximum l1-norm of challs). */
#define NONZERO     36
/* Security level to attain. */
#define LEVEL       128
/* The \infty-norm bound of certain elements. */
#define BETA 	    1
/* Width k of the comming matrix. */
#define WIDTH 	    4
/* Height of the commitment matrix. */
#define HEIGHT 	    1
/* Dimension of the committed messages. */
#ifndef SIZE
#define SIZE        2
#endif
/* Large modulus. */
#define PRIMEQ      "302231454903657293688833"
/* Small modulus. */
#define PRIMEP      3
/* Degree of the irreducible polynomial. */
#define DEGREE      4096
/* Sigma for the commitment gaussian distribution. */
#define SIGMA_C     (1u << 12)
/* Sigma for the approximate amortized proof. */
#define SIGMA_ANEX  (1e66l)
/* Parties that run the distributed decryption protocol. */
#define PARTIES     4
/* Security level for Distributed Decryption. */
#define BGVSEC      40
/* Bound for Distributed Decryption = 2^BGVSEC * q/(2 * p * PARTIES). */
#define BDISTD      "750837175903336127688539820910095018"
/* Norm bound for PI_ANEX */
#define BOUND_ANEX   145

namespace params {
    using poly_p = nfl::poly_from_modulus<uint32_t, DEGREE, 30>;
    using poly_q = nfl::poly_from_modulus<uint64_t, DEGREE, 124>;
    using poly_big = nfl::poly_from_modulus<uint64_t, 4 * DEGREE, 124>;
}

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

/* Class that represents a commitment key pair. */
class commitkey_t {
    public:
       params::poly_q A1[HEIGHT][WIDTH - HEIGHT];
       params::poly_q A2[SIZE][WIDTH];
};

/* Class that represents a commitment in CRT representation. */
class commit_t {
    public:
      params::poly_q c1;
      vector<params::poly_q> c2;
};

/* Class that represents a BGV key pair. */
class bgvkey_t {
    public:
       params::poly_q a;
       params::poly_q b;
};

class cipher_t {
    public:
       params::poly_q u;
       params::poly_q v;
};

#include "util.hpp"

void commit_sample(vector<params::poly_q>& r);
void commit_sample_chall(params::poly_q& f);
bool commit_test_norm(params::poly_q r, uint64_t sigma_sqr);
void commit_doit(commit_t& com, vector<params::poly_q> m, commitkey_t& key, vector<params::poly_q> r);
int commit_open(commit_t& com, vector<params::poly_q> m, commitkey_t& key, vector<params::poly_q> r, params::poly_q& f);
void commit_keygen(commitkey_t& key);

void encrypt_sample_message(params::poly_p& r);
void encrypt_sample_short(params::poly_q& r);
void encrypt_keygen(bgvkey_t& pk, params::poly_q& sk);
void encrypt_doit(cipher_t &c, bgvkey_t& pk, params::poly_p& m);
void encrypt_undo(params::poly_p& m, cipher_t& c, params::poly_q & sk);
