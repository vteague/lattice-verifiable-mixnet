#include "sha.h"
#include "test.h"
#include "bench.h"
#include "common.h"
#include "gaussian.h"

#define R       (HEIGHT+1)
#define V       (HEIGHT+3)
#define TAU     1000
#define NTI     130

static void pianex_hash(uint8_t h[SHA256HashSize], params::poly_q A[R][V], params::poly_q t[TAU][V], params::poly_q W[R][NTI]) {
	SHA256Context sha;
    std::array<mpz_t, params::poly_q::degree> coeffs;
	params::poly_q tmp;

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

	SHA256Reset(&sha);
	/* Hash public key. */
	for (size_t i = 0; i < R; i++) {
		for (int j = 0; j < V; j++) {
            A[i][j].poly2mpz(coeffs);
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                SHA256Input(&sha, (uint8_t *)coeffs[k]->_mp_d, coeffs[k]->_mp_size * sizeof(uint64_t));
            }
		}
	}
    for (size_t i = 0; i < TAU; i++) {
		for (int j = 0; j < V; j++) {
            t[i][j].poly2mpz(coeffs);
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                SHA256Input(&sha, (uint8_t *)coeffs[k]->_mp_d, coeffs[k]->_mp_size * sizeof(uint64_t));
            }
		}
	}
    for (size_t i = 0; i < R; i++) {
		for (int j = 0; j < NTI; j++) {
            W[i][j].poly2mpz(coeffs);
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                SHA256Input(&sha, (uint8_t *)coeffs[k]->_mp_d, coeffs[k]->_mp_size * sizeof(uint64_t));
            }
		}
	}

	SHA256Result(&sha, h);

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

/**
 * Test if the l2-norm is within bounds (BOUND_ANEX).
 *
 * @param[in] r 			- the polynomial to compute the l2-norm.
 * @return the computed norm.
 */
bool pianex_test_norm(params::poly_q r) {
    array<mpz_t, params::poly_q::degree> coeffs;
    mpz_t norm, qDivBy2, tmp, q;

    /// Constructors
    mpz_inits(norm, qDivBy2, tmp, q, nullptr);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    r.poly2mpz(coeffs);
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
    mpz_set_ui(norm, 0);
    mpz_set_str(q, PRIMEQ, 10);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_mod(coeffs[i], coeffs[i], q);
        mpz_mul(tmp, coeffs[i], coeffs[i]);
        mpz_add(norm, norm, tmp);
    }

    mpz_set_ui(tmp, 1);
    mpz_mul_2exp(tmp, tmp, BOUND_ANEX);
	int result = mpz_cmp(norm, tmp) <= 0;

    mpz_clears(norm, qDivBy2, tmp, q, nullptr);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }

    return result;
}


static void pianex_prover(uint8_t h[SHA256HashSize], params::poly_q Z[V][NTI],
        params::poly_q A[R][V], params::poly_q t[TAU][V], params::poly_q s[TAU][V]) {
    params::poly_q Y[V][NTI], W[R][NTI], C[TAU][NTI];
    std::array<mpz_t, params::poly_q::degree> coeffs;
    mpz_t qDivBy2;

    mpz_init(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    /* Prover samples Y from Gaussian. */
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < NTI; j++) {
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                int64_t coeff = discrete_gaussian(0.0);
                mpz_set_si(coeffs[k], coeff);
            }
            Y[i][j].mpz2poly(coeffs);
            Y[i][j].ntt_pow_phi();
        }
    }

    /* Prover computes W = AY. */
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < NTI; j++) {
            W[i][j] = 0;
            for (int k = 0; k < V; k++) {
                W[i][j] = W[i][j] + A[i][k] * Y[k][j];
            }
        }
    }

    pianex_hash(h, A, t, W);

    /* Sample challenge from RNG seeded with hash. */
    nfl::fastrandombytes_seed(h, SHA256HashSize);

    /* Verifier samples challenge matrix C. */
    for (int i = 0; i < TAU; i++) {
        for (int j = 0; j < NTI; j++) {
            C[i][j] = nfl::ZO_dist();
            C[i][j].ntt_pow_phi();
        }
    }

    nfl::fastrandombytes_reseed();

    /* Prover computes Z = Y + SC and performs rejection sampling. */
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < NTI; j++) {
            Z[i][j] = Y[i][j];
            for (int k = 0; k < TAU; k++) {
                Z[i][j] = Z[i][j] + s[k][i] * C[k][j];
            }
        }
    }

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
    mpz_clear(qDivBy2);
}

int pianex_verifier(uint8_t h1[SHA256HashSize], params::poly_q Z[V][NTI],
        params::poly_q A[R][V], params::poly_q t[TAU][V]) {
    params::poly_q W[R][NTI], C[TAU][NTI];
    uint8_t h2[SHA256HashSize];
    int result;

    /* Sample challenge from RNG seeded with hash. */
    nfl::fastrandombytes_seed(h1, SHA256HashSize);

    /* Verifier samples challenge matrix C. */
    for (int i = 0; i < TAU; i++) {
        for (int j = 0; j < NTI; j++) {
            C[i][j] = nfl::ZO_dist();
            C[i][j].ntt_pow_phi();
        }
    }

    /* Verifier checks that W = AZ - TC. */
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < NTI; j++) {
            W[i][j] = 0;
            for (int k = 0; k < V; k++) {
                W[i][j] = W[i][j] + A[i][k] * Z[k][j];
            }
        }
    }
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < NTI; j++) {
            for (int k = 0; k < TAU; k++) {
                W[i][j] = W[i][j] - t[k][i] * C[k][j];
            }
        }
    }

    pianex_hash(h2, A, t, W);

    result = memcmp(h1, h2, SHA256HashSize) == 0;
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < V; j++) {
            Z[i][j].invntt_pow_invphi();
            result &= pianex_test_norm(Z[i][j]);
        }
    }

    return result;
}

static void test() {
    params::poly_q A[R][V], s[TAU][V], t[TAU][V], Z[V][NTI];
    uint8_t h1[SHA256HashSize];
    std::array<mpz_t, params::poly_q::degree> coeffs;
    gmp_randstate_t prng;
	mpz_t q;

    mpz_init(q);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    /* Create instances. */
    gmp_randinit_default(prng);
    mpz_set_str(q, PRIMEQ, 10);
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < V; j++) {
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                mpz_urandomb(coeffs[k], prng, LEVEL);
                mpz_mod(coeffs[k], coeffs[k], q);
            }
            A[i][j].mpz2poly(coeffs);
        }
    }

    /* Create a total of TAU relations t_i = A * s_i */
    for (int i = 0; i < TAU; i++) {
        for (int j = 0; j < V; j++) {
            s[i][j] = nfl::ZO_dist();
            s[i][j].ntt_pow_phi();
        }
        for (int j = 0; j < R; j++) {
            t[i][j] = 0;
            for (int k = 0; k < V; k++) {
                t[i][j] = t[i][j] + A[j][k] * s[i][k];
            }
        }
    }

    TEST_ONCE("ANEX proof is consistent") {
        pianex_prover(h1, Z, A, t, s);
		TEST_ASSERT(pianex_verifier(h1, Z, A, t) == 1, end);
	} TEST_END;

end:

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
    mpz_clear(q);
    gmp_randclear(prng);
    return;
}

static void bench() {
    params::poly_q A[R][V], s[TAU][V], t[TAU][V], Z[V][NTI];
    uint8_t h1[SHA256HashSize];

    /* Create instances. */
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < V; j++) {
            A[i][j] = nfl::uniform();
        }
    }

    /* Create a total of TAU relations t_i = A * s_i */
    for (int i = 0; i < TAU; i++) {
        for (int j = 0; j < V; j++) {
            s[i][j] = nfl::ZO_dist();
            s[i][j].ntt_pow_phi();
        }
        for (int j = 0; j < R; j++) {
            t[i][j] = 0;
            for (int k = 0; k < V; k++) {
                t[i][j] = t[i][j] + A[j][k] * s[i][k];
            }
        }
    }

	BENCH_SMALL("ANEX prover (N relations)", pianex_prover(h1, Z, A, t, s));
	BENCH_SMALL("ANEX verifier (N relations)", pianex_verifier(h1, Z, A, t));
}

int main(int argc, char *argv[]) {
	printf("\n** Tests for lattice-based ANEX proof:\n\n");
	test();

	printf("\n** Benchmarks for lattice-based ANEX proof:\n\n");
	bench();
	printf("\nMultiply prover by 3 due to rejection sampling.");
}
