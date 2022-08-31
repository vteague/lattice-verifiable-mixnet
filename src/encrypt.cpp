#include "test.h"
#include "bench.h"
#include "common.h"
#include <sys/random.h>

namespace params {
    using gauss_struct = nfl::gaussian<uint16_t, uint64_t, 2>;
    using gauss_t = nfl::FastGaussianNoise<uint16_t, uint64_t, 2>;
    gauss_t fg_prng_enc(8.0, LEVEL, DEGREE);
}

void encrypt_sample_message(params::poly_p& r) {
    // Sample a short polynomial.
    std::array<mpz_t, DEGREE> coeffs;
    uint64_t buf;

    size_t bits_in_moduli_product = params::poly_p::bits_in_moduli_product();
    for (size_t i = 0; i < params::poly_p::degree; i++) {
        mpz_init2(coeffs[i], bits_in_moduli_product << 2);
    }

	for (size_t j = 0; j < params::poly_p::degree; j += 32) {
		getrandom(&buf, sizeof(buf), 0);
		for (size_t k = 0; k < 64; k += 2) {
            mpz_set_ui(coeffs[j+k/2], (buf >> k) % 3);
		}
	}
    r.mpz2poly(coeffs);

    for (size_t i = 0; i < params::poly_p::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

void encrypt_sample_short(params::poly_q& r) {
    r = nfl::ZO_dist();
}

// Generate a key pair.
void encrypt_keygen(bgvkey_t& pk, params::poly_q& sk) {
    params::poly_q e = {params::gauss_struct(&params::fg_prng_enc)};
    pk.a = nfl::uniform();

    encrypt_sample_short(sk);
    sk.ntt_pow_phi();
    e.ntt_pow_phi();

    pk.b = pk.a * sk + (e + e + e);
}

void encrypt_keyshare(params::poly_q s[], size_t shares, params::poly_q& sk) {
    params::poly_q t = sk;
    for (size_t i = 1; i < shares; i++) {
        s[i] = nfl::uniform();
        s[i].ntt_pow_phi();
        t = t - s[i];
    }
    s[0] = t;
}

void encrypt_doit(cipher_t &c, bgvkey_t& pk, params::poly_p& m) {
    std::array<mpz_t, DEGREE> coeffs;
    params::poly_q e1{params::gauss_struct(&params::fg_prng_enc)};
    params::poly_q e2{params::gauss_struct(&params::fg_prng_enc)};
    params::poly_q r;

    e1.ntt_pow_phi();
    e2.ntt_pow_phi();
    encrypt_sample_short(r);
    r.ntt_pow_phi();

    c.u = pk.a * r + (e1 + e1 + e1);
    c.v = pk.b * r + (e2 + e2 + e2);

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], params::poly_q::bits_in_moduli_product() << 2);
    }

    m.poly2mpz(coeffs);
    r.mpz2poly(coeffs);
    r.ntt_pow_phi();
    c.v = c.v + r;

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

void encrypt_undo(params::poly_p& m, cipher_t& c, params::poly_q & sk) {
    std::array<mpz_t, DEGREE> coeffs;
    params::poly_q t = c.v - sk * c.u;
    mpz_t qDivBy2;

    mpz_init(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], params::poly_q::bits_in_moduli_product() << 2);
    }
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    // Reduce the coefficients
    t.invntt_pow_invphi();
    t.poly2mpz(coeffs);
	for (size_t i = 0; i < params::poly_q::degree; i ++) {
        util::center(coeffs[i], coeffs[i], params::poly_q::moduli_product(), qDivBy2);
        mpz_mod_ui(coeffs[i], coeffs[i], 3);
	}
    m.mpz2poly(coeffs);

    mpz_clear(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

void encrypt_ddec(params::poly_q& tj, cipher_t& c, params::poly_q& sj) {
    std::array<mpz_t, DEGREE> coeffs;
    params::poly_q mj, Ej;
    mpz_t qDivBy2, bound;

    mpz_init(qDivBy2);
    mpz_init(bound);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], params::poly_q::bits_in_moduli_product() << 2);
    }
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
    mpz_set_str(bound, PRIMEQ, 10);
    mpz_mul_2exp(bound, bound, BGVSEC);
    mpz_fdiv_q_ui(bound, bound, PRIMEP);
    mpz_fdiv_q_ui(bound, bound, PARTIES);

    Ej = nfl::uniform();
    Ej.poly2mpz(coeffs);
    for (size_t i = 0; i < params::poly_q::degree; i ++) {
        util::center(coeffs[i], coeffs[i], params::poly_q::moduli_product(), qDivBy2);
        mpz_mod(coeffs[i], coeffs[i], bound);
	}
    Ej.mpz2poly(coeffs);
    Ej.ntt_pow_phi();
    mj = sj * c.u;
    tj = mj + (Ej + Ej + Ej);
}

void encrypt_comb(params::poly_p& m, cipher_t& c, params::poly_q t[], size_t shares) {
    std::array<mpz_t, DEGREE> coeffs;
    params::poly_q v;
    mpz_t qDivBy2;

    mpz_init(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], params::poly_q::bits_in_moduli_product() << 2);
    }
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    v = c.v - t[0];
    for (size_t i = 1; i < shares; i++) {
        v = v - t[i];
    }
    v.invntt_pow_invphi();
    // Reduce the coefficients modulo p
    v.poly2mpz(coeffs);
	for (size_t i = 0; i < params::poly_q::degree; i ++) {
        util::center(coeffs[i], coeffs[i], params::poly_q::moduli_product(), qDivBy2);
        mpz_mod_ui(coeffs[i], coeffs[i], 3);
	}
    m.mpz2poly(coeffs);

    mpz_clear(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

#ifdef MAIN
static void test() {
    bgvkey_t pk;
    params::poly_q sk, s[PARTIES], t[PARTIES], acc;
    params::poly_p m, _m;
    cipher_t c;

    encrypt_keygen(pk, sk);
    encrypt_sample_message(m);

    TEST_BEGIN("BGV encryption is consistent") {
        encrypt_doit(c, pk, m);
        encrypt_undo(_m, c, sk);
        TEST_ASSERT(m - _m == 0, end);
        encrypt_sample_message(m);
        encrypt_undo(_m, c, sk);
        TEST_ASSERT(m - _m != 0, end);
        encrypt_doit(c, pk, m);
        encrypt_keygen(pk, sk);
        encrypt_undo(_m, c, sk);
        TEST_ASSERT(m - _m != 0, end);
	} TEST_END;

    TEST_BEGIN("BGV distributed decryption is consistent") {
        encrypt_keygen(pk, sk);
        encrypt_doit(c, pk, m);
        encrypt_keyshare(s, PARTIES, sk);
        acc = s[0];
        for (size_t j = 1; j < PARTIES; j++) {
            acc = acc + s[j];
        }
        TEST_ASSERT(sk - acc == 0, end);
        for (size_t j = 0; j < PARTIES; j++) {
            encrypt_ddec(t[j], c, s[j]);
        }
        encrypt_comb(_m, c, t, PARTIES);
        TEST_ASSERT(m - _m == 0, end);
    } TEST_END;

end:
    return;
}

static void bench() {
    bgvkey_t pk;
    params::poly_q sk, s[PARTIES], t[PARTIES], acc;
    params::poly_p m, _m;
    cipher_t c;

    encrypt_keygen(pk, sk);

    BENCH_BEGIN("encrypt_sample_message") {
		BENCH_ADD(encrypt_sample_message(m));
	} BENCH_END;

	BENCH_BEGIN("encrypt_doit") {
		BENCH_ADD(encrypt_doit(c, pk, m));
	} BENCH_END;

	BENCH_BEGIN("encrypt_undo") {
        encrypt_doit(c, pk, m);
		BENCH_ADD(encrypt_undo(_m, c, sk));
	} BENCH_END;

    encrypt_keyshare(t, PARTIES, sk);
    BENCH_BEGIN("encrypt_ddec") {
        encrypt_doit(c, pk, m);
		BENCH_ADD(encrypt_ddec(t[0], c, s[0]));
	} BENCH_END;

    for (size_t i = 1; i < PARTIES; i++) {
        encrypt_ddec(t[i], c, s[i]);
    }

    BENCH_BEGIN("encrypt_comb") {
        BENCH_ADD(encrypt_comb(_m, c, t, PARTIES));
    } BENCH_END;
}

int main(int argc, char *arv[]) {
	printf("\n** Tests for BGV encryption:\n\n");
	test();

	printf("\n** Benchmarks for BGV encryption:\n\n");
	bench();

    return 0;
}
#endif
