#include <bits/stdc++.h>
#ifdef _WIN32
#include <fcntl.h>
#endif

#include <immintrin.h>

#define BUFFER_BLOCK_SIZE 0x80000

uint32_t buffer[BUFFER_BLOCK_SIZE >> 2];

inline void endian_from_little_to_big(uint64_t* a) {
    for (int i = 0; i < 16; i++) {
        std::swap(a[i], a[31 - i]);
        a[i] = __builtin_bswap64(a[i]);
        a[31 - i] = __builtin_bswap64(a[31 - i]);
    }
}

uint32_t count_leading_zero(uint64_t* a) {
    uint32_t ans = 0;

    for (int i = 31; i >= 0; i--) {
        if (a[i] == 0) {
            ans += 64;
        }
        else {
            ans += __builtin_clzll(a[i]);
            return ans;
        }
    }

    return ans;
}

uint64_t compute_neg_m_inverse(uint64_t* m) {
    uint64_t cur = m[0];
    uint64_t ans = 1;

    for (uint64_t i = 1; i < 64; i++) {
        uint64_t tmp = m[0] << i;
        if ((1ULL << i) & cur) {
            cur += tmp;
            ans += (1ULL << i);
        }
    }

    return -ans;
}

constexpr uint8_t L[32] = {
    0xE3, 0xB0, 0xC4, 0x42, 0x98, 0xFC, 0x1C, 0x14, 0x9A, 0xFB, 0xF4, 0xC8, 0x99, 0x6F, 0xB9, 0x24, 
    0x27, 0xAE, 0x41, 0xE4, 0x64, 0x9B, 0x93, 0x4C, 0xA4, 0x95, 0x99, 0x1B, 0x78, 0x52, 0xB8, 0x55
};

static const union {
    uint32_t dw[64];
    __m128i x[16];
} K =
{
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
};

const __m128i byteswapindex = _mm_set_epi8(12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3);

inline void mult_by_2(uint64_t* a) {
    /* this function aimed to implement the compute of R^2 mod m in montgomery modulo */
    bool ok = false;
    bool tmp = false;
    for (int i = 0; i < 32; i++) {
        if (a[i] & (1ULL << 63)) {
            tmp = true;
        }
        else {
            tmp = false;
        }

        a[i] = (a[i] << 1) | (ok ? 1 : 0);

        ok = tmp;
    }
    a[32] |= ok ? 1 : 0;
}

bool cmp_2049bits(const uint64_t* a, const uint64_t* b) {
    for (int i = 32; i >= 0; i--) {
        if (a[i] != b[i]) {
            return a[i] >= b[i] ? true : false;
        }
    }

    return true;
}

inline void trivial_multiple_precision_subtraction_v2(uint64_t* a, uint64_t* b, uint64_t* ans) {
    /* a should be bigger than b */
    uint64_t c = 0;

    for (int i = 0; i < 33; i++) {
        if (a[i] == 0 && b[i] == 0 && c == 0) continue;

        uint64_t tmp = a[i] - c - b[i];

        if ((c == 0 && tmp > a[i]) || (c == 1 && tmp >= a[i])) {
            c = 1;
        }
        else {
            c = 0;
        }

        ans[i] = tmp;
    }

    ans[32] &= 1ULL;
}

inline void montgomery_init_module(uint64_t* a, uint64_t* m) {
    /* warning: this function compute R^2 mod m, in order to implement the montgomery multiplication */
    for (int i = 0; i < 4096; i++) {
        mult_by_2(a);

        if (cmp_2049bits(a, m)) {
            trivial_multiple_precision_subtraction_v2(a, m, a);
        }
    }
    return;
}

inline void montgomery_init_module_v2(uint64_t* a, uint64_t* m) {
    /* warning: this function compute R mod m, in order to implement the montgomery exponentiation */
    for (int i = 0; i < 2048; i++) {
        mult_by_2(a);

        if (cmp_2049bits(a, m)) {
            trivial_multiple_precision_subtraction_v2(a, m, a);
        }
    }
}

bool cmp_2048bits(const uint64_t* a, const uint64_t* b) {
    /* true when a >= b */
    for (int i = 31; i >= 0; i--) {
        if (a[i] != b[i]) {
            return a[i] >= b[i] ? true : false;
        }
    }

    return true;
}

inline void trivial_multiple_precision_subtraction(uint64_t* a, uint64_t* b, uint64_t* ans) {
    /* a should be bigger than b */
    uint64_t c = 0;

    for (int i = 0; i < 32; i++) {
        if (a[i] == 0 && b[i] == 0 && c == 0) {
            ans[i] = 0ULL;
            continue;
        }

        uint64_t tmp = a[i] - c - b[i];

        if ((c == 0 && tmp > a[i]) || (c == 1 && tmp >= a[i])) {
            c = 1;
        }
        else {
            c = 0;
        }

        ans[i] = tmp;
    }
}

inline void montgomery_multiplication_v2(uint64_t* m, const uint64_t* x, const uint64_t* y, uint64_t* ans, uint64_t neg_m_inverse) {
    uint64_t T[34] = {0}; 

    for (int i = 0; i < 32; i++) {
        __uint128_t carry = 0;

        for (int j = 0; j < 32; j++) {
            __uint128_t tmp = (__uint128_t) x[i] * (__uint128_t) y[j] + (__uint128_t) T[j] + carry;
            T[j] = (uint64_t) tmp;
            carry = tmp >> 64;
        }
        __uint128_t tmp2 = (__uint128_t) T[32] + carry;
        T[32] = (uint64_t) tmp2;
        T[33] = tmp2 >> 64;

        uint64_t u = T[0] * neg_m_inverse;

        __uint128_t tmp3 = (__uint128_t) u * (__uint128_t) m[0] + (__uint128_t) T[0];
        carry = tmp3 >> 64;
        
        for (int j = 1; j < 32; j++) {
            __uint128_t tmp4 = (__uint128_t) u * (__uint128_t) m[j] + (__uint128_t) T[j] + (__uint128_t) carry;
            T[j - 1] = (uint64_t) tmp4;
            carry = tmp4 >> 64;
        }

        __uint128_t tmp5 = (__uint128_t) T[32] + carry;
        T[31] = (uint64_t)tmp5;
        carry = tmp5 >> 64;
        
        tmp5 = (__uint128_t) T[33] + carry;
        T[32] = (uint64_t) tmp5;
        T[33] = 0;
    }

    if (T[32] || cmp_2048bits(T, m)) {
        trivial_multiple_precision_subtraction(T, m, ans);
    }
    else {
        memcpy(ans, T, 256);
    }
}

__attribute__((target("sse2,ssse3,sha")))
void sha256(uint32_t* w, __m128i &h0145, __m128i &h2367) {
    const __m128i* msgx = (const __m128i*)w;
    __m128i cw0 = _mm_shuffle_epi8(_mm_loadu_si128(msgx), byteswapindex);
    __m128i cw1 = _mm_shuffle_epi8(_mm_loadu_si128(msgx + 1), byteswapindex);
    __m128i cw2 = _mm_shuffle_epi8(_mm_loadu_si128(msgx + 2), byteswapindex);
    __m128i cw3 = _mm_shuffle_epi8(_mm_loadu_si128(msgx + 3), byteswapindex);

    #define CYCLE_W(CW0, CW1, CW2, CW3) \
        CW0 = _mm_sha256msg1_epu32(CW0, CW1); \
        CW0 = _mm_add_epi32(CW0, _mm_alignr_epi8(CW3, CW2, 4)); \
        CW0 = _mm_sha256msg2_epu32(CW0, CW3);

    __m128i state1 = h0145;
    __m128i state2 = h2367;
    __m128i tmp;

    #define SHA256_ROUNDS_4(cwN, n) \
    tmp = _mm_add_epi32(cwN, K.x[n]); /* w3+K3 : w2+K2 : w1+K1 : w0+K0 */                 \
    state2 = _mm_sha256rnds2_epu32(state2, state1, tmp);/* state2 = a':b':e':f' / state1 = c':d':g':h' */   \
    tmp = _mm_unpackhi_epi64(tmp, tmp);                 /* - : - : w3+K3 : w2+K2 */                         \
    state1 = _mm_sha256rnds2_epu32(state1, state2, tmp);/* state1 = a':b':e':f' / state2 = c':d':g':h' */

    /* w0 - w3 */
    SHA256_ROUNDS_4(cw0, 0);        
    /* w4 - w7 */
    SHA256_ROUNDS_4(cw1, 1);
    /* w8 - w11 */
    SHA256_ROUNDS_4(cw2, 2);
    /* w12 - w15 */
    SHA256_ROUNDS_4(cw3, 3);
    /* w16 - w19 */                                                     
    CYCLE_W(cw0, cw1, cw2, cw3);    /* cw0 = w19 : w18 : w17 : w16 */   
    SHA256_ROUNDS_4(cw0, 4);                                            
    /* w20 - w23 */                                                 
    CYCLE_W(cw1, cw2, cw3, cw0);    /* cw1 = w23 : w22 : w21 : w20 */   
    SHA256_ROUNDS_4(cw1, 5);                                        
    /* w24 - w27 */                                                     
    CYCLE_W(cw2, cw3, cw0, cw1);    /* cw2 = w27 : w26 : w25 : w24 */   
    SHA256_ROUNDS_4(cw2, 6);                                        
    /* w28 - w31 */                                                     
    CYCLE_W(cw3, cw0, cw1, cw2);    /* cw3 = w31 : w30 : w29 : w28 */   
    SHA256_ROUNDS_4(cw3, 7);
    /* w32 - w35 */
    CYCLE_W(cw0, cw1, cw2, cw3);    /* cw0 = w35 : w34 : w33 : w32 */
    SHA256_ROUNDS_4(cw0, 8);
    /* w36 - w39 */
    CYCLE_W(cw1, cw2, cw3, cw0);    /* cw1 = w39 : w38 : w37 : w36 */
    SHA256_ROUNDS_4(cw1, 9);
    /* w40 - w43 */
    CYCLE_W(cw2, cw3, cw0, cw1);    /* cw2 = w43 : w42 : w41 : w40 */
    SHA256_ROUNDS_4(cw2, 10);
    /* w44 - w47 */
    CYCLE_W(cw3, cw0, cw1, cw2);    /* cw3 = w47 : w46 : w45 : w44 */
    SHA256_ROUNDS_4(cw3, 11);
    /* w48 - w51 */
    CYCLE_W(cw0, cw1, cw2, cw3);    /* cw0 = w51 : w50 : w49 : w48 */
    SHA256_ROUNDS_4(cw0, 12);
    /* w52 - w55 */
    CYCLE_W(cw1, cw2, cw3, cw0);    /* cw1 = w55 : w54 : w53 : w52 */
    SHA256_ROUNDS_4(cw1, 13);
    /* w56 - w59 */
    CYCLE_W(cw2, cw3, cw0, cw1);    /* cw2 = w59 : w58 : w57 : w56 */
    SHA256_ROUNDS_4(cw2, 14);
    /* w60 - w63 */
    CYCLE_W(cw3, cw0, cw1, cw2);    /* cw3 = w63 : w62 : w61 : w60 */
    SHA256_ROUNDS_4(cw3, 15);

    h0145 = _mm_add_epi32(state1, h0145);
    h2367 = _mm_add_epi32(state2, h2367);
}

__attribute__((target("sse2,ssse3,sha")))
void MGF_32to223(uint8_t DBMask[223], uint8_t* seed) {
    __m128i h0145 = _mm_set_epi32(0x6a09e667, 0xbb67ae85, 0x510e527f, 0x9b05688c);
    __m128i h2367 = _mm_set_epi32(0x3c6ef372, 0xa54ff53a, 0x1f83d9ab, 0x5be0cd19);
    __m128i h0123, h4567;
    uint32_t w[16] = {0};
    memcpy(w, seed, 32);
    w[9] = 1 << 7;
    w[15] = __builtin_bswap32((uint32_t) (36 << 3));

    for (uint8_t i = 0; i < 6; i++) {
        w[8] = __builtin_bswap32(i);
        sha256(w, h0145, h2367);
        h0123 = _mm_unpackhi_epi64(h2367, h0145);
        h4567 = _mm_unpacklo_epi64(h2367, h0145);
        h0123 = _mm_shuffle_epi8(h0123, _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
        h4567 = _mm_shuffle_epi8(h4567, _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
        memcpy(&DBMask[0 + (i << 5)], &h0123, sizeof(__m128i));
        memcpy(&DBMask[16 + (i << 5)], &h4567, sizeof(__m128i));
        h0145 = _mm_set_epi32(0x6a09e667, 0xbb67ae85, 0x510e527f, 0x9b05688c);
        h2367 = _mm_set_epi32(0x3c6ef372, 0xa54ff53a, 0x1f83d9ab, 0x5be0cd19);
    }

    w[8] = __builtin_bswap32(6);
    sha256(w, h0145, h2367);
    h0123 = _mm_unpackhi_epi64(h2367, h0145);
    h4567 = _mm_unpacklo_epi64(h2367, h0145);
    h0123 = _mm_shuffle_epi8(h0123, _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
    h4567 = _mm_shuffle_epi8(h4567, _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
    memcpy(&DBMask[192], &h0123, sizeof(__m128i));
    memcpy(&DBMask[208], &h4567, 15);
}

__attribute__((target("sse2,ssse3,sha")))
void MGF_223to32(uint8_t SeedMask[32], uint8_t* MaskedDB) {
    __m128i h0145 = _mm_set_epi32(0x6a09e667, 0xbb67ae85, 0x510e527f, 0x9b05688c);
    __m128i h2367 = _mm_set_epi32(0x3c6ef372, 0xa54ff53a, 0x1f83d9ab, 0x5be0cd19);
    __m128i h0123, h4567;
    uint32_t w[16] = {0};

    for (uint8_t i = 0; i < 3; i++) {
        memcpy(w, &MaskedDB[i << 6], 64);
        sha256(w, h0145, h2367);
    }

    memset(w, 0, 64);
    memcpy(w, &MaskedDB[192], 31);
    w[8] |= (1 << 31);
    w[15] = __builtin_bswap32(227 << 3);
    sha256(w, h0145, h2367);
    h0123 = _mm_unpackhi_epi64(h2367, h0145);
    h4567 = _mm_unpacklo_epi64(h2367, h0145);
    h0123 = _mm_shuffle_epi8(h0123, _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
    h4567 = _mm_shuffle_epi8(h4567, _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
    memcpy(SeedMask, &h0123, 16);
    memcpy(&SeedMask[16], &h4567, 16);
}

__attribute__((target("sse2,ssse3,sha")))
int main() {
    #ifdef _WIN32
    setmode(fileno(stdin), O_BINARY);
    setmode(fileno(stdout), O_BINARY);
    #endif
    #ifndef ONLINE_JUDGE
    FILE* in = fopen("52input.bin", "rb");
    FILE* out = fopen("52output.bin", "wb");
    #else
    #define in stdin
    #define out stdout
    #endif

    /* initialize */
    uint64_t p[16], q[16];
    uint64_t n[32], d[32];
    uint64_t d_mod_psubone[16], d_mod_qsubone[16];
    uint64_t q_inv_modp[16];
    uint64_t input[32];

    fread(p, sizeof(uint64_t), 16, in);
    fread(q, sizeof(uint64_t), 16, in);
    fread(n, sizeof(uint64_t), 32, in);
    fread(d, sizeof(uint64_t), 32, in);
    fread(d_mod_psubone, sizeof(uint64_t), 16, in);
    fread(d_mod_qsubone, sizeof(uint64_t), 16, in);
    fread(q_inv_modp, sizeof(uint64_t), 16, in);
    fread(input, sizeof(uint64_t), 32, in);

    endian_from_little_to_big(input);
    endian_from_little_to_big(n);
    endian_from_little_to_big(d);
    uint64_t n_extend[33] = {0};
    memcpy(n_extend, n, 256);

    uint64_t R[33] = {0};
    R[0] = 1ULL;
    uint64_t R2[33] = {0};
    R2[0] = 1ULL;
    uint64_t one[32] = {0};
    one[0] = 1ULL;
    montgomery_init_module_v2(R, n_extend);
    montgomery_init_module(R2, n_extend);
    uint64_t neg_n_inverse = compute_neg_m_inverse(n);

    uint64_t ans_pow[32] = {0};
    ans_pow[0] = 1ULL;
    memcpy(ans_pow, R, 256);

    uint64_t iteration_new[32] = {0};
    uint64_t aR[32] = {0};
    montgomery_multiplication_v2(n, input, R2, aR, neg_n_inverse);
    memcpy(iteration_new, aR, 256);

    for (int i = 0; i < 2048; i++, montgomery_multiplication_v2(n, iteration_new, iteration_new, iteration_new, neg_n_inverse)) {
        if ((d[i >> 6] >> (i % 64)) & 1ULL) {
            montgomery_multiplication_v2(n, ans_pow, iteration_new, ans_pow, neg_n_inverse);
        }
    }

    montgomery_multiplication_v2(n, ans_pow, one, ans_pow, neg_n_inverse);
    endian_from_little_to_big(ans_pow);
    
    /* CRT */

    uint8_t EM[256];
    memcpy(EM, ans_pow, 256);
    if (EM[0] != 0x00) {
        std::cerr << "EM should be started with 0x00" << '\n';
        return 0;
    }

    uint8_t MaskedSeed[32];
    memcpy(MaskedSeed, &EM[1], 32);
    uint8_t MaskedDB[223];
    memcpy(MaskedDB, &EM[33], 223);

    uint8_t SeedMask[32] = {0};
    MGF_223to32(SeedMask, MaskedDB);
    uint8_t seed[32];
    for (int i = 0; i < 32; i++) {
        seed[i] = MaskedSeed[i] ^ SeedMask[i];
    }
    uint8_t DBMask[223] = {0};
    MGF_32to223(DBMask, seed);

    uint8_t DB[223];
    for (int i = 0; i < 223; i++) {
        DB[i] = DBMask[i] ^ MaskedDB[i];
    }

    for (int i = 0; i < 32; i++) {
        if (DB[i] != L[i]) {
            std::cerr << "hashl incorrect" << '\n';
            return 0;
        }
    }

    int ptr = -1;
    for (int i = 32; i < 223; i++) {
        if (DB[i] == 0x00) continue;
        else if (DB[i] == 0x01) {
            if (i < 222) {
                ptr = i + 1;
                break;
            }
            else {
                std::cerr << "plaintext is null" << '\n';
                return 0;
            }
            
        }
        else {
            std::cerr << "no 0x01 expected" << '\n';
            return 0;
        }
    }

    if (ptr != -1) {
        fwrite(&DB[ptr], sizeof(uint8_t), 223 - ptr, out);
    }
    return 0;
}