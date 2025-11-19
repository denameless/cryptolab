#include <bits/stdc++.h>
#ifdef _WIN32
#include <fcntl.h>
#endif

#include <immintrin.h>

#define BUFFER_BLOCK_SIZE 0x80000

uint32_t buffer[BUFFER_BLOCK_SIZE >> 2];

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

/* https://www.officedaytime.com/simd512e/simdimg/sha256.html */
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
int main() {
    #ifdef _WIN32
    setmode(fileno(stdin), O_BINARY);
    setmode(fileno(stdout), O_BINARY);
    #endif
    #ifndef ONLINE_JUDGE
    FILE* in = fopen("dump.bin", "rb");
    FILE* out = fopen("41output.bin", "wb");
    #else
    #define in stdin
    #define out stdout
    #endif

    /* initialize */

    uint32_t H[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    };

    __m128i h0145 = _mm_set_epi32(H[0], H[1], H[4], H[5]);
    __m128i h2367 = _mm_set_epi32(H[2], H[3], H[6], H[7]);

    uint64_t len = 0;
    uint32_t w[64] = {0};
    uint64_t curlen = 0;
    
    while (1) {
        curlen = fread(buffer, sizeof(uint8_t), BUFFER_BLOCK_SIZE, in);

        if (curlen == BUFFER_BLOCK_SIZE) {
            len += BUFFER_BLOCK_SIZE;
            for (uint32_t i = 0; i < BUFFER_BLOCK_SIZE >> 6; i++) {
                // memcpy(w, &buffer[i << 4], 64);
                sha256(&buffer[i << 4], h0145, h2367);
            }

        }
        else {
            break;
        }
    }

    for (uint32_t i = 0; i < curlen >> 6; i++) {
        // memcpy(w, &buffer[i << 4], 64);
        sha256(&buffer[i << 4], h0145, h2367);
    }

    memset(&reinterpret_cast<uint8_t*>(buffer)[curlen], 0, BUFFER_BLOCK_SIZE - curlen);

    // memcpy(w, &buffer[(curlen >> 6) << 4], 64);

    uint32_t base = (curlen >> 6) << 4;

    len += (curlen >> 6) << 6;
    curlen %= 64;

    buffer[base + (curlen >> 2)] |= 1 << (7 + ((curlen % 4) << 3));
    len += curlen;

    if ((len % 64 < 56) && (len % 64 > 0)) {
        buffer[base + 14] = __builtin_bswap32((uint32_t) (len >> 29));
        buffer[base + 15] = __builtin_bswap32((uint32_t) (len << 3) & 0xffffffff);
        sha256(&buffer[base], h0145, h2367);
    }
    else {
        if (len % 64 != 0) {
            sha256(&buffer[base], h0145, h2367);
            memset(&buffer[base], 0, 62);
        }
        else {
            memset(w, 0, 62);
            buffer[base + (len >> 2) % 16] |= 1 << 7;
        }

        buffer[base + 14] = __builtin_bswap32((uint32_t) (len >> 29));
        buffer[base + 15] = __builtin_bswap32((uint32_t) (len << 3) & 0xffffffff);
        sha256(&buffer[base], h0145, h2367);
    }

    __m128i h0123 = _mm_unpackhi_epi64(h2367, h0145);
    __m128i h4567 = _mm_unpacklo_epi64(h2367, h0145);
    h0123 = _mm_shuffle_epi8(h0123, _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
    h4567 = _mm_shuffle_epi8(h4567, _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));

    fwrite(&h0123, sizeof(__m128i), 1, out);
    fwrite(&h4567, sizeof(__m128i), 1, out);
    return 0;
}