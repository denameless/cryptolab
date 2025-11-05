#include <bits/stdc++.h>
#ifdef _WIN32
#include <fcntl.h>
#endif

#include <wmmintrin.h>

using namespace std;

#define key_expanding(i, rcon) \
    tmp = _mm_aeskeygenassist_si128(expandedkey[i - 1], rcon); \
    expandedkey[i] = _mm_xor_si128(expandedkey[i - 1], _mm_shuffle_epi32(tmp, 0xff)); \
    tmp = _mm_slli_si128(expandedkey[i - 1], 0x04); \
    expandedkey[i] = _mm_xor_si128(expandedkey[i], tmp); \
    tmp = _mm_slli_si128(expandedkey[i - 1], 0x08); \
    expandedkey[i] = _mm_xor_si128(expandedkey[i], tmp); \
    tmp = _mm_slli_si128(expandedkey[i - 1], 0x0c); \
    expandedkey[i] = _mm_xor_si128(expandedkey[i], tmp);

constexpr uint8_t rcon[10] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36
};

__attribute__((target("aes,pclmul,sse2")))
void AESencrypt(__m128i &txt, const __m128i* expandedkey) {
    txt = _mm_xor_si128(txt, expandedkey[0]);

    for (uint8_t i = 1; i <= 9; i++) {
        txt = _mm_aesenc_si128(txt, expandedkey[i]);
    }

    txt = _mm_aesenclast_si128(txt, expandedkey[10]);
}

__attribute__((target("aes,pclmul,sse2")))
void AESdecrypt(__m128i &txt, const __m128i* expandedkey) {
    txt = _mm_xor_si128(txt, expandedkey[10]);

    for (uint8_t i = 9; i >= 1; i--) {
        txt = _mm_aesdec_si128(txt, expandedkey[i]);
    }

    txt = _mm_aesdeclast_si128(txt, expandedkey[0]);
}

void CBCencrypt(__m128i &txt, const __m128i* expandedkey, __m128i &y) {
    txt = _mm_xor_si128(txt, y);
    AESencrypt(txt, expandedkey);
    y = txt;
}

void CBCdecrypt(__m128i &txt, const __m128i* expandedkey, __m128i &y) {
    __m128i tmp = txt;
    AESdecrypt(txt, expandedkey);
    txt = _mm_xor_si128(txt, y);
    y = tmp;
}

__attribute__((target("aes,pclmul,sse2")))
int main() {
    #ifdef _WIN32
    setmode(fileno(stdin), O_BINARY);
    setmode(fileno(stdout), O_BINARY);
    #endif
    #ifndef ONLINE_JUDGE
    FILE* in = fopen("31input.bin", "rb");
    FILE* out = fopen("31output.bin", "wb");
    #else
    #define in stdin
    #define out stdout
    #endif

    uint8_t mod;
    uint8_t key[16], IV[16];
    fread(&mod, sizeof(uint8_t), 1, in), fread(key, sizeof(uint8_t), 16, in), fread(IV, sizeof(uint8_t), 16, in);

    uint32_t len;
    fread(&len, sizeof(uint32_t), 1, in);

    __m128i expandedkey[11];
    expandedkey[0] = _mm_loadu_si128(reinterpret_cast<const __m128i *>(key));

    __m128i y, tmp;
    y = _mm_loadu_si128(reinterpret_cast<const __m128i *>(IV));

    key_expanding(1, 0x01);

    key_expanding(2, 0x02);

    key_expanding(3, 0x04);

    key_expanding(4, 0x08);

    key_expanding(5, 0x10);

    key_expanding(6, 0x20);

    key_expanding(7, 0x40);

    key_expanding(8, 0x80);

    key_expanding(9, 0x1b);

    key_expanding(10, 0x36);

    const uint32_t extend = ((len % 16) == 0 ? 16 : 16 - len % 16);
    const uint32_t block = (len >> 4) - 1;
    const uint32_t extendedblock = ((len + extend) >> 4) - 1;

    uint8_t txt[16];
    __m128i block_i128;
    __m128i a;

    switch(mod) {
        case 0x01 :
            for (int i = 0; i < extendedblock; i++) {
                fread(&a, sizeof(__m128i), 1, in);
                CBCencrypt(a, expandedkey, y);
                fwrite(&a, sizeof(__m128i), 1, out);
            }

            fread(txt, sizeof(uint8_t), 16 - extend, in);
            for (int j = 16 - extend; j < 16; j++) {
                txt[j] = extend;
            }

            block_i128 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(txt));
            CBCencrypt(block_i128, expandedkey, y);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(txt), block_i128); 
            fwrite(&txt, sizeof(uint8_t), 16, out);
            break;

        case 0x81 :
            for (uint8_t i = 1; i < 10; i++) {
                expandedkey[i] = _mm_aesimc_si128(expandedkey[i]);
            }
            for (int i = 0; i < block; i++) {
                fread(&a, sizeof(__m128i), 1, in);
                CBCdecrypt(a, expandedkey, y);
                fwrite(&a, sizeof(__m128i), 1, out);
            }

            fread(&a, sizeof(__m128i), 1, in);
            CBCdecrypt(a, expandedkey, y);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(txt), a);

            uint32_t expected_extend = txt[15];
            for (long long i = 15; i > 16 - expected_extend; i--) {
                if (txt[i] != txt[i - 1]) {
                    fwrite(txt, sizeof(uint8_t), 16, out);
                    return 0;
                }
            }

            fwrite(txt, sizeof(uint8_t), 16 - expected_extend, out);
            break;
    }

    return 0;
}