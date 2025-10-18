#include <bits/stdc++.h>
#ifdef _WIN32
#include <fcntl.h>
#endif

#include <immintrin.h>
using namespace std;

void gf_add(vector<uint64_t> &ans, const vector<uint64_t> &a, const vector<uint64_t> &b) {

    for (uint32_t i = 0; i < 3; i++) {
        ans[i] = a[i] ^ b[i];
    }

}

vector<uint64_t> gf_mod(vector<uint64_t> &a) {
    vector<uint64_t> ans(3);
    uint64_t t;
    for (uint8_t i = 4; i >= 3; i--) {
        t = a[i];
        a[i - 3] ^= (t << 61) ^ (t << 62) ^ (t << 63);
        a[i - 2] ^= (t << 10) ^ (t >> 1) ^ (t >> 2) ^ (t >> 3);
        a[i - 1] ^= (t >> 54);
    }

    t = a[2] & 0xfffffffffffffff8;

    ans[0] = a[0] ^ (t << 10) ^ (t >> 1) ^ (t >> 2) ^ (t >> 3);
    ans[1] = a[1] ^ (t >> 54);
    ans[2] = a[2] & 0x7;

    return ans;
}

__attribute__((target("pclmul,sse2")))
void gf_mul(vector<uint64_t> &ans, const vector<uint64_t> &a, const vector<uint64_t> &b) {
    vector<__m128i> t(9);
    vector<__m128i> m_a(3);
    vector<__m128i> m_b(3);
    for (uint8_t i = 0; i < 3; i++) {
        m_a[i] = _mm_set_epi64x(0, a[i]);
        m_b[i] = _mm_set_epi64x(0, b[i]);
    }

    t[0] = _mm_clmulepi64_si128(_mm_xor_si128(m_a[1], m_a[2]), _mm_xor_si128(m_b[1], m_b[2]), 0x00);
    t[1] = _mm_clmulepi64_si128(_mm_xor_si128(m_a[0], m_a[2]), _mm_xor_si128(m_b[0], m_b[2]), 0x00);
    t[2] = _mm_clmulepi64_si128(_mm_xor_si128(m_a[0], m_a[1]), _mm_xor_si128(m_b[0], m_b[1]), 0x00);
    t[3] = _mm_clmulepi64_si128(m_a[0], m_b[0], 0x00);
    t[4] = _mm_clmulepi64_si128(m_a[1], m_b[1], 0x00);
    t[5] = _mm_clmulepi64_si128(m_a[2], m_b[2], 0x00);
    t[6] = _mm_xor_si128(_mm_xor_si128(t[0], t[4]), t[5]);
    t[7] = _mm_xor_si128(_mm_xor_si128(t[2], t[4]), t[3]);
    t[8] = _mm_xor_si128(_mm_xor_si128(t[1], t[3]), _mm_xor_si128(t[4], t[5]));

    uint64_t t3_low = _mm_cvtsi128_si64(t[3]);
    uint64_t t3_high = _mm_cvtsi128_si64(_mm_srli_si128(t[3], 8)); 
    uint64_t t5_low = _mm_cvtsi128_si64(t[5]);
    uint64_t t6_low = _mm_cvtsi128_si64(t[6]);
    uint64_t t6_high = _mm_cvtsi128_si64(_mm_srli_si128(t[6], 8));
    uint64_t t7_low = _mm_cvtsi128_si64(t[7]);
    uint64_t t7_high = _mm_cvtsi128_si64(_mm_srli_si128(t[7], 8));
    uint64_t t8_low = _mm_cvtsi128_si64(t[8]);
    uint64_t t8_high = _mm_cvtsi128_si64(_mm_srli_si128(t[8], 8));

    vector<uint64_t> mod_input(5);
    mod_input[0] = t3_low;
    mod_input[1] = t7_low ^ t3_high;
    mod_input[2] = t8_low ^ t7_high;
    mod_input[3] = t6_low ^ t8_high;
    mod_input[4] = t5_low ^ t6_high;

    ans = gf_mod(mod_input);
}

// void gf_pow2(vector<uint64_t> &ans, const vector<uint64_t> &a) {
//     vector<uint64_t> mod_input(5);
//     vector<uint64_t> tmp = a;
//     uint64_t pattern = 1;

//     for (uint16_t i = 0; i < 32; i++) {
//         if (tmp[0] & 1) {
//             mod_input[0] += pattern;
//         } 
//         tmp[0] >>= 1;
//         pattern <<= 2;
//     }
    
//     pattern = 1;

//     for (uint16_t i = 0; i < 32; i++) {
//         if (tmp[0] & 1) {
//             mod_input[1] += pattern;
//         }
//         tmp[0] >>= 1;
//         pattern <<= 2;
//     }

//     pattern = 1;

//     for (uint16_t i = 0; i < 32; i++) {
//         if (tmp[1] & 1) {
//             mod_input[2] += pattern;
//         }
//         tmp[1] >>= 1;
//         pattern <<= 2;
//     }

//     pattern = 1;

//     for (uint16_t i = 0; i < 32; i++) {
//         if (tmp[1] & 1) {
//             mod_input[3] += pattern;
//         }
//         tmp[1] >>= 1;
//         pattern <<= 2;
//     }

//     for (uint16_t i = 0; i < 3; i++) {
//         if (tmp[2] & 1) {
//             mod_input[4] += ((uint64_t) 1 << (i * 2));
//         }
//         tmp[2] >>= 1;
//     }

//     ans = gf_mod(mod_input);
// }

void gf_inv(vector<uint64_t> &ans, const vector<uint64_t> &y) {
    vector<uint16_t> n = {1, 2, 4, 8, 16, 32, 65, 130};
    vector<vector<uint64_t>> x(8, vector<uint64_t>(3, 1));
    x[0] = y;
    uint16_t a = 0x82;

    for (uint8_t i = 0; i < 7; i++) {
        a <<= 1;
        vector<uint64_t> tmp = x[i];

        if (i == 5) {
            gf_mul(tmp, x[i], x[i]);
            x[i + 1] = tmp;

            for (uint16_t j = 0; j < n[i]; j++) {
                vector<uint64_t> tmp2 = tmp;
                gf_mul(tmp2, tmp, tmp);
                tmp = tmp2;
            }

            vector<uint64_t> mul_result(3);
            gf_mul(mul_result, x[i + 1], tmp);
            gf_mul(x[i + 1], mul_result, y);
        }
        else {
            x[i + 1] = x[i];

            for (uint16_t j = 0; j < n[i]; j++) {
                vector<uint64_t> tmp2 = tmp;
                gf_mul(tmp2, tmp, tmp);
                tmp = tmp2;
            }

            vector<uint64_t> mul_result(3);
            gf_mul(mul_result, x[i + 1], tmp);
            x[i + 1] = mul_result;
        }
    }

    gf_mul(ans, x[7], x[7]);
}

int main() {
    #ifdef _WIN32
    setmode(fileno(stdin), O_BINARY);
    setmode(fileno(stdout), O_BINARY);
    #endif
    #ifndef ONLINE_JUDGE
    FILE* in = fopen("21input.bin", "rb");
    FILE* out = fopen("21output.bin", "wb");
    #else
    #define in stdin
    #define out stdout
    #endif

    uint32_t operatetimes;
    uint8_t type;
    fread(&operatetimes, sizeof(uint32_t), 1, in);

    for (uint32_t i = 0; i < operatetimes; i++) {
        fread(&type, sizeof(uint8_t), 1, in);
        vector<vector<uint64_t>> num(2, vector<uint64_t>(3));
        fread(num[0].data(), sizeof(uint64_t), 3, in);
        fread(num[1].data(), sizeof(uint64_t), 3, in);
        vector<uint64_t> ans(3);

        switch(type) {
            case 0 :
                gf_add(ans, num[0], num[1]);
                break;
            
            case 1 :
                gf_mul(ans, num[0], num[1]);
                break;

            case 2 : 
                gf_mul(ans, num[0], num[0]);
                // gf_pow2(ans, num[0]);
                break;
            
            case 3 : 
                gf_inv(ans, num[0]);
                break;
        }
        
        fwrite(ans.data(), sizeof(uint64_t), 3, out);
    }

    return 0;
}