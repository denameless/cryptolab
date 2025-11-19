#include <bits/stdc++.h>
#ifdef _WIN32
#include <fcntl.h>
#endif

#define BUFFER_BLOCK_SIZE 0x80000

uint32_t buffer[BUFFER_BLOCK_SIZE >> 2];

uint32_t k[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

uint32_t t1 = 0, t2 = 0;

void sha256(uint32_t H[8], uint32_t h[8], uint32_t w[64]) {
    memcpy(h, H, 32);

    for (uint32_t i = 16; i < 64; i++) {
        w[i] = (((w[i - 2] >> 17) | (w[i - 2] << 15)) ^ ((w[i - 2] >> 19) | (w[i - 2] << 13)) ^ (w[i - 2] >> 10)) + w[i - 7] + (((w[i - 15] >> 7) | (w[i - 15] << 25)) ^ ((w[i - 15] >> 18) | (w[i - 15] << 14)) ^ (w[i - 15] >> 3)) + w[i - 16];
    }

    for (uint8_t i = 0; i < 64; i++) {
        t1 = h[7] + (((h[4] >> 6) | (h[4] << 26)) ^ ((h[4] >> 11) | (h[4] << 21)) ^ ((h[4] >> 25) | (h[4] << 7))) + ((h[4] & h[5]) ^ (~h[4] & h[6])) + k[i] + w[i];
        t2 = (((h[0] >> 2) | (h[0] << 30)) ^ ((h[0] >> 13) | (h[0] << 19)) ^ ((h[0] >> 22) | (h[0] << 10))) + ((h[0] & h[1]) ^ (h[0] & h[2]) ^ (h[1] & h[2]));
        h[7] = h[6];
        h[6] = h[5];
        h[5] = h[4];
        h[4] = h[3] + t1;
        h[3] = h[2];
        h[2] = h[1];
        h[1] = h[0];
        h[0] = t1 + t2;
    }

    for (uint8_t i = 0; i < 8; i++) {
        H[i] += h[i];
    }
}

int main() {
    #ifdef _WIN32
    setmode(fileno(stdin), O_BINARY);
    setmode(fileno(stdout), O_BINARY);
    #endif
    #ifndef ONLINE_JUDGE
    FILE* in = fopen("41standard.bin", "rb");
    FILE* out = fopen("41output.bin", "wb");
    #else
    #define in stdin
    #define out stdout
    #endif

    uint32_t H[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    };
    uint32_t h[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    };

    uint64_t len = 0;
    uint32_t w[64] = {0};
    uint64_t curlen = 0;

    while (1) {
        curlen = fread(buffer, sizeof(uint8_t), BUFFER_BLOCK_SIZE, in);

        if (curlen == BUFFER_BLOCK_SIZE) {
            len += BUFFER_BLOCK_SIZE;
            for (uint32_t i = 0; i < BUFFER_BLOCK_SIZE >> 6; i++) {
                for (uint32_t j = 0; j < 16; j++) {
                    buffer[(i << 4) + j] = __builtin_bswap32(buffer[(i << 4) + j]);
                }

                memcpy(w, &buffer[i << 4], 64);
                sha256(H, h, w);
            }

        }
        else {
            break;
        }
    }

    for (uint32_t i = 0; i < curlen >> 6; i++) {
        for (uint32_t j = 0; j < 16; j++) {
            buffer[(i << 4) + j] = __builtin_bswap32(buffer[(i << 4) + j]);
        }

        memcpy(w, &buffer[i << 4], 64);
        sha256(H, h, w);
    }

    memset(&reinterpret_cast<uint8_t*>(buffer)[curlen], 0, BUFFER_BLOCK_SIZE - curlen);

    for (uint32_t i = 0; i < 16; i++) {
        buffer[((curlen >> 6) << 4) + i] = __builtin_bswap32(buffer[((curlen >> 6) << 4) + i]);
    }

    memcpy(w, &buffer[(curlen >> 6) << 4], 64);

    len += (curlen >> 6) << 6;
    curlen %= 64;

    w[curlen >> 2] |= 1 << (31 - ((curlen % 4) << 3));
    len += curlen;

    if ((len % 64 < 56) && (len % 64 > 0)) {
        w[14] = (uint32_t) (len >> 29);
        w[15] = (uint32_t) (len << 3) & 0xffffffff;
        sha256(H, h, w);
    }
    else {
        if (len % 64 != 0) {
            sha256(H, h, w);
            memset(w, 0, 62);
        }
        else {
            memset(w, 0, 62);
            w[(len >> 2) % 16] |= 1 << 31;
        }

        w[14] = (uint32_t) (len >> 29);
        w[15] = (uint32_t) (len << 3) & 0xffffffff;
        sha256(H, h, w);
    }

    for (int i = 0; i < 8; i++) {
        uint32_t tmp = __builtin_bswap32(H[i]);
        fwrite(&tmp, sizeof(uint32_t), 1, out);
    }

    return 0;
}