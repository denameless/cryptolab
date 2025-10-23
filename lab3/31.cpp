#include <bits/stdc++.h>
#ifdef _WIN32
#include <fcntl.h>
#endif

using namespace std;

uint8_t txt[0x1000];

const uint8_t Sbox[256] = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
};

const uint8_t Sbox_inv[256] = {
    0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
    0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
    0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
    0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
    0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
    0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
    0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
    0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
    0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
    0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
    0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
    0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
    0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
    0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
    0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
    0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d
};

void addroundkey(uint8_t txt[16], uint8_t key[16]) {
    for (uint8_t i = 0; i < 16; i++) {
        txt[i] ^= key[i];
    }
}

void subbytesin128(uint8_t txt[16]) {
    for (int i = 0; i < 16; i++) {
        txt[i] = Sbox[txt[i]];
    }
}

void subbytesin128_inv(uint8_t txt[16]) {
    for (int i = 0; i < 16; i++) {
        txt[i] = Sbox_inv[txt[i]];
    }
}

void shiftrows(uint8_t txt[16]) {
    for (uint8_t i = 0; i < 4; i++) {
        for (uint8_t j = 0; j < i; j++) {
            for (uint8_t k = 0; k < 3; k++) {
                swap(txt[k * 4 + i], txt[(k + 1) * 4 + i]);
            }
        }
    }
}

void shiftrows_inv(uint8_t txt[16]) {
    for (uint8_t i = 0; i < 4; i++) {
        for (uint8_t j = 0; j < i; j++) {
            for (int k = 2; k >= 0; k--) {
                swap(txt[k * 4 + i], txt[(k + 1) * 4 + i]);
            }
        }
    }
}

void fieldmol(uint16_t &x) {
    uint16_t ruler = x >> 8;
    x &= 0xff;
    x ^= ruler;
    x ^= (ruler << 1);
    x ^= (ruler << 3);
    x ^= (ruler << 4);
    x ^= (ruler << 8);
    x &= 0xff;
}

uint8_t fieldmult(uint8_t x, uint8_t y) {
    uint16_t ans = 0;
    
    for (uint8_t i = 0; i < 8; i++) {
        for (uint8_t j = 0; j < 8; j++) {
            ans ^= (x >> i) & (y >> j) & 0x01 ? (1 << (i + j)) : 0; 
        }
    }

    fieldmol(ans);
    return ans;
} 

uint8_t fieldmultby2(uint8_t x) {
    if (x < 128) {
        return x << 1;
    }
    else {
        return (x << 1) ^ 27;
    }
}

uint8_t fieldmultby3(uint8_t x) {
    if (x < 128) {
        return x ^ (x << 1);
    }
    else {
        return x ^ (x << 1) ^ 27;
    }
}

uint8_t fieldmultby14(uint8_t x) {
    /*0x0e*/
    uint16_t ans = 0;
    ans ^= (x << 1) ^ (x << 2) ^ (x << 3);

    if ((x >> 5) & 0x01) {
        ans ^= 27;
    }
    if ((x >> 6) & 0x01) {
        ans ^= 27 ^ 54;
    }
    if ((x >> 7) & 0x01) {
        ans ^= 27 ^ 54 ^ 108;
    }

    return ans;
}

uint8_t fieldmultby11(uint8_t x) {
    /*0x0b*/
    uint16_t ans = x;
    ans ^= (x << 1) ^ (x << 3);

    if ((x >> 5) & 0x01) {
        ans ^= 27;
    }
    if ((x >> 6) & 0x01) {
        ans ^= 54;
    }
    if ((x >> 7) & 0x01) {
        ans ^= 27 ^ 108;
    }

    return ans;
}

uint8_t fieldmultby13(uint8_t x) {
    /*0x0d*/
    uint16_t ans = x;
    ans ^= (x << 2) ^ (x << 3);

    if ((x >> 5) & 0x01) {
        ans ^= 27;
    }
    if ((x >> 6) & 0x01) {
        ans ^= 27 ^ 54;
    }
    if ((x >> 7) & 0x01) {
        ans ^= 54 ^ 108;
    }

    return ans;
}

uint8_t fieldmultby9(uint8_t x) {
    /*0x09*/
    uint16_t ans = x;
    ans ^= (x << 3);

    if ((x >> 5) & 0x01) {
        ans ^= 27;
    }
    if ((x >> 6) & 0x01) {
        ans ^= 54;
    }
    if ((x >> 7) & 0x01) {
        ans ^= 108;
    }

    return ans;
}

void mixcolumns(uint8_t txt[16]) {
    uint8_t x = 0x02;

    for (uint8_t i = 0; i < 4; i++) {
        uint8_t u[] = {0, 0, 0, 0};
        u[0] = fieldmultby2(txt[i * 4]) ^ fieldmultby3(txt[i * 4 + 1]) ^ txt[i * 4 + 2] ^ txt[i * 4 + 3];
        u[1] = fieldmultby2(txt[i * 4 + 1]) ^ fieldmultby3(txt[i * 4 + 2]) ^ txt[i * 4 + 3] ^ txt[i * 4];
        u[2] = fieldmultby2(txt[i * 4 + 2]) ^ fieldmultby3(txt[i * 4 + 3]) ^ txt[i * 4] ^ txt[i * 4 + 1];
        u[3] = fieldmultby2(txt[i * 4 + 3]) ^ fieldmultby3(txt[i * 4]) ^ txt[i * 4 + 1] ^ txt[i * 4 + 2];

        txt[i * 4] = u[0];
        txt[i * 4 + 1] = u[1];
        txt[i * 4 + 2] = u[2];
        txt[i * 4 + 3] = u[3];
    }
} 

void mixcolumns_inv(uint8_t txt[16]) {

    for (uint8_t i = 0; i < 4; i++) {
        uint8_t u[] = {0, 0, 0, 0};
        u[0] = fieldmultby14(txt[i * 4]) ^ fieldmultby11(txt[i * 4 + 1]) ^ fieldmultby13(txt[i * 4 + 2]) ^ fieldmultby9(txt[i * 4 + 3]);
        u[1] = fieldmultby14(txt[i * 4 + 1]) ^ fieldmultby11(txt[i * 4 + 2]) ^ fieldmultby13(txt[i * 4 + 3]) ^ fieldmultby9(txt[i * 4]);
        u[2] = fieldmultby14(txt[i * 4 + 2]) ^ fieldmultby11(txt[i * 4 + 3]) ^ fieldmultby13(txt[i * 4]) ^ fieldmultby9(txt[i * 4 + 1]);
        u[3] = fieldmultby14(txt[i * 4 + 3]) ^ fieldmultby11(txt[i * 4]) ^ fieldmultby13(txt[i * 4 + 1]) ^ fieldmultby9(txt[i * 4 + 2]);

        txt[i * 4] = u[0];
        txt[i * 4 + 1] = u[1];
        txt[i * 4 + 2] = u[2];
        txt[i * 4 + 3] = u[3];
    }
} 

void rotword(uint8_t* txt) {
    uint8_t tmp = txt[0];
    txt[0] = txt[1];
    txt[1] = txt[2];
    txt[2] = txt[3];
    txt[3] = tmp;
}

void subword(uint8_t* txt) {
    txt[0] = Sbox[txt[0]];
    txt[1] = Sbox[txt[1]];
    txt[2] = Sbox[txt[2]];
    txt[3] = Sbox[txt[3]];
}

void keyexpansion(uint8_t key[16], uint8_t* w) {
    uint8_t rcon[11][4] = {
        {0, 0, 0, 0}, 
        {0x01, 0x00, 0x00, 0x00},
        {0x02, 0x00, 0x00, 0x00},
        {0x04, 0x00, 0x00, 0x00}, 
        {0x08, 0x00, 0x00, 0x00}, 
        {0x10, 0x00, 0x00, 0x00}, 
        {0x20, 0x00, 0x00, 0x00}, 
        {0x40, 0x00, 0x00, 0x00}, 
        {0x80, 0x00, 0x00, 0x00}, 
        {0x1b, 0x00, 0x00, 0x00}, 
        {0x36, 0x00, 0x00, 0x00}
    };

    for (uint8_t i = 0; i < 4; i++) {
        uint8_t tmp[] = {key[4 * i], key[4 * i + 1], key[4 * i + 2], key[4 * i + 3]};
        memcpy(&w[i * 4], tmp, 4);
    }

    for (uint8_t i = 4; i < 44; i++) {
        uint8_t tmp[4];
        memcpy(&tmp, &w[i * 4 - 4], 4);

        if (i % 4 == 0) {
            rotword(tmp);
            subword(tmp);
            for (uint8_t j = 0; j < 4; j++) {
                tmp[j] ^= rcon[i / 4][j];
            }
        } 

        for (uint8_t j = 0; j < 4; j++) {
            w[i * 4 + j] = w[(i - 4) * 4 + j] ^ tmp[j];
        }

    } 

}

void AESencrypt(uint8_t txt[16], uint8_t key[16]) {
    uint8_t expandedkey[176];
    keyexpansion(key, expandedkey);
    addroundkey(txt, &expandedkey[0]);

    for (uint8_t i = 1; i <= 9; i++) {
        subbytesin128(txt);
        shiftrows(txt);
        mixcolumns(txt);
        addroundkey(txt, &expandedkey[i * 16]);
    }

    subbytesin128(txt);
    shiftrows(txt);
    addroundkey(txt, &expandedkey[160]);
}

void AESdecrypt(uint8_t txt[16], uint8_t key[16]) {
    uint8_t expandedkey[176];
    keyexpansion(key, expandedkey);

    addroundkey(txt, &expandedkey[160]);
    shiftrows_inv(txt);
    subbytesin128_inv(txt);

    for (int i = 9; i >= 1; i--) {
        addroundkey(txt, &expandedkey[i * 16]);
        mixcolumns_inv(txt);
        shiftrows_inv(txt);
        subbytesin128_inv(txt);
    }

    addroundkey(txt, &expandedkey[0]);
}

void encrypt(uint8_t* txt, uint8_t* key, uint8_t* IV, uint32_t block, uint8_t y[16]) {
    for (long long i = 0; i < block; i++) {
        for (uint8_t j = 0; j < 16; j++) {
            txt[i * 16 + j] ^= y[j];
        }
        AESencrypt(&txt[i * 16], key);
        memcpy(y, &txt[i * 16], 16);
    }
}

void decrypt(uint8_t* txt, uint8_t* key, uint8_t* IV, uint32_t block, uint8_t y[16]) {
    uint8_t tmp[16];

    for (long long i = 0; i < block; i++) {
        memcpy(tmp, &txt[i * 16], 16);
        AESdecrypt(&txt[i * 16], key);
        for (uint8_t j = 0; j < 16; j++) {
            txt[i * 16 + j] ^= y[j];
        }

        memcpy(y, tmp, 16);
    }
}

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

    uint8_t mod, key[16], IV[16];
    fread(&mod, sizeof(uint8_t), 1, in), fread(&key, 1, 16, in), fread(&IV, sizeof(uint8_t), 16, in);

    uint32_t len;
    fread(&len, sizeof(uint32_t), 1, in);
    uint32_t extend = (len % 16) == 0 ? 16 : 16 - len % 16;

    uint8_t y[16];
    memcpy(y, IV, 16);

    switch(mod) {
        case 0x01 :
            for (long long i = 0; i < (len + extend) / 16; i++) {
                if (i != (len + extend) / 16 - 1) {
                    fread(txt, sizeof(uint8_t), 16, in);
                }

                else {
                    fread(txt, sizeof(uint8_t), 16 - extend, in);
                    for (int j = 16 - extend; j < 16; j++) {
                        txt[j] = extend;
                    }
                }
                encrypt(txt, key, IV, 1, y);
                fwrite(txt, sizeof(uint8_t), 16, out);
            }
            break;

        case 0x81 :
            for (long long i = 0; i < len / 16; i++) {
                fread(txt, sizeof(uint8_t), 16, in);
                decrypt(txt, key, IV, 1, y);
                if (i != len / 16 - 1) {
                    fwrite(txt, sizeof(uint8_t), 16, out);
                }
                else {
                    uint32_t expected_extended = txt[15];
                    for (long long i = 15; i > 16 - expected_extended; i--) {
                        if (txt[i] != txt[i - 1]) {
                            fwrite(txt, sizeof(uint8_t), 16, out);
                            return 0;
                        }
                    }

                    fwrite(txt, sizeof(uint8_t), 16 - expected_extended, out);
                }
            }
            break;
    }

    return 0;
}