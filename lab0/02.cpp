#include <bits/stdc++.h>

#ifdef _WIN32
#include <fcntl.h>
#endif

using namespace std;

int main() {
    #ifdef _WIN32
    setmode(fileno(stdin), O_BINARY);
    setmode(fileno(stdout), O_BINARY);
    #endif
    #ifndef ONLINE_JUDGE
    FILE* in = fopen("newinput.bin", "rb");
    FILE* out = fopen("newoutput.bin", "wb");
    #else
    #define in stdin
    #define out stdout
    #endif

    ios::sync_with_stdio(false), cin.tie(0), cout.tie(0);

    uint8_t a;
    uint32_t len;
    fread(&a, sizeof(uint8_t), 1, in);
    fread(&len, sizeof(uint32_t), 1, in);
    uint8_t cur;
    
    for (uint32_t i = 0; i < len; i++) {
        fread(&cur, sizeof(uint8_t), 1, in);
        cur ^= a;
        fwrite(&cur, sizeof(uint8_t), 1, out);
    }
    return 0;
}