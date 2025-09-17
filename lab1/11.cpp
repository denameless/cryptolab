#include <bits/stdc++.h>
#ifdef _WIN32
#include <fcntl.h>
#endif

using namespace std;

long double expectedvalue = 0.065601;

vector<long double> refer = {
    0.082, 0.015, 0.028, 0.043, 0.127, 0.022, 0.020, 0.061, 0.070, 0.002, 0.008, 0.040,
    0.024, 0.067, 0.075, 0.019, 0.001, 0.060, 0.063, 0.091, 0.028, 0.010, 0.023, 0.001,
    0.020, 0.001
};

long double analyze(string a) {
    int len = a.length();
    vector<int> hashtable(26);
    for (auto ch : a) {
        hashtable[(int)(ch - 'a')]++;

    }
    long double ans = 0;
    for (int cur : hashtable) {
        ans += cur * (cur - 1);
    }
    return ans / (len * (len - 1));
}

long double newanalyze(string a) {
    long double len = a.length();
    vector<int> hashtable(26);
    for (auto ch : a) {
        hashtable[(int)(ch - 'a')]++;

    }
    long double ans = 0;
    for (int i = 0; i < 26; i++) {
        ans += abs(((long double)hashtable[i] / len) - refer[i]);
    }
    return ans;
}

string shiftstring(string a, int shift) {
    string ans = a;

    for (char &ch : ans) {
        int offset = ch - 'a';
        int newoffset = (offset - shift + 26) % 26;
        ch = 'a' + newoffset;
    }

    return ans;
}
int main() {
    string secret;
    string alphasecret;
    #ifdef ONLINE_JUDGE
    string line;
    while (getline(cin, line)) {
        secret += line;
        secret += '\n';
    }
    #else
    ifstream infile("11input.txt");
    string line;
    if (infile.is_open()) {
        while(getline(infile, line)) {
            secret += line;
            secret += '\n';
        }
    }
    infile.close();
    #endif

    for (auto ch : secret) {
        if (isalpha(ch)) {
            alphasecret += tolower(ch);
        }
    }

    int len = secret.length() / 80 + 1;
    int keylen = -1;

    long double minvalue = INT_MAX;
    for (int i = 1; i <= len; i++) {
        int ptr = 0;
        vector<string> str(i);
        for (auto ch : alphasecret) {
            str[ptr] += ch;
            ptr = (ptr + 1) % i; 
        }
        long double curanalysis = 0;
        for (string cur : str) {
            curanalysis += abs(analyze(cur) - expectedvalue);
        }
        curanalysis /= i;
        if (curanalysis < minvalue) {
            minvalue = curanalysis;
            keylen = i;
        }
    }

    vector<string> str(keylen);
    int ptr = 0;
    for (auto ch : alphasecret) {
        str[ptr] += ch;
        ptr = (ptr + 1) % keylen; 
    }

    string crypto;
    for (int i = 0; i < keylen; i++) {
        int ansshift = -1;
        long double curmin = INT_MAX;

        for (int j = 0; j < 26; j++) {
            string temp = str[i];
            long double curanalysis = newanalyze(shiftstring(temp, j));
            if (curanalysis < curmin) {
                curmin = curanalysis;
                ansshift = j;
            }
        }
        crypto += (char)('A' + ansshift);
    }

    cout << crypto << '\n';

    int newptr = 0;

    for (char &ch : secret) {
        if (!isalpha(ch)) {
            cout << ch;
        }
        else {
            if (islower(ch)) {
                ch = 'a' + (26 - (crypto[newptr] - 'A') + (ch - 'a')) % 26;
            }
            else {
                ch = 'A' + (26 - (crypto[newptr] - 'A') + (ch - 'A')) % 26;
            }
            cout << ch;
            newptr = (newptr + 1) % keylen;
        }
    }
    return 0;
}