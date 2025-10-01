#include <bits/stdc++.h>
#ifdef _WIN32
#include <fcntl.h>
#endif

using namespace std;

vector<int> mulinv = {
    -1, 1, -1, 9,
    -1, 21, -1, 15,
    -1, 3, -1, 19,
    -1, -1, -1, 7,
    -1, 23, -1, 11,
    -1, 5, -1, 17,
    -1, 25
};

int det(vector<vector<int>> a) {
    int ans = 0;
    int n = a.size();
    if (n == 1) return a[0][0] % 26;
    
    for (int i = 0; i < n; i++) {
        vector<vector<int>> b(n - 1, vector<int>(n - 1));
        int tmp = 0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            for (int k = 1; k < n; k++) {
                b[tmp][k - 1] = a[j][k];
            }
            tmp++;
        }
        ans = i % 2 == 0 ? (ans + (a[i][0] * det(b) % 26)) % 26 : (ans + (26 - det(b)) * a[i][0] % 26) % 26;
    }

    return ans % 26;
}

vector<vector<int>> add(vector<vector<int>> a, vector<vector<int>> b) {
    int m = a.size(), n = a[0].size();
    vector<vector<int>> ans(m, vector<int>(n));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ans[i][j] = (a[i][j] + b[i][j]) % 26;
        }
    }

    return ans;
}

vector<vector<int>> sub(vector<vector<int>> a, vector<vector<int>> b) {
    int m = a.size(), n = a[0].size();
    vector<vector<int>> ans(m, vector<int>(n));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ans[i][j] = (a[i][j] + (26 - b[i][j])) % 26;
        }
    }

    return ans;
}

vector<int> subvec(vector<int> a, vector<int> b) {
    int n = a.size();
    vector<int> ans(n);

    for (int i = 0; i < n; i++) {
        ans[i] = (a[i] + (26 - b[i])) % 26;
    }

    return ans;
}

vector<vector<int>> multiply(vector<vector<int>> a, vector<vector<int>> b) {
    int m = a.size(), n = b[0].size(), p = a[0].size();
    vector<vector<int>> ans(m, vector<int>(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < p; k++) {
                ans[i][j] += a[i][k] * b[k][j]; 
            }
            ans[i][j] %= 26;
        }
    }

    return ans;
}

vector<vector<int>> mul(vector<vector<int>> a, int x) {
    int m = a.size(), n = a[0].size();

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = a[i][j] * x % 26;
        }
    }

    return a;
}

vector<vector<int>> trans(vector<vector<int>> a) {
    int m = a.size(), n = a[0].size();
    vector<vector<int>> ans(n, vector<int>(m));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ans[i][j] = a[j][i];
        }
    }

    return ans;
}

int cofactor(vector<vector<int>> a, int x, int y) {
    int m = a.size(), n = a[0].size();
    vector<vector<int>> ans(m - 1, vector<int>(n - 1));
    
    int tempi = 0;
    for (int i = 0; i < m; i++) {
        if (i == x) continue;
        int tempj = 0;
        for (int j = 0; j < n; j++) {
            if (j == y) continue;
            ans[tempi][tempj++] = a[i][j]; 
        }
        tempi++;
    }

    return (x + y) % 2 == 0 ? det(ans) : 26 - det(ans);
}

vector<vector<int>> adj(vector<vector<int>> a) {
    int m = a.size(), n = a[0].size();
    vector<vector<int>> ans(m, vector<int>(n));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ans[i][j] = cofactor(a, j, i);
        }
    }

    return ans;
}

int main() {
    string plain, cipher;
    int m;
    #ifdef ONLINE_JUDGE
    cin >> m >> plain >> cipher;
    #else
    cin >> m >> plain >> cipher;
    #endif

    vector<vector<int>> L(m, vector<int>(m));
    vector<int> b(m);
    int len = plain.length();
    vector<vector<int>> plain_to_dig;
    vector<vector<int>> cipher_to_dig;
    for (int i = 0; i < len / m; i++) {
        vector<int> x, y;
        for (int j = 0; j < m; j++) {
            x.push_back((int)(plain[m * i + j] - 'A'));
            y.push_back((int)(cipher[m * i + j] - 'A'));
        }
        plain_to_dig.push_back(x);
        cipher_to_dig.push_back(y);
    }

    vector<vector<int>> tmpplain, tmpcipher;
    bool found = false;

    auto dfs = [&](this auto&& dfs, int i) -> void {
        int d = m + 1 - tmpplain.size(); 
        if (d == 0) { 
            vector<vector<int>> difplain, difcipher;
            for (int i = 0; i < m; i++) {
                difplain.push_back(subvec(tmpplain[i], tmpplain[i + 1]));
                difcipher.push_back(subvec(tmpcipher[i], tmpcipher[i + 1]));
            }
            int det_difplain = det(difplain);
            if (det_difplain != 0 && gcd(det_difplain, 26) == 1) {
                L = multiply(mul(adj(difplain), mulinv[det_difplain]), difcipher);
                vector<vector<int>> x = {plain_to_dig[0]};
                vector<vector<int>> y = {cipher_to_dig[0]};
                b = sub(y, multiply(x, L))[0];
                found = true;
                return;
            }
            else {
                return;
            }
        }

        for (int j = i; j >= d; j--) {
            tmpplain.push_back(plain_to_dig[j]);
            tmpcipher.push_back(cipher_to_dig[j]);
            dfs(j - 1);
            if (found) return;
            tmpplain.pop_back();
            tmpcipher.pop_back();
        }
    };
    
    dfs(len / m - 1);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            cout << L[i][j] << ' ';
        }
        cout << '\n';
    }

    for (int i = 0; i < m; i++) {
        cout << b[i] << ' ';
    }

    return 0;
}