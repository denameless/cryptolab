#include <bits/stdc++.h>

constexpr uint64_t N = 2;

long long cnt = 0;

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

uint32_t count_leading_zero(uint64_t* a) {
    uint32_t ans = 0;

    for (int i = N - 1; i >= 0; i--) {
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

inline void endian_from_little_to_big(uint64_t* a) {
    for (int i = 0; i < N / 2; i++) {
        std::swap(a[i], a[N - 1 - i]);
        a[i] = __builtin_bswap64(a[i]);
        a[N - 1 - i] = __builtin_bswap64(a[N - 1 - i]);
    }
}

inline void cvrt_str_to_array(std::string str, uint64_t* ans) {
    /* convert decimal string to uint64_t array. the converted number's length is 2048 bits*/
    int len = str.length();
    int ptr = 0;
    
    for (uint32_t i = 0; i < N; i++) {
        for (uint32_t j = 0; j < 64; j++) {
            int borrow = 0;
            bool ok = false;

            for (int k = ptr; k < str.size(); k++) {
                if (k == ptr && str[k] == '1') ok = true;

                int cur = str[k] - '0' + 0;
                char tmp = (char) ('0' + (borrow * 10 + cur) / 2);
                str[k] = tmp;

                if (cur & 1) {
                    borrow = 1;
                }
                else {
                    borrow = 0;
                }
            }

            if (borrow) {
                ans[i] |= (1ULL << j);
            }

            if (ok) ptr++;
        }
    }

    return;
}

inline void cvrt_array_to_string(std::string &str, uint64_t* arr) {
    uint64_t q[N] = {0};
    int msb = ((N << 6) - 1 - (int) count_leading_zero(arr)) / 64;

    while (1) {
        bool greaterthanzero = false;
        uint64_t cur_remainder = 0;

        for (int i = msb; i >= 0; i--) {
            if (arr[i]) {
                greaterthanzero = true;
                break;
            }
            else {
                msb--;
            }
        }

        if (!greaterthanzero) {
            break;
        }

        for (int i = msb; i >= 0; i--) {
            __uint128_t cur = ((__uint128_t) cur_remainder << 64) + arr[i];
            cur_remainder = cur % 10;
            arr[i] = (uint64_t) (cur / (__uint128_t) 10);
        }

        char ch = '0' + cur_remainder;
        str.insert(0, 1, ch);
    }
}

inline void trivial_multiple_precision_subtraction(uint64_t* a, uint64_t* b, uint64_t* ans) {
    /* a should be bigger than b */
    uint64_t c = 0;

    for (int i = 0; i < N; i++) {
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

inline void trivial_multiple_precision_subtraction_v2(uint64_t* a, uint64_t* b, uint64_t* ans) {
    /* a should be bigger than b */
    uint64_t c = 0;

    for (int i = 0; i < N + 1; i++) {
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

    ans[N] &= 1ULL;
}

inline bool trivial_multiple_precision_addition(uint64_t* a, uint64_t* b, uint64_t* ans) {
    /* result should be 2048 bits */
    uint64_t c = 0;

    for (int i = 0; i < N - 1; i++) {
        uint64_t tmp = a[i] + b[i] + c;

        if (a[i] == 0 && b[i] == 0 && c == 0) continue;

        if ((c == 1 && (tmp <= a[i] || tmp <= b[i])) || (c == 0 && tmp < a[i] && tmp < b[i])) {
            c = 1;
        }
        else {
            c = 0;
        }

        ans[i] = tmp;
    }

    if (a[N - 1] == 0 && b[N - 1] == 0 && c == 0) return false;

    uint64_t tmp = a[N - 1] + b[N - 1] + c;

    if (tmp <= a[N - 1] || tmp <= b[N - 1]) {
        ans[N - 1] = tmp;
        return true;
    }
    else {
        ans[N - 1] = tmp;
        return false;
    }
}

inline bool trivial_multiple_precision_addition_v2(uint64_t* a, uint64_t* b, uint64_t* ans) {
    /* this function implement the montgomery multiplication. the only difference from previous function is that operands can be 2049 bits */
    uint64_t c = 0;

    for (int i = 0; i < N; i++) {
        uint64_t tmp = a[i] + b[i] + c;

        if (a[i] == 0 && b[i] == 0 && c == 0) continue;

        if ((c == 1 && (tmp <= a[i] || tmp <= b[i])) || (c == 0 && tmp < a[i] && tmp < b[i])) {
            c = 1;
        }
        else {
            c = 0;
        }

        ans[i] = tmp;
    }

    if (a[N] == 0 && b[N] == 0 && c == 0) return false;

    uint64_t tmp = a[N] + b[N] + c;

    if ((c == 1 && (tmp <= a[N] || tmp <= b[N])) || (c == 0 && tmp < a[N] && tmp < b[N])) {
        ans[N] = tmp;
        return true;
    }
    else {
        ans[N] = tmp;
        return false;
    }
}

inline bool trivial_multiple_precision_addition_v3(uint64_t* a, uint64_t* b, uint64_t* ans) {
    /* this function implement the montgomery multiplication. the only difference from previous function is that operands can be 2049 bits */
    uint64_t c = 0;

    for (int i = 0; i < N; i++) {
        uint64_t tmp = a[i] + b[i] + c;

        if (a[i] == 0 && b[i] == 0 && c == 0) continue;

        if ((c == 1 && (tmp <= a[i] || tmp <= b[i])) || (c == 0 && tmp < a[i] && tmp < b[i])) {
            c = 1;
        }
        else {
            c = 0;
        }

        ans[i] = tmp;
    }

    if (a[N] == 0 && b[N] == 0 && c == 0) return false;

    uint64_t tmp = a[N] + b[N] + c;

    if ((c == 1 && (tmp <= a[N] || tmp <= b[N])) || (c == 0 && tmp < a[N] && tmp < b[N])) {
        ans[N] = tmp & 1ULL;
        return true;
    }
    else {
        ans[N] = tmp & 1ULL;
        return false;
    }
    // tmp &= 1ULL;
    // if (!a[32] && !b[32]) tmp &= 0ULL;
    // return true;
}

inline void trivial_multiple_precision_addition_v4(uint64_t* a, uint64_t* b, uint64_t* ans) {
    uint64_t c = 0;

    for (int i = 0; i < N + 1; i++) {
        uint64_t tmp = a[i] + b[i] + c;

        if (a[i] == 0 && b[i] == 0 && c == 0) continue;

        if ((c == 1 && (tmp <= a[i] || tmp <= b[i])) || (c == 0 && tmp < a[i] && tmp < b[i])) {
            c = 1;
        }
        else {
            c = 0;
        }

        ans[i] = tmp;
    }

    ans[N + 1] += c;

    return;
}

inline void multi_mult_single(const uint64_t* a, const uint64_t b, uint64_t* ans) {
    uint64_t hi_save = 0;

    for (int i = 0; i < N; i++) {
        __uint128_t tmp = (__uint128_t) a[i] * (__uint128_t) b + (__uint128_t) hi_save;
        ans[i] = (uint64_t) tmp;
        hi_save = tmp >> 64;
    }

    ans[N] = hi_save;
}

inline void rlli(uint64_t* a) {
    for (int i = 0; i < N + 1; i++) {
        a[i] = a[i + 1];
    }
    
    a[N + 1] = 0;
}

inline void mult_by_2(uint64_t* a) {
    /* this function aimed to implement the compute of R^2 mod m in montgomery modulo */
    bool ok = false;
    bool tmp = false;
    for (int i = 0; i < N; i++) {
        if (a[i] & (1ULL << 63)) {
            tmp = true;
        }
        else {
            tmp = false;
        }

        a[i] = (a[i] << 1) | (ok ? 1 : 0);

        ok = tmp;
    }
    a[N] |= ok ? 1 : 0;
}

bool cmp_2048bits(const uint64_t* a, const uint64_t* b) {
    /* true when a >= b */
    for (int i = N - 1; i >= 0; i--) {
        if (a[i] != b[i]) {
            return a[i] >= b[i] ? true : false;
        }
    }

    return true;
}

bool cmp_2049bits(const uint64_t* a, const uint64_t* b) {
    for (int i = N; i >= 0; i--) {
        if (a[i] != b[i]) {
            return a[i] >= b[i] ? true : false;
        }
    }

    return true;
}

bool cmp_2048bits_v2(const uint64_t* a, const uint64_t* b) {
    /* true when $a > $b */
    for (int i = N - 1; i >= 0; i--) {
        if (a[i] != b[i]) {
            return a[i] > b[i] ? true : false;
        }
    }

    return true;
}

bool cmp_2049bits_v2(const uint64_t* a, const uint64_t* b) {
    /* true when $a > $b */
    for (int i = N; i >= 0; i--) {
        if (a[i] != b[i]) {
            return a[i] > b[i] ? true : false;
        }
    }

    return true;
}

inline void montgomery_init_module(uint64_t* a, uint64_t* m) {
    /* warning: this function compute R^2 mod m, in order to implement the montgomery multiplication */
    for (int i = 0; i < (N << 7); i++) {
        mult_by_2(a);

        if (cmp_2049bits(a, m)) {
            trivial_multiple_precision_subtraction_v2(a, m, a);
        }
    }
    return;
}

inline void montgomery_init_module_v2(uint64_t* a, uint64_t* m) {
    /* warning: this function compute R mod m, in order to implement the montgomery exponentiation */
    for (int i = 0; i < (N << 6); i++) {
        mult_by_2(a);

        if (cmp_2049bits(a, m)) {
            trivial_multiple_precision_subtraction_v2(a, m, a);
        }
    }
}

// inline void montgomery_multiplication_v2(uint64_t* m, const uint64_t* x, const uint64_t* y, uint64_t* ans, uint64_t neg_m_inverse) {
//     /* OUTPUT: (x * y * (R^-1)) mod m */
//     uint64_t A[34] = {0};

//     for (int i = 0; i < 32; i++) {
//         uint64_t u_i = (uint64_t) (((__uint128_t) A[0] + (__uint128_t) x[i] * (__uint128_t) y[0]) * (__uint128_t) neg_m_inverse);
//         uint64_t xi_y[33] = {0};
//         multi_mult_single(y, x[i], xi_y);
//         uint64_t ui_m[33] = {0};
//         multi_mult_single(m, u_i, ui_m);

//         trivial_multiple_precision_addition_v4(A, xi_y, A);
//         trivial_multiple_precision_addition_v4(A, ui_m, A);
//         rlli(A);
//     }

//     if (A[32] > 0 || cmp_2048bits(A, m)) {
//         trivial_multiple_precision_subtraction(A, m, A);
//     }

//     memcpy(ans, A, 256);
// }

inline void montgomery_multiplication_v2(uint64_t* m, const uint64_t* x, const uint64_t* y, uint64_t* ans, uint64_t neg_m_inverse) {
    uint64_t T[N + 2] = {0}; 

    for (int i = 0; i < N; i++) {
        __uint128_t carry = 0;

        for (int j = 0; j < N; j++) {
            __uint128_t tmp = (__uint128_t) x[i] * (__uint128_t) y[j] + (__uint128_t) T[j] + carry;
            T[j] = (uint64_t) tmp;
            carry = tmp >> 64;
        }
        __uint128_t tmp2 = (__uint128_t) T[N] + carry;
        T[N] = (uint64_t) tmp2;
        T[N + 1] = tmp2 >> 64;

        uint64_t u = T[0] * neg_m_inverse;

        __uint128_t tmp3 = (__uint128_t) u * (__uint128_t) m[0] + (__uint128_t) T[0];
        carry = tmp3 >> 64;
        
        for (int j = 1; j < N; j++) {
            __uint128_t tmp4 = (__uint128_t) u * (__uint128_t) m[j] + (__uint128_t) T[j] + (__uint128_t) carry;
            T[j - 1] = (uint64_t) tmp4;
            carry = tmp4 >> 64;
        }

        __uint128_t tmp5 = (__uint128_t) T[N] + carry;
        T[N - 1] = (uint64_t)tmp5;
        carry = tmp5 >> 64;
        
        tmp5 = (__uint128_t) T[N + 1] + carry;
        T[N] = (uint64_t) tmp5;
        T[N + 1] = 0;
    }

    if (T[N] || cmp_2048bits(T, m)) {
        trivial_multiple_precision_subtraction(T, m, ans);
    }
    else {
        memcpy(ans, T, (N << 3));
    }
}

inline void rlli_33words_byonebit(uint64_t* a) {
    uint64_t remain = 0;
    for (int i = N; i >= 0; i--) {
        uint64_t tmp = (a[i] >> 1) | remain;
        remain = a[i] << 63;
        a[i] = tmp;
    }

    if (a[N - 1] & (1ULL << 63)) a[N] |= 1ULL;
    return;
}

inline void rlli_33words_byonebit_v2(uint64_t* a) {
    uint64_t remain = 0;
    for (int i = N; i >= 0; i--) {
        uint64_t tmp = (a[i] >> 1) | remain;
        remain = a[i] << 63;
        a[i] = tmp;
    }

    return;
}

inline void rlli_32words(uint64_t* a, uint64_t shift) {
    if (shift >= (N << 6)) {
        memset(a, 0, (N << 3));
        return;
    }

    uint64_t word_shift = shift >> 6;
    uint64_t bit_shift = shift & 0x3f;

    if (word_shift > 0) {
        memmove(a, a + word_shift, (N - word_shift) << 3);
        memset(a + N - word_shift, 0, word_shift << 3);
    }

    if (bit_shift > 0) {
        uint64_t remain = 0;
        for (int i = N - 1 - word_shift; i >= 0; i--) {
            uint64_t tmp = (a[i] >> bit_shift) | remain;
            remain = a[i] << (64 - bit_shift);
            a[i] = tmp;
        }
    }

    return;
}

// inline void binary_extended_gcd(uint64_t* a, uint64_t* b, uint64_t* x, uint64_t* y) {
//     uint64_t c[N + 1] = {0};
//     uint64_t d[N + 1] = {0};
//     uint64_t u[N + 1], v[N + 1];
//     memcpy(u, x, ((N + 1) << 3));
//     memcpy(v, y, ((N + 1) << 3));

//     d[0] = 1ULL;

//     while (1) {
//         bool isequal = true;
//         for (int i = N; i >= 0; i--) {
//             if (u[i] != v[i]) {
//                 isequal = false;
//                 break;
//             }
//         }
//         if (isequal) return;

//         if (u[0] % 2 == 0) {
//             rlli_32words(u, 1);

//             if (a[0] % 2 == 0 && b[0] % 2 == 0) {
//                 rlli_33words_byonebit(a);
//                 rlli_33words_byonebit(b);
//             }
//             else {
//                 trivial_multiple_precision_addition_v3(a, y, a);
//                 rlli_33words_byonebit(a);

//                 trivial_multiple_precision_subtraction_v2(b, x, b);
//                 rlli_33words_byonebit(b);
//             }
//         }
//         else if (v[0] % 2 == 0) {
//             rlli_32words(v, 1);

//             if (c[0] % 2 == 0 && d[0] % 2 == 0) {
//                 rlli_33words_byonebit(c);
//                 rlli_33words_byonebit(d);
//             }
//             else {
//                 trivial_multiple_precision_addition_v3(c, y, c);
//                 rlli_33words_byonebit(c);

//                 trivial_multiple_precision_subtraction_v2(d, x, d);
//                 rlli_33words_byonebit(d);
//             }
//         }
//         else if (cmp_2049bits_v2(u, v)) {
//             trivial_multiple_precision_subtraction_v2(u, v, u);
//             trivial_multiple_precision_subtraction_v2(a, c, a);
//             trivial_multiple_precision_subtraction_v2(b, d, b);
//         }
//         else {
//             trivial_multiple_precision_subtraction_v2(v, u, v);
//             trivial_multiple_precision_subtraction_v2(c, a, c);
//             trivial_multiple_precision_subtraction_v2(d, b, d);
//         }
//     }

//     return;
// }

inline void binary_extended_gcd_v2(uint64_t* a, uint64_t* b, uint64_t* x, uint64_t* y) {
    uint64_t c[N + 1] = {0};
    uint64_t d[N + 1] = {0};
    uint64_t u[N + 1], v[N + 1];
    memcpy(u, x, ((N + 1) << 3));
    memcpy(v, y, ((N + 1) << 3));

    d[0] = 1ULL;

    while (1) {
        bool isequal = true;
        for (int i = N; i >= 0; i--) {
            if (u[i] != v[i]) {
                isequal = false;
                break;
            }
        }
        if (isequal) return;

        while (u[0] % 2 == 0) {
            rlli_32words(u, 1);

            if (a[0] % 2 == 0) {
                rlli_33words_byonebit_v2(a);
            }
            else {
                trivial_multiple_precision_addition_v2(a, y, a);
                rlli_33words_byonebit_v2(a);
            }
        }
        while (v[0] % 2 == 0) {
            rlli_32words(v, 1);

            if (b[0] % 2 == 0) {
                rlli_33words_byonebit_v2(b);
            }
            else {
                trivial_multiple_precision_addition_v2(b, y, b);
                rlli_33words_byonebit_v2(b);
            }
        }

        bool ok = true;
        for (int i = N; i >= 0; i--) {
            if (u[i] != v[i]) {
                ok = false;
                break;
            }
        }
        if (ok) return;

        if (cmp_2048bits(u, v)) {
            trivial_multiple_precision_subtraction(u, v, u);
            if (!cmp_2048bits(a, b)) trivial_multiple_precision_addition_v2(a, y, a);
            trivial_multiple_precision_subtraction_v2(a, b, a);
        }
        else {
            trivial_multiple_precision_subtraction(v, u, v);
            if (!cmp_2048bits(b, a)) trivial_multiple_precision_addition_v2(b, y, b);
            trivial_multiple_precision_subtraction_v2(b, a, b);
        }
    }

    return;
}

/* new */

// __uint128_t ex_gcd(__uint128_t a, __uint128_t b, __uint128_t& x, __uint128_t& y) {
//   if (!b) {
//     x = 1;
//     y = 0;
//     return a;
//   } else {
//     __uint128_t d = ex_gcd(b, a % b, y, x);
//     y -= a / b * x;
//     return d;
//   }
// }

// __uint128_t solve_linear_congruence_equation(__uint128_t a, __uint128_t b, __uint128_t n) {
//   __uint128_t x, y;
//   __uint128_t d = ex_gcd(a, n, x, y);
//   if (b % d) return -1;
//   n /= d;
//   return (x * (b / d) % n + n) % n;
// }

int main() {
    std::string ps, ns, alphas, betas;
    std::cin >> ps >> ns >> alphas >> betas;

    uint64_t p[N] = {0}, n[N] = {0}, alpha[N] = {0}, beta[N] = {0}, p_extend[N + 1] = {0}, n_extend[N + 1] = {0};
    cvrt_str_to_array(ps, p);
    cvrt_str_to_array(ns, n);
    cvrt_str_to_array(alphas, alpha);
    cvrt_str_to_array(betas, beta);
    memcpy(p_extend, p, (N << 3));
    memcpy(n_extend, n, (N << 3));

    uint64_t x[N] = {0}, a[N] = {0}, b[N] = {0};
    memcpy(x, beta, (N << 3));
    b[0] = 1ULL;
    uint64_t x_[N] = {0}, a_[N] = {0}, b_[N] = {0};

    uint64_t R_mod_p[N + 1] = {0}; R_mod_p[0] = 1ULL;
    uint64_t R_mod_n[N + 1] = {0}; R_mod_n[0] = 1ULL;
    uint64_t R2_mod_p[N + 1] = {0}; R2_mod_p[0] = 1ULL;
    uint64_t R2_mod_n[N + 1] = {0}; R2_mod_n[0] = 1ULL;

    uint64_t one[N] = {0}; one[0] = 1ULL;
    uint64_t two[N] = {0}; two[0] = 2ULL;

    montgomery_init_module_v2(R_mod_p, p_extend);
    montgomery_init_module_v2(R_mod_n, n_extend);
    montgomery_init_module(R2_mod_p, p_extend);
    montgomery_init_module(R2_mod_n, n_extend);

    uint64_t twoR_mod_n[N]; memcpy(twoR_mod_n, R_mod_n, (N << 3));
    mult_by_2(twoR_mod_n);

    if (cmp_2048bits(twoR_mod_n, n)) {
        trivial_multiple_precision_subtraction(twoR_mod_n, n, twoR_mod_n);
    }

    uint64_t neg_p_inverse = compute_neg_m_inverse(p);
    uint64_t neg_n_inverse = compute_neg_m_inverse(n);

    montgomery_multiplication_v2(p, x, R2_mod_p, x, neg_p_inverse);
    montgomery_multiplication_v2(n, b, R2_mod_n, b, neg_n_inverse);
    memcpy(x_, x, (N << 3)), memcpy(a_, a, (N << 3)), memcpy(b_, b, (N << 3));

    uint64_t alpha_mult_R_mod_p[N] = {0};
    memcpy(alpha_mult_R_mod_p, alpha, (N << 3));
    montgomery_multiplication_v2(p, alpha_mult_R_mod_p, R2_mod_p, alpha_mult_R_mod_p, neg_p_inverse);
    uint64_t beta_mult_R_mod_p[N] = {0};
    memcpy(beta_mult_R_mod_p, beta, (N << 3));
    montgomery_multiplication_v2(p, beta_mult_R_mod_p, R2_mod_p, beta_mult_R_mod_p, neg_p_inverse);

    uint64_t tmpx[N] = {0}, tmpx_[N] = {0};
    memcpy(tmpx, beta, (N << 3));
    memcpy(tmpx_, beta, (N << 3));

    if (tmpx_[0] % 3 == 0) {
        montgomery_multiplication_v2(p, x_, x_, x_, neg_p_inverse);
        montgomery_multiplication_v2(n, a_, twoR_mod_n, a_, neg_n_inverse);
        montgomery_multiplication_v2(n, b_, twoR_mod_n, b_, neg_n_inverse);
        montgomery_multiplication_v2(p, x_, one, tmpx_, neg_p_inverse);
    }
    else if (tmpx_[0] % 3 == 1) {
        montgomery_multiplication_v2(p, x_, beta_mult_R_mod_p, x_, neg_p_inverse);
        trivial_multiple_precision_addition(b_, R_mod_n, b_);
        if (cmp_2048bits(b_, n)) trivial_multiple_precision_subtraction(b_, n, b_);

        montgomery_multiplication_v2(p, x_, one, tmpx_, neg_p_inverse);
    }
    else {
        montgomery_multiplication_v2(p, x_, alpha_mult_R_mod_p, x_, neg_p_inverse);
        trivial_multiple_precision_addition(a_, R_mod_n, a_);
        if (cmp_2048bits(a_, n)) trivial_multiple_precision_subtraction(a_, n, a_);

        montgomery_multiplication_v2(p, x_, one, tmpx_, neg_p_inverse);
    }

    

    while (1) {
        // if (cnt % 10000 == 0) {
        //     std::cout << cnt << '\n';
        // }
        // cnt++;
        
        bool isequal = true;
        for (int i = N - 1; i >= 0; i--) {
            if (x[i] != x_[i]) {
                isequal = false;
                break;
            }
        }
        if (isequal) break;

        if (tmpx[0] % 3 == 0) {
            montgomery_multiplication_v2(p, x, x, x, neg_p_inverse);
            montgomery_multiplication_v2(n, a, twoR_mod_n, a, neg_n_inverse);
            montgomery_multiplication_v2(n, b, twoR_mod_n, b, neg_n_inverse);
            montgomery_multiplication_v2(p, x, one, tmpx, neg_p_inverse);
        }
        else if (tmpx[0] % 3 == 1) {
            montgomery_multiplication_v2(p, x, beta_mult_R_mod_p, x, neg_p_inverse);
            trivial_multiple_precision_addition(b, R_mod_n, b);
            if (cmp_2048bits(b, n)) trivial_multiple_precision_subtraction(b, n, b);

            montgomery_multiplication_v2(p, x, one, tmpx, neg_p_inverse);
        }
        else {
            montgomery_multiplication_v2(p, x, alpha_mult_R_mod_p, x, neg_p_inverse);
            trivial_multiple_precision_addition(a, R_mod_n, a);
            if (cmp_2048bits(a, n)) trivial_multiple_precision_subtraction(a, n, a);

            montgomery_multiplication_v2(p, x, one, tmpx, neg_p_inverse);
        }

        if (tmpx_[0] % 3 == 0) {
            montgomery_multiplication_v2(p, x_, x_, x_, neg_p_inverse);
            montgomery_multiplication_v2(n, a_, twoR_mod_n, a_, neg_n_inverse);
            montgomery_multiplication_v2(n, b_, twoR_mod_n, b_, neg_n_inverse);
            montgomery_multiplication_v2(p, x_, one, tmpx_, neg_p_inverse);
        }
        else if (tmpx_[0] % 3 == 1) {
            montgomery_multiplication_v2(p, x_, beta_mult_R_mod_p, x_, neg_p_inverse);
            trivial_multiple_precision_addition(b_, R_mod_n, b_);
            if (cmp_2048bits(b_, n)) trivial_multiple_precision_subtraction(b_, n, b_);

            montgomery_multiplication_v2(p, x_, one, tmpx_, neg_p_inverse);
        }
        else {
            montgomery_multiplication_v2(p, x_, alpha_mult_R_mod_p, x_, neg_p_inverse);
            trivial_multiple_precision_addition(a_, R_mod_n, a_);
            if (cmp_2048bits(a_, n)) trivial_multiple_precision_subtraction(a_, n, a_);

            montgomery_multiplication_v2(p, x_, one, tmpx_, neg_p_inverse);
        }

        if (tmpx_[0] % 3 == 0) {
            montgomery_multiplication_v2(p, x_, x_, x_, neg_p_inverse);
            montgomery_multiplication_v2(n, a_, twoR_mod_n, a_, neg_n_inverse);
            montgomery_multiplication_v2(n, b_, twoR_mod_n, b_, neg_n_inverse);
            montgomery_multiplication_v2(p, x_, one, tmpx_, neg_p_inverse);
        }
        else if (tmpx_[0] % 3 == 1) {
            montgomery_multiplication_v2(p, x_, beta_mult_R_mod_p, x_, neg_p_inverse);
            trivial_multiple_precision_addition(b_, R_mod_n, b_);
            if (cmp_2048bits(b_, n)) trivial_multiple_precision_subtraction(b_, n, b_);

            montgomery_multiplication_v2(p, x_, one, tmpx_, neg_p_inverse);
        }
        else {
            montgomery_multiplication_v2(p, x_, alpha_mult_R_mod_p, x_, neg_p_inverse);
            trivial_multiple_precision_addition(a_, R_mod_n, a_);
            if (cmp_2048bits(a_, n)) trivial_multiple_precision_subtraction(a_, n, a_);

            montgomery_multiplication_v2(p, x_, one, tmpx_, neg_p_inverse);
        }
    }

    montgomery_multiplication_v2(n, a, one, a, neg_n_inverse);
    montgomery_multiplication_v2(n, b, one, b, neg_n_inverse);
    montgomery_multiplication_v2(n, a_, one, a_, neg_n_inverse);
    montgomery_multiplication_v2(n, b_, one, b_, neg_n_inverse);

    uint64_t A[N] = {0};
    if (cmp_2048bits(b_, b)) {
        trivial_multiple_precision_subtraction(b_, b, A);
    }
    else {
        trivial_multiple_precision_subtraction(n, b, A);
        trivial_multiple_precision_addition(A, b_, A);
    }

    uint64_t B[N] = {0};
    if (cmp_2048bits(a, a_)) {
        trivial_multiple_precision_subtraction(a, a_, B);
    }
    else {
        trivial_multiple_precision_subtraction(n, a_, B);
        trivial_multiple_precision_addition(B, a, B);
    }

    // if (std::gcd((b_ >= b ? (b_ -  b) : (n - b + b_)), n) == 1) {
    //     __uint128_t inv = 1;
    //     __uint128_t cur = (b_ >= b ? (b_ -  b) : (n - b + b_));
    //     for (__uint128_t i = 0; i < n - 2; i++) {
    //         inv = inv * cur % n;
    //     }

    //     __uint128_t tmp = (a >= a_ ? (a - a_) * inv % n : (n - a_ + a) * inv % n);

    //     std::string ans;
    //     cvrt_array_to_string(ans, tmp);
    //     std::cout << ans << '\n';
    //     return 0;
    // }
    // else {
    //     __uint128_t A = (b_ >= b ? (b_ -  b) : (n - b + b_)), B = (a >= a_ ? (a - a_) : (n - a_ + a));
    //     __uint128_t tmp = solve_linear_congruence_equation(A, B, n);

    //     std::string ans;
    //     cvrt_array_to_string(ans, tmp);
    //     std::cout << ans << '\n';
    //     return 0;
    // }
    
    uint64_t A_mult_R_mod_n[N] = {0};
    montgomery_multiplication_v2(n, A, R2_mod_n, A_mult_R_mod_n, neg_n_inverse);
    uint64_t iteration_new[N] = {0};
    memcpy(iteration_new, A_mult_R_mod_n, (N << 3));
    uint64_t B_mult_R_mod_n[N] = {0};
    montgomery_multiplication_v2(n, B, R2_mod_n, B_mult_R_mod_n, neg_n_inverse);

    uint64_t n_sub_two[N] = {0};
    trivial_multiple_precision_subtraction(n, two, n_sub_two);

    uint64_t ans_pow[N] = {0};
    memcpy(ans_pow, R_mod_n, (N << 3));

    for (int i = 0; i < (N << 6); i++, montgomery_multiplication_v2(n, iteration_new, iteration_new, iteration_new, neg_n_inverse)) {
        if ((n_sub_two[i >> 6] >> (i % 64)) & 1ULL) {
            montgomery_multiplication_v2(n, ans_pow, iteration_new, ans_pow, neg_n_inverse);
        }
    }

    uint64_t ans[N] = {0};
    montgomery_multiplication_v2(n, B_mult_R_mod_n, ans_pow, ans, neg_n_inverse);
    montgomery_multiplication_v2(n, ans, one, ans, neg_n_inverse);

    std::string ans_string;
    cvrt_array_to_string(ans_string, ans);
    std::cout << ans_string << '\n';
    return 0;
}