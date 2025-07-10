#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#define DEBUG 0

using namespace std;

typedef vector<char> BigInt;
typedef vector<int64_t> Vector;
typedef vector<Vector> Matrix;

// s contains only digits [0-9]
BigInt toBigInt(const string& s) {
    BigInt ret;
    ret.reserve(s.size());
    auto c = s.begin();
    if (c != s.end() && *c == '-') {
        c++;
    }
    for(; c != s.end(); c++) {
        ret.push_back(*c - '0');
    }
    if (!ret.size()) {
        ret.push_back(0);
    }
    return ret;
}

BigInt toBigInt(uint64_t x) {
    BigInt ret;
    if (x == 0) {
        ret.push_back(0);
        return ret;
    }
    while (x > 0) {
        ret.push_back(static_cast<char>(x % 10));
        x /= 10;
    }
    reverse(ret.begin(), ret.end());
    return ret;
}

string toString(const BigInt& u) {
    string ret; ret.reserve(u.size() + 1);
    for (auto c : u) {
        ret += c + '0';
    }
    return ret;
}

// u > v => 1
// u = v => 0
// u < v => -1
int cmp(const BigInt& u, const BigInt& v) {
    int64_t i, j;
    if (u.size() > v.size()) {
        i = 0, j = v.size() - u.size();
    } else {
        i = u.size() - v.size(), j = 0;
    }
    for (; i < u.size(); i++, j++) {
        char ui = i < 0 ? 0 : u[i];
        char vj = j < 0 ? 0 : v[j];
        if (ui > vj) {
            return 1;
        } else if (ui < vj) {
            return -1;
        }
    }
    return 0;
}

BigInt add(const BigInt& u, const BigInt& v) {
    BigInt ret; ret.reserve(max(u.size(), v.size()));
    bool carry = 0;
    for (auto cu = u.rbegin(), cv = v.rbegin(); cu != u.rend() || cv != v.rend();) {
        char sum = (cu == u.rend() ? 0 : *cu) + (cv == v.rend() ? 0 : *cv) + carry;
        ret.push_back(sum % 10);
        carry = sum / 10;
        if (cu != u.rend()) cu++;
        if (cv != v.rend()) cv++;
    }
    if (carry) {
        ret.push_back(carry);
    }
    reverse(ret.begin(), ret.end());
    return ret;
}

BigInt sub(const BigInt& u, const BigInt& v) {
    BigInt ret; ret.reserve(max(u.size(), v.size()));
    bool borrow = 0;
    for (auto cu = u.rbegin(), cv = v.rbegin(); cu != u.rend() || cv != v.rend();) {
        char minuend = (cu == u.rend() ? 0 : *cu);
        char subtrahend = (cv == v.rend() ? 0 : *cv) + borrow;
        if (minuend < subtrahend) {
            minuend += 10;
            borrow = 1;
        } else {
            borrow = 0;
        }
        ret.push_back(minuend - subtrahend);
        if (cu != u.rend()) cu++;
        if (cv != v.rend()) cv++;
    }
    reverse(ret.begin(), ret.end());
    return ret;
}

// u => u / 2
BigInt div2(const BigInt& u) {
    BigInt ret; ret.reserve(u.size());
    int remainder = 0;
    for (auto c = u.begin(); c != u.end(); c++) {
        remainder = remainder * 10 + *c;
        ret.push_back(remainder / 2);
        remainder %= 2;
    }
    return ret;
}

void printMatrix(const Matrix& A) {
    for (auto row : A) {
        for (auto col : row) {
            cout << col << "\t";
        }
        cout << endl;
    }
}

Matrix matMod(const Matrix& A, int d) {
    Matrix ret(A.size(), Vector(A[0].size(), 0));
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[0].size(); j++) {
            ret[i][j] = A[i][j] % d;
        }
    }
    return ret;
}

Matrix matMul(const Matrix& A, const Matrix& B, int d) {
    Matrix ret(A.size(), Vector(B[0].size(), 0));
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < B[0].size(); j++) {
            for (size_t k = 0; k < A[0].size(); k++) {
                ret[i][j] = (ret[i][j] + (A[i][k] * B[k][j]) % d) % d;
            }
        }
    }
    return ret;
}

Matrix identityMatrix(int n) {
    Matrix ret(n, Vector(n, 0));
    for (size_t i = 0; i < n; i++) {
        ret[i][i] = 1;
    }
    return ret;
}

// A^b % d
Matrix matPow(const Matrix& A, BigInt b, int d) {
    // cout << "matPow(" << endl;
    // printMatrix(A);
    // cout << ", " << toString(b) << ", " << d << ")" << endl;
    Matrix lut[10] = {identityMatrix(A.size()), matMod(A, d)}; // {A^0, A^1, A^2, ..., A^9}
    auto ret = lut[0];
    for (int i = 2; i < 10; i++) {
        lut[i] = matMul(lut[1], lut[i - 1], d);
    }
    // for (int i = 0; i < 10; i++) {
    //     cout << "A^" << i << ":" << endl;
    //     printMatrix(lut[i]);
    //     cout << endl;
    // }
    for (char r : b) {
        auto ret2 = matMul(ret, ret, d);
        auto ret4 = matMul(ret2, ret2, d);
        auto ret8 = matMul(ret4, ret4, d);
        auto ret10 = matMul(ret8, ret2, d);
        ret = matMul(ret10, lut[r], d);
        // cout << "ret2:\n"; printMatrix(ret2); cout << endl;
        // cout << "ret4:\n"; printMatrix(ret4); cout << endl;
        // cout << "ret8:\n"; printMatrix(ret8); cout << endl;
        // cout << "ret10:\n"; printMatrix(ret10); cout << endl;
        // cout << "ret:\n"; printMatrix(ret); cout << endl;
    }
    return ret;
}

// t: target point to reach
// s: step = 2 * r - 1
// d: modulo divisor
int64_t pathCount(const BigInt& t, int s, int d) {
    if (DEBUG) cout << "pathCount(" << toString(t) << ", " << s << ", " << d << ")" << endl;
    if (cmp(t, toBigInt(0)) == 0) return 1;
    // cout << "cmp(" << toString(t) << ", " << toString(toBigInt(s)) << ") = " << cmp(t, toBigInt(s)) << endl;
    if (cmp(t, toBigInt(s)) == -1) return matPow({{2}}, sub(t, toBigInt(1)), d)[0][0];

    Matrix A(s, Vector(s, 0));
    A[s - 1][0] = 1;
    for (int i = 1; i < s; i++) {
        A[s - 1][i] = 1;
        A[i - 1][i] = 1;
    }
    Matrix x(s, Vector(1, 1));
    for (int i = 2; i < s; i++) {
        x[i][0] = (2 * x[i - 1][0]) % d;
    }
    // cout << "x:" << endl;
    // printMatrix(x);
    // cout << endl;
    auto A_to_t_s_1 = matPow(A, sub(t, toBigInt(s - 1)), d);
    if (DEBUG) cout << "A^{t-s+1}:" << endl;
    if (DEBUG) printMatrix(A_to_t_s_1);
    if (DEBUG) cout << endl;
    if (DEBUG) cout << "A^{t-s+1}x:" << endl;
    if (DEBUG) printMatrix(matMul(A_to_t_s_1, x, d));
    if (DEBUG) cout << endl;
    return matMul(A_to_t_s_1, x, d)[s - 1][0];
}

int main() {
    cout << endl;
    int r;
    while (cin >> r) {
        string su, sv; cin >> su >> sv;
        BigInt u = toBigInt(su);
        BigInt v = toBigInt(sv);
        int d; cin >> d;
        if (DEBUG) cout << "########################################" << endl;
        if (DEBUG) cout << r << " " << toString(u) << " " << toString(v) << " " << d << endl;
        if (DEBUG) cout << "########################################" << endl;
        if (r == 0) {
            if (cmp(u, toBigInt(0)) != 0 || cmp(v, toBigInt(0)) != 0) {
                cout << 0 << endl;
            } else {
                cout << 1 << endl;
            }
            continue;
        }
        if (cmp(u, v) == -1 || (u.back() - v.back()) % 2 != 0) {
            cout << 0 << endl;
            continue;
        }
        // cout << toString(u) << " + " << toString(v) << " = " << toString(add(u, v)) << endl;
        // cout << toString(u) << " - " << toString(v) << " = " << toString(sub(u, v)) << endl;
        // cout << "(" << toString(u) << " + " << toString(v) << ") / 2 = " << toString(div2(add(u, v))) << endl;
        // cout << "(" << toString(u) << " - " << toString(v) << ") / 2 = " << toString(div2(sub(u, v))) << endl;


        cout << (pathCount(div2(add(u, v)), 2 * r - 1, d) * 
                 pathCount(div2(sub(u, v)), 2 * r - 1, d)) % d << endl;
    }
    return 0;
}