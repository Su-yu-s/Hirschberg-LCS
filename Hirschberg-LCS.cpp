#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/* 安全分配内存 */
int* allocIntArray(int size) {
    int *arr = (int*)calloc(size, sizeof(int));
    if (arr == NULL) {
        printf("内存分配失败！\n");
        exit(1);
    }
    return arr;
}

/* 计算 X[0..m-1] 与 Y 的各前缀 LCS 长度 */
void prefixLengths(char *X, int m, char *Y, int n, int *left) {
    for (int j = 0; j <= n; j++) left[j] = 0;

    for (int i = 0; i < m; i++) {
        int prev = 0;
        for (int j = 0; j < n; j++) {
            int temp = left[j + 1];
            if (X[i] == Y[j]) {
                left[j + 1] = prev + 1;
            } else {
                left[j + 1] = MAX(left[j + 1], left[j]);
            }
            prev = temp;
        }
    }
}

/* 计算 X[0..m-1] 与 Y 的各后缀 LCS 长度 */
void suffixLengths(char *X, int m, char *Y, int n, int *right) {
    for (int j = 0; j <= n; j++) right[j] = 0;

    for (int i = m - 1; i >= 0; i--) {
        int prev = 0;
        for (int j = n - 1; j >= 0; j--) {
            int temp = right[j];
            if (X[i] == Y[j]) {
                right[j] = prev + 1;
            } else {
                right[j] = MAX(right[j], right[j + 1]);
            }
            prev = temp;
        }
    }
}

/* Hirschberg 分治算法：还原一条 LCS */
void hirschberg(char *X, int m, char *Y, int n, char *lcs, int *pos) {
    if (m == 0 || n == 0) return;

    if (m == 1) {
        for (int j = 0; j < n; j++) {
            if (X[0] == Y[j]) {
                lcs[(*pos)++] = X[0];
                break;
            }
        }
        return;
    }

    int mid = m / 2;

    int *left = allocIntArray(n + 1);
    int *right = allocIntArray(n + 1);

    prefixLengths(X, mid, Y, n, left);
    suffixLengths(X + mid, m - mid, Y, n, right);

    int best = 0;
    int maxLen = -1;

    for (int j = 0; j <= n; j++) {
        int curLen = left[j] + right[j];
        if (curLen > maxLen) {
            maxLen = curLen;
            best = j;
        }
    }

    free(left);
    free(right);

    hirschberg(X, mid, Y, best, lcs, pos);
    hirschberg(X + mid, m - mid, Y + best, n - best, lcs, pos);
}

/* 对外调用函数 */
void getLCS_Hirschberg(char X[], char Y[], char lcs[]) {
    int pos = 0;
    int m = strlen(X);
    int n = strlen(Y);

    hirschberg(X, m, Y, n, lcs, &pos);
    lcs[pos] = '\0';
}

/* 使用滚动数组计算 LCS 长度 */
int lcsLength(char X[], char Y[]) {
    int m = strlen(X);
    int n = strlen(Y);

    int *dp = allocIntArray(n + 1);

    for (int i = 1; i <= m; i++) {
        int prev = 0;
        for (int j = 1; j <= n; j++) {
            int temp = dp[j];
            if (X[i - 1] == Y[j - 1]) {
                dp[j] = prev + 1;
            } else {
                dp[j] = MAX(dp[j], dp[j - 1]);
            }
            prev = temp;
        }
    }

    int result = dp[n];
    free(dp);
    return result;
}

int main() {
    char *X[] = {
        "ABCBDAB",
        "AAAAA",
        "ABCDEF",
        "ABCDEF",
        "AGGTAB",
        "",
        "A",
        "XMJYAUZ"
    };

    char *Y[] = {
        "BDCABA",
        "AAA",
        "GHIJKL",
        "ABCDEF",
        "GXTXAYB",
        "ABC",
        "BCA",
        "MZJAWXU"
    };

    int testCount = 8;

    for (int t = 0; t < testCount; t++) {
        int maxLcsLen = MIN(strlen(X[t]), strlen(Y[t]));
        char *lcs = (char*)malloc((maxLcsLen + 1) * sizeof(char));

        if (lcs == NULL) {
            printf("内存分配失败！\n");
            return 1;
        }

        int len = lcsLength(X[t], Y[t]);
        getLCS_Hirschberg(X[t], Y[t], lcs);

        printf("========== 测试用例 %d ==========\n", t + 1);
        printf("序列X: %s\n", X[t]);
        printf("序列Y: %s\n", Y[t]);
        printf("LCS长度: %d\n", len);
        printf("LCS序列: \"%s\"\n", lcs);
        printf("\n");

        free(lcs);
    }

    return 0;
}
