# Longest Common Subsequence (LCS) Solver & Sequence Reconstruction

Dynamic programming with Hirschberg's divide-and-conquer algorithm

## Table of Contents

1. [Introduction](#1-introduction)
2. [Problem Analysis](#2-problem-analysis)
3. [Mathematical Model](#3-mathematical-model)
4. [Algorithm Design](#4-algorithm-design)
5. [Code Walkthrough](#5-code-walkthrough)
6. [Test Cases & Results](#6-test-cases--results)
7. [Complexity Analysis](#7-complexity-analysis)
8. [Comparison with Traditional DP](#8-comparison-with-traditional-dp)
9. [Future Improvements](#9-future-improvements)
10. [Build & Run](#10-build--run)
11. [References](#11-references)

---

## 1. Introduction

The **Longest Common Subsequence (LCS)** problem is a classic dynamic programming challenge: given two sequences, find the longest subsequence common to both (elements need not be contiguous but must maintain relative order).

This project implements the LCS solver using **dynamic programming combined with Hirschberg's divide-and-conquer algorithm**, reconstructing an actual LCS sequence alongside computing its length. By leveraging rolling arrays and recursive partitioning, the space complexity is reduced from O(mn) to **O(n)** — essential for handling longer sequences within memory constraints.

---

## 2. Problem Analysis

### 2.1 Background

LCS has wide applications across multiple domains:

| Domain | Application |
|--------|-------------|
| **Bioinformatics** | DNA sequence alignment, gene similarity assessment |
| **Version Control** | Core algorithm behind `git diff` for computing file differences |
| **Text Processing** | Document similarity measurement, spell checking |
| **GIS** | Line feature similarity metrics |

### 2.2 Why Dynamic Programming?

A brute-force approach would enumerate all possible subsequences, yielding **O(2^min(m,n))** time complexity — infeasible for anything beyond trivial inputs.

LCS possesses two key properties that make DP applicable:

1. **Optimal substructure**: The LCS of two sequences contains within it the LCS of their prefixes
2. **Overlapping subproblems**: Many prefix-pair LCS computations are reused across different branches

---

## 3. Mathematical Model

### 3.1 State Definition

Let X = ⟨x₁, x₂, …, xₘ⟩ and Y = ⟨y₁, y₂, …, yₙ⟩.

Define `dp[i][j]` as the length of the LCS between **the first i elements of X** and **the first j elements of Y**.

### 3.2 State Transition Equation

```
                    ⎧ 0                                    if i = 0 or j = 0
dp[i][j] =          ⎨ dp[i-1][j-1] + 1                     if xᵢ = yⱼ
                    ⎩ max(dp[i-1][j], dp[i][j-1])          if xᵢ ≠ yⱼ
```

### 3.3 Optimal Substructure Proof

Let Z = ⟨z₁, z₂, …, zₖ⟩ be any LCS of X and Y:

1. **If xₘ = yₙ**: then zₖ = xₘ = yₙ, and the prefix Z[1..k-1] is an LCS of X[1..m-1] and Y[1..n-1]
2. **If xₘ ≠ yₙ and zₖ ≠ xₘ**: then Z is an LCS of X[1..m-1] and Y
3. **If xₘ ≠ yₙ and zₖ ≠ yₙ**: then Z is an LCS of X and Y[1..n-1]

---

## 4. Algorithm Design

### 4.1 Limitations of Traditional DP

The traditional approach requires a full `(m+1) × (n+1)` 2D DP table plus a direction-tracking table for backtracking — impractical for long sequences.

**Two improvements:**
- Computing LCS **length**: Use a rolling 1D array → space O(mn) → O(n)
- Reconstructing the LCS **sequence**: Hirschberg's divide-and-conquer → space O(n)

### 4.2 Hirschberg's Divide-and-Conquer Algorithm

#### Core Idea

Split sequence X at its midpoint into X_left and X_right. Compute:
- `left[j]`: LCS length of X_left with Y[0..j-1] (forward pass)
- `right[j]`: LCS length of X_right with Y[j..n-1] (backward pass)

Find the split point `best` that maximizes `left[j] + right[j]`, divide Y accordingly, and recurse on both halves.

#### Correctness Argument

For any split point j, `left[j] + right[j]` represents the total LCS length when Y is divided at position j. The maximum identifies the optimal split, ensuring the concatenated recursive results form a global optimum.

#### Pseudocode

```
function Hirschberg(X[0..m-1], Y[0..n-1]):
    if m == 0 or n == 0: return ""
    if m == 1:
        find first occurrence of X[0] in Y
        return that character or ""

    mid = m / 2
    left  = computePrefixLCS(X[0..mid-1], Y)   // O(mid * n)
    right = computeSuffixLCS(X[mid..m-1], Y)   // O((m-mid) * n)

    best = argmax_j(left[j] + right[j])

    return Hirschberg(X[0..mid-1], Y[0..best-1])
         + Hirschberg(X[mid..m-1], Y[best..n-1])
```

### 4.3 Rolling Array Optimization

Both prefix and suffix LCS computations use a single 1D array with in-place updates:

```
prev stores dp[i-1][j-1] (the diagonal value from the 2D table)
dp[j]  currently holds dp[i][j]
```

Each update uses `prev + 1` (match) or `max(dp[j], dp[j-1])` (no match), eliminating the need for a full 2D table.

---

## 5. Code Walkthrough

### 5.1 Module Overview

```
Hirschberg-LCS.cpp
│
├── allocIntArray(size)          # Safe memory allocation, exits on failure
├── prefixLengths(X, m, Y, n)   # Computes prefix LCS lengths
├── suffixLengths(X, m, Y, n)   # Computes suffix LCS lengths
├── hirschberg(X, m, Y, n)      # Hirschberg divide-and-conquer
├── getLCS_Hirschberg(X, Y)     # Public interface
└── lcsLength(X, Y)             # Rolling array length computation
```

### 5.2 Key Functions

#### `allocIntArray(int size)`

```c
int* allocIntArray(int size) {
    int *arr = (int*)calloc(size, sizeof(int));
    if (arr == NULL) { exit(1); }
    return arr;
}
```
Uses `calloc` to allocate and zero-initialize. Terminates immediately on failure to prevent null-pointer crashes downstream.

#### `prefixLengths(char *X, int m, char *Y, int n, int *left)`

Performs a forward滚动 update of the `left` array. After completion, `left[j]` holds the LCS length between X[0..m-1] and Y[0..j-1]. The `prev` variable captures the top-left diagonal value from the conceptual 2D table.

#### `suffixLengths(char *X, int m, char *Y, int n, int *right)`

Symmetric to `prefixLengths` but traverses backwards. `right[j]` holds the LCS length between X[0..m-1] and Y[j..n-1], used by Hirschberg to assess the right-half matching potential.

#### `hirschberg(char *X, int m, char *Y, int n, char *lcs, int *pos)`

The divide-and-conquer core. Base cases:
- `m == 0 || n == 0`: empty sequence, return immediately
- `m == 1`: search for X[0] in Y, append if found

Recursive step:
1. `mid = m / 2`, split X into two halves
2. Compute `left[]` and `right[]`
3. Find optimal split `best = argmax(left[j] + right[j])`
4. Recurse on left and right subproblems, append results to `lcs`

#### `lcsLength(char X[], char Y[])`

Standalone length computation using a rolling array. Returns `dp[n]` directly — no direction table needed.

---

## 6. Test Cases & Results

8 built-in test cases covering the following scenarios:

| # | Sequence X | Sequence Y | LCS Length | LCS Sequence | Scenario |
|---|------------|------------|------------|--------------|----------|
| 1 | `ABCBDAB` | `BDCABA` | 4 | `BDAB` | Classic non-contiguous match; multiple optima exist (BCBA also valid) |
| 2 | `AAAAA` | `AAA` | 3 | `AAA` | Repeated characters; optimal selection under order constraint |
| 3 | `ABCDEF` | `GHIJKL` | 0 | `(empty)` | No common characters; verifies empty-result handling |
| 4 | `ABCDEF` | `ABCDEF` | 6 | `ABCDEF` | Identical sequences; verifies identity case |
| 5 | `AGGTAB` | `GXTXAYB` | 4 | `GTAB` | Dispersed character matching |
| 6 | `(empty)` | `ABC` | 0 | `(empty)` | Empty sequence boundary condition |
| 7 | `A` | `BCA` | 1 | `A` | Base case where X has length 1 |
| 8 | `XMJYAUZ` | `MZJAWXU` | 4 | `MJAU` | Complex non-contiguous matching |

### Sample Output

```
========== 测试用例 1 ==========
序列X: ABCBDAB
序列Y: BDCABA
LCS长度: 4
LCS序列: "BDAB"

========== 测试用例 2 ==========
序列X: AAAAA
序列Y: AAA
LCS长度: 3
LCS序列: "AAA"
...
```

---

## 7. Complexity Analysis

### Time Complexity

| Step | Cost | Notes |
|------|------|-------|
| `lcsLength` rolling array | O(mn) | Iterates all character pairs |
| `prefixLengths` | O(mid × n) | Forward pass, left half |
| `suffixLengths` | O((m-mid) × n) | Backward pass, right half |
| `hirschberg` per level | O(n) | Finding optimal split point |
| **Total recursion** | **O(mn)** | Each level totals O(mn); depth O(log m) but work halves each level, sum remains O(mn) |

**Conclusion: Overall time complexity O(mn), identical to traditional DP.**

### Space Complexity

| Data Structure | Size | Notes |
|----------------|------|-------|
| `dp` (rolling array) | O(n) | Used during length computation |
| `left` / `right` | O(n) × 2 | Used during Hirschberg split |
| `lcs` (result) | O(min(m,n)) | Stores the output sequence |
| **Recursion stack** | O(log m) | Binary partition depth |

**Conclusion: Overall space complexity O(n), a significant improvement over the traditional O(mn).**

---

## 8. Comparison with Traditional DP

| Aspect | Traditional DP | This Implementation (Hirschberg) |
|--------|---------------|----------------------------------|
| Time complexity | O(mn) | O(mn) |
| Space complexity | O(mn) | O(n) |
| Sequence reconstruction | Requires backtracking | Direct recursive construction |
| Direction table needed | Yes | No |
| Best suited for | Short sequences | Long sequences, memory-constrained |

---

## 9. Future Improvements

1. **Enumerate all LCS**: Currently returns one optimal solution; can be extended via DFS/backtracking to list all
2. **Input flexibility**: Replace hardcoded sequences with file I/O or interactive input
3. **Large-scale optimization**: Add pruning strategies or parallel computing for very long sequences
4. **Multi-sequence extension**: Generalize from two sequences to multiple sequence alignment

---

## 10. Build & Run

### Environment

- OS: Windows
- Compiler: Any C++23-compliant compiler (g++, clang++, MSVC)
- IDE: Dev-C++, Visual Studio 2022

### Compile

```bash
g++ -std=c++23 Hirschberg-LCS.cpp -o hirschberg-lcs
```

### Run

```bash
./hirschberg-lcs
```

---

## 11. References

1. Hirschberg D S. A linear space algorithm for computing maximal common subsequences[J]. Communications of the ACM, 1975, 18(6): 341-343.
2. 吴东根, 周小安. 基于最长公共子序列的 DNA 序列相似性分析[J]. 智能计算机与应用, 2018, 8(06): 22-26.
3. 郭文月, 刘海砚, 孙群, 等. 利用最长公共子序列度量线要素相似性的方法[J]. 测绘科学技术学报, 2018, 35(05): 518-523.
4. 郑子君, 王洪, 余成. 求解最长循环公共子序列问题的两个算法[J]. 计算机应用研究, 2020, 37(11): 3334-3337.
5. 王防修, 周康. 基于最长公共子序列的随机路径选择算法设计[J]. 计算机工程与设计, 2014, 35(06): 2170-2173.
