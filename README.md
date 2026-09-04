# 最长公共子序列（LCS）求解与序列还原

基于动态规划与 Hirschberg 分治算法

## 一、项目简介

**最长公共子序列（Longest Common Subsequence，LCS）** 是动态规划领域的经典问题：给定两个序列，找出它们共有的、长度最长的子序列（元素不必连续，但必须保持相对顺序）。

本项目采用 **动态规划 + Hirschberg 分治算法**，在求出 LCS 长度的同时，完整还原出一条具体的 LCS 序列。相比传统 DP 需要保存完整的二维状态表（空间 O(mn)），本实现通过滚动数组和分治策略将空间复杂度降至 **O(n)**，适合处理较长序列的场景。

---

## 二、问题分析

### 2.1 问题背景

LCS 问题在多个领域有重要应用：

| 领域 | 应用场景 |
|------|---------|
| **生物信息学** | DNA 序列比对，判断基因相似度 |
| **版本控制** | `git diff` 计算文件差异的核心算法之一 |
| **文本处理** | 计算两篇文章的相似度、拼写检查 |
| **地理信息系统** | 线要素相似性度量 |

### 2.2 暴力枚举 vs 动态规划

对于长度为 m 和 n 的两个序列，暴力枚举需要检查所有子集组合，时间复杂度为 **O(2^min(m,n))**，随着序列增长急剧恶化。

LCS 问题具有两个关键性质，使其适合动态规划：

1. **最优子结构**：两个序列的 LCS 的子序列，必然是其前缀的 LCS
2. **重叠子问题**：不同前缀对的 LCS 计算存在大量重复

---

## 三、数学模型

### 3.1 状态定义

设序列 X = ⟨x₁, x₂, …, xₘ⟩，Y = ⟨y₁, y₂, …, yₙ⟩

定义 `dp[i][j]` 为 **X 的前 i 个元素**与**Y 的前 j 个元素**的最长公共子序列长度。

### 3.2 状态转移方程

```
                    ⎧ 0                                    if i = 0 or j = 0
dp[i][j] =          ⎨ dp[i-1][j-1] + 1                     if xᵢ = yⱼ
                    ⎩ max(dp[i-1][j], dp[i][j-1])          if xᵢ ≠ yⱼ
```

### 3.3 最优子结构证明

设 Z = ⟨z₁, z₂, …, zₖ⟩ 是 X 和 Y 的任意一个 LCS：

1. **若 xₘ = yₙ**：则 zₖ = xₘ = yₙ，且 Z 的前 k-1 项是 Xₘ₋₁ 和 Yₙ₋₁ 的 LCS
2. **若 xₘ ≠ yₙ 且 zₖ ≠ xₘ**：则 Z 是 Xₘ₋₁ 和 Y 的 LCS
3. **若 xₘ ≠ yₙ 且 zₖ ≠ yₙ**：则 Z 是 X 和 Yₙ₋₁ 的 LCS

---

## 四、算法设计

### 4.1 传统 DP 的局限

传统方法需要创建一个 `(m+1) × (n+1)` 的二维 DP 表，并用方向数组记录每个状态的前驱路径以回溯还原序列。当序列较长时，内存消耗巨大。

**改进方向**：
- 求 LCS **长度**：可用滚动数组，空间从 O(mn) → O(n)
- 还原 LCS **序列**：Hirschberg 分治算法，空间同样控制在 O(n)

### 4.2 Hirschberg 分治算法

#### 核心思想

将序列 X 从中间分成左右两部分 X_left 和 X_right，分别计算：
- `left[j]`：X_left 与 Y[0..j-1] 的 LCS 长度（前缀 DP）
- `right[j]`：X_right 与 Y[j..n-1] 的 LCS 长度（后缀 DP）

寻找使 `left[j] + right[j]` 最大的分割点 `best`，将 Y 也从该点分割，递归处理左右两个子问题。

#### 算法正确性

对于任意分割点 j，`left[j] + right[j]` 表示「X 左半部分与 Y 前 j 个字符的 LCS」加上「X 右半部分与 Y 后 n-j 个字符的 LCS」的总长度。取最大值即找到最优的 Y 分割位置，保证递归结果拼接后仍是全局最优解。

#### 伪代码

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

### 4.3 滚动数组优化

前缀/后缀 LCS 长度计算均使用一维数组 + 滚动更新：

```
prev 记录 dp[i-1][j-1]（左上角的旧值）
dp[j] 当前记录 dp[i][j]
```

每次更新 `dp[j]` 时用 `prev + 1`（匹配）或 `max(dp[j], dp[j-1])`（不匹配），避免分配完整的二维表。

---

## 五、代码解析

### 5.1 模块总览

```
Hirschberg-LCS.cpp
│
├── allocIntArray(size)          # 安全内存分配，失败时直接终止
├── prefixLengths(X, m, Y, n)   # 计算 X 前缀与 Y 各前缀的 LCS 长度
├── suffixLengths(X, m, Y, n)   # 计算 X 后缀与 Y 各后缀的 LCS 长度
├── hirschberg(X, m, Y, n)      # Hirschberg 分治，还原 LCS 序列
├── getLCS_Hirschberg(X, Y)     # 对外接口，封装调用
└── lcsLength(X, Y)             # 滚动数组求 LCS 长度
```

### 5.2 关键函数详解

#### `allocIntArray(int size)`

```c
int* allocIntArray(int size) {
    int *arr = (int*)calloc(size, sizeof(int));
    if (arr == NULL) { exit(1); }
    return arr;
}
```
用 `calloc` 分配并初始化为 0，防止未初始化内存导致错误结果；分配失败直接退出，避免后续空指针崩溃。

#### `prefixLengths(char *X, int m, char *Y, int n, int *left)`

按行方向滚动更新 `left` 数组。`left[j]` 最终表示 `X[0..m-1]` 与 `Y[0..j-1]` 的 LCS 长度。内层循环中 `prev` 变量保存上一行、前一列的值（即二维 DP 中的 `dp[i-1][j-1]`）。

#### `suffixLengths(char *X, int m, char *Y, int n, int *right)`

与 `prefixLengths` 对称，但遍历方向相反（从末尾向前）。`right[j]` 表示 `X[0..m-1]` 与 `Y[j..n-1]` 的 LCS 长度，用于 Hirschberg 分治时计算右半部分的匹配能力。

#### `hirschberg(char *X, int m, char *Y, int n, char *lcs, int *pos)`

分治核心。递归基线：
- `m == 0 || n == 0`：空序列，直接返回
- `m == 1`：在 Y 中查找 X[0]，找到则写入结果

递归步骤：
1. 取 `mid = m / 2`，将 X 分为左右两半
2. 计算 `left[]` 和 `right[]`
3. 找最优分割点 `best = argmax(left[j] + right[j])`
4. 递归求解左右子问题，结果依次写入 `lcs`

#### `lcsLength(char X[], char Y[])`

独立的长度计算函数，使用滚动数组，无需额外方向表，直接返回 `dp[n]` 即为 LCS 长度。

---

## 六、测试用例与结果

程序内置 8 组测试用例，覆盖以下场景：

| # | 序列 X | 序列 Y | LCS 长度 | LCS 序列 | 测试场景 |
|---|--------|--------|----------|----------|----------|
| 1 | `ABCBDAB` | `BDCABA` | 4 | `BDAB` | 经典非连续匹配，存在多个最优解（BCBA 亦可） |
| 2 | `AAAAA` | `AAA` | 3 | `AAA` | 重复字符，顺序约束下的最优选择 |
| 3 | `ABCDEF` | `GHIJKL` | 0 | `(空)` | 无公共字符，验证空结果处理 |
| 4 | `ABCDEF` | `ABCDEF` | 6 | `ABCDEF` | 完全相同序列，验证恒等情形 |
| 5 | `AGGTAB` | `GXTXAYB` | 4 | `GTAB` | 分散字符匹配 |
| 6 | `(空)` | `ABC` | 0 | `(空)` | 空序列边界情况 |
| 7 | `A` | `BCA` | 1 | `A` | X 长度为 1 的递归基线 |
| 8 | `XMJYAUZ` | `MZJAWXU` | 4 | `MJAU` | 复杂非连续匹配 |

### 运行输出示例

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

## 七、复杂度分析

### 时间复杂度

| 步骤 | 计算量 | 说明 |
|------|--------|------|
| `lcsLength` 滚动数组 | O(mn) | 遍历所有字符对 |
| `prefixLengths` | O(mid × n) | 前半部分前缀计算 |
| `suffixLengths` | O((m-mid) × n) | 后半部分后缀计算 |
| `hirschberg` 单层 | O(n) | 寻找最优分割点 |
| **整体递归** | **O(mn)** | 每层合计 O(mn)，递归深度 O(log m)，但各层计算量逐层减半，总和仍为 O(mn) |

**结论：整体时间复杂度 O(mn)，与传统 DP 相同。**

### 空间复杂度

| 数据结构 | 大小 | 说明 |
|----------|------|------|
| `dp`（滚动数组） | O(n) | 计算长度时使用 |
| `left` / `right` | O(n) × 2 | Hirschberg 分治时使用 |
| `lcs`（结果） | O(min(m,n)) | 存储 LCS 序列 |
| **递归栈** | O(log m) | 二分递归深度 |

**结论：整体空间复杂度 O(n)，相比传统 DP 的 O(mn) 显著优化。**

---

## 八、与传统 DP 对比

| 特性 | 传统 DP | 本实现（Hirschberg） |
|------|---------|---------------------|
| 时间复杂度 | O(mn) | O(mn) |
| 空间复杂度 | O(mn) | O(n) |
| 能否还原序列 | 需要额外回溯 | 直接递归构建 |
| 是否需要方向表 | 是 | 否 |
| 适合场景 | 序列较短 | 序列较长，内存受限 |

---

## 九、改进方向

1. **输出全部 LCS**：当前只返回一条，可通过回溯枚举所有最优路径输出全部解
2. **输入方式扩展**：目前序列硬编码，可改为文件读取或交互式输入
3. **大规模数据适配**：引入剪枝策略或并行计算，进一步优化超长序列的比对效率
4. **多序列扩展**：当前为双序列，可扩展为多序列比对

---

## 十、编译与运行

### 环境要求

- 操作系统：Windows
- 编译器：支持 ISO C++23 标准的 g++ / clang++ / MSVC
- 开发工具：Dev-C++、Visual Studio 2022

### 编译命令

```bash
g++ -std=c++23 Hirschberg-LCS.cpp -o hirschberg-lcs
```

### 运行

```bash
./hirschberg-lcs
```

---

## 十一、参考文献

1. 吴东根, 周小安. 基于最长公共子序列的 DNA 序列相似性分析[J]. 智能计算机与应用, 2018, 8(06): 22-26.
2. 郭文月, 刘海砚, 孙群, 等. 利用最长公共子序列度量线要素相似性的方法[J]. 测绘科学技术学报, 2018, 35(05): 518-523.
3. 郑子君, 王洪, 余成. 求解最长循环公共子序列问题的两个算法[J]. 计算机应用研究, 2020, 37(11): 3334-3337.
4. 王防修, 周康. 基于最长公共子序列的随机路径选择算法设计[J]. 计算机工程与设计, 2014, 35(06): 2170-2173.
5. Hirschberg D S. A linear space algorithm for computing maximal common subsequences[J]. Communications of the ACM, 1975, 18(6): 341-343.
