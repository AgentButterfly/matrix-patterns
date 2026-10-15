# HB/bcsstk03 — structure and detected patterns

[SuiteSparse Matrix Collection — HB/bcsstk03](https://sparse.tamu.edu/HB/bcsstk03)

MATLAB indices, starting at 1. Unless stated otherwise, indices follow the original file ordering. The identities below were verified on the stored coefficients, without thresholding or approximation.

## General properties

| Property | Result |
|---|---|
| Size | $112\times112$ |
| Type | Real; stiffness matrix of a small test structure |
| Symmetry | Exactly symmetric : $A=A^\top$ |
| Nonzero entries | 640 |
| Nonnegative entries — are all entries nonnegative? | **No**: 340 strictly positive entries and 300 strictly negative entries |
| Zero entries | 11,904 |
| Density — fraction of nonzero entries | $640/112^2=5/98\approx5.102\,\%$ |
| Sparsity — fraction of zero entries | $1-640/112^2=93/98\approx94.898\,\%$ |
| Connected components of the support graph | **2**, with 56 vertices each |
| Half-bandwidth $\max_{A_{ij}\ne0}\lvert i-j\rvert$ | 7 |

Entry counts include both triangles of the matrix. Here, sparsity means the fraction of zero entries; density is also provided to avoid ambiguity.

## Support plot

![Original support, block close-up, and grouped connected components](images/bcsstk03.png)

Each point represents a nonzero entry, regardless of its value or sign. The last panel uses the permutation defined in the connected-components section.

## 1. Exact block-tridiagonal structure with blocks of size 4

Define, for $i=1,\ldots,28$,

$$
I_i=\{4i-3,4i-2,4i-1,4i\},\qquad D_i=A(I_i,I_i),
$$

and, for $i=1,\ldots,27$, $T_i=A(I_i,I_{i+1})$. Then

$$
A=\begin{pmatrix}
D_1&T_1&0&\cdots&0\\
T_1^\top&D_2&T_2&\ddots&\vdots\\
0&T_2^\top&D_3&\ddots&0\\
\vdots&\ddots&\ddots&\ddots&T_{27}\\
0&\cdots&0&T_{27}^\top&D_{28}
\end{pmatrix}.
$$

This decomposition is exact and requires no permutation. Each block couples only to itself and its immediate neighbours.

## 2. Common form of the 27 upper off-diagonal blocks

All upper off-diagonal blocks follow the same exact three-parameter template

$$
T_i=\begin{pmatrix}
a_i&0&0&b_i\\
0&a_i&-b_i&0\\
0&b_i&c_i&0\\
-b_i&0&0&c_i
\end{pmatrix}.
$$

The parameters are defined by

$$
a_i=A_{4i-3,4i+1},\qquad b_i=A_{4i-3,4i+4},\qquad c_i=A_{4i-1,4i+3}.
$$

- Local entries $(1,1)$ and $(2,2)$ are equal.
- Entries $(3,3)$ and $(4,4)$ are equal.
- Entries $(1,4)$ and $(3,2)$ equal $b_i$; $(2,3)$ and $(4,1)$ equal $-b_i$.
- All other entries are zero.
- Each $T_i$ has exactly 8 nonzero entries; all 27 blocks therefore have the same support.
- The lower off-diagonal blocks are the transposes $T_i^\top$; the $T_i$ themselves are not symmetric.

## 3. Diagonal-block template and exceptions

For all indices $i=1,\ldots,28$,

$$
D_i=\begin{pmatrix}
u_i&0&0&v_i\\
0&u_i&-v_i&0\\
0&-v_i&w_i&0\\
v_i&0&0&w_i
\end{pmatrix}+\delta_i E_{22},
$$

where $E_{22}$ is the $4\times4$ matrix with a 1 at position $(2,2)$ and zeros elsewhere, and

$$
u_i=A_{4i-3,4i-3},\quad v_i=A_{4i-3,4i},\quad
w_i=A_{4i-1,4i-1},\quad \delta_i=A_{4i-2,4i-2}-u_i.
$$

**25 of the 28 blocks** follow the template without correction. The only exceptions are

$$
\delta_i=0\quad\text{for }i\notin\{21,25,27\}.
$$

| Exceptional block | Modified local entry | Global entry |
|---|---|---|
| $D_{21}$ | $(2,2)$ | $(82,82)$ |
| $D_{25}$ | $(2,2)$ | $(98,98)$ |
| $D_{27}$ | $(2,2)$ | $(106,106)$ |

Within each block, entries $(3,3)$ and $(4,4)$ remain equal. Entries $(1,4)$ and $(4,1)$ equal $v_i$, while $(2,3)$ and $(3,2)$ equal $-v_i$.

### Two diagonal-block support patterns

$$
v_i=0\quad\Longleftrightarrow\quad i\in\{4,5,6,18\}.
$$

These four blocks are diagonal and have 4 nonzero entries. The other 24 blocks have 8 nonzero entries and share the same support. The three corrections $\delta_i$ do not change the support.

This recovers the total number of nonzero entries:

$$
\operatorname{nnz}(A)=\underbrace{4\times4+24\times8}_{208\text{ in the }D_i}
+2\underbrace{(27\times8)}_{216\text{ in the }T_i}=640.
$$

## 4. Exact repetitions between blocks

| Family | Exactly identical blocks |
|---|---|
| Diagonal blocks | $D_4=D_5=D_6$ |
| Upper off-diagonal blocks | $T_3=T_4=T_5=T_6=T_7$ |
| Upper off-diagonal blocks | $T_{17}=T_{18}$ |

These are the only exact repetition groups within each family: there are **26 distinct diagonal blocks** and **22 distinct upper off-diagonal blocks**. A shared symbolic template therefore does not imply identical block values.

## 5. Scalar-diagonal description

The only occupied offsets are $0,\pm1,\pm3,\pm4,\pm5,\pm7$. The formulas below completely describe the upper triangle; the lower triangle follows by symmetry. Unlisted entries are zero.

| Offset $k$ | Position $(j,j+k)$ | Value | Index range |
|---|---|---|---|
| 0 | $(4i-3,4i-3)$ | $u_i$ | $1\le i\le28$ |
| 0 | $(4i-2,4i-2)$ | $u_i+\delta_i$ | $1\le i\le28$ |
| 0 | $(4i-1,4i-1)$ and $(4i,4i)$ | $w_i$ | $1\le i\le28$ |
| 1 | $(4i-2,4i-1)$ | $-v_i$ | $1\le i\le28$ |
| 1 | $(4i,4i+1)$ | $-b_i$ | $1\le i\le27$ |
| 3 | $(4i-3,4i)$ | $v_i$ | $1\le i\le28$ |
| 3 | $(4i-1,4i+2)$ | $b_i$ | $1\le i\le27$ |
| 4 | $(4i-3,4i+1)$ and $(4i-2,4i+2)$ | $a_i$ | $1\le i\le27$ |
| 4 | $(4i-1,4i+3)$ and $(4i,4i+4)$ | $c_i$ | $1\le i\le27$ |
| 5 | $(4i-2,4i+3)$ | $-b_i$ | $1\le i\le27$ |
| 7 | $(4i-3,4i+4)$ | $b_i$ | $1\le i\le27$ |

The zeros of $v_i$ at indices $4,5,6,18$ create gaps along the diagonals with offsets 1 and 3.

## 6. Two completely decoupled connected components

The two components are

$$
G_1=\bigcup_{i=1}^{28}\{4i-3,4i\},\qquad
G_2=\bigcup_{i=1}^{28}\{4i-2,4i-1\}.
$$

Thus $G_1=(1,4,5,8,\ldots,109,112)$ and $G_2=(2,3,6,7,\ldots,110,111)$, in the indicated order. No nonzero entry connects the two groups. Self-loops corresponding to diagonal entries do not affect the connected components.

Setting $p=(G_1,G_2)$ and $P=[e_{p_1},\ldots,e_{p_{112}}]$,

$$
P^\top AP=\begin{pmatrix}B&0\\0&C\end{pmatrix},\qquad B,C\in\mathbb R^{56\times56}.
$$

### Internal structure and relationship between the two subsystems

In the ordering above, $B$ and $C$ are block-tridiagonal with $2\times2$ blocks. Their diagonal and upper off-diagonal blocks are

$$
B_i=\begin{pmatrix}u_i&v_i\\v_i&w_i\end{pmatrix},\qquad
C_i=\begin{pmatrix}u_i+\delta_i&-v_i\\-v_i&w_i\end{pmatrix},
$$

$$
U_i=\begin{pmatrix}a_i&b_i\\-b_i&c_i\end{pmatrix},\qquad
V_i=\begin{pmatrix}a_i&-b_i\\b_i&c_i\end{pmatrix}.
$$

With $J=\operatorname{diag}(1,-1)$ and $Q=I_{28}\otimes J$, we obtain

$$
C=QBQ+\Delta,\qquad
\Delta=\sum_{i\in\{21,25,27\}}\delta_i e_{2i-1}e_{2i-1}^\top.
$$

The two subsystems are therefore related by a sign change and three diagonal corrections, at positions 41, 49 and 53 of $C$.

## Sources and verification

- [SuiteSparse entry](https://sparse.tamu.edu/HB/bcsstk03).
- [Matrix Market data, NIST mirror](https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcsstruc1/bcsstk03.mtx.gz).
- Block templates, repetitions, signs and connected components were verified on the Matrix Market file on 15 October 2026. No threshold was applied.
