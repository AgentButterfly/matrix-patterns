# HB/bcsstk01 — structure and detected patterns

[SuiteSparse Matrix Collection — HB/bcsstk01](https://sparse.tamu.edu/HB/bcsstk01)

All indices start at 1. Unless explicitly stated otherwise, the original matrix ordering is used. A block support is the set of its nonzero positions. Exact identities below refer to the stored coefficients; approximate scaling relations are labelled separately.

## General properties

| Property | Result |
|---|---|
| Size | $48\times48$ |
| Type | Real stiffness matrix |
| Symmetry | Exactly symmetric: $A=A^\top$ |
| Nonzero entries | 400, counting both triangles |
| Nonnegative entries? | **No:** 200 strictly positive and 200 strictly negative entries |
| Row sparsity | Maximum of **12** nonzero entries per row |
| Column sparsity | Maximum of **12** nonzero entries per column |
| $d$-sparsity | **$d=12$** |
| Connected components of the support graph | **1**, containing all 48 indices |

Here sparsity means

$$
d=\max\left\{\max_i\#\{j:A_{ij}\ne0\},\ \max_j\#\{i:A_{ij}\ne0\}\right\}.
$$

The original matrix cannot be split by a simultaneous row/column permutation into nonempty independent diagonal subsystems. The two chains described below arise by removing specified couplings, not by finding connected components of the original matrix.

## Support plots

![Original support, block occupancy, and the two-chain decomposition](images/bcsstk01.png)

The first panel shows the original support with a $6\times6$ block grid. The second shows which blocks are nonzero. In the third, blue entries belong to the retained two-chain matrix, while orange entries are the removed couplings.

## 1. Partition into eight blocks of size 6

Define $K_{ij}\in\mathbb R^{6\times6}$ entrywise by

$$
(K_{ij})_{rs}=A_{6(i-1)+r,\,6(j-1)+s},
\qquad 1\le i,j\le8,\quad1\le r,s\le6.
$$

Symmetry implies $K_{ji}=K_{ij}^\top$.

All eight diagonal blocks are nonzero. The only strictly upper nonzero blocks are

$$
K_{12},K_{14},K_{15},K_{23},K_{26},K_{34},
K_{37},K_{38},K_{48},K_{56},K_{67},K_{78}.
$$

Together with their transposes, this gives **32 nonzero blocks out of 64**. Every other block is identically zero.

## 2. Exact block repetitions and transpose relationships

The following identities hold exactly:

$$
K_{15}=K_{26}=K_{37}=K_{48}=C,
\qquad K_{14}=K_{23}=Y,
$$

$$
K_{34}=K_{12}^\top,\qquad K_{78}=K_{56}^\top.
$$

The four copies of $C$ follow the fixed placement rule $j=i+4$, for $i=1,\ldots,4$.

The transpose relationships compare different blocks within the upper triangle; they are additional structure beyond the global identity $K_{ji}=K_{ij}^\top$.

Set

$$
D_i=K_{ii},\quad X=K_{12},\quad Y=K_{23},\quad
\widehat X=K_{56},\quad\widehat Y=K_{67},\quad F=K_{38}.
$$

The full matrix is then

$$
A=\begin{pmatrix}
D_1&X&0&Y&C&0&0&0\\
X^\top&D_2&Y&0&0&C&0&0\\
0&Y^\top&D_3&X^\top&0&0&C&F\\
Y^\top&0&X&D_4&0&0&0&C\\
C^\top&0&0&0&D_5&\widehat X&0&0\\
0&C^\top&0&0&\widehat X^\top&D_6&\widehat Y&0\\
0&0&C^\top&0&0&\widehat Y^\top&D_7&\widehat X^\top\\
0&0&F^\top&C^\top&0&0&\widehat X&D_8
\end{pmatrix}.
$$

The eight diagonal blocks are numerically distinct.

## 3. Templates for the ordinary off-diagonal blocks

Write $E_{rs}^{(n)}$ for the $n\times n$ matrix with a single 1 at $(r,s)$. All remaining entries are zero. A shared template specifies arrangement and sign relationships, not equality of all parameter values between blocks.

### 3.1 Template for $X$ and $\widehat X$

Both have the form

$$
\mathcal X(a,b,c,d,e,f,p,q)=
\begin{pmatrix}
a&0&0&0&p&0\\
0&b&0&q&0&0\\
0&0&c&0&0&0\\
0&-q&0&d&0&0\\
-p&0&0&0&e&0\\
0&0&0&0&0&f
\end{pmatrix}.
$$

Equivalently,

$$
\mathcal X=\operatorname{diag}(a,b,c,d,e,f)
+p(E_{15}^{(6)}-E_{51}^{(6)})
+q(E_{24}^{(6)}-E_{42}^{(6)}).
$$

Thus the off-diagonal part is skew-symmetric. Transposition preserves the diagonal and changes the signs of $p,q$. This accounts for $K_{34}=X^\top$ and $K_{78}=\widehat X^\top$.

### 3.2 Template for $Y$ and $\widehat Y$

Both have the form

$$
\mathcal Y(a,b,c,d,e,f,p,q)=
\begin{pmatrix}
a&0&0&0&0&0\\
0&b&0&0&0&p\\
0&0&c&0&q&0\\
0&0&0&d&0&0\\
0&0&-q&0&e&0\\
0&-p&0&0&0&f
\end{pmatrix},
$$

or

$$
\mathcal Y=\operatorname{diag}(a,b,c,d,e,f)
+p(E_{26}^{(6)}-E_{62}^{(6)})
+q(E_{35}^{(6)}-E_{53}^{(6)}).
$$

The symbols in $\mathcal X$ and $\mathcal Y$ are independent parameters. The same construction is reused with different coupled index pairs.

### 3.3 Template for the repeated block $C$

The four identical blocks $C$ have the form

$$
C=\begin{pmatrix}
\alpha&0&0&0&0&p\\
0&\beta&0&0&0&0\\
0&0&\alpha&-p&0&0\\
0&0&p&\gamma&0&0\\
0&0&0&0&\eta&0\\
-p&0&0&0&0&\gamma
\end{pmatrix}.
$$

This is a third support template: the coupled pairs are $\{1,6\}$ and $\{3,4\}$. Its repeated diagonal entries and the signs of the four off-diagonal entries are part of the template.

### 3.4 Support-template inventory

| Family in the strictly upper block triangle | Blocks | Nonzeros per block |
|---|---|---:|
| $\mathcal X$ support | $K_{12},K_{34},K_{56},K_{78}$ | 10 |
| $\mathcal Y$ support | $K_{14},K_{23},K_{67}$ | 10 |
| $C$ support | $K_{15},K_{26},K_{37},K_{48}$ | 10 |
| Exceptional support | $K_{38}=F$ | 18 |

There are **three support templates among the eleven 10-entry blocks**, plus one exceptional block. Transposition does not change these three supports, although it changes coefficients.

## 4. Common template for the diagonal blocks

Each $D_i=K_{ii}$ has the exact form

$$
D_i=\begin{pmatrix}
d_{i,1}&0&0&0&p_i&q_i\\
0&d_{i,2}&h_i&r_i&0&s_i\\
0&h_i&d_{i,3}&t_i&u_i&0\\
0&r_i&t_i&d_{i,4}&0&0\\
p_i&0&u_i&0&d_{i,5}&\ell_i\\
q_i&s_i&0&0&\ell_i&d_{i,6}
\end{pmatrix}.
$$

The parameters are the entries at the displayed positions; parameters from different blocks are not implicitly equal.

The exceptional couplings obey

$$
h_i=\ell_i=0\quad\text{for }i\notin\{3,8\},
$$

and the missing couplings obey

$$
s_5=u_5=s_8=u_8=0.
$$

| Blocks | Support description | Nonzeros per block |
|---|---|---:|
| $D_1,D_2,D_4,D_6,D_7$ | Common support with $h_i=\ell_i=0$ | 18 |
| $D_3$ | Common support plus symmetric pairs $(2,3)$ and $(5,6)$ | 22 |
| $D_5$ | Common support without symmetric pairs $(2,6)$ and $(3,5)$ | 14 |
| $D_8$ | Missing pairs as in $D_5$, plus the exceptional pairs of $D_3$ | 18 |

Each pair includes its transposed position. The diagonal blocks contain 144 nonzeros altogether, giving

$$
\operatorname{nnz}(A)=144+2(11\times10+18)=400.
$$

## 5. Compact form of the upper-right region

Split the matrix after index 24:

$$
A=\begin{pmatrix}L&U\\U^\top&R\end{pmatrix},\qquad L,U,R\in\mathbb R^{24\times24}.
$$

The complete upper-right region is

$$
U=\begin{pmatrix}
K_{15}&K_{16}&K_{17}&K_{18}\\
K_{25}&K_{26}&K_{27}&K_{28}\\
K_{35}&K_{36}&K_{37}&K_{38}\\
K_{45}&K_{46}&K_{47}&K_{48}
\end{pmatrix}
=\begin{pmatrix}
C&0&0&0\\0&C&0&0\\0&0&C&F\\0&0&0&C
\end{pmatrix}.
$$

Therefore

$$
\boxed{U=I_4\otimes C+E_{34}^{(4)}\otimes F.}
$$

All unlisted blocks in this region are zero. Global block $K_{38}$ occupies local position $(3,4)$ within $U$.

## 6. Two block-tridiagonal chains plus structured couplings

Set $U,U^\top,K_{14},K_{41}$ to zero, retaining all other entries. The resulting matrix is

$$
\widetilde A=\operatorname{diag}(B,D),
$$

where

$$
B=\begin{pmatrix}
D_1&X&0&0\\
X^\top&D_2&Y&0\\
0&Y^\top&D_3&X^\top\\
0&0&X&D_4
\end{pmatrix},\qquad
D=\begin{pmatrix}
D_5&\widehat X&0&0\\
\widehat X^\top&D_6&\widehat Y&0\\
0&\widehat Y^\top&D_7&\widehat X^\top\\
0&0&\widehat X&D_8
\end{pmatrix}.
$$

Here $B,D$ are each $24\times24$; $D_i$ denotes a $6\times6$ diagonal block. They are two independent four-block chains in $\widetilde A$, with a zero coupling between global blocks 4 and 5.

Writing

$$
H=E_{14}^{(4)}\otimes Y+E_{41}^{(4)}\otimes Y^\top,
$$

the original matrix has the exact decomposition

$$
\boxed{A=\widetilde A+\mathcal R,\qquad
\mathcal R=\begin{pmatrix}H&U\\U^\top&0\end{pmatrix}.}
$$

The symbol $\mathcal R$ denotes the removed couplings, whereas $R$ above denotes the original lower-right region. This is an exact additive decomposition, not a claim that discarding $\mathcal R$ gives a sufficiently accurate approximation.

## 7. Relationships between the two chains

### 7.1 The same arrangement of connecting blocks

Their upper connecting blocks follow the same template:

$$
(X,Y,X^\top)\qquad\text{and}\qquad
(\widehat X,\widehat Y,\widehat X^\top).
$$

Corresponding connecting blocks have exactly the same supports. Therefore, outside their diagonal $6\times6$ blocks, the two chain supports coincide exactly.

### 7.2 Localized differences between diagonal-block supports

All coordinates below are local to a $6\times6$ block. The table lists the symmetric difference of the supports, not differences of coefficient values.

| Corresponding blocks | Differing support positions |
|---|---|
| $D_1,D_5$ | $(2,6),(6,2),(3,5),(5,3)$ |
| $D_2,D_6$ | None |
| $D_3,D_7$ | $(2,3),(3,2),(5,6),(6,5)$ |
| $D_4,D_8$ | All eight positions listed above |

The full $24\times24$ supports differ at only **16 positions**.

### 7.3 Approximate scaling relationships with single-entry corrections

The coefficient patterns additionally satisfy

$$
\widehat X=\frac54X-255\,000E_{33}^{(6)}+\varepsilon_X,
$$

$$
\widehat Y=\frac54Y-425\,000E_{11}^{(6)}+\varepsilon_Y.
$$

These are **not exact scaling identities** for the stored data. The maximum absolute residual entries are approximately

$$
\max_{r,s}|(\varepsilon_X)_{rs}|=7.5\times10^{-9},\qquad
\max_{r,s}|(\varepsilon_Y)_{rs}|=7.5\times10^{-4}.
$$

This records a common factor $5/4$ with one distinguished diagonal correction per connecting template, up to the stated small residuals. It does not imply that the entire second chain is $5/4$ times the first.

## Sources and verification

- [SuiteSparse Matrix Collection](https://sparse.tamu.edu/HB/bcsstk01).
- [Matrix Market data, NIST mirror](https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcsstruc1/bcsstk01.mtx.gz).
- Compiled and checked against the stored Matrix Market coefficients on 19 October 2026. Support tests use literal nonzeros, without thresholding. Approximate relationships are isolated in Section 7.3.
