# HB/bcsstk01 — structure and detected patterns

[SuiteSparse Matrix Collection — HB/bcsstk01](https://sparse.tamu.edu/HB/bcsstk01)

MATLAB indices, starting at 1. Unless stated otherwise, indices follow the original file ordering. The studied matrix will be denoted A.

## General properties

| Property | Result |
|---|---|
| Size | $48\times48$ |
| Type | Real stiffness matrix |
| Symmetry | symmetric |
| Nonzero entries | 400 |
| Sparsity coefficient | 12 |
| Connected components of the support graph | **1**, containing all 48 indices |

## Support plots

![Original support, block occupancy, and the two-chain decomposition](images/bcsstk01.png)

The first panel shows the original support with a $6\times6$ block grid. The second shows which blocks are nonzero. In the third, blue entries belong to the retained two-chain matrix, while orange entries are the removed couplings.

## 1. Partition into eight blocks of size 6

Define $K_{ij}\in\mathbb R^{6\times6}$ by

$$
(K_{ij})_{rs}=A_{6(i-1)+r, 6(j-1)+s},
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


## 5. Compact form of the upper-right region

Split the matrix after index 24:

$$
A=\begin{pmatrix}L&U\\U^\top&R\end{pmatrix},\qquad L,U,R\in\mathbb R^{24\times24}.
$$

The complete upper-right region is

$$
U=\begin{pmatrix}
C&0&0&0\\0&C&0&0\\0&0&C&F\\0&0&0&C
\end{pmatrix}.
$$

Therefore

$$
\boxed{U=I_4\otimes C+E_{34}^{(4)}\otimes F.}
$$

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
