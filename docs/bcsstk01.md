# HB/bcsstk01 

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
(K_{ij})_{rs}=A_{6(i-1)+r, 6(j-1)+s}
\qquad 1\le i,j\le8,\quad1\le r,s\le6.
$$

Symmetry implies $K_{ji}=K_{ij}^\top$.

All eight diagonal blocks are nonzero. The only strictly upper nonzero blocks are

$$
K_{12},K_{14},K_{15},K_{23},K_{26},K_{34},
K_{37},K_{38},K_{48},K_{56},K_{67},K_{78}.
$$

Together with their transposes, this gives **32 nonzero blocks out of 64**.

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

Write $E_{rs}^{(n)}$ for the $n\times n$ matrix with a single 1 at $(r,s)$. All remaining entries are zero. 

### 3.1 Template for $X$ and $\widehat X$

Both have the form

$$
\mathcal X=\mathrm{diag}(a,b,c,d,e,f)
+p(E_{15}^{(6)}-E_{51}^{(6)})
+q(E_{24}^{(6)}-E_{42}^{(6)}).
$$

Thus the off-diagonal part is skew-symmetric. 

### 3.2 Template for $Y$ and $\widehat Y$

Both have the form

$$
\mathcal Y=\mathrm{diag}(a,b,c,d,e,f)
+p(E_{26}^{(6)}-E_{62}^{(6)})
+q(E_{35}^{(6)}-E_{53}^{(6)}).
$$

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
A=\begin{pmatrix}
L&U\\
U^\top&R
\end{pmatrix},\qquad L,U,R\in\mathbb R^{24\times24}.
$$

The complete upper-right region is

$$
U=\begin{pmatrix}
C&0&0&0\\
0&C&0&0\\
0&0&C&F\\
0&0&0&C
\end{pmatrix} =I_4\otimes C+E_{34}^{(4)}\otimes F
$$

## 6. Two block-tridiagonal chains plus structured couplings

Set $U,U^\top,K_{14},K_{41}$ to zero, retaining all other entries. The resulting matrix is

$$
\widetilde A=\mathrm{diag}(B,D),
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

Here $B,D$ are each $24\times24$ and block tridiagonal. Their upper connecting blocks follow the same template:

$$
(X,Y,X^\top)\qquad\text{and}\qquad
(\widehat X,\widehat Y,\widehat X^\top).
$$

Writing

$$
H=E_{14}^{(4)}\otimes Y+E_{41}^{(4)}\otimes Y^\top,
$$

the original matrix has the exact decomposition

$$
\boxed{A=\widetilde A+\mathcal R,\qquad
\mathcal R=\begin{pmatrix}H&U\\ U^\top&0 \end{pmatrix}.}
$$

The symbol $\mathcal R$ denotes the removed couplings, whereas $R$ above denotes the original lower-right region. This is an exact additive decomposition, not a claim that discarding $\mathcal R$ gives a sufficiently accurate approximation.

## 8. Row representations derived from the block structure

For a matrix $M$, let $R_M(i)$ denote the list of pairs $(j,M_{ij})$
corresponding to the nonzero entries of row $i$.

We use $\Vert$ for list concatenation and define the column-shift operation

```math
\mathcal S_t(L)=\bigl[(s+t,v):(s,v)\in L\bigr].
```

### 8.1 Additive decomposition of A

Using the two block-tridiagonal chains $B,D$ defined above,

```math
A=
\mathrm{diag}(B,D)
+E_{14}^{(8)}\otimes Y
+E_{41}^{(8)}\otimes Y^\top
+
\begin{pmatrix}
0 & U \\
U^\top & 0
\end{pmatrix},
```

where

```math
U=I_4\otimes C+E_{34}^{(4)}\otimes F.
```

Here $E_{ab}^{(n)}$ denotes the $n\times n$ matrix with a single
1 at position $(a,b)$ and zeros elsewhere.

The chains $B,D$ have size $24\times24$, and the blocks
$Y,C,F$ have size $6\times6$. In particular, $B$ does not
include the couplings $K_{14}=Y$ and $K_{41}=Y^\top$:
these are supplied by the two explicit Kronecker terms.

The four terms in this decomposition have pairwise disjoint
supports. Their row lists can therefore be concatenated
without combining contributions at the same column.

### 8.2 First representation: composition of the larger regions

For the first half of the matrix, start with the corresponding
row of $B$, add the row of $U$ with a column shift of 24, and
include the additional $Y$ or $Y^\top$ coupling where applicable.

For the second half, combine the corresponding row of $U^\top$
with the shifted row of $D$.

```math
R_A^{(1)}(i)=
\begin{cases}
R_B(i)
\Vert \mathcal S_{24}\!\left(R_U(i)\right)
\Vert \mathcal S_{18}\!\left(R_Y(i)\right),
& 1\le i\le6, \\[2mm]
R_B(i)
\Vert \mathcal S_{24}\!\left(R_U(i)\right),
& 7\le i\le18, \\[2mm]
R_B(i)
\Vert \mathcal S_{24}\!\left(R_U(i)\right)
\Vert R_{Y^\top}(i-18),
& 19\le i\le24, \\[2mm]
R_{U^\top}(i-24)
\Vert \mathcal S_{24}\!\left(R_D(i-24)\right),
& 25\le i\le48.
\end{cases}
```

The arguments of each local representation are within its valid
row range. No out-of-range indices or undefined blocks are introduced.

#### Row representations of U and its transpose

For a local row $t\in\{1,\ldots,24\}$, define

```math
q=\left\lfloor\frac{t-1}{6}\right\rfloor+1,
\qquad r=t-6(q-1).
```

Then

```math
R_U(t)=
\begin{cases}
\mathcal S_{6(q-1)}\!\left(R_C(r)\right),
& q\ne3, \\[2mm]
\mathcal S_{12}\!\left(R_C(r)\right)
\Vert \mathcal S_{18}\!\left(R_F(r)\right),
& q=3,
\end{cases}
```

and

```math
R_{U^\top}(t)=
\begin{cases}
\mathcal S_{6(q-1)}\!\left(R_{C^\top}(r)\right),
& q\ne4, \\[2mm]
\mathcal S_{12}\!\left(R_{F^\top}(r)\right)
\Vert \mathcal S_{18}\!\left(R_{C^\top}(r)\right),
& q=4.
\end{cases}
```

The four copies of $C$ are handled by one placement rule.
The block $F$ contributes only to block row 3 of $U$, and
$F^\top$ only to block row 4 of $U^\top$.

### 8.3 Second representation: direct composition of 6 × 6 blocks

For a global row $i\in\{1,\ldots,48\}$, define

```math
q=\left\lfloor\frac{i-1}{6}\right\rfloor+1,
\qquad r=i-6(q-1).
```

Here $q$ is the global block-row index and $r$ is the row
index within a $6\times6$ block.

Let $\mathcal L_q$ be the list of pairs $(j,M)$ specifying
the nonzero blocks in block row $q$: $j$ is the block-column
index and $M=K_{qj}$.

| Block row $q$ | List $\mathcal L_q$ |
|---|---|
| 1 | $(1,D_1),(2,X),(4,Y),(5,C)$ |
| 2 | $(1,X^\top),(2,D_2),(3,Y),(6,C)$ |
| 3 | $(2,Y^\top),(3,D_3),(4,X^\top),(7,C),(8,F)$ |
| 4 | $(1,Y^\top),(3,X),(4,D_4),(8,C)$ |
| 5 | $(1,C^\top),(5,D_5),(6,\widehat X)$ |
| 6 | $(2,C^\top),(5,\widehat X^\top),(6,D_6),(7,\widehat Y)$ |
| 7 | $(3,C^\top),(6,\widehat Y^\top),(7,D_7),(8,\widehat X^\top)$ |
| 8 | $(3,F^\top),(4,C^\top),(7,\widehat X),(8,D_8)$ |

The complete row representation is

```math
R_A^{(2)}(i)=
\underset{(j,M)\in\mathcal L_q}{\Vert}
\mathcal S_{6(j-1)}\!\left(R_M(r)\right).
```

The concatenation follows the order of the table. Equivalently,
each pair $(s,v)\in R_M(r)$ becomes the global pair

```math
\bigl(6(j-1)+s,\ v\bigr).
```

This formulation directly reuses the row representations of
$D_q,X,Y,\widehat X,\widehat Y,C,F$ and their transposes.

### 8.4 Relationship between the two representations

Both constructions list exactly the nonzero entries of the
same row of $A$, with identical coefficients and no duplicates.

Their list order can differ. The second construction gives
increasing global column indices if every local row list is
itself sorted. The first construction need not be sorted.

Both contain at most 12 pairs per row, consistently with
the matrix sparsity $d=12$.

The first representation explicitly uses the additive
decomposition and repeated placement rules. The second
uses one block-selection table and one local-to-global
column conversion. These are two descriptions of the same
matrix; neither description alone establishes a smaller
circuit gate count.
