# HB/bcsstk03 — structure and detected patterns

[SuiteSparse Matrix Collection — HB/bcsstk03](https://sparse.tamu.edu/HB/bcsstk03)

MATLAB indices, starting at 1. Unless stated otherwise, indices follow the original file ordering. 

## General properties

| Property | Result |
|---|---|
| Size | $112\times112$ |
| Type | Real; stiffness matrix |
| Symmetry | symmetric |
| Nonzero entries | 640 |
| Sparsity coefficient | 6 |
| Connected components of the support graph | **2**, with 56 vertices each |
| Half-bandwidth $\max_{A_{ij}\ne0}\lvert i-j\rvert$ | 7 |


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

## 2. Form of the upper off-diagonal blocks

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

- Each $T_i$ has exactly 8 nonzero entries; all 27 blocks therefore have the same support.
- The lower off-diagonal blocks are the transposes $T_i^\top$; the $T_i$ themselves are not symmetric.

## 3. Form of the Diagonal blocks

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

Only 3 blocks have $$\delta_i \ne 0$$:

$$
\delta_i=0\quad\text{for }i\notin\{21,25,27\}.
$$

### Diagonal blocks support patterns

$$
v_i=0\quad\Longleftrightarrow\quad i\in\{4,5,6,18\}.
$$

These four blocks are diagonal and have 4 nonzero entries. The other 24 blocks have 8 nonzero entries and share the same support. The three corrections $\delta_i$ do not change the support.

## 4. Exact repetitions amongst blocks

| Family | Exactly identical blocks |
|---|---|
| Diagonal blocks | $D_4=D_5=D_6$ |
| Upper off-diagonal blocks | $T_3=T_4=T_5=T_6=T_7$ |
| Upper off-diagonal blocks | $T_{17}=T_{18}$ |

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

## 6. Connected components

The two components are

$$
G_1=\bigcup_{i=1}^{28}\{4i-3,4i\},\qquad
G_2=\bigcup_{i=1}^{28}\{4i-2,4i-1\}.
$$
