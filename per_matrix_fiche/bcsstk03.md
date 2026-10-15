# HB/bcsstk03

[SuiteSparse Matrix Collection — HB/bcsstk03](https://sparse.tamu.edu/HB/bcsstk03)

MATLAB indices, starting at 1. Unless stated otherwise, indices follow the original file ordering. The matrix in HB/bcsstk03 will be denoted A.

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


## Support visualization

![Original support, block close-up, and grouped connected components](images/bcsstk03.png)

Each point represents a nonzero entry, regardless of its value or sign. The last panel uses the permutation defined in the connected-components section.

## 1. Exact structure of A

Using MATLAB submatrix notation, define for $i=1,\ldots,28$,

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

Only 3 blocks have $$\delta_i \ne 0$$ :

$$
\delta_i=0\quad\text{for }i\notin\{21,25,27\}.
$$

### Diagonal blocks support patterns

$$
v_i=0 for i\in\{4,5,6,18\}.
$$

These four blocks are diagonal and have 4 nonzero entries. The other 24 blocks have 8 nonzero entries and share the same support. The three corrections $\delta_i$ do not change the support.

## 4. Exact repetitions amongst blocks

| Family | Identical blocks |
|---|---|
| Diagonal blocks | $D_4=D_5=D_6$ |
| Upper off-diagonal blocks | $T_3=T_4=T_5=T_6=T_7$ |
| Upper off-diagonal blocks | $T_{17}=T_{18}$ |

## 5. Row representation $R_A$

For a row index $i\in\{1,\ldots,112\}$, define its block index by

$$
k=\left\lfloor\frac{i-1}{4}\right\rfloor+1.
$$

We describe $R_A(i)$ as a list of pairs $(j,A_{ij})$,
where $j$ is a column index. The parameters are those defined
in the block templates above.

### Case 1: $i\equiv1\pmod4$

$$
R_A(i)=\bigl[
(i-4,a_{k-1}),
(i-1,-b_{k-1}),
(i,u_k),
(i+3,v_k),
(i+4,a_k),
(i+7,b_k)
\bigr].
$$

### Case 2: $i\equiv2\pmod4$

$$
R_A(i)=\bigl[
(i-4,a_{k-1}),
(i-3,b_{k-1}),
(i,u_k+\delta_k),
(i+1,-v_k),
(i+4,a_k),
(i+5,-b_k)
\bigr].
$$

### Case 3: $i\equiv3\pmod4$

$$
R_A(i)=\bigl[
(i-5,-b_{k-1}),
(i-4,c_{k-1}),
(i-1,-v_k),
(i,w_k),
(i+3,b_k),
(i+4,c_k)
\bigr].
$$

### Case 4: $i\equiv0\pmod4$

$$
R_A(i)=\bigl[
(i-7,b_{k-1}),
(i-4,c_{k-1}),
(i-3,v_k),
(i,w_k),
(i+1,-b_k),
(i+4,c_k)
\bigr].
$$

### Boundary and zero-entry conventions

The lists above are interpreted with the following rules:

- For $k=1$, omit the first two pairs: there is no previous
  coupling block $T_0$.
- For $k=28$, omit the last two pairs: there is no next
  coupling block $T_{28}$.
- Omit any pair whose coefficient is zero. In particular,
  $v_k=0$ for $k\in\{4,5,6,18\}$.

Boundary pairs are omitted before evaluating their coefficients,
so undefined parameters such as $a_0$ or $b_{28}$ are never used.

With these conventions, each nonzero entry of row $i$ appears
exactly once in $R_A(i)$. The column indices are distinct and
listed in increasing order.

Each list contains at most six pairs, and this maximum is
attained. Thus the matrix has row sparsity 6 and, by symmetry,
column sparsity 6:

$$
d=6.
$$
## 6. Connected components

The two components are

$$
G_1=\bigcup_{i=1}^{28}\{4i-3,4i\},\qquad
G_2=\bigcup_{i=1}^{28}\{4i-2,4i-1\}.
$$
