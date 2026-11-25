# Research handoff for Codex: periodic pattern mining in sparse block matrices

## Context and objective

This is a research prototype for identifying structural patterns in sparse matrices after partitioning them into square blocks of size `b x b`. The broader research examines which block sizes reveal useful structure; examples previously discussed include SuiteSparse `bcsstk` matrices and block sizes 2, 3, and intermediate candidate sizes 14 through 24 in another matrix. Block half-bandwidth `h_b`, including the case `h_b = 1` (block tridiagonal), is an existing structural descriptor. Do **not** assume we have finalized a criterion for the *best* block size; retain candidates and characterize them rather than prematurely ranking them. The role of `R_a` in the broader research has not yet been specified enough to implement here.

The research idea is to adapt Galbrun et al., **"Mining Periodic Patterns with a MDL Criterion"** (arXiv:1807.01706, supplied separately as a PDF), to the positions of blocks in a sparse matrix.

Proposed correspondence:

- event `alpha` -> a **block type** (not a unique block occurrence)
- timestamp `t` -> integer block position within a chosen one-dimensional scan
- cycle period `p` -> integer distance between repeated block occurrences
- start `tau` -> first block position in the pattern
- repetitions `r` -> number of occurrences
- corrections `E` -> position spacing deviations from the ideal cycle
- optional pattern trees -> nested or multi-event repetitions (later milestone)

We can extract a 1D event sequence along any fixed block row (`t=j`), fixed block column (`t=i`), or fixed diagonal parallel to the main diagonal (`j-i=d`, `t=i`). Treat these as **separate analyses** initially. A diagonal identifier and line origin/length must be retained in metadata so an encoded position can be mapped back to `(i,j)`. Do not assume different scan directions yield independent coverage of the same matrix: a block can appear in all three scans.

## Important research choices still open

1. **Event definition:** The user has proposed `event = block`. Distinguish (a) binary block *presence* and (b) exact binary support signature of a block. Build the type-extraction interface so both are possible, but avoid assuming one is scientifically preferable. Numerical values should not matter for sparsity-pattern equivalence, unless explicitly requested later.
2. **Cost model:** Initially test the paper's original one-dimensional MDL code lengths; do not invent a replacement cost or assume it measures computational advantages of a block structure. A relative code length may measure positional compressibility without measuring quality for sparse computation.
3. **Comparison across block sizes:** Do not compare partial pattern code lengths across `b` as though they were complete matrix description lengths. The number of blocks, distinct block signatures, domain sizes, and dictionary costs change with `b`. A principled cross-size MDL comparison would require a consistent *lossless* encoding of the same original matrix, including block dictionary/contents, positional/domain metadata, patterns, and residuals.
4. **Shift deviations:** The paper treats small timestamp offsets as compressible noise. A block displacement may have structural significance, so expose `allow_corrections` and examine results with `E=0` before enabling imperfect cycles.
5. **Tree patterns and 2D periodicity:** Do NOT implement these in the first iteration. They are research extensions to examine after the 1D foundation is tested.

## Relevant source: formulas for the original single-event cycle model

Read Sections 3–4 and 6 of the attached paper (especially pp. 4–10, 10–13, 27–30). Base-2 logarithms. For a 1D event sequence `S`, let `N=|S|`, `Delta(S) = t_max-t_min`, and `n_alpha = |S^(alpha)|`. With the paper's conditioning and separately transmitted sequence metadata, the residual cost is

`L(o=(t, alpha)) = log2(Delta(S)+1) - log2(n_alpha/N)`.

A cycle is `C=(alpha, r, p, tau, E)`, `r>1`, `p>0`, `tau>=0`, `E=[e_1,...,e_(r-1)]`, where reconstructed positions obey `t_1=tau` and `t_(k+1)=t_k+p+e_k`. Let `sigma(E)=sum(E)`. The paper derives the simplified cycle code length

`L(C)=log2(N) + log2(floor((Delta(S)-sigma(E))/(r-1))) + log2(Delta(S)-sigma(E)-(r-1)*p+1) + 2*(r-1) + sum(abs(e) for e in E)`.

This combines the event and cycle-length costs (`L(alpha)+L(r)=log2(N)`). Verify that the constraints required for all log arguments and the reconstructed cover are valid; infeasible cycles should be rejected, not assigned misleading scores. The source conditions on fixed sequence bounds and event frequencies, whose metadata is omitted for within-sequence comparisons. If implementing a *fully decodable* compressor, account for those metadata bits explicitly. In particular, scan length and origin need encoding if they are not known from outside.

For a list of occurrence positions of one event, the paper chooses the median of adjacent inter-occurrence distances as the integer period minimizing `sum |e|`. For the initial prototype, exact (`E=0`) cycles are easier to validate.

The MDL objective for a collection `C` is

`L(C,S) = sum_{cycle in C} L(cycle) + sum_{o in S not covered by C} L(o)`.

Use the *same one-dimensional sequence S* for every term in a given score. Avoid counting one occurrence twice as a residual. Clearly document conventions when candidate cycles overlap. Simple cycles and periodic trees are distinct: do not treat the full tree encoding as identical to the simple-cycle encoding; the paper discusses a delimiter overhead.

## Desired implementation approach, in stages

**Stage 1 (do this first):** A minimal, testable Python prototype:

1. Load a scipy sparse matrix (e.g. Matrix Market `.mtx`) or accept a `scipy.sparse.spmatrix` as input. Avoid dense materialization of large matrices.
2. Partition it into `b x b` blocks in block-index coordinates. Specify handling for dimensions not divisible by `b` (padding versus partial blocks) and make behavior testable.
3. Extract **one** chosen family of 1D sequences at a time: rows, columns, or parallel main-diagonals. Preserve line identifiers, block coordinates, and scan extent.
4. Implement configurable event typing: `presence` and `binary_support_signature`. Define whether zero blocks are excluded (suggest initial default: exclude them, with explicit information that the scan domain includes empty positions).
5. Find exact single-event arithmetic progression cycles (`r >= 3` by default, `p > 0`, `E=0`); output every candidate's `alpha`, `r`, `p`, `tau`, covered positions, source line, and cost. Start with small, transparent synthetic matrices and correctness, not optimized mining.
6. Implement the paper's `L(o)`, `L(C)`, residual cost, and objective for a **fixed block size and scan sequence**, with well-defined domain bounds and edge-case handling; report whether a cycle is cheaper than residual encoding of the same occurrences.
7. Unit tests: an exact row repeat, an exact column repeat, an offset diagonal repeat, mixed block signatures, no repetition, partial boundary blocks, and all-zero lines. Include a test where an apparently periodic pattern is *not* cost-effective.
8. Present a small demonstration plus plain-text or tabular outputs; run tests and share the outcomes.

**Stage 2 (only after discussion):** Candidate pruning and greedy cover selection; compare recurrence detection and MDL selection.

**Stage 3 (only after discussion):** Imperfect cycles with accumulated corrections, and then horizontal/vertical pattern-tree combinations using the paper's specific encodings.

**Stage 4 (separate research question):** Define a complete encoding for comparing *different* block sizes and explore whether compression is predictive of structural benefits (e.g. block half-bandwidth or later-defined `R_a`).

## Coding priorities

- Python 3, `numpy`, `scipy`, `pytest`, standard packaging/typing where practical.
- Small functions and explicit data structures (`BlockOccurrence`, `LineSequence`, `CycleCandidate`, `MDLScore` or similar); no complex framework.
- Keep exact block coordinates and a reversible mapping from line positions to `(block_row, block_col)`.
- No unverified claim of equivalence to the full Galbrun et al. algorithm: the first version is intentionally restricted to exact simple cycles.
- Explain theoretical assumptions and implementation limitations rather than disguising them.
- First inspect any existing code/repository and reuse relevant utilities instead of overwriting work.

## First interaction requested from Codex

Please first inspect the repository and the attached paper, summarize a minimal Stage 1 interface and implementation plan, and then implement Stage 1 with unit tests. If a necessary input or repository detail is unavailable, state the assumption you chose and keep it configurable. Work one step at a time; do not expand into tree-pattern mining, global optimization, or 2D periodicity without a subsequent request.

## Context that needs to be supplied separately

- Attach `1807.01706v1.pdf` (the Galbrun et al. paper) to the Codex workspace or place it in the repository.
- Link or open the matrix-patterns code repository discussed previously: `https://github.com/AgentButterfly/matrix-patterns/tree/main/scripts` (if it is the repository to extend).
- Add representative `.mtx` test matrices or existing pattern output files, if available. Examples previously mentioned in discussion: `bcsstk15` with block size 2, `bcsstk16` with block size 3, `bcsstk14` with block size 2. These are research leads, *not validated expected outputs* for this new algorithm.
