# 4D Tucker reconstruction

Computes X\_(1) = A G\_(1) (D kron C kron B)^T.

## Usage

``` r
reconstruct_array4(A, B_mat, C, D_mat, G_mat)
```

## Arguments

- A:

  Mode-1 factor matrix (age).

- B_mat:

  Mode-2 factor matrix (time).

- C:

  Mode-3 factor matrix (group).

- D_mat:

  Mode-4 factor matrix (cause).

- G_mat:

  Mode-1 unfolded core tensor.

## Value

4D array \[n_a x n_t x n_g x n_d\].
