# Reconstruct a 3D array from Tucker3 factors

Computes `X = A %*% G %*% kronecker(C, B)^T` where G is the unfolded
core tensor, then reshapes to `[n_a x n_t x n_g]`.

## Usage

``` r
reconstruct_array(A, B, C, G_mat)
```

## Arguments

- A:

  Mode-1 factor matrix (ages).

- B:

  Mode-2 factor matrix (years).

- C:

  Mode-3 factor matrix (groups).

- G_mat:

  Mode-1 unfolded core tensor.

## Value

3D array `[nrow(A) x nrow(B) x nrow(C)]`.
