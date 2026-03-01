# Modified Cholesky–Banachiewicz Solver for Tridiagonal Systems (MATLAB)

This project presents a MATLAB implementation of a memory-optimized Cholesky–Banachiewicz decomposition for solving linear systems with symmetric positive-definite tridiagonal matrices.

The algorithm stores only the main and subdiagonal of the matrix, performs decomposition, and solves the linear system using forward and backward substitution. Numerical properties such as determinant accuracy, conditioning, stability, and relative error are analyzed and compared with MATLAB’s built-in solver.

## Features

- Cholesky–Banachiewicz decomposition for tridiagonal matrices
- Memory-optimized storage (two diagonals only)
- Linear system solver Ax = b
- Determinant computation
- Numerical stability and conditioning analysis
- Comparison with MATLAB built-in solver (A\b)
- Multiple test cases with varying matrix properties

## Project Structure
- src/ MATLAB source code
- report/ Project report (PDF)


## Main Components

- `algorytm.m` – Cholesky–Banachiewicz decomposition  
- `rozwiazanie_ukladu_rownan.m` – linear system solver  
- `main.m` – test script and numerical analysis  


## Numerical Analysis

The implementation is evaluated using multiple tridiagonal matrices with different magnitudes and conditioning properties. Results include:

- determinant error
- relative solution error
- condition number
- decomposition error
- stability and correctness coefficients

## Report

Detailed mathematical description and experimental results:

👉 [Cholesky–Banachiewicz Report (PDF)](report/Sprawozdanie_projekt2.pdf)
