# AMRVW

Implementation of core-chasing algorithms for finding eigenvalues of
factored matrices.

This repository provides a Julia package implementing the methods, as
applied to the problem of finding the roots of polynomials through the
computation of the eigenvalues of a sparse factorization of the
companion matrix in:

*Fast and Backward Stable Computation of Roots of Polynomials.* Jared
L. Aurentz, Thomas Mach, Raf Vandebril, and David S. Watkins SIAM
J. Matrix Anal. Appl., 36(3), 942–973. (2015)
[https://doi.org/10.1137/140983434](https://doi.org/10.1137/140983434)

*Fast and backward stable comptation of roots of polynomials, Part II:
backward error analysis; companion matrix and companion pencil.* By
Jared L. Aurentz, Thomas Mach, utation of roots of polynomials, Part
II: backward error analysis; companion matrix and companion pencil. By
Jared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, David
S. Watkins; arXiv:1611.02435

The methods are summarized in monograph format:

*Core-Chasing Algorithms for the Eigenvalue Problem.* by Jared
L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, and David
S. Watkins;
[https://doi.org/10.1137/1.9781611975345](https://doi.org/10.1137/1.9781611975345)

As well, the twisted algorithm from *A generalization of the multishift QR algorithm* by Raf Vandebril and David S. Watkins; [https://doi.org/10.1137/11085219X](https://doi.org/10.1137/11085219X) is implemented here.

The core-chasing algorithms utilize Francis's QR algorithm on sparse
factorizations of the respected companion matrix. For polynomials with
real coefficients, the storage requirements are ``O(n)`` and the
algorithm requires ``O(n)`` flops per iteration, or ``O(n^2)`` flops
overall. The basic QR algorithm applied to the full companion matrix
would require ``O(n^2)`` storage and ``O(n^3)`` flops overall.

## Reference


```@autodocs
Modules = [AMRVW]
```
