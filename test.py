"""
A Python implementation of NNLS algorithm
References:
[1]  Lawson, C.L. and R.J. Hanson, Solving Least-Squares Problems, Prentice-Hall, Chapter 23, p. 161, 1974.
Contributed by Klaus Schuch (schuch@igi.tugraz.at)
based on MATLAB's lsqnonneg function
https://gist.github.com/jdelafon/b7fdc7a0bae42af56366fc7786cc5d54
"""

import numpy as np

def lsqnonneg(C, d, x0=None, tol=None, itmax_factor=3):
    '''Linear least squares with nonnegativity constraints.
    (x, resnorm, residual) = lsqnonneg(C,d) returns the vector x that minimizes norm(d-C*x)
    subject to x >= 0, C and d must be real

    A Python implementation of NNLS algorithm
    References:
    [1]  Lawson, C.L. and R.J. Hanson, Solving Least-Squares Problems, Prentice-Hall, Chapter 23, p. 161, 1974.
    Contributed by Klaus Schuch (schuch@igi.tugraz.at)
    based on MATLAB's lsqnonneg function
    https://gist.github.com/jdelafon/b7fdc7a0bae42af56366fc7786cc5d54
    '''

    eps = 2.22e-16    # from matlab
    def norm1(x):
        return abs(x).sum().max()

    def msize(x, dim):
        s = x.shape
        if dim >= len(s):
            return 1
        else:
            return s[dim]

    if tol is None: tol = 10*eps*norm1(C)*(max(C.shape)+1)

    C = np.asarray(C)

    (m,n) = C.shape
    P = np.zeros(n)
    Z = np.arange(1, n+1)

    if x0 is None: x=P
    else:
        if any(x0 < 0): x=P
        else: x=x0

    ZZ = Z

    resid = (d - np.dot(C, x))
    #resid = (d - np.dot(C, x)) / d
    w = np.dot(C.T, resid)

    outeriter=0; it=0
    itmax=itmax_factor*n
    exitflag=1

    # outer loop to put variables into set to hold positive coefficients
    while np.any(Z) and np.any(w[ZZ-1] > tol):
        outeriter += 1

        t = w[ZZ-1].argmax()
        t = ZZ[t]

        P[t-1]=t
        Z[t-1]=0

        PP = np.where(P != 0)[0]+1
        ZZ = np.where(Z != 0)[0]+1

        CP = np.zeros(C.shape)

        CP[:, PP-1] = C[:, PP-1]
        CP[:, ZZ-1] = np.zeros((m, msize(ZZ, 1)))

        #z=np.dot(np.linalg.pinv(CP), d)
        D = np.diag(1/d)
        D_square = np.square(D)
        z = np.linalg.pinv(CP.T @ D_square @ CP) @ CP.T @ D_square @ d

        z[ZZ-1] = np.zeros((msize(ZZ,1), msize(ZZ,0)))

        # inner loop to remove elements from the positve set which no longer belong
        while np.any(z[PP-1] <= tol):
            it += 1

            if it > itmax:
                max_error = z[PP-1].max()
                #return (z, sum(resid*resid), resid)
                raise Exception(f'Exiting: Iteration count {it} exceeded\n  Try raising the tolerance tol. {max_error}')

            QQ = np.where((z <= tol) & (P != 0))[0]
            alpha = min(x[QQ]/(x[QQ] - z[QQ]))
            x = x + alpha*(z-x)

            ij = np.where((abs(x) < tol) & (P != 0))[0]+1
            Z[ij-1] = ij
            P[ij-1] = np.zeros(max(ij.shape))
            PP = np.where(P != 0)[0]+1
            ZZ = np.where(Z != 0)[0]+1

            CP[:, PP-1] = C[:, PP-1]
            CP[:, ZZ-1] = np.zeros((m, msize(ZZ, 1)))

            z=np.dot(np.linalg.pinv(CP), d)
            D = np.diag(1/d)
            D_square = np.square(D)
            #z = np.linalg.pinv(CP.T @ D_square @ CP) @ CP.T @ D_square @ d
            z[ZZ-1] = np.zeros((msize(ZZ,1), msize(ZZ,0)))

        x = z
        #resid = (d - np.dot(C, x))
        resid = (d - np.dot(C, x)) / d
        w = np.dot(C.T, resid)

    return (x, sum(resid * resid), resid)


# Unittest
if __name__ == '__main__':

    from scipy.optimize import nnls
    A = np.array([[1, 0], [1, 0], [0, 1]])
    b = np.array([2, 1, 1])
    x, res = nnls(A, b)
    print(f'Scipy NNLS: {x, res**2}')
    [x, resnorm, residual] = lsqnonneg(A, b)
    print(f'Current: {x, residual}')
    x, res, _, _ = np.linalg.lstsq(A, b, rcond=None)
    print(f'Numpy LLS: {x, res}')
    D = np.diag(1/b)
    D_square = np.square(D)
    x = np.linalg.pinv(A.T @ D_square @ A) @ A.T @ D_square @ b
    print(f'% min: {x}')

    print()
    b = np.array([-1, -1, -1])

    x, res = nnls(A, b)
    print(f'Scipy NNLS: {x, res**2}')
    [x, resnorm, residual] = lsqnonneg(A, b)
    print(f'Current: {x, residual}')
    x, res, _, _ = np.linalg.lstsq(A, b, rcond=None)
    print(f'Numpy LLS: {x, res}')

    D = np.diag(1/b)
    D_square = np.square(D)
    x = np.linalg.pinv(A.T @ D_square @ A) @ A.T @ D_square @ b
    print(f'% min: {x}')
