function invM = chol2invchol(M)
%Function that computes the inverse of a matrix using its Cholesky decomposition 
L = chol(M);
invM = inv(L)*inv(L'); 