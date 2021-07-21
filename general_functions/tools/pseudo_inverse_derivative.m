function dpAdx = pseudo_inverse_derivative(A, pA, dAdx)
%dpAdx= -pA*dAdx*pA + pA*pA'*dAdx'*(1-A*pA)+(1-pA*A)*dAdx'*pA'*pA;
% see : https://mathoverflow.net/questions/25778/analytical-formula-for-numerical-derivative-of-the-matrix-pseudo-inverse

invAtA = inv(A'*A);
dpAdx = -invAtA*dAdx'*A*pA-pA*dAdx*pA+inv(A'*A)*dAdx';

