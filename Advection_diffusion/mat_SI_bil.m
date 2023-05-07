function [Xnext] = mat_SI_bil(matpde,ic,u)
% Computation of the solution for the next level with control "u"

A11 = matpde.A11;
A22 = matpde.A22;
hh=matpde.h;

rhs= ic*(1 + hh*u); % bilinear control
Xnext = lyap(A11, A22, -rhs);

end