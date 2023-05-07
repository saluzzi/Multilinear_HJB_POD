function [Xnext] = mat_SI_pod_transport(matpde,ic,u)
% Computation of the solution for the next level with control "u" for the
% reduced dynamics

VA = matpde.VA; VB = matpde.VB;
VAi = matpde.VAi; VBi = matpde.VBi;
L = matpde.L;
hh=matpde.h;
rhs= ic*(1 + hh*u); %bilinear
Xnext = VA*( L.*(VAi*(rhs)*VB) )*VBi;

end