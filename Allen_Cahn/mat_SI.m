function [Xnext,Fcomp] = mat_SI(matpde,Xi,Fcomp,u)

% Computation of the solution for the next level with control "u"

X0 = matpde.X0;
h = matpde.h;

nn1 = size(matpde.redvA,1);
nn2 = size(matpde.redvB,2);
L = eye(nn1,nn2);

%% Calculate X1

dd1 = diag(matpde.redDA);
dd2 = diag(matpde.redDB);
VAI = matpde.redvAi;
VBI = matpde.redvBi;

stpsiz = h;
for i = 1:nn1
    for j = 1:nn2
        L(i,j) = 1/((1-stpsiz*dd1(i)) - (stpsiz*dd2(j)));
    end
end

rhs = VAI*(Xi + h*Fcomp+ h*X0*u)*matpde.redvB;
rhs2 = L.*rhs;
Xnext = matpde.redvA*rhs2*VBI;

end
