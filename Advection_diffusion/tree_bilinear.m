function [nodes,lengths,matrix_adiacenza]=tree_bilinear(PDE)
% Construction of the tree for a vectorized bilinear dynamics

x0 = PDE.ic;
nx=length(x0);
control = PDE.control_snap;
dt = PDE.dt_snap;
nt = PDE.nt_snap;
na = PDE.na_snap;
nodes = zeros(nx,nt*(nt+1)/2);
nodes(:,1) = x0;
vold = x0;
matrix_adiacenza = [];
lengths = 1:nt;
new_adiacenza = zeros(na,1);
contnodes = 1;
A2 = PDE.A2;
A3 = PDE.A3;
I = sparse(eye(size(A3)));
dtcontrol = dt*control;

for time = 1:nt-1

    v = zeros(nx,time+1);
    cont = 0;

    for j = 1:time

        new = A2\((A3+I*dtcontrol(1))*vold(:,j));
        cont = cont+1;
        v(:,cont) = new(:,end);
        contnodes = contnodes+1;
        new_adiacenza = [contnodes;contnodes+1];

    end

    new = (A2)\((A3+I*dtcontrol(2))*vold(:,time));
    matrix_adiacenza = [matrix_adiacenza new_adiacenza];
    cont = cont+1;
    v(:,cont) = new(:,end);
    contnodes = contnodes+1;
    vold = v;
    nodes(:,contnodes-cont+1:contnodes) = v;
    
end