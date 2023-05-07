function [nodes,lengths,matrix_adiacenza] = tree_pod_pruning(PDE)
% Construction of the tree with the geometrical pruning criterion

nodes = PDE.icred;
control = PDE.control_mor;
dt = PDE.dt_mor;
nt = PDE.nt_mor;
na = PDE.na_mor;
tol = PDE.tol;

vold = nodes;
matrix_adiacenza = [];
number_nodes = 1;
s = 1;
v = [];
lengths = zeros(1,nt);
lengths(1) = 1;
tolinv = 1/tol;
flo = floor((nodes(1))*tolinv);
M = containers.Map(flo,[1 1]);
tol2 = tol*tol;
new_adiacenza = zeros(na,1);
contnodes = 1;
A = PDE.Apod;
A2 = eye(size(A))-dt*A;
A2 = sparse(A2);
A3 = sparse(eye(size(A)));
I = sparse(eye(size(A3)));
dtcontrol = dt*control;


for time = 2:nt
    for j = 1:s
        for k = 1:na

            new = A2\((A3+I*dtcontrol(k))*vold(:,j));
            flo = floor(sum(new.^2)*tolinv);

            nodestocheck = [];
            flag_flo = 0;
            if isKey(M,flo)
                nodestocheck = M(flo);
                flag_flo = 1;
            end
            if isKey(M,flo-1)
                nodestocheck = [nodestocheck M(flo-1)];
            end
            if isKey(M,flo+1)
                nodestocheck = [nodestocheck M(flo+1)];
            end

            if ~isempty(nodestocheck)
                [flag,indice] = check(new,nodes(:,nodestocheck),tol2);
            else
                flag = 0;
            end
           
                
                if flag
                    new_adiacenza(k) = indice;
                else
                    v = [v new(:,end)];
                    nodes = [nodes new(:,end)];
                    contnodes = contnodes+1;
                    new_adiacenza(k) = contnodes;
                    if flag_flo
                        M(flo) = [M(flo) contnodes];
                    else
                        M(flo) = contnodes;
                    end
                end


        end
        matrix_adiacenza = [matrix_adiacenza new_adiacenza];
    end
    vold = v;
    s = size(v,2);
    v = [];
    lengths(time) = s;
    number_nodes = number_nodes+s;
end
end