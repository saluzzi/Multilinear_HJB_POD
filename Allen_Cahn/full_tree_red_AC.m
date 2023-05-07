function [norme,lengths,matrix_adiacenza,nodes]=full_tree_red_AC(matpde,F)
% Computation of the tree structure driven by the reduced dynamics

X0 = matpde.X0;
control = matpde.control;
h = matpde.h;
tf = matpde.T;
VII = matpde.VII;
WJJ = matpde.WJJ;
Dleft = matpde.Dleft;
Dright = matpde.Dright;
nt = length(0:h:tf);
na = length(control);
[nx,ny] = size(X0);

nodes(:,:,1) = X0;
vold = X0;
norme(1) = norm(X0,'fro')^2;
number_nodes = 1;
s = 1;
lengths = zeros(1,nt);
lengths(1) = 1;
new_adiacenza = zeros(na,1);
contnodes = 1;

cont_ad = 1;
fprintf('Time:\n')
for time = 1:nt-1
    s1 = na^(time-1);
    v= zeros(nx,ny,s);
    cont = 0;
    fprintf('..%d..',time)

    for j=1:s1

        Feval = Dleft*F(VII*vold(:,:,j)*WJJ')*Dright';

        for k=1:na

            new_mat = mat_SI(matpde, vold(:,:,j),Feval,control(k));
            cont = cont+1;
            v(:,:,cont) = new_mat;
            nodes(:,:,end+1) = new_mat;

            contnodes = contnodes+1;
            new_adiacenza(k) = contnodes;
            norme(contnodes) = norm(new_mat,'fro')^2;

        end

        matrix_adiacenza(:,cont_ad) = new_adiacenza;
        cont_ad=cont_ad+1;

    end

    vold = v;
    lengths(time+1) = cont;
    number_nodes = number_nodes+cont;

end
end
