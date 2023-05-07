function [nodes,norme,lengths,matrix_adiacenza]=full_tree_pruned_transport(matpde)
% Construction of the tree with the geometrical pruning criterion for the
% matrix version

X0 = matpde.X0;
control = matpde.control;
h = matpde.h;
tf = matpde.T;
nt = length(0:h:tf);

tol = h^2;
tolinv = 1/tol;
[~,s1,~] = svds(X0,1);
flo = floor(s1*tolinv);
M = containers.Map(flo,[1,1]);

na = length(control);
nodes(:,:,1)=X0;
checknorms(:,1) = s1;
vold = X0;
norme(1) = norm(X0,'fro')^2;

lengths(1) = 1;
new_adiacenza = zeros(na,1);
contnodes = 1;

cont_ad = 1;
prunes = 0;
for time = 1:nt-1

    cont = 0;
    s1 = size(vold,3);

    for j = 1:s1

        for k = 1:na

            new_mat = mat_SI_pod_transport(matpde, vold(:,:,j),control(k));
            s1 = norm(new_mat,'fro');
            flo = floor(s1*tolinv);

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
                [flag,indice] = check_mat(new_mat,nodes(:,:,nodestocheck),s1,checknorms(:,nodestocheck),10*tol);
            else
                flag = 0;
            end


            if flag
                new_adiacenza(k) = indice;
                prunes = prunes+1;
            else
                cont = cont+1;
                v(:,:,cont) = new_mat;
                contnodes = contnodes+1;
                nodes(:,:,contnodes) = new_mat;
                checknorms(:,contnodes) = s1;
                new_adiacenza(k) = contnodes;
                norme(contnodes) = norm(new_mat,'fro')^2;

                if flag_flo
                    M(flo) = [M(flo) contnodes];
                else
                    M(flo) = contnodes;
                end
            end

        end

        matrix_adiacenza(:,cont_ad) = new_adiacenza;
        cont_ad = cont_ad+1;

    end

    vold = v;
    lengths(time+1) = cont;

end
end
