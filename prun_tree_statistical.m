function [norme,lengths,matrix_adiacenza,nodes] = prun_tree_statistical(matpde,F,nodes_flag,lengths_flag,tol_flag)
% Construction of the tree with the statical pruning criteria

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
nodes(:,:,1) = X0;

vold = X0;

norme(1) = norm(X0,'fro')^2;

s1 = 1;
lengths = zeros(1,nt);
lengths(1) = 1;
new_adiacenza = zeros(na,1);
contnodes = 1;
cont_ad = 1;

for time = 1:nt-1
    cont = 0;
    ind = sum(lengths_flag(1:time))+1:sum(lengths_flag(1:time+1));
    actual_nodes = nodes_flag(:,:,ind);
    Xmin = min(actual_nodes,[],3)-tol_flag;
    Xmax = max(actual_nodes,[],3)+tol_flag;
    v = [];
    while isempty(v)
    for j = 1:s1
        Feval = Dleft*F(VII*vold(:,:,j)*WJJ')*Dright';
        for k = 1:na
            new_mat = mat_SI(matpde, vold(:,:,j),Feval,control(k));
            if time>=4
                if new_mat>=Xmin & new_mat<=Xmax
                    flag = 1;
                else
                    flag = 0;
                end
            else
                flag = 1;
            end
            if flag
                cont = cont+1;
                v(:,:,cont) = new_mat;
                nodes(:,:,end+1)=new_mat;
                
                contnodes=contnodes+1;
                new_adiacenza(k)=contnodes;
                norme(contnodes) = norm(new_mat,'fro')^2;
            else
                new_adiacenza(k)=1;
            end
            
        end
        
        matrix_adiacenza(:,cont_ad)=new_adiacenza;
        cont_ad=cont_ad+1;
        
    end
    end
    
    vold=v;
    s1=size(v,3);
    lengths(time+1)=s1;
    
end
end