function [nodes_flag,lengths_flag]=refiniment_tree(PDE,nodes,lengths,V,rho,tol_flag,cost_fun)

% Construction of the refined tree with ratio "rho" retaining the nodes with lowest 
% values for the value function

N = 3; % starting level
nodes_flag = nodes(:,:,1:sum(lengths(1:N-1)));
lengths_flag = lengths(1:N-1);

for k = N:PDE.nt
    len = sum(lengths(1:k));
    len1 = sum(lengths(1:k-1))+1;
    node = nodes(:,:,len1:len);
    appo = V(k);
    appo = appo(len1:len);
    number = floor(length(appo)*rho);
    [~,index] = mink(appo,number);

    if max(appo(index))<cost_fun(k)
        new_node = node(:,:,(appo<=cost_fun(k)+tol_flag));
        lengths_flag(k) = size(new_node,3);
        nodes_flag(:,:,end+1:end+lengths_flag(k)) = new_node;
    else
        nodes_flag(:,:,end+1:end+number) = node(:,:,index);
        lengths_flag(k) = number;
    end

end

