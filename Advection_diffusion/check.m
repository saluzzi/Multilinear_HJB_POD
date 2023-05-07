function [flag,indice] = check(u,nodes,tol)
% Check for the vector pruning criterion

v = sum((nodes-u).^2);
flag = logical(v<tol);

if flag
    [~,indice] = min(v);
else
    indice=0;
end

end