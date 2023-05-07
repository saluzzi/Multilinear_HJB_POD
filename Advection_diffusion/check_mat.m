function [flag,indice] = check_mat(u,nodes,snew,sold,tol)
% Check for the matrix pruning criterion

smallnorm = abs(snew -  sold);
[~,jj] = min(smallnorm);
nrm = norm(u - nodes(:,:,jj));

flag=logical(nrm<tol);

if flag

    indice = jj;
else

    indice=0;
end

end
