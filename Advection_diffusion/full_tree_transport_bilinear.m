function [Uset,Vset] = full_tree_transport_bilinear(matpde,tau,kappa)
% Computation of the tree structure driven by the reduced dynamics


X0 = matpde.X0;
control = matpde.control;
tf = matpde.T;
nt = length(0:matpde.h:tf);

[uu,ss,vv] = svd(X0);
gamma = sum(diag(ss)/ss(1,1)>tau);
vleftold{1}(:,:,1) = uu(:,1:gamma);
vrightold{1}(:,:,1) = vv(:,1:gamma);
vsingsold{1} = diag(ss(1:gamma,1:gamma));
uset2 = vleftold{1}(:,:,1);
vset2 = vrightold{1}(:,:,1);
sset2 = vsingsold{1}';

number_nodes = 1;
lengths = zeros(1,nt);
lengths(1) = 1;
contnodes = 1;
nodecount = 0;
ii = 0;

for time = 1:nt-1

    s1 = lengths(time);
    cont = 0;

    for j = 1:s1

        ii = ii+1;
        Xi = vleftold{1}(:,:,j)*diag(vsingsold{1}(:,j))*vrightold{1}(:,:,j)';
        [new_mat]= mat_SI_bil(matpde, Xi,control(1));

        [uu1,ss1,vv1] = svd(new_mat);

        cont = cont+1;
        gamma2 = 10;
        vleft(:,:,cont) = uu1(:,1:gamma2);
        vright(:,:,cont) = vv1(:,1:gamma2);
        vsings(:,cont) = diag(ss1(1:gamma2,1:gamma2));

        nodecount = nodecount+1;
        uset2 = [uset2,uu1];
        vset2 = [vset2,vv1];
        sset2 = [sset2, diag(ss1)'];

        if size(uset2,2) > kappa
            [ssort, id] = sort(sset2, 'descend');
            uset2 = uset2(:,id(1:kappa));
            vset2 = vset2(:,id(1:kappa));
            sset2 = ssort(1:kappa);
        end

        contnodes = contnodes+1;
    end

    ii = ii +1;
    Xi = vleftold{1}(:,:,s1)*diag(vsingsold{1}(:,s1))*vrightold{1}(:,:,s1)';
    [new_mat]= mat_SI_bil(matpde, Xi,control(2));

    cont = cont+1;
    [uu1,ss1,vv1] = svd(new_mat);


    gamma2 = 10;
    vleft(:,:,cont) = uu1(:,1:gamma2);
    vright(:,:,cont) = vv1(:,1:gamma2);
    vsings(:,cont) = diag(ss1(1:gamma2,1:gamma2));

    nodecount = nodecount+1;
    uset2 = [uset2,uu1];
    vset2 = [vset2,vv1];
    sset2 = [sset2, diag(ss1)'];

    if size(uset2,2) > kappa
        [ssort, id] = sort(sset2, 'descend');
        uset2 = uset2(:,id(1:kappa));
        vset2 = vset2(:,id(1:kappa));
        sset2 = ssort(1:kappa);
    end

    contnodes = contnodes+1;

    vleftold{1} = vleft;
    vrightold{1} = vright;
    vsingsold{1} = vsings;

    number_nodes = number_nodes+cont;
    lengths(time+1) = cont;
end


[Uset, su, ~] = svd(uset2,0);
[Vset, sv, ~] = svd(vset2,0);

ssetu = diag(su);
ssetv = diag(sv);
tauh = tau/sqrt(kappa);

k1 = length(ssetu)-sum((cumsum(ssetu(end:-1:1))/sum(ssetu))<tauh);
k2 = length(ssetv)-sum((cumsum(ssetv(end:-1:1))/sum(ssetv))<tauh);

Uset = Uset(:,1:k1);
Vset = Vset(:,1:k2);

end