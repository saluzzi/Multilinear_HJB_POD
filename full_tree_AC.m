function [norme,Uset,Vset,usetnl,vsetnl,lengths,matrix_adiacenza] = full_tree_AC(matpde,F,tau,kappa)
% Construction of the tree structure and computation of the HO-POD-DEIM
% basis

h = matpde.h;
hi = h;
X0 = matpde.X0;
control = matpde.control;
tf = matpde.T;
nt = length(0:hi:tf);
na = length(control);
tauh = tau/sqrt(kappa);
Fnew = F(X0);
[uu,ss,vv] = svd(X0,'econ');
[uunl,ssnl,vvnl] = svd(Fnew,'econ');
norme = zeros(1,(na^(nt)-1)/(na-1));
norme(1) = sum(diag(ss).^2);
gamma = sum(diag(ss)./ss(1,1)>tau);
vleftold{1}(:,:,1)=uu(:,1:gamma);
vrightold{1}(:,:,1)= vv(:,1:gamma);
vsingsold{1} = diag(ss(1:gamma,1:gamma));
uset2 = vleftold{1}(:,:,1);
vset2 = vrightold{1}(:,:,1);
uset = uset2;
vset = vset2;
sset2 = vsingsold{1}';
gammanl = sum(diag(ssnl)./ssnl(1,1)>tau);
usetnl2 = uunl(:,1:gammanl);
vsetnl2 = vvnl(:,1:gammanl);
ssetnl2 = diag(ssnl(1:gammanl, 1:gammanl))';
matrix_adiacenza=zeros(na,(na^(nt-1)-1)/(na-1));
number_nodes=1;
lengths=zeros(1,nt);
lengths(1)=1;
new_adiacenza=zeros(na,1);
contnodes=1;
cont_ad=1;
nodecount=0;
matpde.h = hi;

fprintf('Time:')
for time=1:nt-1
    fprintf('..%d..',time)
    s1=lengths(time);
    cont=0;

    for j=1:s1
        Xi = vleftold{1}(:,:,j)*diag(vsingsold{1}(:,j))*vrightold{1}(:,:,j)';
        Fnew = F(Xi);
        for k=1:na

            [new_mat,Fnew]= mat_SI(matpde,Xi,Fnew,control(k));

            cont=cont+1;
            [uu1,ss1,vv1] = svd(new_mat,'econ');
            vleft(:,:,cont)=uu1(:,1:kappa);
            vright(:,:,cont)=vv1(:,1:kappa);
            vsings(:,cont) = diag(ss1(1:kappa,1:kappa));

            e(contnodes) = norm(new_mat - uset*(uset'*new_mat*vset)*vset', 'fro')/norm(new_mat,'fro');

            if e(contnodes) > tauh

                nodecount = nodecount+1;
                gamma = sum(diag(ss1)./ss1(1,1) > tau);
                uu = uu1(:,1:gamma); vv = vv1(:,1:gamma);ss = ss1(1:gamma, 1:gamma);
                uset2 = [uset2,uu];
                vset2 = [vset2,vv];
                sset2 = [sset2, diag(ss)'];

                if size(uset2,2) > kappa
                    [ssort, id] = sort(sset2, 'descend');
                    uset2 = uset2(:,id(1:kappa));
                    vset2 = vset2(:,id(1:kappa));
                    sset2 = ssort(1:kappa);
                end

                [uset,~] = qr(uset2,0);
                [vset,~] = qr(vset2,0);

                [uunl,ssnl,vvnl] = svd(Fnew,'econ');

                gamma = sum(diag(ssnl)./ssnl(1,1) > tau);
                uunl = uunl(:,1:gamma); vvnl = vvnl(:,1:gamma);ssnl = ssnl(1:gamma, 1:gamma);

                usetnl2 = [usetnl2,uunl];
                vsetnl2 = [vsetnl2,vvnl];
                ssetnl2 = [ssetnl2, diag(ssnl)'];

                if size(usetnl2,2) > kappa
                    [ssortnl, id] = sort(ssetnl2, 'descend');
                    usetnl2 = usetnl2(:,id(1:kappa));
                    vsetnl2 = vsetnl2(:,id(1:kappa));
                    ssetnl2 = ssortnl(1:kappa);
                end

            end

            contnodes=contnodes+1;
            norme(contnodes) = sum(diag(ss).^2);
            new_adiacenza(k)=contnodes;

        end

        matrix_adiacenza(:,cont_ad)=new_adiacenza;
        cont_ad=cont_ad+1;

    end


    vleftold{1}=vleft;
    vrightold{1} = vright;
    vsingsold{1} = vsings;

    number_nodes=number_nodes+cont;
    lengths(time+1)=cont;
end


[Uset, su, ~] = svd(uset2,0);
[Vset, sv, ~] = svd(vset2,0);

ssetu = diag(su);
ssetv = diag(sv);

k1=length(ssetu)-sum((cumsum(ssetu(end:-1:1))/sum(ssetu))<tauh);
k2=length(ssetv)-sum((cumsum(ssetv(end:-1:1))/sum(ssetv))<tauh);

Uset = Uset(:,1:k1);
Vset = Vset(:,1:k2);

[usetnl, sunl, ~] = svd(usetnl2,0);
[vsetnl, svnl, ~] = svd(vsetnl2,0);

ssetu = diag(sunl); ssetv = diag(svnl);

p1=length(ssetu)-sum((cumsum(ssetu(end:-1:1))/sum(ssetu))<tauh);
p2=length(ssetv)-sum((cumsum(ssetv(end:-1:1))/sum(ssetv))<tauh);

usetnl = usetnl(:,1:p1);
vsetnl = vsetnl(:,1:p2);

clear vleftold
clear vrightold
clear vleft
clear vright


end