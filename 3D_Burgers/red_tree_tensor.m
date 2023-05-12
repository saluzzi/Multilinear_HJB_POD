function [norme,lengths,matrix_adiacenza]=red_tree_tensor(matpde,dim)

redode = matpde;
Deim = redode.Deim;
DAsmall = redode.DAsmall;
B = redode.B; Vbasis = redode.Vbasis; 
FVbasis =  redode.FVbasis; II = redode.II;
X0 = matpde.Y0;
control = matpde.control;
h = matpde.h;
tf = matpde.T;
nt = length(0:h:tf);

na = length(control);

for k = 1:3
    vold{k,1}=X0{k};
end
norme(1) = normT(X0{1})^2 + normT(X0{2})^2 + normT(X0{3})^2;
number_nodes=1;

s=1;
lengths=zeros(1,nt);
lengths(1)=1;
new_adiacenza=zeros(na,1);
contnodes=1;

t_grid = 0;
cont_ad=1;

for time=1:nt-1
    time
    t_grid = [t_grid(end), t_grid(end) + h];
    s1=na^(time-1);
    cont=0;
    
    for j=1:s1
     
        
        %%% Evaluate Nonlinear term %%%
        
        IC = vold(:,j);
        
        for k = 1:3
            ddAsmall{k,3} = diag(DAsmall{k,3});
            Ynew{k} = IC{k};
        end
        
        
        for q = 1:3
            for k = 1:3
                Yrowred{k,q} = Ynew{k};
                for r = 1:3
                    Yrowred{k,q} = ttm(Yrowred{k,q},Vbasis{k,r}(II{q,r},:),r);
                end
            end
        end
        
        
        %%% Evaluate low-dimensional nonlinear term
        
        for k = 1:3
            
            
            Fsmall1{k} = Ynew{k};
            Fsmall2{k} = Ynew{k};
            Fsmall3{k} = Ynew{k};
            
            Fsmall1{k} = ttm(Fsmall1{k}, B{k,1}(II{k,1},:)*Vbasis{k,1},1);
            Fsmall1{k} = ttm(Fsmall1{k}, Vbasis{k,2}(II{k,2},:),2);
            Fsmall1{k} = ttm(Fsmall1{k}, Vbasis{k,3}(II{k,3},:),3);
            Fsmall1{k} = Fsmall1{k}.*Yrowred{1,k};
            
            Fsmall2{k} = ttm(Fsmall2{k}, B{k,2}(II{k,2},:)*Vbasis{k,2},2);
            Fsmall2{k} = ttm(Fsmall2{k}, Vbasis{k,1}(II{k,1},:),1);
            Fsmall2{k} = ttm(Fsmall2{k}, Vbasis{k,3}(II{k,3},:),3);
            Fsmall2{k} = Fsmall2{k}.*Yrowred{2,k};
            
            Fsmall3{k} = ttm(Fsmall3{k}, Vbasis{k,1}(II{k,1},:),1);
            Fsmall3{k} = ttm(Fsmall3{k}, Vbasis{k,2}(II{k,2},:),2);
            Fsmall3{k} = ttm(Fsmall3{k}, B{k,3}(II{k,3},:)*Vbasis{k,3},3);
            
            Fsmall3{k} = Fsmall3{k}.*Yrowred{3,k};
            
            Fnewsmall{k} =  Fsmall1{k} + Fsmall2{k} + Fsmall3{k};
            
        end
        
        
        
        
        for k = 1:3
            for r = 1:3
                Fnewsmall{k} = ttm(Fnewsmall{k} ,Deim{k,r},r);
            end
        end
        
        for kk=1:na
            
            if dim
                new_ten= ten_SI(matpde, vold(:,j),control(kk));
            else
                new_ten= ten_SI_poddeim(matpde, IC,control(kk));
            end
            
            
            cont=cont+1;
            
            for k = 1:3
                v{k,cont}=new_ten{k};
            end
            
            
            contnodes=contnodes+1;
            new_adiacenza(kk)=contnodes;
            
            norme(contnodes) = normT(new_ten{1})^2 + normT(new_ten{2})^2 + normT(new_ten{3})^2;
            
            
        end
        
        matrix_adiacenza(:,cont_ad)=new_adiacenza;
        cont_ad=cont_ad+1;
        
    end
    
    vold=v;
    %     nodes(:,:,(na^(time)-1)/(na-1)+1:(na^(time+1)-1)/(na-1))=v;
    lengths(time+1)=cont;
    number_nodes=number_nodes+cont;
    
end
end