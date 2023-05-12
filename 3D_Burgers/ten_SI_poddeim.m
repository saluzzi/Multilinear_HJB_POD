function [Ynew]= ten_SI_poddeim(redode, IC,u);


Ak = redode.Ak;
Deim = redode.Deim;
VAismall = redode.VAismall;
VAsmall = redode.VAsmall;  DAsmall = redode.DAsmall;
B = redode.B; Vbasis = redode.Vbasis; h = redode.h; tau = redode.tau;
FVbasis =  redode.FVbasis; II = redode.II; p = redode.p;


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




for k = 1:3
    
    
    rhsmall{k} = (1+h*u)*Ynew{k} + h*Fnewsmall{k};
    
    romtimer2=tic;
    rhsmall{k} = double(tenmat(rhsmall{k},3,'t'));
    rhstransmall{k} = rhsmall{k}*VAsmall{k,3};
    
    for l = 1:p{k,3} % Solve a sequence of nn sylvester equations to determine the solution
        
        RHSmall{k} = reshape(rhstransmall{k}(:,l),p{k,1},p{k,2});
        
        A11small{k} = (1-h*ddAsmall{k,3}(l))*eye(size(Ak{k,1})) - h*Ak{k,1};
        A21small{k} = -h*Ak{k,2}';
        YYsmall{k} = lyap(A11small{k},A21small{k},-RHSmall{k});
        XXsmall{k}(:,l) = YYsmall{k}(:);
        
    end
    
    XXsmall{k} = XXsmall{k}*VAismall{k,3};
    
    Ynew{k} = tensor(XXsmall{k},[p{k,1},p{k,2},p{k,3}]);
    
end

