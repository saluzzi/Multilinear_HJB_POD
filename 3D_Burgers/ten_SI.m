function [Xnew, Fnew]= ten_SI(ode, IC,u);


A = ode.A;

for k = 1:3
    [VA3{k},DA3{k}] = eig(full(A{k,3})');
    VAi3{k} = inv(VA3{k});
    ddA3{k} = diag(DA3{k});
end

B = ode.B;
h = ode.h;
nn = ode.nn+1;

for k = 1:3
    
    Xnew{k} = IC{k};
    
end

for k = 1:3
    
    for r = 1:3
        Fnew{k,r} = ttm(Xnew{k},B{r},r).*Xnew{r};
    end
    
    
    
    Fnew{k} = Fnew{k,1}+Fnew{k,2}+Fnew{k,3};
    
    rhs = h*Fnew{k} + (1+h*u)*Xnew{k};
    rhs = double(tenmat(rhs,3,'t'));
    rhstrans = rhs*VA3{k};
    
    for l = 1:nn % Solve a sequence of nn sylvester equations to determine the solution
        
        RHS = reshape(rhstrans(:,l),nn,nn);
        
        A11 = (1-h*ddA3{k}(l))*eye(size(A{k,1})) - h*A{k,1};
        A21 = -h*A{k,2}';
        YY = lyap(A11,A21,-RHS);
        
        XX(:,l) = YY(:);
        
    end
    
    X3{k} = XX*VAi3{k};
    Xnew{k} = tensor(X3{k},[nn,nn,nn]);
    
end




for k = 1:3
    
    for r = 1:3
        Fnew{k,r} = ttm(Xnew{k},B{r},r).*Xnew{r};
    end
    
    
    
    Fnew{k} = Fnew{k,1}+Fnew{k,2}+Fnew{k,3};
end


end
