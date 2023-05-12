function [norme,Vbasis,FVbasis,lengths,matrix_adiacenza]=integrate_and_basis_tensor(ode, kappa)


X0 = ode.X0;
A = ode.A;
B = ode.B;
tau = ode.tau;

tauh = tau;
control = ode.control;
h = ode.h;
tf = ode.T;
nt = length(0:h:tf);

na = length(control);

nn = ode.nn+1;

addpath TENSORTOOLS



%%% initiate tree storing and first basis vectors %%%

for k = 1:3
    
    for r = 1:3
        Fnew{k,r} = ttm(X0{k},B{r},r).*X0{r};
    end
    
    
    
    Fnew{k} = Fnew{k,1}+Fnew{k,2}+Fnew{k,3};
end


for k = 1:3
    
    [u{k}, pset, s{k}] = sthosvd(X0{k},kappa);
    [Fu{k}, Fpset, Fs{k}] = sthosvd(Fnew{k},kappa);
    coretreeold{k,1} = u{k}.core;
    
    for r = 1:3
        
        utreeold{k,r,1} = u{k}.U{r};
        
        N{k,r}=sum(s{k}(r,:)/s{k}(r,1)>tau);
        S{k,r} = s{k}(r,1:N{k,r});
        
        Su{k,r} = S{k,r}; Suold{k,r} = Su{k,r};
        U{k,r} = u{k}.U{r}(:,1:N{k,r}); Uold{k,r} = U{k,r}; p{k,r} = size(U{k,r},2); ku{k,r} = p{k,r};
        
        Vbasis{k,r} = U{k,r}(:,1:p{k,r});
        
        FN{k,r}=sum(Fs{k}(r,:)/Fs{k}(r,1)>tau);
        FS{k,r} = Fs{k}(r,1:FN{k,r});
        
        FSu{k,r} = FS{k,r}; FSuold{k,r} = FSu{k,r};
        FU{k,r} = Fu{k}.U{r}(:,1:FN{k,r}); FUold{k,r} = FU{k,r}; Fp{k,r} = size(FU{k,r},2); Fku{k,r} = Fp{k,r};
        
        FVbasis{k,r} = FU{k,r}(:,1:Fp{k,r});
        
        
    end
    
    norme{k,1} = norm(coretreeold{k,1}) ;
    
    
    
end

norme{1} = normT(coretreeold{1,1})^2 + normT(coretreeold{2,1})^2 + normT(coretreeold{3,1})^2;


number_nodes=1;

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
        
        for kk=1:na
            
      
            [Xnew, Fnew]= sylvsolve_tensor(ode, utreeold(:,:,j),coretreeold(:,j),control(kk));
       
            
            cont=cont+1;
            
            
            %%% Update the tree storing and  calculate the new basis vec
            %%% tors
            
            for k = 1:3
                
                [u{k}, pset, s{k}] = sthosvd(Xnew{k},kappa);
                [Fu{k}, Fpset, Fs{k}] = sthosvd(Fnew{k},kappa);
                
                coretree{k,cont} = u{k}.core;
                
                for r = 1:3
                    
                    
                    utree{k,r,cont} = u{k}.U{r};
                    
                    N{k,r}=sum(s{k}(r,:)/s{k}(r,1)>tau);
                    S{k,r} = s{k}(r,1:N{k,r});
                    
                    FN{k,r}=sum(Fs{k}(r,:)/Fs{k}(r,1)>tau);
                    FS{k,r} = Fs{k}(r,1:FN{k,r});
                    
                end
                
                
            end
            
            
            
            for k = 1:3 % Update the solution bases for all three equations
                
                
                
                for r = 1:3 %Update each of the three spacial directions
                    
                    CU{k,r} = u{k}.U{r}(:,1:N{k,r});
                    
                    Uold{k,r} = Vbasis{k,r};
                    Suold{k,r} = Suold{k,r}(1:p{k,r});
                    
                    if ku{k,r} >= nn
                        U{k,r} = Uold{k,r};
                        Su{k,r} = Suold{k,r};
                    else
                        [U{k,r},su1{k,r},~]=svd([Uold{k,r}*diag(Suold{k,r}),CU{k,r}*diag(S{k,r})],0); strunc{k,r}=diag(su1{k,r});
                        p{k,r}=size([Uold{k,r},CU{k,r}],2)-sum((cumsum(strunc{k,r}(end:-1:1))/sum(strunc{k,r}))<tauh);
                        Suold{k,r} = strunc{k,r};
                        ku{k,r} = size([Uold{k,r},CU{k,r}],2);
                    end
                    
                    
                    
                    Vbasis{k,r} = U{k,r}(:,1:p{k,r});
                    
                end
            end
            
            
            
            for k = 1:3 % Update the nonlinear bases for all three equations
                
                for r = 1:3 %Update each of the three spacial directions
                    
                    FCU{k,r} = Fu{k}.U{r}(:,1:FN{k,r});
                    
                    if isempty(FUold{k,r})
                        FSu{k,r} = FS{k,r}; FSuold{k,r} = FSu{k,r};
                        FU{k,r} = FCU{k,r}; FUold{k,r} = FU{k,r}; Fp{k,r} = size(FU{k,r},2); Fku{k,r} = Fp{k,r};
                    else
                        
                        FUold{k,r} = FVbasis{k,r};
                        FSuold{k,r} = FSuold{k,r}(1:Fp{k,r});
                        
                        if Fku{k,r} >= nn
                            FU{k,r} = FUold{k,r};
                            FSu{k,r} = FSuold{k,r};
                        else
                            [FU{k,r},Fsu1{k,r},~]=svd([FUold{k,r}*diag(FSuold{k,r}),FCU{k,r}*diag(FS{k,r})],0); Fstrunc{k,r}=diag(Fsu1{k,r});
                            Fp{k,r}=size([FUold{k,r},FCU{k,r}],2)-sum((cumsum(Fstrunc{k,r}(end:-1:1))/sum(Fstrunc{k,r}))<tauh);
                            FSuold{k,r} = Fstrunc{k,r};
                            Fku{k,r} = size([FUold{k,r},FCU{k,r}],2);
                        end
                        
                    end
                    
                    FVbasis{k,r} = FU{k,r}(:,1:Fp{k,r});
                    
                end
            end
            
            
            
            contnodes=contnodes+1;
            
            coretree;
            
            norme{contnodes} = normT(u{1}.core)^2 + normT(u{2}.core)^2 + normT(u{3}.core)^2;
            
            new_adiacenza(kk)=contnodes;
            
        end
        
        matrix_adiacenza(:,cont_ad)=new_adiacenza;
        cont_ad=cont_ad+1;
        
    end
    
    
   utree
    
    utreeold=utree;
    coretreeold = coretree;
    
    lengths(time+1)=cont;
    number_nodes=number_nodes+cont;
    


end


clear utreeold
clear utree

end