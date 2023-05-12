function [norme,snaps,lengths,matrix_adiacenza]=integrate_and_basis_vector(ode)


X0 = ode.X0;
AAbig = ode.AAbig;

B = ode.B;
tau = ode.tau;
addpath '/Users/gk/Dropbox/Sylvester_Image_Project/DEIM/Progress_Deim/DEIM_FOR_SYSTEMS/Systems_Matlab/TENSORTOOLS'

tauh = tau;
control = ode.control;
h = ode.h;
tf = ode.T;
nt = length(0:h:tf);

na = length(control);

nn = ode.nn+1;

Adiscrete = speye(size(AAbig)) - h*AAbig;

%%% initiate tree storing and first basis vectors %%%

for k = 1:3
    
    for r = 1:3
        Fnew{k,r} = ttm(X0{k},B{r},r).*X0{r};
    end
    
    
    
    Fnew{k} = Fnew{k,1}+Fnew{k,2}+Fnew{k,3};
end


for k = 1:3
    
    snaps{k}=X0{k}(:);

end

norme{1} = norm(snaps{1},2)^2 + norm(snaps{2},2)^2 + norm(snaps{3},2)^2;


% matrix_adiacenza=zeros(na,(na^(nt-1)-1)/(na-1));
number_nodes=1;

lengths=zeros(1,nt);
lengths(1)=1;
new_adiacenza=zeros(na,1);
contnodes=1;

t_grid = 0;
cont_ad=1;

opts.type='ict'; opts.droptol=1e-4;
L=ichol(Adiscrete,opts);

for time=1:nt-1
    time

    t_grid = [t_grid(end), t_grid(end) + h];
    s1=na^(time-1);
    
    cont=0;
    
    for j=1:s1
        
        for kk=1:na
            
            
            for k = 1:3
                
                yold = snaps{k}(:,na^(time-1) + j - 1);
                rhs = (1+h*control(kk))*yold + h*Fnew{k}(:);
                [ynew,flag] = pcg(Adiscrete,rhs,1e-4,100,L,L');
                snaps{k}(:,contnodes) = ynew;
                

                Xnew{k} = tensor(reshape(ynew,nn,nn,nn));
                
            end
            
            for k = 1:3
                
                for r = 1:3
                    Fnew{k,r} = ttm(Xnew{k},B{r},r).*Xnew{r};
                end
                
                
                
                Fnew{k} = Fnew{k,1}+Fnew{k,2}+Fnew{k,3};
            end
            
            norme{contnodes+1} = norm(snaps{1}(:,contnodes),2)^2 + norm(snaps{2}(:,contnodes),2)^2 + norm(snaps{3}(:,contnodes),2)^2;
            
            

            
            
            
            contnodes=contnodes+1;
            new_adiacenza(kk)=contnodes;
            
        end
        
        matrix_adiacenza(:,cont_ad)=new_adiacenza;
        cont_ad=cont_ad+1;
        
    end
    
   
    
    lengths(time+1)=cont;
    number_nodes=number_nodes+cont;
    
end
whos

end