
clear all
close all
clc


fulltime = [];
n = 50;  % Discrete dimension in each spatial direction. N = n^d
ts = 10; % Number of timesteps for constructing the HO-POD-DEIM bases
r = 100; % Reynold's number
array  = 'tensor';

%velocity parameters

alfa1=1; 
alfa2=1;

% Grid Setup

xx = linspace(0,1,n+1); 
yy = linspace(0,1,n+1);
zz = xx;

[x,y,z] = meshgrid(xx,yy,zz);


nh=  n+1;

% Matrix Setup

Nx = nh; dx = 1/(Nx);
e = ones(Nx, 1);
A_1D = spdiags([e, -2*e, e], -1:1, Nx, Nx);
AA = A_1D;
sigma = 1/r;


AA = 1/dx^2*sigma*(AA);

spa = speye(size(AA));
AAbig = sparse(kron(kron(AA,spa), spa) + kron(kron(spa,spa), AA) + kron(kron(spa,AA), spa));
ode.AAbig = AAbig;

for k = 1:3
    for r = 1:3
        A{k,r} = AA;
    end
end
ode.A = A;

BB = eye(Nx) - diag(ones(Nx-1,1),-1);

CC = spdiags([-1*ones(Nx,1), ones(Nx,1)], [-1,1], Nx, Nx);
% CC(1,2) = 0 ; CC(end,end-1) = 0;
for k = 1:3
    for r = 1:3
        B{k,r} = (1/(2*dx))*CC;
    end
end
ode.B = B;



%%% Initial condition

X0{1} = (1/10)*sin(2.*pi.*x).*sin(2.*pi.*y).*cos(2.*pi.*z);
X0{2} = (1/10)*sin(2.*pi.*x).*cos(2.*pi.*y).*sin(2.*pi.*z);
X0{3} = (1/10)*cos(2.*pi.*x).*cos(2.*pi.*y).*sin(2.*pi.*z);

for k = 1:3
    X0{k} = tensor(X0{k},[n+1,n+1,n+1]);
end


ode.X0 = X0;

burger = 1;
tf = 1;
hh = tf/ts;

ode.mu = 5;

ode.T= tf;
ode.h = hh;
ode.control = [-2,0];
ode.dx = dx;
ode.gamma = 0.1;
kappa = 20;
ode.tau = 1e-2;
ode.nn = n;

%%
fulltimer = tic;
switch array
    case 'tensor'
[norme,Vbasis,FVbasis,lengths_snap2,matrix_adiacenza_snap2]=integrate_and_basis_tensor(ode, kappa);
    case 'vector'
[norme,snaps,lengths_snap2,matrix_adiacenza_snap2]=integrate_and_basis_vector(ode);
end
fulltime = [fulltime,toc(fulltimer)];


%%
for k = 1:3
       
    Y0{k} = X0{k};

    for r = 1:3
        
        Y0{k} = ttm(Y0{k},Vbasis{k,r},r,'t');
    Fp{k,r} = size(FVbasis{k,r},2)
        p{k,r} = size(Vbasis{k,r},2)
        [~,~,II{k,r}]=qr(FVbasis{k,r}','vector');  II{k,r}=II{k,r}(1:Fp{k,r})';
        Deimapprox{k,r}=FVbasis{k,r}/FVbasis{k,r}(II{k,r},:);
        
        % Form the ROM by projection
        
        Ak{k,1} = Vbasis{k,1}'*A{k,1}*Vbasis{k,1};
        Ak{k,2} = Vbasis{k,2}'*A{k,2}*Vbasis{k,2};
        Ak{k,3} = Vbasis{k,3}'*A{k,3}*Vbasis{k,3};
        
        
        Deim{k,r} = Vbasis{k,r}'*Deimapprox{k,r};
        
        
    end
    
    
    [VAsmall{k,3},DAsmall{k,3}] = eig(full(Ak{k,3})'); %eigenvalue decomposition of symmetric A
    ddAsmall{k,3} = diag(DAsmall{k,3});
    VAismall{k,3} = inv(VAsmall{k,3});
    
end


%% SET UP Reduced Model



nt = 10;
h = 1/nt;
redode.Y0 =Y0; redode.Deim = Deim;
redode.Ak = Ak; redode.VAismall = VAismall;
redode.VAsmall = VAsmall; redode.DAsmall = DAsmall;
redode.B = ode.B;
redode.T= 1;
redode.h = 1/nt;
redode.control = linspace(-2,0,2);   
redode.mu = 5;
redode.Vbasis= Vbasis;
redode.FVbasis= FVbasis;
redode.II = II;
redode.tau = 1e-4;
redode.p = p;

tic
[norme,lengths,matrix_adiacenza]=red_tree_tensor(redode,0);
redtime = toc

ode.h = redode.h; %stepSize

new=ode.dx^2*norme;
nt = length(0:ode.h:ode.T);
V=containers.Map(nt,new);
index_control=containers.Map(nt,zeros(1,2));
control2dtRuu=ode.gamma*ode.h*redode.control.^2';

len=length(norme);
deltadtdx=ode.h*ode.dx^2;
tic
norme=deltadtdx*norme;
for i=nt-1:-1:1
    len=len-lengths(i+1);
    new2=zeros(1,len);
    new3=zeros(1,len);
    for j=1:len     %Qui sto scorrendo tutto il sotto albero, si potrebbe cambiare mettendo solo il livello i
        [app, new3(j)]=min(new(matrix_adiacenza(:,j))+control2dtRuu');
        new2(j)=app+norme(j);
    end
    index_control(i)=new3;
    V(i)=new2;
    new=new2;
end
toc;

%% COMPUTATION OPTIMAL CONTROL AND TRAJECTORY

alfaoptimal_euler=zeros(1,nt-1);
dynamics=zeros(1,nt);
dynamics(1)=1;
s=1;

X0 = ode.X0;

whos X0
optimal_dyn{1,:} = X0;
optimal_dyn_red{1,:} = Y0;
figure
[XX,YY] = meshgrid(xx,y);
fullvec = [];
redvec= [];

for j=1:nt-1


    vett=index_control(j);
    alfaoptimal_euler(j)=redode.control(vett(s));
    
    
    optimal_dyn{j+1,:} = ten_SI(ode,optimal_dyn{j},alfaoptimal_euler(j));
    optimal_dyn_red{j+1,:} = ten_SI_poddeim(redode,optimal_dyn_red{j},alfaoptimal_euler(j));
    
    
    
    for k = 1:3
        Xcell = optimal_dyn{j+1};
        Xapprox{k} = Xcell{k};
        
        for r = 1:3
            
            Xapprox{k} = ttm(Xapprox{k},(Vbasis{k,r}*Vbasis{k,r}'),r);
            
        end
    end
    
    for k = 1:3
        
        if j == 1 || j== 5 || j == nt-1
        
        figure(j)
        Xrealcell = Xapprox;
        TT =  double(tenmat(Xrealcell{2},1));
        mesh(abs(TT))
        axis([0 1000 0 30 0 0.7])

        end
        
        err(j,k)=  normT(Xrealcell{k} - Xapprox{k})/normT(Xrealcell{k});
        
        
    end
        
    %keyboard
    s=matrix_adiacenza(vett(s),s);
    dynamics(j+1)=s;
end




costfunctional_euler=zeros(1,nt);
costfunctional_euler(1)=V(1);
costfunctional=zeros(1,nt);
Xlargecell = optimal_dyn_red{end};
costfunctional(end)= ode.dx^2 * (normT(Xlargecell{1})^2 + normT(Xlargecell{2})^2 + normT(Xlargecell{3})^2);
for time=2:nt
    
    app=V(time);
    costfunctional_euler(time)=app(dynamics(time));
    Xlargecell = optimal_dyn_red{nt-time+1};
    costfunctional(nt-time+1) = deltadtdx*(normT(Xlargecell{1}).^2 + normT(Xlargecell{2}).^2 + normT(Xlargecell{3}).^2)+costfunctional(nt-time+2) +ode.gamma*ode.h*alfaoptimal_euler(nt-time+1).^2;
    
end
%% Compare the cost functions obtained by the reduced model and the full model


figure
plot(linspace(0,1,length(costfunctional_euler)),costfunctional_euler,'b-*')
hold on
plot(linspace(0,1,length(costfunctional_euler)),costfunctional,'r--')












