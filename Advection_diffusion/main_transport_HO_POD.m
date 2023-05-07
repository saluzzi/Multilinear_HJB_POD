%% OPTIMAL CONTROL OF THE CONVECTION DIFFUSION EQUATION VIA TSA-HO-POD APPROACH

close all;
clear

a = -5; b = 5; % domain [a,b]
nn = 161; % number of grid points per dimension
hh = 0.05;
c1 = 2; % convection paramters
c2 = 1;
sigma = 1; % diffusion parameter
%hh = 0.05; % time step for the snapshots
redpde.h = 0.05; % time step for the online phase
matpde.T = 1; % Final time
ts = floor(matpde.T/hh);
umin = -3; umax = -1; % control set [-umin,umax]
x = linspace(a,b,nn);
dx = x(2)-x(1);
sigmadx2=  sigma/dx^2;
[X,Y] = meshgrid(x,x);
n = nn*nn;
tau = 1.e-4; % tolerance for the projection error
kappa = 50; % maximum number of reduced basis
matpde.delta2 = 1; % coefficient for the cost functional
matpde.delta1 = 1;
matpde.gamma = 0.1;
M = 3; % number of controls for the online phase

ic = @(x,y) max(2-(x).^2-(y).^2,0); % initial condition
matpde.X0 = ic(X,Y);

tic
% Construction of the matrices for the approximation of the dynamical
e = ones(nn, 1);
A_1D = spdiags([e, -2*e, e], -1:1, nn, nn);
sigma = 0;
A = (1/dx^2)*sigma*(A_1D);
BB = eye(nn) - diag(ones(nn-1,1),-1);
BB = (1/(dx))*BB;
matpde.B = BB;
A1 = A - diag(c2*ones(nn,1))*BB;
A2 = A - diag(c1*ones(nn,1))*BB;

matpde.h = hh;
matpde.control = linspace(umin,umax,2);
matpde.dx = dx;
matpde.A1 = A1;
matpde.A2 = A2;
Xold = matpde.X0;
matpde.A11=eye(size(A1))-(hh)*A1;
matpde.A22=-(hh)*A2';


%% TREE CREATION FULL USING MATRICES WITH LR


[uset,vset] = full_tree_transport_bilinear(matpde,tau,kappa);

%% CONSTRUCTION OF THE POD BASIS AND REDUCED MATRICES

VL = uset;
WR = vset;

disp(['Number of the reduced basis: [',num2str(size(VL,2)),', ',num2str(size(WR,2)), ']'])

% Create reduced matrices
redA = VL'*A1*VL;
redB = WR'*A2*WR;
redpde.X0 = VL'*matpde.X0*WR;
offline_H0_POD = toc;
disp(['Offline phase CPU time: ',num2str(offline_H0_POD)])


%% ONLINE PHASE



redpde.T= matpde.T;
redpde.control = linspace(umin,umax,M);
redpde.t_grid=0:redpde.h:redpde.T;
redpde.redA = redA;
redpde.redB = redB;
redpde.VL = VL;
redpde.WR = WR;
redpde.redA1=eye(size(redA))-redpde.h*redA;
redpde.redB1=- redpde.h*redB';

[redpde.VA, redpde.DA] = eig(redA); dd1 = diag(redpde.DA);
[redpde.VB, redpde.DB] = eig(redB'); dd2 = diag(redpde.DB);

redpde.VAi = inv(redpde.VA);
redpde.VBi = inv(redpde.VB);

L = zeros(size(redA,2),size(redB,2));
for ii = 1:size(redA,2)
    for jj = 1:size(redB,2)
        L(ii,jj) = 1/((1-hh*dd1(ii)) - (hh*dd2(jj)));
    end
end

redpde.L = L;

tic
[nodes,norme1,lengths1,matrix_adiacenza] = full_tree_pruned_transport(redpde);
lengths = lengths1;
norme = norme1;
online_2s = toc;
len = length(norme);
cardinality_s = len;
disp(['The cardinality of the pruned tree: ',num2str(len)])

%% RESOLUTION OF THE HJB EQUATION

new = matpde.delta2*matpde.dx^2*norme;
nt = length(0:redpde.h:matpde.T);
V = containers.Map(nt,new);
index_control = containers.Map(nt,zeros(1,2));
control2dtRuu = matpde.gamma*matpde.h*redpde.control.^2';

len = length(norme);
deltadtdx = matpde.h*matpde.dx^2;
norme = matpde.delta1*deltadtdx*norme;
for i = nt-1:-1:1
    len = len-lengths(i+1);
    new2 = zeros(1,len);
    new3 = zeros(1,len);
    for j=1:len
        [app, new3(j)] = min(new(matrix_adiacenza(:,j))+control2dtRuu');
        new2(j) = app+norme(j);
    end
    index_control(i) = new3;
    V(i) = new2;
    new = new2;
end
online_H0_POD = toc;
disp(['Online phase CPU time: ',num2str(online_H0_POD)])

%% COMPUTATION OPTIMAL CONTROL AND TRAJECTORY

alfaoptimal_euler = zeros(1,nt-1);
dynamics = zeros(1,nt);
dynamics(1) = 1;
s = 1;
[nx,ny] = size(matpde.X0);
optimal_dyn = zeros(nx,ny,nt);
optimal_dyn(:,:,1) = matpde.X0;

optimal_dyn_red(:,:,1) = VL'*matpde.X0*WR;
fullvec = [];
redvec = [];
matpde.h = hh;
for j = 1:nt-1
    vett = index_control(j);
    alfaoptimal_euler(j)=redpde.control(vett(s));

    xnew = mat_SI_bil(matpde,optimal_dyn(:,:,j),alfaoptimal_euler(j));
    optimal_dyn(:,:,j+1) = xnew;

    s = matrix_adiacenza(vett(s),s);
    dynamics(j+1) = s;
end
%% PLOT OF THE RESULTS

figure

minimum = min(min(min(optimal_dyn)));
maximum = max(max(max(optimal_dyn)));
for i = 1:floor(nt/20):nt
    mesh(X,Y,optimal_dyn(:,:,i))
    zlim([minimum maximum])
    drawnow
    pause(.3)
end

