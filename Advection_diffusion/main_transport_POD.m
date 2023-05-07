%% OPTIMAL CONTROL OF THE CONVECTION DIFFUSION EQUATION VIA TSA-POD APPROACH

clear
close all

a = -5; b = 5; % domain [a,b]
nx = 161; % number of grid points per dimension
c1 = 1; % convection paramters
c2 = 0;
sigma = 0; % diffusion parameter
PDE.dt_snap = 0.05; % time step for the snapshots 
PDE.dt_mor = 0.05; % time step for the online phase
T = 1; % Final time
umin = -3; umax = -1; % control set [-umin,umax]
t = 0:PDE.dt_snap:T;
PDE.nt_snap = length(t);
x = linspace(a,b,nx);
dx = x(2)-x(1);
sigmadx2=  sigma/dx^2;
[X,Y] = meshgrid(x,x);
nx2 = nx*nx;
tau = 1.e-4; % tolerance for the projection error
PDE.delta2 = 1; % coefficient for the cost functional
PDE.delta = 1;
PDE.gamma = 0.1;
M = 3; % number of controls for the online phase

ic = @(x,y) max(2-(x).^2-(y).^2,0); % initial condition
ic2 = ic(X,Y);
ic3 = reshape(ic2,nx2,1);
PDE.ic=ic3;

tic
% Construction of the matrices for the approximation of the dynamical
% system
v1 = ones(nx2,1)*(-(c1+c2)/dx-4*sigmadx2);
v2 = (c2/dx+sigmadx2)*ones(nx2-1,1);
v3 = (sigmadx2)*ones(nx2-1,1);
v4 = (sigmadx2+c1/dx)*ones(nx2-nx,1);
v5 = (sigmadx2)*ones(nx2-nx,1);
A = diag(v1)+diag(v2,-1)+diag(v3,1)+diag(v4,-nx)+diag(v5,nx);
PDE.A = A;
PDE.A2 = sparse(eye(size(A))-PDE.dt_snap*PDE.A);
PDE.A3 = sparse(eye(size(A)));

%% CONSTRUCTION OF THE SNAPSNOT MATRIX

PDE.control_snap = [umin umax];
PDE.na_snap=length(PDE.control_snap);
[nodes_snap,lengths_snap2,matrix_adiacenza_snap2] = tree_bilinear(PDE);

%% CONSTRUCTION OF THE POD BASIS AND REDUCED MATRICES

[Psi2,U,~] = svd(nodes_snap,'econ');
r = sum(diag(U)/U(1,1) > tau);
disp(['Number of the reduced basis: ',num2str(r)])
Psi = Psi2(:,1:r);
PDE.Apod = Psi'*A*Psi;
PDE.icred = Psi'*PDE.ic;
offline_POD = toc;
disp(['Offline phase CPU time: ',num2str(offline_POD)])

%% ONLINE PHASE

PDE.na_mor = M; % number of discrete controls for the online phase
PDE.control_mor = linspace(umin,umax,PDE.na_mor);
PDE.nt_mor = length(0:PDE.dt_mor:T);
PDE.tol = PDE.dt_mor^2/M^2; % tolerance for the pruning

tic
[nodes_POD2,lengths_POD2,matrix_adiacenza_POD2] = tree_pod_pruning(PDE);
num = size(nodes_POD2,2);
cardinality_tree = num;
disp(['The cardinality of the pruned tree: ',num2str(num)])

norme = zeros(1,num);
for i = 1:num
    norme(i) = sum((Psi*nodes_POD2(:,i)).^2);
end

%% RESOLUTION OF THE HJB EQUATION

lengths = lengths_POD2;
matrix_adiacenza = matrix_adiacenza_POD2;

new = dx^2*PDE.delta2*norme;
V = containers.Map(PDE.nt_mor,new);
index_control = containers.Map(PDE.nt_mor,zeros(1,2));
control2dtRuu = PDE.gamma*PDE.dt_mor*PDE.control_mor.^2';
len = num;
deltadtdx = PDE.delta*PDE.dt_mor*dx^2;
norme = deltadtdx*norme;
for i = PDE.nt_mor-1:-1:1
    len = len-lengths(i+1);
    new2 = zeros(1,len);
    new3 = zeros(1,len);
    for j = 1:len
        [app, new3(j)] = min(new(matrix_adiacenza(:,j))+control2dtRuu');
        new2(j) = app+norme(j);
    end
    index_control(i) = new3;
    V(i) = new2;
    new = new2;
end
online_POD = toc;
disp(['Online phase CPU time: ',num2str(online_POD)])

%% COMPUTATION OPTIMAL CONTROL AND TRAJECTORY

alfaoptimal_euler = zeros(1,PDE.nt_mor-1);
dynamics = zeros(1,PDE.nt_mor);
dynamics(1) = 1;
s = 1;
for j = 1:PDE.nt_mor-1
    vett = index_control(j);
    alfaoptimal_euler(j) = PDE.control_mor(vett(s));
    s = matrix_adiacenza(vett(s),s);
    dynamics(j+1) = s;
end

opt_dyn = nodes_POD2(:,dynamics(:));
opt_dyn_euler = Psi*opt_dyn;

costfunctional_euler = zeros(1,PDE.nt_mor);
costfunctional_euler(1) = V(1);

for time = 2:PDE.nt_mor
    app = V(time);
    costfunctional_euler(time) = app(dynamics(time));
end

%% PLOT OF THE RESULTS

figure

minimum = min(min(opt_dyn_euler));
maximum = max(max(opt_dyn_euler));
for i = 1:floor(PDE.nt_mor/20):PDE.nt_mor
    mesh(X,Y,reshape(opt_dyn_euler(:,i),nx,nx))
    zlim([minimum maximum])
    drawnow
    pause(.3)
end