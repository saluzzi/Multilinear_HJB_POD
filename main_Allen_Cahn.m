%% OPTIMAL CONTROL OF THE ALLEN CAHN EQUATION VIA TSA-HO-POD-DEIM APPROACH

close all;
clear

isLoadOffline = 1; % fixing 1 it skips the offline stage loading the precomputed basis
isLoadOnline = 1; % fixing 1 it gives directly the result

if ~isLoadOnline
    
tau = 1e-3; % tolerance for HO-POD-DEIM
nn = 601; % number of grid points per dimension
ts = 5; % number of time steps
umin=-2; umax=0; % minimum and maxumum control
kappa = 10; % number of singular values considered per snapshot
nmax = ts;
xx =linspace(-1,1,nn+1)'; 
y = linspace(-1,1,nn+1);
[XX,YY] = meshgrid(xx,y);
nh =  nn+1;
n = nh^2; % This will be the size of the matrx solution

% CONSTRUCTION OF THE DIFFUSION MATRIX
x = xx(2:end-1);
Nx = nh; 
dx = 1/(Nx);
e = ones(Nx, 1);
A_1D = spdiags([e, -2*e, e], -1:1, Nx, Nx);
sigma = 0.1;
B = zeros(nh); B(1,1) = 2; B(1,2) = -1/2; B(end,end-1) = -1/2; B(end,end) = 2;
B = (2/3)*B;
B = sparse(B);
A = 1/dx^2*sigma*(A_1D+B);
hh = 1/ts;
[redvA, redDA] = eig(full(A));
[redvB, redDB] = eig(full(A'));
matpde.Aleft = speye(nn+1) - hh*A   ;
matpde.Aright = -hh*A;

[X,Y] = meshgrid(xx,y);
matpde.X0 = 2+cos(2*pi*X).*cos(2*pi*Y); % initial condition

matpde.T= 1; % final time
matpde.h = hh;
M = 2; % we start with 2 discrete controls
matpde.control = linspace(umin,umax,2); 
matpde.dx = dx;
matpde.gamma = 0.01;

matpde.redvA = redvA;
matpde.redvB = redvB;
matpde.redvAi = inv(redvA);
matpde.redvBi = inv(redvB);
matpde.redDA = redDA;
matpde.redDB = redDB;
matpde.nmax = nmax;
% dd1 = diag(matpde.redDA);
% dd2 = diag(matpde.redDB);
% nn1 = size(matpde.redvA,1);
% nn2 = size(matpde.redvB,2);
% L = zeros(nn1,nn2);
% for i = 1:nn1
%     for j = 1:nn2
%         L(i,j) = 1/((1-hh*dd1(i)) - (hh*dd2(j)));
%     end
% end
% matpde.L = L;

F = @(u) u - u.^3; % nonlinearity

%% OFFLINE STAGE

if ~isLoadOffline

%%% CONSTRUCTION OF THE TREE FOR THE OFFLINE STAGE
    
[~,uset,vset,usetnl,vsetnl]=full_tree_AC(matpde,F,tau,kappa);

%%% CREATE BASIS MATRICES FOR LINEAR AND NONLINEAR PARTS

[Vleftdeim,Sleftdeim,~] = svd(usetnl,'econ');
[Wrightdeim,Srightdeim,~] = svd(vsetnl,'econ');

[Vleft,Sleft,~] = svd(uset,'econ');
[Wright,Sright,~] = svd(vset,'econ');


%%% TRUNCATE BASIS MATRICES
b1 = sum(diag(Sleft)/Sleft(1,1) > tau);
b2 = sum(diag(Sright)/Sright(1,1) > tau);
VL = Vleft(:,1:b1);
WR = Wright(:,1:b2);

%%% CREATE DEIM PROJECTION

p1 = sum(diag(Sleftdeim)/Sleftdeim(1,1) > tau);
p2 = sum(diag(Srightdeim)/Srightdeim(1,1) > tau);


% SELECTION VECTOR

[~,~,II]=qr(Vleftdeim(:,1:p1)','vector'); II=II(1:p1)';
[~,~,JJ]=qr(Wrightdeim(:,1:p2)','vector'); JJ=JJ(1:p2)';

%  COMPUTE LEFT AND RIGHT BASIS

VLD=Vleftdeim(:,1:p1)/Vleftdeim(II,1:p1);
WRD=Wrightdeim(:,1:p2)/Wrightdeim(JJ,1:p2);

% CREATE REDUCED MATRICES

redA = VL'*A*VL;
redB = WR'*A*WR;

[redvA, redDA] = eig(redA);
[redvB, redDB] = eig(redB');

redpde.Dleft = VL'*VLD;
redpde.Dright = WR'*WRD;
X0 = matpde.X0;
redpde.X0 = VL'*X0*WR;

else
    load offline_AC.mat
    redpde.A = redA;
    redpde.B = redB';
end

%% ONLINE STAGE

redpde.T= matpde.T;
redpde.h = hh;
Mnew = 3; % we can consider more controls online
redpde.control = linspace(umin,umax,Mnew);
redpde.VII = VL(II,:);
redpde.WJJ = WR(JJ,:);
redpde.redvA = redvA;
redpde.redvB = redvB;
redpde.redvAi = inv(redvA);
redpde.redvBi = inv(redvB);
redpde.redDA = redDA;
redpde.redDB = redDB;


dd1 = diag(redpde.redDA); dd2 = diag(redpde.redDB);
nn1 = size(redpde.redDA,1);
nn2 = size(redpde.redDB,2);
% Lred = zeros(nn1,nn2);
% for ii = 1:size(redpde.X0,1)
%     for jj = 1:size(redpde.X0,2)
%         Lred(ii,jj) = 1/((redpde.h*dd1(ii)) + (redpde.h*dd2(jj)));
%     end
% end
% redpde.L = Lred;

[norme1,lengths1,matrix_adiacenza,nodes]=full_tree_red_AC(redpde,F);
lengths = lengths1;
norme = norme1;

[V,alfaoptimal_euler,optimal_dyn,optimal_dyn_red,costfunctional_euler,costfunctional] = value_function(matpde,redpde,norme,matrix_adiacenza,F,lengths,VL,WR);
V0 = V(1);
%% APPLICATION OF THE STATISTICAL PRUNING

redpde.nt = length(0:redpde.h:matpde.T);
rho = 0.2; % pruning ratio
M_vec = [5 9 17 33]; % vector of the number of controls considered in the iteration
len_vec = []; % vector of the cardinality of the trees
V_vec = []; % vector of the total costs
N = length(M_vec);
fprintf('\n Iteration: \n')
 
for i=1:N
   
    fprintf('..%d..',i)
    M = M_vec(i);
    redpde.control = linspace(umin,umax,M);
    tol_flag = 0;
    [nodes_flag,lengths_flag] = refiniment_tree(redpde,nodes,lengths,V,rho,tol_flag,costfunctional_euler);
    [norme_label,lengths_label,matrix_adiacenza_label,nodes_label] = prun_tree_statistical(redpde,F,nodes_flag,lengths_flag,tol_flag);
    new = matpde.dx^2*norme_label;
    nt = length(0:redpde.h:matpde.T);
    len = length(norme_label);
    full_size = (M^(nt+1)-1)/(M-1);
    [V_label,alfaoptimal_euler_label,optimal_dyn_label,optimal_dyn_red_label,costfunctional_euler_label,costfunctional_label] = value_function(matpde,redpde,norme_label,matrix_adiacenza_label,F,lengths_label,VL,WR);
    len_vec = [len_vec length(norme_label)];
    V_vec = [V_vec V_label(1)];
    nodes = nodes_label;
    lengths = lengths_label;
    V = V_label;
    costfunctional_euler = costfunctional_euler_label;
   
end
else

    load final_results_AC.mat
    
end

%% PLOTS OF THE RESULTS

semilogy( M_vec(1:N), len_vec(1:N),'-s','LineWidth',1)
 xlabel('# controls')
title('Cardinality of the tree')
grid

figure
plot([2 M_vec(1:N)], [V0 V_vec(1:N)],'-*','LineWidth',1)
 xlabel('# controls')
title('Total cost')
grid

figure
t_grid = 0:redpde.h:matpde.T;
plot(0:redpde.h:matpde.T-redpde.h,alfaoptimal_euler_label,'LineWidth',1)
grid
xlabel('t')
title('Optimal control')
for i = 1:floor(size(optimal_dyn_label,3)/2):size(optimal_dyn_label,3)
    figure
    mesh(XX,YY,abs(real(optimal_dyn_label(:,:,i))))
    axis([-1 1 -1 1 0 3])
    colorbar
end


