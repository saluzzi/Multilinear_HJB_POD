clear all
close all
clc
nn =100; % Problem size (n in our notes)
nsnaps = 200; % Number of snapshots (m in our notes)


nsing = 5;

x3 = linspace(0.1,0.9,nn); x2 = linspace(0.1,0.9,nn); x1 = linspace(0.1,0.9,nn);

[X,Y,Z] = meshgrid(x1,x2,x3);

tf = 0.5;

% Nonlinear function (f(t_i) will be Y_i that we have in our notes)

f = @(t) 1./sqrt( 0.001*(X) + 0.001*(Y) + (Z-t).^2 + 0.1.^2 );
ftxy = @(X,Y,Z,t)  1./sqrt( 0.001*(X) + 0.001*(Y) + (Z-t).^2 + 0.1.^2 );

time = 0;
h = tf/nsnaps;

U1 = [];
U2 = [];
U3 = [];

while time <=  tf
    
    fshot = f(time);
    
        for i =1:3
        ff(:,:,i) = double(tenmat(fshot,i));
        [u(:,:,i),~,~] = svds(ff(:,:,i),'econ');
        end
        
        U1 = [U1,u(:,1:2,1)];
        U2 = [U2,u(:,1:2,2)];
        U3 = [U3,u(:,1:2,3)];
       
        
        time = time + h
    
    
    
end

[u1,s1,v1] = svd(U1);
[u2,s2,v2] = svd(U2);
[u3,s3,v3] = svd(U3);

figure
semilogy(s1,'x')
hold on
semilogy(s2,'x')
semilogy(s3,'x')



% Generate the bases

tau = 1e-2;

p1=sum(diag(s1)/s1(1,1)>tau)
p2=sum(diag(s2)/s2(1,1)>tau)
p3=sum(diag(s3)/s3(1,1)>tau)


% Truncate the bases and perform DEIM in each dimension

[~,~,II]=qr(u1(:,1:p1)','vector'); II=II(1:p1)';
[~,~,JJ]=qr(u2(:,1:p2)','vector'); JJ=JJ(1:p2)';
[~,~,KK]=qr(u3(:,1:p3)','vector'); KK=KK(1:p3)';

Uapp1=u1(:,1:p1)/u1(II,1:p1);
Uapp2=u2(:,1:p2)/u2(JJ,1:p2);
Uapp3=u3(:,1:p3)/u3(KK,1:p3);

nt = 25

[X1,Y1,Z1] = meshgrid(x1(II),x2(JJ),x3(KK));

%%
figure
jj = 0;
for tt = linspace(0,tf,nt) 
    
    jj = jj+1;  
    
    FF=ftxy(X1,Y1,Z1,tt);
    FF = tensor(FF);
    Fapprox = ttm(FF, Uapp1, 2);
size(Fapprox)
    Fapprox = ttm(Fapprox, Uapp2, 1); 
    Fapprox = ttm(Fapprox, Uapp3, 3);
    
    ff = f(tt);
%     err(jj) = norm(double(tenmat(ff(:,:,1),1)) - double(tenmat(Fapprox(:,:,1),1)))/norm(double(tenmat(ff(:,:,1),1)))
    
    T = ff - Fapprox;
   
    subplot(1,2,1)
    mesh(double(squeeze(ff(:,10,:))))
    subplot(1,2,2)
    mesh(double(squeeze(Fapprox(:,10,:))))
    pause(0.1)
    T = tensor(T);
    n = normT(T);
    
    Tf = tensor(ff);
    vf = reshape(Tf.data, numel(Tf.data), 1);
    nf = norm(vf);
    
    err(jj) = n/nf;
   
      
end
    

error  = sum(err)/nt
