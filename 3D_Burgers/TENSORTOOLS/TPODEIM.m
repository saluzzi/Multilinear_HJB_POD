clear all
close all
clc
nnn  = [20,50,100,200];
svdtime = [];
for nn = 200 % Problem size (n in our notes)
    
    clear U1
    clear U2
    clear U3
    clear ff
    clear u
    clear s
    clear s1
    clear s2
    clear s3
    
nsnaps = 100; % Number of snapshots (m in our notes)


notrips = 5;

x3 = linspace(0.1,0.9,nn); x2 = linspace(0.1,0.9,nn); x1 = linspace(0.1,0.9,nn);

[X,Y,Z] = meshgrid(x1,x2,x3);

tf = 0.5;

% Nonlinear function (f(t_i) will be Y_i that we have in our notes)

set = '3';
    
    switch set
        
         case '1'
            
f = @(t) 1./sqrt( (X-t).^2 + (Y-t).^2 + (Z-t).^2 + 0.1.^2 );
ftxy = @(X,Y,Z,t)   1./sqrt( (X-t).^2 + (Y-t).^2 + (Z-t).^2 + 0.1.^2 );
        
         case '2'
            
f = @(t) 1./sqrt(0.001*X + (Y-3*t).^2 + (Z-t).^2 + 0.1.^2 );
ftxy = @(X,Y,Z,t)  1./sqrt(0.001*X + (Y-3*t).^2 + (Z-t).^2 + 0.1.^2 );
        
         case '3'
            
f = @(t) 1./sqrt((X-3*t).^2 + (Y-3*t).^2 + 0.001*Z + 0.1.^2 );
ftxy = @(X,Y,Z,t) 1./sqrt((X-3*t).^2 + (Y-3*t).^2 + 0.001*Z + 0.1.^2 );

        case '4'
            
f = @(t) 1./sqrt((X-t).^2 + 0.001*(Y) + (Z-t).^2 + 0.1.^2 );
ftxy = @(X,Y,Z,t)  1./sqrt(  (X-t).^2 + 0.001*(Y) + (Z-t).^2 + 0.1.^2 );

    end

time = 0;
h = tf/nsnaps;

U1 = [];
U2 = [];
U3 = [];

s1 = [];
s2 = [];
s3 = [];

kmax = 2*nsnaps;

tic
while time <=  tf
    
    fshot = f(time);
    fshot = tensor(fshot);
    
    [Uset,p,sset] = sthosvd(fshot,[notrips,notrips,notrips]);
    
    u1 = Uset.U{1};
    u2 = Uset.U{2};
    u3 = Uset.U{3};
    
%         for i =1:3
%         ff(:,:,i) = double(tenmat(fshot,i));
%         [u(:,:,i),s(:,:,i),~] = svd(ff(:,:,i),'econ');
%         end
%         
         U1 = [U1,u1];
         U2 = [U2,u2];
         U3 = [U3,u3];
%         
         s1 = [s1, sset(1,:)];
         s2 = [s2, sset(2,:)];
         s3 = [s3, sset(3,:)];
       
        
        time = time + h
        
        
        
     if size(s1,2) > kmax %% When the collected singular triplets become more than kmax, truncate and keep the top kmax
        
        
       [sdecrease1,id1] = sort(s1,'descend');
       s1 = sdecrease1(1:kmax);
       
       [sdecrease2,id2] = sort(s2,'descend');
       s2 = sdecrease2(1:kmax);
       
       [sdecrease3,id3] = sort(s3,'descend');
       s3 = sdecrease3(1:kmax);
       
       U1 = U1(:,id1);
       U1 = U1(:,1:kmax);
       
       U2 = U2(:,id2);
       U2 = U2(:,1:kmax);
       
       U3 = U3(:,id3);
       U3 = U3(:,1:kmax);
       
    end
    
    
    
end

[u1,ss1,v1] = svd(U1);
[u2,ss2,v2] = svd(U2);
[u3,ss3,v3] = svd(U3);

svdtime = [svdtime,toc]

end


figure
loglog(nnn,svdtime,'bo')
pause
figure
semilogy(diag(ss2),'bx')
hold on
semilogy(diag(ss1),'ro')
semilogy(diag(ss3),'ks')


leg = legend('Entries of $\overline{\bf S}_{1, {\cal F}}$','Entries of $\overline{\bf S}_{2, {\cal F}}$', 'Entries of $\overline{\bf S}_{3, {\cal F}}$')
set(leg,'Interpreter','latex')
set(gca,'FontSize',20)
pause

%%

% Generate the bases 
for tau = [1e-2,1e-4,1e-6];

p1=sum(diag(ss1)/ss1(1,1)>tau)
p2=sum(diag(ss2)/ss2(1,1)>tau)
p3=sum(diag(ss3)/ss3(1,1)>tau)


% Truncate the bases and perform DEIM in each dimension

[~,~,II]=qr(u1(:,1:p1)','vector'); II=II(1:p1)';
[~,~,JJ]=qr(u2(:,1:p2)','vector'); JJ=JJ(1:p2)';
[~,~,KK]=qr(u3(:,1:p3)','vector'); KK=KK(1:p3)';

Uapp1=u1(:,1:p1)/u1(II,1:p1);
Uapp2=u2(:,1:p2)/u2(JJ,1:p2);
Uapp3=u3(:,1:p3)/u3(KK,1:p3);

nt = 25

[X1,Y1,Z1] = meshgrid(x2(JJ),x1(II),x3(KK));

%%
figure
jj = 0;
for tt = linspace(0,tf,nt) 
    
    jj = jj+1;  
    
    FF=ftxy(X1,Y1,Z1,tt);
    FF = tensor(FF);
    Fapprox = ttm(FF, Uapp1, 1);
    [p1,p2,p3]
    Fapprox = ttm(Fapprox, Uapp2, 2); 
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
  

error  = sum(err)/nt;

[p1,p2,p3,error]
pause

end
