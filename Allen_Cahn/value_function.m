function [V,alfaoptimal_euler,optimal_dyn,optimal_dyn_red,costfunctional_euler,costfunctional] = value_function(matpde,redpde,norme,matrix_adiacenza,F,lengths,VL,WR)

% Given the tree structure, this routine computes the value function on the
% tree, the optimal trajectory and the optimal control

new = matpde.dx^2*norme;
nt = length(0:redpde.h:matpde.T);
V = containers.Map(nt,new);
index_control = containers.Map(nt,zeros(1,2));
control2dtRuu = matpde.gamma*redpde.h*redpde.control.^2';

len = length(norme);
deltadtdx = redpde.h*matpde.dx^2;
norme = deltadtdx*norme;
for i = nt-1:-1:1
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

alfaoptimal_euler = zeros(1,nt-1);
dynamics = zeros(1,nt);
dynamics(1) = 1;
s = 1;

[nx,ny] = size(matpde.X0);
optimal_dyn = zeros(nx,ny,nt);
optimal_dyn(:,:,1) = matpde.X0;
optimal_dyn_red(:,:,1) = VL'*matpde.X0*WR;

matpde.h2 = redpde.h;
for j = 1:nt-1
    vett = index_control(j);
    alfaoptimal_euler(j) = redpde.control(vett(s));
    
    F1 = F(optimal_dyn(:,:,j));
    F2 = F(optimal_dyn_red(:,:,j));
    optimal_dyn(:,:,j+1) = mat_SI(matpde,optimal_dyn(:,:,j),F1,alfaoptimal_euler(j));
    optimal_dyn_red(:,:,j+1) = mat_SI(redpde,optimal_dyn_red(:,:,j),F2,alfaoptimal_euler(j));
  
    s = matrix_adiacenza(vett(s),s);
    dynamics(j+1) = s;
end


costfunctional_euler = zeros(1,nt);
costfunctional_euler(1) = V(1);
costfunctional = zeros(1,nt);
costfunctional(end) = matpde.dx^2 * norm(optimal_dyn(:,:,end),'fro')^2;
for time = 2:nt
    app = V(time);
    costfunctional_euler(time) = app(dynamics(time));
    costfunctional(nt-time+1) = deltadtdx* norm(optimal_dyn(:,:,nt-time+1),'fro').^2+costfunctional(nt-time+2) +matpde.gamma*redpde.h*alfaoptimal_euler(nt-time+1).^2;
end

end
