% se calcoli errore, commenta le prime 3 righe

% clc
% clear all
% close all
% 
% area_max = 0.002;

global epsilon  geom  u_d  b_tot beta gamma  U_gdl u phi_1 phi_2 phi_3 phi_4 phi_5 phi_6    
global grad_phi_1 grad_phi_2 grad_phi_3 grad_phi_4 grad_phi_5 grad_phi_6 
global hess_xx hess_yx hess_xy hess_yy A_D_tot  A_tot_lumping

Sample_Square_Dirichlet 
P2 

% assicurarsi che questi coefficienti rispettino la coercività della forma
% bilineare
epsilon = @ (x,y) 10^(0); 
beta = @ (x,y) [0;0] ;
gamma = @ (x,y) 10000 ;

% funzioni test, gradienti ed hessiane
phi_1 = @(x,y) 2 * x * (x - 1/2);
phi_2 = @(x,y) 2 * y * (y - 1/2); 
phi_3 = @(x,y) 2 * (1 - x - y) * (1- x - y - 1/2);
phi_4 = @(x,y) 4 * x * y;
phi_5 = @(x,y) 4 * y * (1 - x - y);
phi_6 = @(x,y) 4 * x * (1 - x - y);

grad_phi_1 = @(x,y) [4 * x - 1;0]; 
grad_phi_2 = @(x,y) [0; 4 * y - 1];
grad_phi_3 = @(x,y) [4 * x + 4 * y - 3;4 * x + 4 * y - 3];
grad_phi_4 = @(x,y) [4 * y; 4 * x];
grad_phi_5 = @(x,y) [- 4 * y; 4 * (1 - x - 2 * y)];
grad_phi_6 = @(x,y) [4 * (1 - 2 * x - y); - 4 * x];

phi_1_xx = 4;
phi_1_xy = 0;
phi_1_yx = 0;
phi_1_yy = 0;

phi_2_xx = 0;
phi_2_xy = 0;
phi_2_yx = 0;
phi_2_yy = 4;

phi_3_xx = 4;
phi_3_xy = 4;
phi_3_yx = 4;
phi_3_yy = 4;

phi_4_xx = 0;
phi_4_xy = 4;
phi_4_yx = 4;
phi_4_yy = 0;

phi_5_xx = 0;
phi_5_xy = -4;
phi_5_yx = -4;
phi_5_yy = -8;

phi_6_xx = -8;
phi_6_xy = -4;
phi_6_yx = -4;
phi_6_yy = 0;

hess_xx = [phi_1_xx; phi_2_xx; phi_3_xx; phi_4_xx; phi_5_xx; phi_6_xx];
hess_yx = [phi_1_yx; phi_2_yx; phi_3_yx; phi_4_yx; phi_5_yx; phi_6_yx];
hess_xy = [phi_1_xy; phi_2_xy; phi_3_xy; phi_4_xy; phi_5_xy; phi_6_xy];
hess_yy = [phi_1_yy; phi_2_yy; phi_3_yy; phi_4_yy; phi_5_yy; phi_6_yy];


assembla_matrice ()

% Np = geom.nelements.nVertexes;

% U_gdl = pcg(A_tot_lumping, b_tot - A_D_tot * u_d',1.0E-5,100);
 U_gdl = A_tot_lumping \ (b_tot - A_D_tot * u_d');

u = zeros (nnode,1);


for j = 1 : nnode
    jj = geom.pivot.pivot(j);
    if jj > 0
        u (j) = U_gdl(jj);
    elseif jj<0
        u(j) = u_d(-jj);
    end
end

u_ex = sol_esatta_bis (geom.elements.coordinates(:,1), geom.elements.coordinates(:,2)); % sol esatta oppure sol esatta bis

if r < lung 
    close(figure(1))
end

if r == lung

    % plot soluzione approssimata
    figure(2)
    plotta_soluz(geom.elements.triangles(:,1:3), geom.elements.coordinates(:,1), geom.elements.coordinates(:,2), u )
    hold on
    title ('soluzione approssimata')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold off

    % plot soluzione esatta

    figure(3)
    plotta_soluz(geom.elements.triangles(:,1:3), geom.elements.coordinates(:,1), geom.elements.coordinates(:,2), u_ex )
    hold on
    title ('soluzione esatta')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold off

end



% L_Inf =  max ( abs ( u - u_ex ) )





