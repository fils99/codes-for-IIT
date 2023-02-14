close all %P1 parabolico
clear all
clc

global Nt L_q zita csi eta omega geom N_dof U_gdl  grad_phi  u phi_1 phi_2 phi_3 T_fin istanti
global D C R 
lung = 7;
lista_aree=zeros(lung,1);
for g = 1:lung
    lista_aree (g) = 0.02 * 1 / (2^(g-1));
end

lista_lati_max = lista_aree.^(1/2);
err_L2 = zeros(lung,1);
err_H1 = zeros(lung,1);
err_L_inf = zeros(lung,1);

numero_condiz_stiff = zeros(lung,1);


for r=1:lung
    area_max = lista_aree(r);

main


tot = zeros(Nt, 1);
tot_grad = zeros(Nt, 1);
for e=1:Nt

    % trovo coordinate vertici triangoli e i delta
    x1 =  geom.elements.coordinates (geom.elements.triangles (e, 1) , 1);
    x2 =  geom.elements.coordinates (geom.elements.triangles (e, 2) , 1);
    x3 =  geom.elements.coordinates (geom.elements.triangles (e, 3) , 1);
    y1 =  geom.elements.coordinates (geom.elements.triangles (e, 1) , 2);
    y2 =  geom.elements.coordinates (geom.elements.triangles (e, 2) , 2);
    y3 =  geom.elements.coordinates (geom.elements.triangles (e, 3) , 2);
    dx1 = x3 - x2;
    dx2 = x1 - x3;
    dx3 = x2 - x1;
    dy1 = y2 - y3;
    dy2 = y3 - y1;
    dy3 = y1 - y2;

    area = geom.support.TInfo(e).Area;

    % matrice e termine noto della trasformazione affine
    B = [dx2, -dx1; -dy2, dy1];
    B_inv = 1 / (2 * area) * [dy1, dx1; dy2, dx2];
    B_inv_tras = 1 / (2 * area) * [dy1, dy2; dx1, dx2];
    vertice_3 = [x3;y3];

    % nodi e pesi di quadratura
    [zita,csi,eta,omega] = int_nodes_weights(5,1);


    differ = 0;
    differ_grad = 0;
    for q=1: L_q % ciclo sui nodi di quadratura

        nodo_q_locale = [csi(q);eta(q)];

        % valuto le funzioni test nei nodi di quadratura

        phi_val_1 = phi_1 (csi(q),eta(q));
        phi_val_2 = phi_2 (csi(q),eta(q));
        phi_val_3 = phi_3 (csi(q),eta(q));
        phi = [phi_val_1; phi_val_2; phi_val_3];


        % tramite trasformaz affine, passo da triang di rif a triang fisico
        nodo_q_globale = vertice_3 + B * (nodo_q_locale);


        sol_ex_valutata = sol_esatta_bis(nodo_q_globale(1),  nodo_q_globale(2),T_fin); % sol esatta o sol esatta bis
        grad_sol_ex_valutato = grad_soluz_esatta_bis( nodo_q_globale(1),  nodo_q_globale(2),T_fin); % grad esatto o grad esatto bis

        

        sum_approx = 0;
        grad_approx = 0;
        for k = 1:3
            kk =  geom.elements.triangles(e,k);
                sum_approx = sum_approx + u(kk,istanti) * phi(k);
                grad_approx = grad_approx +   u(kk,istanti) * grad_phi(:,k);
        end
        differ = differ + omega(q) * ( sol_ex_valutata - sum_approx)^2;
        differ_grad = differ_grad + omega(q) * (grad_sol_ex_valutato - B_inv_tras * grad_approx)' * (grad_sol_ex_valutato - B_inv_tras * grad_approx) ;

    end % for q

    tot(e) = differ * 2 * area;
    tot_grad(e) = differ_grad * 2 * area;

end % for e

err_L2(r) = sqrt(sum(tot));
err_H1(r) = sqrt(sum(tot) +  sum(tot_grad));
err_L_inf(r) = max ( abs ( u(:,istanti) - u_ex(:,istanti) ) );

numero_condiz_stiff (r) = cond(D + C + R,2);


 
%  figure (4)
%  plot (log(lista_aree), (P_L2));
%  hold off
%  
%  figure(5)
%  plot (log(lista_aree),log(err_L2));
%  hold off

end

p_L2 = polyfit (log(lista_lati_max), log (err_L2), 1);
Poli_L2 = polyval(p_L2,log(lista_lati_max));
figure (4)
subplot(1,3,1)
plot(log(lista_lati_max), Poli_L2 , 'LineWidth', 2  )
hold on
scatter(log(lista_lati_max), log(err_L2) , 'LineWidth', 2  )
legend('Retta di regressione', 'Andamento reale')
title('Errore norma $L^2$','FontSize',18,'interpreter','latex')
xlabel('log(h)','FontSize',18,'interpreter','latex')
ylabel('log(err $L^2$)','FontSize',18,'interpreter','latex')
axis equal
hold off


p_H1 = polyfit (log(lista_lati_max), log (err_H1), 1);
Poli_H1 = polyval(p_H1,log(lista_lati_max));
subplot(1,3,2)
plot(log(lista_lati_max), Poli_H1, 'LineWidth', 2 )
hold on
scatter(log(lista_lati_max), log(err_H1) , 'LineWidth', 2  )
legend('Retta di regressione', 'Andamento reale')
title('Errore norma $H^1$','FontSize',18,'interpreter','latex')
xlabel('log(h)','FontSize',18,'interpreter','latex')
ylabel('log(err $H^1$)','FontSize',18,'interpreter','latex')
axis equal
hold off


p_L_inf = polyfit (log(lista_lati_max), log (err_L_inf), 1);
Poli_L_inf = polyval(p_L_inf,log(lista_lati_max));
subplot(1,3,3)
plot(log(lista_lati_max), Poli_L_inf , 'LineWidth', 2  )
hold on
scatter(log(lista_lati_max), log(err_L_inf) , 'LineWidth', 2  )
legend('Retta di regressione', 'Andamento reale')
title('Errore norma $L^{\infty}$','FontSize',18,'interpreter','latex')
xlabel('log(h)','FontSize',18,'interpreter','latex')
ylabel('log(err $L^{\infty}$)','FontSize',18,'interpreter','latex')
axis equal
hold off


p_cond = polyfit (log(lista_lati_max), log (numero_condiz_stiff), 1);
Poli_cond = polyval(p_cond,log(lista_lati_max));
figure(5)
plot(log(lista_lati_max), Poli_cond , 'LineWidth', 2  )
hold on
scatter(log(lista_lati_max), log(numero_condiz_stiff) , 'LineWidth', 2  )
legend('Retta di regressione', 'Andamento reale')
title('Numero di condizionamento matrice di rigidezza','FontSize',18,'interpreter','latex')
xlabel('log(h)','FontSize',18,'interpreter','latex')
ylabel('$log(cond_2 (A))$','FontSize',18,'interpreter','latex')
axis equal
hold off

salvataggio = round([lista_lati_max, err_L2, err_H1, err_L_inf, numero_condiz_stiff],7);
 