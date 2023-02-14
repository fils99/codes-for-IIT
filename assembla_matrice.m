function assembla_matrice () % P1

global Nt Ndof  A_tot b epsilon geom N_dri A_D u_d b_N b_tot beta gamma L_q grad_phi 
global phi_1 phi_2 phi_3 area R_tot_lumping A_D_tot R_D_tot A_tot_lumping D_tot Pe_h 
global tau_h R_tot_lumping D C R
global istanti T_fin t delta_t M M_D term_noto

Ndof = max(geom.pivot.pivot);
N_dri = -min (geom.pivot.pivot);

D = zeros (Ndof);
C = zeros (Ndof);
R = zeros (Ndof);
M = zeros (Ndof);


D_stab = zeros (Ndof);
C_stab = zeros (Ndof);
R_stab = zeros (Ndof);


D_D = zeros (Ndof, N_dri);
C_D = zeros (Ndof, N_dri);
R_D = zeros (Ndof, N_dri);
M_D = zeros (Ndof, N_dri);

D_D_stab = zeros (Ndof, N_dri);
C_D_stab = zeros (Ndof, N_dri);
R_D_stab = zeros (Ndof, N_dri);


b = zeros (Ndof,istanti);
b_stab = zeros (Ndof,istanti);

m_k = 1/3;

Nt = geom.nelements.nTriangles;  % numero triangoli



for e = 1 : Nt
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


    area = geom.support.TInfo(e).Area; % area triangolo
    % area_bis = 1 / 2 * (dx2 * dy1 - dy2 * dx1);

    h_1 = sqrt (dx1^2 + dy1^2);
    h_2 = sqrt (dx2^2 + dy2^2);
    h_3 = sqrt (dx3^2 + dy3^2);
    h = [h_1; h_2; h_3];

    h_E = max(h); % lato max triangolo


    % coordinate baricentro triangolo
    x_bar = geom.support.TInfo(e).CG(1); 
    y_bar = geom.support.TInfo(e).CG(2);

    beta_bar = beta (x_bar, y_bar);
    epsilon_bar = epsilon (x_bar, y_bar);
    gamma_bar = gamma (x_bar, y_bar);
   

    Pe_h = m_k * norm(beta_bar) * h_E / (2 * epsilon_bar);

    Se_h = gamma_bar * h_E / (6 * epsilon_bar);

    if  (0 <= Pe_h) && (Pe_h <= 1)
        tau_h = m_k * h_E ^2 / (4 * epsilon_bar);
    elseif Pe_h > 1
        tau_h = h_E / (2* norm(beta_bar));
        
    end

    % nodi e pesi di quadratura
    [zita,csi,eta,omega] = int_nodes_weights(5,1);

    % matrice e termine noto della trasformazione affine 
    B = [dx2, -dx1; -dy2, dy1];
    B_inv = 1 / (2 * area) * [dy1, dx1; dy2, dx2];
    B_inv_tras = 1 / (2 * area) * [dy1, dy2; dx1, dx2];
    vertice_3 = [x3;y3];
 
     
    % raccolgo i delta come vettori, per comodità
    dx = [dx1; dx2; dx3];
    dy = [dy1; dy2; dy3];

   
    for j = 1 : 3
        jj = geom.pivot.pivot (geom.elements.triangles(e,j));
        if jj > 0
                                
                                
            for k = 1 : 3
                kk = geom.pivot.pivot (geom.elements.triangles(e,k));
                 

                L_q = length(omega);
                    diff = 0;
                    conv = 0;
                    reaz = 0;
                    term_noto = zeros(1,istanti);

                    stab_diff = 0;
                    stab_conv = 0;
                    stab_reaz = 0;
                    stab_termine_noto = zeros(1,istanti);

                    for q = 1: L_q
                        nodo_q_locale = [csi(q);eta(q)];
                        % valuto le funzioni test nei nodi di quadratura
                        phi_val_1 = phi_1 (csi(q),eta(q));
                        phi_val_2 = phi_2 (csi(q),eta(q));
                        phi_val_3 = phi_3 (csi(q),eta(q));
                        phi = [phi_val_1; phi_val_2; phi_val_3];

                        % tramite trasformaz affine, passo da triang di rif
                        % a triang fisico
                        nodo_q_globale = vertice_3 + B * (nodo_q_locale);
                   
                        % calcolo i contributi di diffusione, convezione,
                        % reazione; essi andranno nella matrice A o A_D
                        diff = diff + omega(q) * epsilon (nodo_q_globale(1), nodo_q_globale(2)) * grad_phi(:,k)' * B_inv * B_inv_tras * grad_phi(:,j) * 2 * area;  
                        conv = conv + omega(q) * beta( nodo_q_globale(1), nodo_q_globale(2))' * B_inv_tras * grad_phi(:,k) * phi(j)  * 2 * area;
                        reaz = reaz + omega(q) * gamma(nodo_q_globale(1), nodo_q_globale(2)) * phi(k) * phi(j) * 2 * area;

                        for n = 1:istanti
                            term_noto(n) = term_noto(n)  + omega(q) * forzante_bis(nodo_q_globale(1), nodo_q_globale(2), t(n)) * phi(j) * 2 * area; % forzante o forzante bis
                        end

                        % stabilizzazione
                        % stab_massa = ?
                        stab_diff = stab_diff - 0; % per i P1 e basta; occhio al -
                        stab_conv = stab_conv + tau_h * omega(q) * beta( nodo_q_globale(1), nodo_q_globale(2))' * B_inv_tras  * grad_phi(:,k) * beta( nodo_q_globale(1), nodo_q_globale(2))' * B_inv_tras  * grad_phi(:,j) * (2* area);
                        stab_reaz = stab_reaz + tau_h * omega(q) * gamma (nodo_q_globale(1), nodo_q_globale(2)) * phi(k) * beta( nodo_q_globale(1), nodo_q_globale(2))' * B_inv_tras  * grad_phi(:,j) * (2* area);

                        for n = 1:istanti
                            stab_termine_noto(n) = stab_termine_noto(n) + tau_h * omega(q) * forzante_bis (nodo_q_globale(1), nodo_q_globale(2),t(n)) ...
                                                 * beta ( nodo_q_globale(1), nodo_q_globale(2))' * B_inv_tras  * grad_phi(:,j) * (2* area);  % forzante o forzante bis
                        end

                    end % ciclo sui nodi di quadratura

                    massa = area / 12 * (1 + (j==k));

                if kk > 0

                    % costruisco le matrici di diffusione, convezione e
                    % reazione
                    D(jj,kk) = D(jj,kk) + diff;
                    C(jj,kk) = C(jj,kk) + conv;
                    R(jj,kk) = R(jj,kk) + reaz;
                    M(jj,kk) = M(jj,kk) + massa;
                    
                 
                    D_stab(jj,kk) = D_stab(jj,kk) + stab_diff;
                    C_stab(jj,kk) = C_stab(jj,kk) + stab_conv;
                    R_stab(jj,kk) = R_stab(jj,kk) + stab_reaz;
                     % stabiliz sulla massa?
 
                    % coeff costanti, parte commentata
%                     A(jj,kk) = A(jj,kk) + ( epsilon / (4 * area) ) * (dy(k) * dy(j) + dx(k) * dx(j)) ...
%                                + (beta(1) * dy(k) + beta(2) * dx(k)) / 6 ...
%                                + gamma * area / 12 * ( 1+ j == k );
                
                else % costruisco la matrice A_D relativa a condizioni di Dirichlet

                    D_D(jj,-kk) = D_D(jj,-kk) + diff;
                    C_D(jj,-kk) = C_D(jj,-kk) + conv;
                    R_D(jj,-kk) = R_D(jj,-kk) + reaz;
                    M_D(jj,-kk) = M_D(jj,-kk) + massa;

                    D_D_stab(jj,-kk) = D_D_stab(jj,-kk) + stab_diff;
                    C_D_stab(jj,-kk) = C_D_stab(jj,-kk) + stab_conv;
                    R_D_stab(jj,-kk) = R_D_stab(jj,-kk) + stab_reaz;
                    % stabiliz sulla massa?


                    % coeff costanti, parte commentata
                    

%                     A_D (jj, -kk) = A_D (jj, -kk) + ( epsilon / (4 * area) ) * (dy(k) * dy(j) + dx(k) * dx(j)) ...
%                                     + (beta(1) * dy(k) + beta(2) * dx(k)) / 6 ...
%                                     + gamma * area / 12 * ( 1+ j == k );

                                 
                end
                
            end
            b (jj,:) = b (jj,:) + term_noto(1,:);
            b_stab (jj,:) = b_stab(jj,:) + stab_termine_noto(1,:);
        end
    end
   
end

D_tot = D;% + D_stab;
C_tot = C;% + C_stab;
R_tot = R;% + R_stab;
somme_righe_R = sum(R_tot);
%R_tot_lamping = diag (somme_righe_R);
R_tot_lumping = R_tot;

A_tot_lumping = D_tot + C_tot + R_tot_lumping;


D_D_tot = D_D;% + D_D_stab;
C_D_tot = C_D;% + C_D_stab;
R_D_tot = R_D;% + R_D_stab;


A_D_tot = D_D_tot + C_D_tot + R_D_tot;

% il b tot lo calcolo dopo Neumann

% soluzione nei nodi di Dirichlet
for n = 1:(istanti)
    for i = 1 : N_dri
        xd = geom.elements.coordinates ( geom.pivot.Di(i,1),1 );
        yd = geom.elements.coordinates ( geom.pivot.Di(i,1),2 );
        marker_dri = geom.pivot.Di(i,2);
        u_d (i,n) = cond_dirichlet (xd, yd, marker_dri, t(n));
    end
end

% aggiungiamo ora anche il contributo delle condizioni di Neumann

assembla_neumann()

b_tot = b + b_N;% + b_stab;


end

