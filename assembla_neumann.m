function assembla_neumann () % P2
global geom Ndof b_N
Le = size(geom.pivot.Ne,1);

b_N = zeros (Ndof,1);
for e = 1 : Le
    l = geom.pivot.Ne(e,1); % indice lato che sta su bordo Neumann
    marker_neu = geom.pivot.Ne(e,2); % marker condizione Neumann

    i_b = geom.elements.borders (l,1); % indice vertice inizio lato
    i_m = geom.elements.borders (l,5); % indice punto medio lato
    i_e = geom.elements.borders (l,2); % indice vertice fine lato

    % coordinate di punti b ed e
    xb =  geom.elements.coordinates (i_b , 1);
    xm =  geom.elements.coordinates (i_m , 1);
    xe =  geom.elements.coordinates (i_e , 1);
    yb =  geom.elements.coordinates (i_b , 2);
    ym =  geom.elements.coordinates (i_m , 2);
    ye =  geom.elements.coordinates (i_e , 2);

    lung_l = sqrt ( (xe - xb)^2 + (ye - yb)^2 );

    gb = gN(xb, yb, marker_neu);
    gm = gN(xm, ym, marker_neu);
    ge = gN(xe, ye, marker_neu);

    
    % gdl associati ai nodi di Neumann
    ii_b = geom.pivot.pivot(i_b); 
    ii_m = geom.pivot.pivot(i_m);
    ii_e = geom.pivot.pivot(i_e);
    
    
    if ii_b > 0
        b_N (ii_b) = b_N(ii_b) + (2/15 * gb + 1/15 * gm - 1/30 * ge) * lung_l;
    end
    if ii_m > 0
        b_N (ii_m) = b_N(ii_m) + (1/15 * gb + 8/15 * gm + 1/15 * ge) * lung_l;
    end
    if ii_e > 0
        b_N (ii_e) = b_N(ii_e) + (-1/30 * gb + 1/15 * gm + 2/15 * ge) * lung_l;
    end
end


end
