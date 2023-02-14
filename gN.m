function z = gN (x,y, marker)
global epsilon

if marker == 2 % (lato sotto)
    z = epsilon(x,y) * grad_soluz_esatta_bis(x,y)' * [0;-1] ;
   
elseif marker == 4 % lato destra
    z = epsilon(x,y) * grad_soluz_esatta_bis(x,y)' * [1;0];
    
   
elseif marker == 6 % lato sopra
    z = epsilon(x,y) *grad_soluz_esatta_bis(x,y)' * [0;1];
   
elseif marker == 8 % lato sx
   z = epsilon(x,y) * grad_soluz_esatta_bis(x,y)' * [-1;0];
   
end
end