function z = gN (x,y, marker,t)

global epsilon 

if marker == 2 % (lato sotto)
    z = grad_soluz_esatta_bis(x,y,t)' *[0;-1]; 
   
elseif marker == 4 % lato destra
    z = grad_soluz_esatta_bis(x,y,t)' *[1;0]; 

elseif marker == 6 % lato sopra
    z = grad_soluz_esatta_bis(x,y,t)' *[0;1]; 

elseif marker == 8 % lato sx
    z = grad_soluz_esatta_bis(x,y,t)' *[-1;0]; 

end
end