function z = cond_dirichlet(x,y,marker)
if marker == 3
     z = sol_esatta_bis(x,y); % caso 1
    %z =  sin(x) .* cos(y); % caso 2
elseif marker == 5
     z = sol_esatta_bis(x,y);
    %z =  sin(x) .* cos(y);
elseif marker == 7
     z = sol_esatta_bis(x,y);
    %z =  sin(x) .* cos(y);
elseif marker == 9
     z = sol_esatta_bis(x,y); 
    %z =  sin(x) .* cos(y);
elseif marker == 1
     z = sol_esatta_bis(x,y);
    %z =  sin(x) .* cos(y);
end
end