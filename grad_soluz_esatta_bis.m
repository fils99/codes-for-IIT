function z = grad_soluz_esatta_bis(x,y)

%% caso 1

z = 4* 16^4 * [(y*(1-y))^4 * (x*(1-x))^3 * (1-2*x); (x*(1-x))^4 * (y*(1-y))^3 * (1-2*y)]; 


%% caso 2
%  a =10;
%  c= 10;
%  z =  [a * cos(a*x)*cos(c*y); - c * sin(a*x)*sin(c*y)]; 

 %% caso 3

% z =[cosh(x); 0];
end