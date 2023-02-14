function z = laplaciano_sol_esatta_bis (x,y)

z = 4 *16^4 * (y*(1-y))^4 * (3 * (x*(1-x))^2 * (1-2*x)^2 -2 * (x*(1-x))^3)...
  + 4 *16^4 * (x*(1-x))^4 * (3 * (y*(1-y))^2 * (1-2*y)^2 -2 * (y*(1-y))^3);

%% caso 2 

% a =5;
% c= 5;
% z = - (a^2+c^2) *  sin(a*x) * cos(c*y) ;

end