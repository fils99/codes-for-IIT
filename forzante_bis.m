function z = forzante_bis (x, y)
global epsilon beta gamma

%z = - 4 * 16^4 * 10^(0) * ( y^4 * (1-y)^4 * (3 * (x - x^2)^2*(1-2*x)^2 -2* (x-x^2)^3)  +  x^4 * (1-x)^4 * (3 * (y - y^2)^2*(1-2*y)^2 -2 *(y-y^2)^3) ) + (1) * (16 * x * (1 - x) * y * (1 - y))^4; %caso 1

%z = - 10^(-6) / (exp(10^6) - exp(-10^6)) * 10^12 * (exp(10^6 * x) - exp(- 10^6 *x)) + (exp(10^6 * x) - exp(- 10^6 * x)) / (exp(10^6) - exp(-10^6));

 z =  - laplaciano_sol_esatta_bis(x,y) * epsilon(x,y) + beta(x,y)' * grad_soluz_esatta_bis(x,y) + gamma(x,y) * sol_esatta_bis(x,y);
end