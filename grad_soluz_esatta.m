function z = grad_soluz_esatta(x,y)
z = [16 * y * (1-y) * (1-2*x) ; 16 * x * (1 - 2 * y) * (1-x)];
end