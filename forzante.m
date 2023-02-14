function z = forzante (x, y, tempo)

% z = 32 * 10 ^(0) * ( x * (1 - x) + y * (1 - y)); % per epsilon = 1

% z = -16 * y * (1-y) *(1-4*x) +32 * (1-x) *x^2; % per epsilon = x

% z = - 32 * ( y * (1-y) * (x - 3*x^2 - y^2 -10) + x * (1-x) * (y - x^2 - 3 * y^2 -10)); %  per epsilon = x^2 + y^2 +1

% z = 32 * y^2 *(1-y) + 16*x *(1 - x) * (3 * y^2 -8 * y^3); %per epsilon=y^3

% z = 32 * y^3 *(1-y) -16 * x * (1-x)*(2 * y - 6 * y^2); % per epsilon = y^2

% z = 32 * y * (y^2 + 10)*(1 - y) - 16 * x * (1-x) * (-6 * y^2 + 2 * y -20);  % per epsilon = y^2 +10

% z = exp(x) * (16 * y * (1-y) * (2 *x + 1) +32 * (1-x) * x);  % per epsilon = exp(x)

% z = 32 * y * (1 - y) * (sin(y)^2 + 10) - 16 * x * (1 - x) * ( -2 * (sin(y)^2 +10) + 2 * (1 - 2*y) * sin(y) * cos(y)); % sin(y)^2 + 10; in caso cambia il 10 con 1 e vedi che non funziona

% z = - 32 * y * (1 - y) * (x /( x^2 + 2) * (1 - 2*x) - log(x^2 +2)) + 32 * x * log(x^2 +2) * (1 - x); % per epsilon = log(x^2 + 2)

% z = 16 * ( (1 - y) * (1 + 4 * x + 2 * y) + (1 - x) * (1 + 2 * x + 4 * y)); % per epsilon = x + y + 1

% z = -32 * y * (1 - y) * ((1 - 2*x) * sin(x) * cos(x) - sin(x)^2 - cos(y)^2 - 10) -32 * x * (1 - x) * ( - (1 - 2*y) * sin(y)*cos(y) - sin(x)^2 - cos(y)^2 - 10); % per epsilon = sin(x)^2 + cos(y)^2 + 10

% z =   16 * (1) * y * (1-y) *(1-2*x) + 16 * (1) * x * (1-x) *(1-2*y) + 10^(-5) * 32 * (y*(1-y) + x*(1-x)) + (0) * 16 * x*(1-x)*y*(1-y); % epsilon cost, beta e sigma come ti pare

z = x^2;

%z = 16 * (y * (1 - y) * (2* 10^(-15) - (1-2*x) * (x^3 + 10*x)) + x * (1 - x) * (2* 10^(-15) - (1 - 2*y)*(y^3 + 10*y))); % per epsilon = 10^(-15), beta = [-x^3 - 10*x; -y^3 - 10 *y]

% z = 16 * (y * (1 - y) * (2 + (1-2*x) * 2) + x * (1 - x) * (2 + (1-2*y)*2)); % per epsilon =1, beta = [2;2]

%z = - 16 * 10^(-15) * y * (1-y) * (-2 * log(x^2 + 2) + (1-2*x) / (x^2 + 2) * 2 * x) + 32 * 10^(-15) * x * (1 -x)* log(x^2+2) -16 * x * y * ((1-y)*(1- 2 * x) + (1-x) * (1-2*y)) + 16 * x * (1-x) * y^3 * (1-y); % epsilon = 10^-5 * log(x^2 + 2), beta = [-x; -y], gamma = y^2
end
