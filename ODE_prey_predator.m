function dy = ODE_prey_predator(y,par)
dy = zeros(3,1);

%% ODE system
% x' = rx - ax^2 + sy1 - by2
% y1'= kxy2 + (1-exp(-zeta*x))b2y2 - gxy1 - m1y1 - exp(-rho*x)Mmy1 - Dy1
% y2' = Dy1 - m2y2 - exp(-rho*x)Mmy2
x = y(1);
y1 = y(2);
y2 = y(3);

dy(1) = (par(1) - par(2)*x + par(6)*y1 - par(4)*y2)*x;

dy(2) = par(3)*x*y2 + (1-exp(-par(13)*x))*par(7)*y2 - par(5)*x*y1 - ...
    par(8)*y1 - exp(-par(11)*x)*par(12)*y1 - par(10)*y1;

dy(3) = par(10)*y1 - par(9)*y2 - exp(-par(11)*x)*par(12)*y2;

%1 r = 0.4; % intrinsic growth rate for prey
%2 a = 0.01; % interspecies competition coefficient for prey
%3 k = 0.3; % birth rate has term k*prey*indicator(maturity age)
%4 bb = 0.4;% effect on the prey by being eaten by adult predator
%5 g = 0.1; % juvenile predator damage coefficient from being eaten by prey
%6 s = 0.2; % prey benefit coefficient from eating juvenile predator
%7 b2 = 8.256614e-02;
%8 m1 = 2.091426e-02;
%9 m2 = 4.147095e-02;
%10 D = 9.250893e+00;
% 
%11 rho = 5;  % 5 an exponential factor in the death rate; 10 = oscillations at g = 1
%12 Mm = 1.0; % 1 a coefficient in the death rate: the death rate in the absense of prey
%13 z = 10;  % an exponential factor  in the birth rate
end