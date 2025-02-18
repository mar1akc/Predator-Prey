function dy = DDE_prey_predator(y,ydelay,par)
dy = zeros(3,1);

%% DDE system
    % x' = x(r - ax + sy_1 - by_2),  
    % y_1' =  kxy_2 + (1-e^{-\zeta x})b_2y_2 - gxy_1 - m_1y_1 - \mu_M e^{-\rho x}y_1 -u(t,\tau^*), 
    % y_2' = - m_2y_2 - \mu_M e^{-\rho x} y_2 + u(t,\tau^*)

% If t \ge \tau^*,
% u(t,\tau^*) = \left(kx(t-\tau^*)y_2(t-\tau^*) + 
% \left(1-e^{-\zeta x(t-\tau^*)}\right)b_2y_2(t-\tau^*)\right)
% \times\thinspace&\exp\left(-\frac{k\tau^*}{2}[x(t-\tau^*) + x(t)] - 
% M_B - \frac{\mu_M\tau^*}{2}\left(e^{-\rho x(t-\tau^*)} + 
% e^{-\rho x(t)}\right)\right)

% If t < tau*,
% u(t,\tau^*) =  u_0(\tau^*-t)\exp\left(-\frac{k\tau^*}{2}[x(\tau^*-t) + x(t)]
% - M_B 
% -\frac{\mu_M\tau^*}{2}\left(e^{-\rho x(\tau^*-t)} + e^{-\rho x(t)}\right)\right)



x = y(1);
x_delay = ydelay(1);
y1 = y(2);
y1_delay = ydelay(2);
y2 = y(3);
y2_delay = ydelay(3);

k = par(3);
z = par(13);
tau = par(14);
MB = par(15);
Mm = par(12);
rho = par(11);

uu1 = k*x_delay*y2_delay + (1 - exp(-z*x_delay))*par(7)*y2_delay;
uu2 = -0.5*k*tau*(x_delay + x) - par(15) - 0.5*Mm*tau*(exp(-rho*x_delay)+exp(-rho*x));
uu = uu1*exp(uu2);

dy(1) = (par(1) - par(2)*x + par(6)*y1 - par(4)*y2)*x;

dy(2) = par(3)*x*y2 + (1-exp(-par(13)*x))*par(7)*y2 - par(5)*x*y1 - ...
    par(8)*y1 - exp(-par(11)*x)*par(12)*y1 - uu;

dy(3) = -par(9)*y2 - exp(-par(11)*x)*par(12)*y2 + uu;

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
%14 tau
%15 MB
end