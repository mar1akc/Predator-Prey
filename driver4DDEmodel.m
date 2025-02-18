function driver4DDEmodel()
L = 30;
r = 0.4; % intrinsic growth rate for prey
a = 0.1; % interspecies competition coefficient for prey
k = 0.3; % birth rate has term k*prey*indicator(maturity age)
bb = 0.8; % effect on the prey by being eaten by adult predator
s = 0.2; % prey benefit coefficient from eating juvenile predator
rho = 5;  % 5 an exponential factor in the death rate; 10 = oscillations at g = 1
Mm = 1.0; % 1 a coefficient in the death rate: the death rate in the absense of prey
z = 10;  % an exponential factor  in the birth rate
birth_par = 0.05; %used to be 0.05
birth_exp_par = 0.1;
death_exp_par = 0.1;
mu_b_par = 8.0*birth_par*exp(-L*death_exp_par); 
save("DDEmodel_parameters.mat","r","a","k","bb","s","z","Mm","rho", "mu_b_par", "birth_par", "birth_exp_par", "death_exp_par");

for tau = 1.18:0.01:3.0
    fprintf("tau = %d\n",tau);
    EqulibriumSolutionsDDE(tau);
end
end