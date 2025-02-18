function [L,h_age,dt,age_vector,AgeNode,par_ifun,r,a,k,bb,s,rho,Mm,...
    z,birth_par,birth_exp_par,death_exp_par,mu_b_par] = assign_parameters(tau)
L = 30; % life spane
h_age = 0.1; % step in age
dt = h_age; % step in time
age_vector = (0:h_age:L)'; % partition of life span
AgeNode = floor(tau/h_age+1);
% parameters for smooth approximations for indicator functions
par_ifun = 100; % parameter in the indicator function 100

% parameters of the model
r = 0.4; % intrinsic growth rate for prey
a = 0.1; % interspecies competition coefficient for prey
k = 0.3; % birth rate has term k*prey*indicator(maturity age)
bb = 0.8; % effect on the prey by being eaten by adult predator % used to be 0.4
s = 0.2; % prey benefit coefficient from eating juvenile predator

rho = 5;  % 5 an exponential factor in the death rate; 10 = oscillations at g = 1
Mm = 1; % 1 a coefficient in the death rate: the death rate in the absense of prey
z = 10;  % an exponential factor  in the birth rate
birth_par = 0.05; 
birth_exp_par = 0.1;
death_exp_par = 0.1;
mu_b_par = 0.4;

end