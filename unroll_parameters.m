function [L,tau,h_age,dt,age_vector,AgeNode,par_ifun,r,a,k,bb,s,rho,z,...
    Mm,mu_b_par,birth_par,birth_exp_par,death_exp_par] ...
    = unroll_parameters(par)


% par = [dat.L, 1
%     dat.tau, 2
%     dat.h_age, 3
%     dat.dt, 4
%     dat.par_ifun, 5
%     dat.r, 6
%     dat.a, 7
%     dat.k,8
%     dat.bb, 9
%     dat.s, 10
%     dat.z, 11
%     dat.Mm, 12
%     dat.rho, 13
%     dat.mu_b_par, 14
%     dat.birth_par,15
%     dat.birth_exp_par, 16
%     dat.death_exp_par]; 17

L = par(1); % life spane
tau = par(2); % maturation age
h_age = par(3); % step in age
dt = h_age; % step in time
age_vector = 0:h_age:L; % partition of life span
%aa = dt/h_age;
AgeNode = floor(tau/h_age+1);
% parameters for smooth approximations for indicator functions
par_ifun = par(5);

% parameters of the model
% %s < g, k < bb
r = par(6); % intrinsic growth rate for prey
a = par(7); % interspecies competition coefficient for prey
k = par(8);
bb = par(9);
s = par(10);

rho = par(13);
z = par(11);
Mm = par(12);
mu_b_par = par(14);
birth_par = par(15);
birth_exp_par = par(16);
death_exp_par = par(17);



end