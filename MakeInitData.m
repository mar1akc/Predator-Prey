function [prey,u] = MakeInitData(g,T,par)
%%%%%%%%%%%%

%%
[L,tau,h_age,dt,age_vector,AgeNode,par_ifun,r,a,k,bb,s,rho,z,...
    Mm,mu_b_par,birth_par,birth_exp_par,death_exp_par] ...
    = unroll_parameters(par);

[I1,I2,mu_base,birth_rate_base] = assign_parameter_functions(par_ifun,...
    tau,L,mu_b_par,death_exp_par,birth_par,birth_exp_par);

t = 0:dt:T; % time vector

N_age = length(age_vector);
Nt = length(t);
u = zeros(N_age,1);

% Boundary conditions
prey = 0.5;
% % Initial age distribution for the predator
initial_age_density = @(x) 0.1*heaviside(tau-x) + 0.05*heaviside(x-tau);
% Initial condition: the initial age distribution
u = initial_age_density(age_vector(:));

h05 = 0.5*h_age;

% juvenile predator at time 0
u1(1) = h05*(u(1)+u(AgeNode))+h_age*sum(u(2:AgeNode-1));
% adult predatir at time 0
u2(1) = h05*(u(AgeNode)+u(end))+h_age*sum(u(AgeNode+1:end-1));

for i=1:Nt-1
    [prey,u,~,~] = PreyPredTimeStep(prey,u,g,par);
end

end