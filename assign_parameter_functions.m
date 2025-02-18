function [I1,I2,mu_base,birth_rate_base] = assign_parameter_functions(par_ifun,...
    tau,L,mu_b_par,death_exp_par,birth_par,birth_exp_par)

% Smooth approximations of [0,tau*),(tau*,\infty)
I1 = @(x) 1./(1+exp(-par_ifun*(tau-x)));
I2 = @(x) 1./(1+exp(-par_ifun*(x-tau)));
% birth- and death-rate base functions
mu_base = @(x) mu_b_par*exp(death_exp_par*(x-L));
birth_rate_base = @(x) heaviside(x-tau).*...
    (birth_par*exp(-birth_exp_par*(x-tau))+birth_par);

end