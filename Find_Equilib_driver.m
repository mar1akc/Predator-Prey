function Find_Equilib_driver(g)

graphix_flag = 0;
%% load parameter file
dat = load("model_parameters.mat");
% save("model_parameters.mat","L","tau","h_age","dt","par_ifun","r","a",
% "k","bb","s","z","Mm","rho", "mu_b_par", "birth_par", "birth_exp_par", "death_exp_par");
par = [dat.L, dat.tau, dat.h_age, dat.dt, dat.par_ifun, dat.r, dat.a, dat.k,...
    dat.bb, dat.s, dat.z, dat.Mm, dat.rho,dat.mu_b_par,dat.birth_par,...
    dat.birth_exp_par,dat.death_exp_par];
%%
fprintf("parameters:\n");
fprintf("L = %d\n",dat.L);
fprintf("tau = %d\n",dat.tau);
fprintf("h_age = %d\n",dat.h_age);
fprintf("dt = %d\n",dat.dt);
fprintf("par_ifun = %d\n",dat.par_ifun);
fprintf("r = %d\n",dat.r);
fprintf("a = %d\n",dat.a);
fprintf("k = %d\n",dat.k);
fprintf("bb = %d\n",dat.bb);
fprintf("s = %d\n",dat.s);
fprintf("z = %d\n",dat.z);
fprintf("Mm = %d\n",dat.Mm);
fprintf("rho = %d\n",dat.rho);
fprintf("mu_b_par = %d\n",dat.mu_b_par);
fprintf("birth_par = %d\n",dat.birth_par);
fprintf("birth_exp_par = %d\n",dat.birth_exp_par);
fprintf("death_exp_par = %d\n",dat.death_exp_par);

% g = 0.88;
tau = dat.tau;
Age_Node = floor(dat.tau/dat.h_age+1);
age_vector = (0:dat.h_age:dat.L)';
[~,~,mu_base,birth_rate_base] = assign_parameter_functions(dat.par_ifun,...
    dat.tau,dat.L,dat.mu_b_par,dat.death_exp_par,dat.birth_par,dat.birth_exp_par);


fname = sprintf('Data/Initial_guess_g%.2f.mat',g);
% save(fname,'prey_end','u','g');
data = load(fname);
x = data.prey_end;
u = data.u;

gmin = 0.0;
gstep = 0.01;
gmax = data.g;
g = data.g;

gvals = linspace(gmin,gmax,round((gmax-gmin)/gstep)+1); 
gvals = gvals(end:-1:1);
Nvals = length(gvals);
xu12vals = zeros(3,Nvals);

save("gvalues.mat","gvals");
%%
% 
% We iterate:
% 
% W_{n+1} = F(W_{n}) where W = [x; u]
% 
% Let W^* be the equilibrium. Then 
% 
% W_{n+1} = F(W_n) = F(W^*) + J(W^*)(W_n - W^*) + o(||W_n - W^*||)
%
% Subtract W^* from both sides. Denote e_n:=W_n-W^*. Since F(W^*) = W^*, we get
%
% e_{n+1} = J(W^*)e_n + o(||e_n||)
%
% W^* is stable if evals of J(W^*) are less or equal 1 in absolute value, 
% and those that are equal to 1 in absolute value have multiplicity 1. 
% 
%
% 
%%

fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
[x,u,u1,u2,evals] = Newton(x,u,g,par);
fprintf('tau = %.4f, g = %.4f, max_abs_eval = %.4e\n',tau,g,max(abs(evals)));
if(max(abs(evals))<= 1 ) 
    fprintf("stable\n");
else
    fprintf("unstable\n");
end
xu12vals(:,1) = [x;u1;u2];
save(fname,'tau','g','x','u','u1','u2','evals');

% for ODE and DDE models
D_ode = zeros(Nvals,1);
B_ode = zeros(Nvals,1);
M1_ode = zeros(Nvals,1);
M2_ode = zeros(Nvals,1);

for j = 1 : Nvals
    g = gvals(j);
    fprintf('tau = %.4f, g = %.4f\n',tau,g);
    [x,u,u1,u2,evals] = Newton(x,u,g,par);
    xu12vals(:,j) = [x;u1;u2];
    fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
    fprintf('tau = %.4f, g = %.4f,max_abs_eval = %.4e\n',tau,g,max(abs(evals)));
    save(fname,'tau','g','x','u','u1','u2','evals');
    if(max(abs(evals))<= 1 ) 
        fprintf("stable\n");
    else
        fprintf("unstable\n");
    end
    % save parameters for ODE and DDE models
    D_ode(j) = u(Age_Node)/(u1);
    B_ode(j) = sum(birth_rate_base(age_vector((Age_Node+1):end)).*u((Age_Node+1):end))/sum(u((Age_Node+1):end));
    M1_ode(j) = sum(mu_base(age_vector(1:Age_Node)).*u(1:Age_Node))/sum(u(1:Age_Node));
    M2_ode(j) = sum(mu_base(age_vector((Age_Node+1):end)).*u((Age_Node+1):end))/sum(u((Age_Node+1):end));

end

fname = sprintf("ODEparams/ODEparams_tau%.2f.mat",tau);
save(fname,"gvals","D_ode","B_ode","M1_ode","M2_ode");

%%
if graphix_flag == 1
    fsz = 20;
    figure(1);clf; hold on;
    plot(gvals,xu12vals(1,:),'Linewidth',2,'DisplayName','prey')
    plot(gvals,xu12vals(2,:),'Linewidth',2,'DisplayName','Juv. predator')
    plot(gvals,xu12vals(3,:),'Linewidth',2,'DisplayName','Adult predator')
    set(gca,'Fontsize',fsz);
    xlabel('g','FontSize',fsz);
    ylabel('Equilibrium counts','Fontsize',fsz);
    legend()
end
end

function [x,u,u1,u2,evals] = Newton(x,u,g,par)
% Newton's method
Nu = length(u);
tol = 1e-12;
h = 1e-6;
[xstep,ustep,~,~] = PreyPredTimeStep(x,u,g,par);
F = [x;u] - [xstep;ustep];
nor = norm(F);
iter = 0;
fprintf('iter = %d, norm(F) = %d\n',iter,nor);
iter_max = 40;
while nor > tol && iter < iter_max
    % compute the Jacobian 
    J = Jacobian(x,u,g,par,h,F);
    xu_new = [x;u] - J\F;
    x = xu_new(1);
    u = xu_new(2:end);
    [xstep,ustep,~,~] = PreyPredTimeStep(x,u,g,par);
    F = [x;u] - [xstep;ustep];
    nor = norm(F);
    iter = iter + 1;
    fprintf('iter = %d, norm(F) = %d\n',iter,nor);
end
[x,~,u1,u2] = PreyPredTimeStep(x,u,g,par);
fprintf('x = %.14e, u1 = %.14e, u2 = %.14e\n',x,u1,u2);
J = Jacobian(x,u,g,par,h,F);
evals = eig(J);
evals = 1-evals;
end

%%
function J = Jacobian(x,u,g,par,h,F)
Nu = length(u);
% compute the Jacobian 
J = zeros(Nu+1);
[xshift_step,ushift_step,~,~] = PreyPredTimeStep(x+h,u,g,par);
Fshift = [x+h;u] - [xshift_step;ushift_step];
J(:,1) = (Fshift - F)/h;
for j = 1 : Nu
    ushift = u;
    ushift(j) = u(j) + h;
    [xshift_step,ushift_step,~,~] = PreyPredTimeStep(x,ushift,g,par);
    Fshift = [x;ushift] - [xshift_step;ushift_step];
    J(:,j+1) = (Fshift - F)/h;
end
end
