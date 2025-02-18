function Find_periodic_solution_driver(gmax)
%% load parameters
fsz = 20;
graphix_flag = 0;
%% load parameter file
dat = load("model_parameters.mat");
% save("model_parameters.mat","L","tau","h_age","dt","par_ifun","r","a",
% "k","bb","s","z","Mm","rho", "mu_b_par", "birth_par", "birth_exp_par", "death_exp_par");
par = [dat.L, dat.tau, dat.h_age, dat.dt, dat.par_ifun, dat.r, dat.a, dat.k,...
    dat.bb, dat.s, dat.z, dat.Mm, dat.rho,dat.mu_b_par,dat.birth_par,...
    dat.birth_exp_par,dat.death_exp_par];
%%

tau = dat.tau;

gmin = 0;
gstep = 0.01;
gvals = (gmin:gstep:gmax)';
Ngvals = length(gvals);

Tmax = 1e3;

u1_eq = zeros(Ngvals,1);
u2_eq = zeros(Ngvals,1);
prey_eq = zeros(Ngvals,1);
u1_loop_min = zeros(Ngvals,1);
u1_loop_max = zeros(Ngvals,1);
u2_loop_min = zeros(Ngvals,1);
u2_loop_max = zeros(Ngvals,1);
prey_loop_min = zeros(Ngvals,1);
prey_loop_max = zeros(Ngvals,1);
period = zeros(Ngvals,1);

iter_max = round(100/dat.dt);
% Poincare section: hyperplane prey = prey_equilib
for j = 1 : Ngvals
    g = gvals(j);
    fprintf("j = %d, tau = %d, g = %.4f\n",j,tau,g);
    feq_name = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
    % save(fname,'k','g','x','u','u1','u2','evals');
    eqdata = load(feq_name);
    prey_equilib = eqdata.x;
    u1_eq(j) = eqdata.u1;
    u2_eq(j) = eqdata.u2;
    prey_eq(j) = eqdata.x;

    % make initial data for the limit cycle
    [prey,u] = MakeInitData(g,Tmax,par);
    % integrate trajectory until the intersection with 
    % the hyperplane prey = prey_equilib
    prey_sign = sign(prey - prey_equilib);
    Nu = length(u);
    if prey_sign == 0
       [prey,u,~,~] = PreyPredTimeStep(prey,u,g,par);
       if prey == prey_equilib
           fprintf("prey_sign = 0\n");
           return
       end
       prey_sign = sign(prey - prey_equilib);
    end
    iter = 0;
    while sign(prey - prey_equilib) == prey_sign && iter < iter_max
        prey_prev = prey;
        u_prev = u;
        [prey,u,~,~] = PreyPredTimeStep(prey,u,g,par);
        iter = iter + 1;
        if iter == iter_max
            prey = NaN;
        end
    end
    if isfinite(prey)
        fprintf("Look for the limit cycle\n");
        % find the point of intersection of the trajectory with the hyperplane
        ustar = u_prev + (u - u_prev)*(prey_equilib - prey_prev)/(prey - prey_prev);
    
        % start looking for the limit cycle
        u = LevenbergMarquardt(prey_equilib,ustar,g,par);
        if isfinite(u)
        %     u = Newton1loop(prey_equilib,u,k,g,par);
            J = Jacobian(prey_equilib,u,g,par) + eye(Nu);
            evals = eig(J);
            max_abs_eval = max(abs(evals));
            min_abs_evals = min(abs(evals));
            fprintf("max abs eval = %d, mn abs eval = %d\n",max_abs_eval,min_abs_evals);
            fname = sprintf("Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g);
            save(fname,"prey_equilib","prey","u","g","evals");
            % find min and max vals along the limit cycle
            [u1_loop_min(j),u1_loop_max(j),u2_loop_min(j),u2_loop_max(j),...
                prey_loop_min(j),prey_loop_max(j),period(j)] = ...
                    find_loop_min_max(prey_equilib,prey,u,g,par);
            fprintf("Results for j = %d, g = %.4f:\n",j,g);    
            fprintf("u1: [%d,%d]\n", u1_loop_min(j),u1_loop_max(j));
            fprintf("u2: [%d,%d]\n", u2_loop_min(j),u2_loop_max(j));
            fprintf("prey: [%d,%d]\n", prey_loop_min(j),prey_loop_max(j));
            fprintf("period = %d\n",period(j));
        else
            u1_loop_min(j) = u1_eq(j);
            u1_loop_max(j) = u1_eq(j);
            u2_loop_min(j) = u2_eq(j);
            u2_loop_max(j) = u2_eq(j);
            prey_loop_min(j) = prey_eq(j);
            prey_loop_max(j) = prey_eq(j);
            period(j) = NaN;
        end
    else
        u1_loop_min(j) = u1_eq(j);
        u1_loop_max(j) = u1_eq(j);
        u2_loop_min(j) = u2_eq(j);
        u2_loop_max(j) = u2_eq(j);
        prey_loop_min(j) = prey_eq(j);
        prey_loop_max(j) = prey_eq(j);
        period(j) = NaN;
    end
end

save("FindCycleOutput.mat","gvals","u1_eq","u2_eq","prey_eq","u1_loop_min",...
    "u1_loop_max","u2_loop_min","u2_loop_max","prey_loop_min","prey_loop_max",...
    "period");
%% graphics
if graphix_flag == 1
    
    fsz = 20;
    
    color_prey = [0,0,1]; % blue
    color_u1 =  [1,0,0]; % red
    color_u2 = [1,0.8,0.2]; % dark yellow
    figure(1);clf; hold on;
    hand = zeros(1,9);
    hand(1) = plot(gvals,prey_eq,'Linewidth',2,'color',color_prey,'DisplayName','prey');
    hand(2) = plot(gvals,prey_loop_min,'Linewidth',2,'color',color_prey,'DisplayName','prey');
    hand(3) = plot(gvals,prey_loop_max,'Linewidth',2,'color',color_prey,'DisplayName','prey');
    hand(4) = plot(gvals,u1_eq,'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
    hand(5) = plot(gvals,u1_loop_min,'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
    hand(6) = plot(gvals,u1_loop_max,'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
    hand(7) = plot(gvals,u2_eq,'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
    hand(8) = plot(gvals,u2_loop_min,'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
    hand(9) = plot(gvals,u2_loop_max,'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
    set(gca,'Fontsize',fsz);
    xlabel('g','FontSize',fsz);
    ylabel('Species counts','Fontsize',fsz);
    legend(hand(1:3:7))
end

end

%% integrate for one lomit cycle
function uloop = integrate1loop(prey_equilib,u,g,par)
% make one step to escape from the Poincare section
[prey,u,~,~] = PreyPredTimeStep(prey_equilib,u,g,par); 
iter = 0;
iter_max = round(100/par(4));

% We needed to detect the crossing of 
% the Poincare section in the same direction
prey_sign = sign(prey - prey_equilib);
while prey_sign == sign(prey - prey_equilib) && iter < iter_max
    prey_prev = prey;
    u_prev = u;
   [prey,u,~,~] = PreyPredTimeStep(prey,u,g,par);
    iter = iter + 1;
    if iter == iter_max
        prey = NaN;
    end
end
if isfinite(prey)
    iter = 0;
    while prey_sign ~= sign(prey - prey_equilib) && iter < iter_max
        prey_prev = prey;
        u_prev = u;
       [prey,u,~,~] = PreyPredTimeStep(prey,u,g,par);
        iter = iter + 1;
        if iter == iter_max
            prey = NaN;
        end
    end   
end
if isfinite(prey)
    % find the point of intersection of the trajectory with the hyperplane
    aux = (prey_equilib - prey_prev)/(prey - prey_prev);
    uloop = u_prev + (u - u_prev)*aux;
else
    uloop = NaN;
end
end
%%
function r = residual(prey_equilib,u,g,par)    
r = integrate1loop(prey_equilib,u,g,par) - u;
end
%%
function J = Jacobian(prey_equilib,u,g,par)
ustar = integrate1loop(prey_equilib,u,g,par);
if isfinite(ustar)
    Nu = length(u);
    J = zeros(Nu,Nu);
    e = zeros(size(u));
    du = 1e-6;
    fprintf("Computing the Jacobian\n");
    for j = 1 : Nu
        e1 = e;
        e1(j) = 1;
        u0 = u + du*e1;
        u1 = integrate1loop(prey_equilib,u0,g,par);
        J(:,j) = (u1 - ustar)/du;
    end
    J = J - eye(Nu);
else
    J = NaN;
end
end
%%
function ustar = Newton1loop(prey_equilib,u,g,par)
ustar = integrate1loop(prey_equilib,u,g,par);
nor = norm(ustar - u);
tol = 1e-12;
iter = 0;
fprintf("iter = %d, res = %d\n",iter,nor);
while nor > tol 
    % find Jacobian matrix for the map u_new = integrate1loop(u)
    Nu = length(u);
    J = zeros(Nu,Nu);
    e = zeros(size(u));
    du = 1e-6;
    fprintf("Computing the Jacobian\n");
    for j = 1 : Nu
        e1 = e;
        e1(j) = 1;
        u0 = u + du*e1;
        u1 = integrate1loop(prey_equilib,u0,g,par);
        J(:,j) = (u1 - ustar)/du;
    end
    J = J - eye(Nu);
    % Newton's step
    u = u - J\ustar;
    ustar = integrate1loop(prey_equilib,u,g,par);
    nor = norm(ustar - u);
    iter = iter + 1;
    fprintf("iter = %d, res = %d\n",iter,nor);
end
ustar = integrate1loop(prey_equilib,u,g,par);
end

%%
%%
function u = LevenbergMarquardt(prey_equilib,u,g_par,par)
tol = 1e-6;
iter_max = 100;
%%
Rmax = 1;
Rmin = 1e-14;
rho_good = 0.75;
rho_bad = 0.25;
eta = 0.01;
% iter_max = 10000;
% tol = 5e-3;
%% setup training mesh
% nt = 5;
%%
r = residual(prey_equilib,u,g_par,par);
J = Jacobian(prey_equilib,u,g_par,par);
if ~isfinite(r)
    fprintf("Failed to find limit cycle\n");
    u = NaN;
    return
end
f = F(r); % the objective function -- see line 164
g = J'*r; % the gradient of the objective function
nor = norm(g);
R = Rmax/5; % initial trust region radius

fprintf('Initially: f = %d, nor(g) = %d\n',f,nor); 
%% The trust region template
tic
rho = 0;
iter = 0;
flag = 1;
I = eye(length(u));
% quadratic model: m(p) = (1/2)||r||^2 + p'*J'*r + (1/2)*p'*J'*J*p;
while nor > tol && iter < iter_max
    B = J'*J + (1e-12)*I;
    pstar = -B\g;
    if norm(pstar) <= R
        p = pstar;
        fprintf('Global min of quad model\n');
    else % solve constrained minimization problem
        isub = 0;
        lam = 1;
        while 1
            B1 = B + lam*I;
            C = chol(B1);
            p = -C\(C'\g); %C'*C = B1
            np = norm(p);
            dd = abs(np - R);
            if dd < 0.01*R
                break
            end
            q = C'\p;
            nq = norm(q);
            lamnew = lam +(np/nq)^2*(np - R)/R;
            if lamnew < 0
                lam = 0.5*lam;
            else
                lam = lamnew;
            end
            isub = isub + 1;
        end
        fprintf('Contraint minimization: %d substeps\n',isub);
    end
    iter = iter + 1;  
    if flag == 0
        break;
    end
    % assess the progress
    unew = u + p;
    rnew = residual(prey_equilib,unew,g_par,par);
    Jnew = Jacobian(prey_equilib,unew,g_par,par);
    mnew = 0.5*(r'*r) + g'*p + 0.5*p'*B*p;
    fnew = F(rnew);
    rho = (f - fnew + 1e-14)/(f - mnew + 1e-14);
    

    % adjust the trust region
    if rho < rho_bad
        R = max([0.25*R,Rmin]);
    else
        if rho > rho_good && abs(norm(p) - R) < tol
            R = min([Rmax,2*R]);
        end
    end
    % accept or reject step
    if rho > eta            
        u = unew;
        r = rnew;
        J = Jnew;
        f = fnew;
        g = J'*r;
        nor = norm(g);        
        fprintf('iter # %d: f = %.14f, |df| = %.4e, rho = %.4e, R = %.4e\n',iter,f,nor,rho,R);
    end
end
fprintf('iter # %d: f = %.14f, |df| = %.4e, rho = %.4e, R = %.4e\n',iter,f,nor,rho,R);
cputime = toc;
fprintf('CPUtime = %d, iter = %d\n',cputime,iter);
end

%%
function p = cauchy_point(B,g,R)
    ng = norm(g);
    ps = -g*R/ng;
    aux = g'*B*g;
    if aux <= 0
        p = ps;
    else
        p = min(ng^3/(R*aux),1);
    end
end
%%
function f = F(r)
    f = 0.5*r'*r;
end

