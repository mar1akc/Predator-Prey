function EquilibriumSolutionsODE()
%% parameters that do not change
dat = load("model_parameters.mat");
% save("model_parameters.mat","L","tau","h_age","dt",
% "par_ifun","r","a","k","bb","s","z","Mm","rho");
% par = [dat.L, dat.tau, dat.h_age, dat.dt, dat.par_ifun, dat.r, dat.a, dat.k,...
%     dat.bb, dat.s, dat.z, dat.Mm, dat.rho, dat.mu_b_par, dat.birth_par, ...
%     dat.birth_exp_par, dat.death_exp_par];
tau = dat.tau;
par_ifun = dat.par_ifun;
% r = 0.4; % intrinsic growth rate for prey
% a = 0.1; % interspecies competition coefficient for prey
% k = 0.3; % birth rate has term k*prey*indicator(maturity age)
% bb = 0.4;% effect on the prey by being eaten by adult predator
% s = 0.2; % prey benefit coefficient from eating juvenile predator
% 
% rho = 5;  % 5 an exponential factor in the death rate; 10 = oscillations at g = 1
% Mm = 1.0; % 1 a coefficient in the death rate: the death rate in the absense of prey
% z = 10;  % an exponential factor  in the birth rate
%% parameters
par = [dat.r,dat.a,dat.k,dat.bb,0,dat.s,0,0,0,0,dat.rho,dat.Mm,dat.z];

%fname = sprintf("ODEparams_tau%.2f.mat",tau);
fname = sprintf("ODEparams/ODEparams_tau%.2f.mat",tau);
% fname = "ODEparams_tau1ppp00.mat";
% fname = "ODEparams_tau1pppp50.mat";
% fname = "test.mat";
%save(fname,"gvals","D_ode","B_ode","M1_ode","M2_ode");
par_data = load(fname);
gvals = par_data.gvals;
%gvals = flip(gvals);
D = par_data.D_ode;
M1 = par_data.M1_ode;
M2 = par_data.M2_ode;
B = par_data.B_ode;

Ngvals = length(gvals);

%% Setting up

y1_equilib = zeros(Ngvals,1);
y2_equilib = zeros(Ngvals,1);
prey_equilib = zeros(Ngvals,1);
stab = zeros(Ngvals,1);
y1_loop_min = zeros(Ngvals,1);
y1_loop_max = zeros(Ngvals,1);
y2_loop_min = zeros(Ngvals,1);
y2_loop_max = zeros(Ngvals,1);
prey_loop_min = zeros(Ngvals,1);
prey_loop_max = zeros(Ngvals,1);
period = zeros(Ngvals,1);

%% solving ODE
prey0 = 0.5;
y10 = 9.875000e-02;
y20 = 1.451250e+00;
x0 = [prey0,y10,y20];
% Find the equilibrium, determine if it is table or unstable, and find the
% stable limit cycle if the equilibrium is unstable

par(5) = gvals(1);
par(7) = B(1);
par(8) = M1(1);
par(9) = M2(1);
par(10) = D(1);

% find initial approximation to search for the equilibrium
odefun = @(t,x) ODE_prey_predator(x,par);    
rtol = 1e-9;
atol = 1e-9;
Tmax = 5e3;
options = odeset('RelTol',rtol,'AbsTol',atol);
[T,Y] = ode45(odefun,[0,Tmax],[1,1,1],options);
yeq_iguess = Y(end,:);
figure(100)
plot(T,Y)
drawnow;
for i = 1 : Ngvals    
    par(5) = gvals(i);
    par(7) = B(i);
    par(8) = M1(i);
    par(9) = M2(i);
    par(10) = D(i);


%% find equilibrium solutions

    fun = @(x) ODE_prey_predator(x,par);
    y = lsqnonlin(fun,yeq_iguess);
    yeq_iguess = y;
    fprintf("tau = %.2f, g = %d, b2 = %d, m1 = %d, m2 = %d, D = %d\n",tau,gvals(i),B(i),M1(i),M2(i),D(i));
    fprintf("prey = %d, y1 = %d, y2 = %d\n",y(1,1),y(1,2),y(1,3));

    prey_equilib(i) = y(1);
    y1_equilib(i) = y(2);
    y2_equilib(i) = y(3);
   
    J = Jacobian_eq(y,fun);
    [~,eval] = eig(J);
    ee = diag(eval);
    [esort,jsort] = sort(real(ee),'descend');
    ee = ee(jsort);
    fprintf("i = %d, ee = (%d+i%d,%d+i%d,%d+i%d)\n",i,...
        real(ee(1)),imag(ee(1)),real(ee(2)),imag(ee(2)),real(ee(3)),imag(ee(3)));
    if esort(1) <= 0
        stab(i) = 1;
        fprintf("g = %d, equilib is stable\n",gvals(i));
    
    else
        stab(i) = 0;
        fprintf("g = %d, equilib is unstable\n",gvals(i));
    end
    %% find limit cycle if the equilibrium is unstable
    if esort(1) > 0 % unstable equilibrium
        % find the initial approximation
        odefun = @(t,x) ODE_prey_predator(x,par);    
        rtol = 1e-9;
        atol = 1e-9;
        Tmax = 5e3;

        fprintf("Start with: (%d,%d,%d)\n",...
            y(1),y(2),y(3));

        options = odeset('RelTol',rtol,'AbsTol',atol);
        [T,Y] = ode45(odefun,[0,Tmax],y.*[1.1,1.1,1.1],options);
        figure(101); hold on;
        plot(T,Y)
        % find a point on the Poincare section

        events_cycle = @(t,y)PoincareSection(y,prey_equilib(i),1);
        options_cycle = odeset('RelTol',rtol,'AbsTol',atol,'Events',events_cycle);
        y0 = Y(end,:);
        fprintf("initial approximation for the limit cycle: (%d,%d,%d)\n",...
            y0(1),y0(2),y0(3));
        % find a point on the Poincare section
        [~,Y] = ode45(odefun,[0,Tmax],y0,options_cycle);
        fprintf("prey_equilib = %d, [prey,y1,y2] = [%d,%d,%d]\n",...
            prey_equilib(i),Y(end,1),Y(end,2),Y(end,3));
        u = Y(end,2:3)';
        % start looking for the limit cycle
        u = LevenbergMarquardt(prey_equilib(i),u,par);
        ycycle = [prey_equilib(i),u(1),u(2)];
        if isfinite(u)
        %     u = Newton1loop(prey_equilib,u,k,g,par);
            J = Jacobian(prey_equilib(i),u,par) + eye(2);
            evals = eig(J);
            max_abs_eval = max(abs(evals));
            min_abs_evals = min(abs(evals));
            fprintf("max abs eval = %d, mn abs eval = %d\n",max_abs_eval,min_abs_evals);
            fname = sprintf("DataODE_DDE/LimitCycle_tau%.4f_g%.4f.mat",tau,gvals(i));
            save(fname,"ycycle","evals");
            % find min and max vals along the limit cycle
            [y1_loop_min(i),y1_loop_max(i),y2_loop_min(i),y2_loop_max(i),...
                prey_loop_min(i),prey_loop_max(i),period(i)] = ...
                    find_loop_min_max(prey_equilib(i),u,par);
            fprintf("Results for j = %d, g = %.4f:\n",i,gvals(i));    
            fprintf("y1: [%d,%d]\n", y1_loop_min(i),y1_loop_max(i));
            fprintf("y2: [%d,%d]\n", y2_loop_min(i),y2_loop_max(i));
            fprintf("prey: [%d,%d]\n", prey_loop_min(i),prey_loop_max(i));
            fprintf("period = %d\n",period(i));
        else
            y1_loop_min(i) = y1_equilib(i);
            y1_loop_max(i) = y1_equilib(i);
            y2_loop_min(i) = y2_equilib(i);
            y2_loop_max(i) = y2_equilib(i);
            prey_loop_min(i) = prey_equilib(i);
            prey_loop_max(i) = prey_equilib(i);
            period(i) = NaN;
        end
    else
        y1_loop_min(i) = y1_equilib(i);
        y1_loop_max(i) = y1_equilib(i);
        y2_loop_min(i) = y2_equilib(i);
        y2_loop_max(i) = y2_equilib(i);
        prey_loop_min(i) = prey_equilib(i);
        prey_loop_max(i) = prey_equilib(i);
        period(i) = NaN;
    end

end

fname = sprintf("DataODE_DDE/ODE_BifurDiag_par_ifun%.1f_tau%.2f.mat",par_ifun,tau);

save(fname,"gvals","prey_equilib","y1_equilib","y2_equilib",...
    "prey_loop_min","prey_loop_max","y1_loop_min","y1_loop_max",...
    "y2_loop_min","y2_loop_max","period","stab");


Istab = find(stab == 1);
Inonstab = find(stab == 0);

    
fsz = 20;

color_prey = [0,0,1]; % blue
color_u1 =  [1,0,0]; % red
color_u2 = [1,0.8,0.2]; % dark yellow
if isempty(Inonstab) == 0 
    figure(1);clf; hold on;
    hand = zeros(1,9);
    hand(1) = plot(gvals(Istab),prey_equilib(Istab),'Linewidth',2,'color',color_prey,'DisplayName','prey');
    hand(2) = plot(gvals(Inonstab),prey_equilib(Inonstab),'--','Linewidth',2,'color',color_prey,'DisplayName','prey');
    hand(3) = plot(gvals,prey_loop_min,'Linewidth',2,'color',color_prey,'DisplayName','prey');
    hand(4) = plot(gvals,prey_loop_max,'Linewidth',2,'color',color_prey,'DisplayName','prey');

    hand(5) = plot(gvals(Istab),y1_equilib(Istab),'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
    hand(6) = plot(gvals(Inonstab),y1_equilib(Inonstab),'--','Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
    hand(7) = plot(gvals,y1_loop_min,'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
    hand(8) = plot(gvals,y1_loop_max,'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');

    hand(9) = plot(gvals(Istab),y2_equilib(Istab),'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
    hand(10) = plot(gvals(Inonstab),y2_equilib(Inonstab),'--','Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
    hand(11) = plot(gvals,y2_loop_min,'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
    hand(12) = plot(gvals,y2_loop_max,'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
    set(gca,'Fontsize',fsz);
    xlabel('Juv. Pred. Damage Parameter, g','FontSize',fsz,'Interpreter','latex');
    ylabel('Population Size','Fontsize',fsz,'Interpreter','latex');
    legend(hand([1,5,9]),'Interpreter','latex')
    title(['$\tau^* = $ ' num2str(tau) ],'Interpreter','latex')
    grid on

else
    figure(1);clf; hold on;
    plot(gvals(Istab),prey_equilib(Istab),'Linewidth',2,'color',color_prey,'DisplayName','prey');
    plot(gvals(Istab),y1_equilib(Istab),'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
    plot(gvals(Istab),y2_equilib(Istab),'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
    set(gca,'Fontsize',fsz)
    set(gcf,'color','white')
    legend('Location','northeast','Interpreter','latex') 
    xlabel('Juv. Pred. Damage Parameter, g','FontSize',fsz,'Interpreter','latex');
    ylabel('Population Size','Fontsize',fsz,'Interpreter','latex');
    title(['$\tau^* = $ ' num2str(tau) ],'Interpreter','latex')
    grid on
end

end
%%
function J = Jacobian_eq(y,fun)
J = zeros(3);
h = 1e-6;
f0 = fun(y);
f1 = fun(y+h*[1;0;0]');
f2 = fun(y+h*[0;1;0]');
f3 = fun(y+h*[0;0;1]');
J(:,1) = (f1 - f0)/h;
J(:,2) = (f2 - f0)/h;
J(:,3) = (f3 - f0)/h;
end


%%
%% integrate for one limit cycle
function uloop = integrate1loop(prey_equilib,u,par)
% step off the Poincare section
rtol = 1e-9;
atol = 1e-9;
Tmax = 5e1;
options = odeset('RelTol',rtol,'AbsTol',atol);
odefun = @(t,x) ODE_prey_predator(x,par);  
[~,Y] = ode45(odefun,[0,1e-3],[prey_equilib;u],options);
fprintf("integrate1loop: prey_equilib = %d,u0 = (%d,%d)\n",prey_equilib,u(1),u(2));
% find a point on the Poincare section
events_cycle = @(t,y)PoincareSection(y,prey_equilib,1);
options_cycle = odeset('RelTol',rtol,'AbsTol',atol,'Events',events_cycle);
y0 = Y(end,:);
fprintf("integrate1loop: p = %d, u1 = (%d,%d)\n",y0(1),y0(2),y0(3));

% find a point on the Poincare section
[~,Y] = ode45(odefun,[0,Tmax],y0,options_cycle);
uloop = Y(end,2:3)';
fprintf("integrate1loop: u_out = (%d,%d)\n",uloop(1),uloop(2));

end
%%
function r = residual(prey_equilib,u,par)    
r = integrate1loop(prey_equilib,u,par) - u;
end
%%
function J = Jacobian(prey_equilib,u,par)
ustar = integrate1loop(prey_equilib,u,par);
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
        u1 = integrate1loop(prey_equilib,u0,par);
        J(:,j) = (u1 - ustar)/du;
    end
    J = J - eye(Nu);
else
    J = NaN;
end
end



%% Levenberg-Marguardt
function u = LevenbergMarquardt(prey_equilib,u,par)
tol = 1e-9;
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
r = residual(prey_equilib,u,par);
J = Jacobian(prey_equilib,u,par);
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
        fprintf('Global min of quad model, norm(p) = %d\n',norm(p));
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
    rnew = residual(prey_equilib,unew,par);
    Jnew = Jacobian(prey_equilib,unew,par);
    mnew = 0.5*(r'*r) + g'*p + 0.5*p'*B*p;
    fnew = F(rnew);
    rho = (f - fnew + 1e-14)/(f - mnew + 1e-14);
    fprintf("rnew = %d, rho = %d\n",rnew,rho);

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
%%

function [u1min,u1max,u2min,u2max,pmin,pmax,period] = ...
    find_loop_min_max(prey_equilib,u,par)
% make one step to escape from the Poincare section
% step off the Poincare section
rtol = 1e-9;
atol = 1e-9;
Tmax = 5e2;
t0 = 1e-3;
options = odeset('RelTol',rtol,'AbsTol',atol);
odefun = @(t,x) ODE_prey_predator(x,par);    
[~,Y1] = ode45(odefun,[0,t0],[prey_equilib;u],options);
% find a point on the Poincare section
events_cycle = @(t,y)PoincareSection(y,prey_equilib,1);
options_cycle = odeset('RelTol',rtol,'AbsTol',atol,'Events',events_cycle);
y0 = Y1(end,:);
% find a point on the Poincare section
[T,Y2] = ode45(odefun,[t0,Tmax],y0,options_cycle);
period = T(end);
prey = [Y1(:,1);Y2(:,1)];
y1 = [Y1(:,2);Y2(:,2)];
y2 = [Y1(:,3);Y2(:,3)];
pmin = min(prey);
pmax = max(prey);
u1min = min(y1);
u1max = max(y1);
u2min = min(y2);
u2max = max(y2);
end
