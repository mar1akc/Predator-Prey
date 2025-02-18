function EqulibriumSolutionsDDE(tau)
par_ifun = 100;
stabcheck_tol = 1e-2;
%% parameters that do not change
dat = load("DDEmodel_parameters.mat");
MB = (dat.mu_b_par/dat.death_exp_par)*(exp(dat.death_exp_par*tau) - 1);
%% parameters
par = [dat.r,dat.a,dat.k,dat.bb,0,dat.s,0,0,0,0,dat.rho,dat.Mm,dat.z,0,MB];

fname = sprintf("ODEparams/ODEparams_tau%.2f.mat",tau);
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

%% solving DDE
% Find the equilibrium, determine if it is table or unstable, and find the
% stable limit cycle if the equilibrium is unstable

par(5) = gvals(1);
par(7) = B(1);
par(8) = M1(1);
par(9) = M2(1);
par(10) = D(1);
par(14) = tau;
% find initial approximation to search for the equilibrium
[x0,u10,u20] = AgeStructuredPreyPredator4DDE(tau,gvals(1));
init_cond = [x0,u10,u20];

ddefun = @(t,x,xdelay) DDE_prey_predator(x,xdelay,par);    
rtol = 1e-9;
atol = 1e-9;
Tmax = 1e4;
delay = tau;
options = ddeset('RelTol',rtol,'AbsTol',atol);

% In dde23, "history" refers to the solution values of the equation 
% for times prior to the initial integration point (t0), 
% essentially providing the necessary past information needed to 
% calculate the solution at subsequent time points where the 
% equation depends on delayed values. 
history = init_cond;

sol = dde23(ddefun,delay,history,[0,Tmax],options);
Y = sol.y;
T = sol.x;
yeq_iguess = Y(:,end);
figure(2); hold on;
plot(T,Y,'LineWidth',2)
drawnow;
for i = 1 : Ngvals-1    
    par(5) = gvals(i);
    par(7) = B(i);
    par(8) = M1(i);
    par(9) = M2(i);
    par(10) = D(i);

%% find equilibrium solutions

    fun = @(x) DDE_prey_predator(x,x,par);
    y = lsqnonlin(fun,yeq_iguess);
    yeq_iguess = y;
    fprintf("tau = %.2f, g = %d, b2 = %d, m1 = %d, m2 = %d, D = %d\n",tau,gvals(i),B(i),M1(i),M2(i),D(i));
    fprintf("prey = %d, y1 = %d, y2 = %d\n",y(1),y(2),y(3));

    prey_equilib(i) = y(1);
    y1_equilib(i) = y(2);
    y2_equilib(i) = y(3);

    % Stability of DDEs
    % https://www.math.fsu.edu/~bertram/lectures/delay.pdf

    % We check stability experimentally
    fac = 1;
    sol = dde23(ddefun,delay,y + fac*[1;1;1],[0,Tmax],options);
    endpt = sol.y(:,end);

    figure(2); clf; hold on
    plot(sol.x,sol.y,'LineWidth',2)
    plot(sol.x,y*ones(size(sol.x)),'LineWidth',2);
    drawnow;

    fprintf("tau = %d, g = %d, norm(y-endpt) = %d\n",tau,gvals(i),norm(y-endpt));
    if norm(y-endpt) < stabcheck_tol
        stab(i) = 1;
        fprintf("g = %d, equilib is stable\n",gvals(i));
    
    else
        stab(i) = 0;
        fprintf("g = %d, equilib is unstable\n",gvals(i));
    end
    %% find limit cycle if the equilibrium is unstable
    if stab(i) == 0 % unstable equilibrium
        % find the initial approximation
        ddefun = @(t,x,xdelay) DDE_prey_predator(x,xdelay,par);    

        % fprintf("Start with: (%d,%d,%d)\n",...
        %     y(1),y(2),y(3));

        sol = dde23(ddefun,delay,y + fac*[1;1;1],[0,Tmax],options);
        % [T,Y] = ode45(odefun,[0,Tmax],y.*[1.1,1.1,1.1],options);
        % figure(101); hold on;
        % plot(sol.x,sol.y)
        % Y = sol.y;
        % find a point on the Poincare section

        % events_cycle = @(t,y,ydelay)PoincareSection(y,prey_equilib(i),1);
        % 
        % % options_cycle = odeset('RelTol',rtol,'AbsTol',atol,'Events',events_cycle);
        % options_cycle = ddeset('RelTol',rtol,'AbsTol',atol,'Events',events_cycle);
        T = sol.x;
        nT = length(T);
        % need a nontrivial history
        % find a point on the Poincare section
        sol = dde23(ddefun,delay,sol,[Tmax,Tmax + 1e2],options);
        % T = sol.x;
        % nT = length(T);
        % sol = dde23(ddefun,delay,sol,[sol.x(end),sol.x(end)+2e-3],options_cycle);
        Y = sol.y(:,nT:end);
        T = sol.x(nT:end);

        figure(3);
        plot3(Y(1,:),Y(2,:),Y(3,:));
        drawnow;
        

        fname = sprintf("DataODE_DDE/LimitCycleDDE_tau%.4f_g%.4f.mat",tau,gvals(i));
        save(fname,"T","Y");
            % find min and max vals along the limit cycle
        y1_loop_min(i) = min(Y(2,:));
        y1_loop_max(i) = max(Y(2,:));
        y2_loop_min(i) = min(Y(3,:));
        y2_loop_max(i) = max(Y(3,:));
        prey_loop_min(i) = min(Y(1,:));
        prey_loop_max(i) = max(Y(1,:));
        period(i) = T(end)-T(1);
        if max([y1_loop_max(i)-y1_loop_min(i),...
                y2_loop_max(i)-y2_loop_min(i),...
                prey_loop_max(i)-prey_loop_min(i)]) < stabcheck_tol
            stab(i) = 1;
            
            fprintf("period = %d, length(T) = %d, loop is degenerate\n",period(i),length(T));
            fprintf("g = %d, equilib is stable\n",gvals(i));
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

fname = sprintf("DataODE_DDE/DDE_BifurDiag_par_ifun%.1f_tau%.2f.mat",par_ifun,tau);

save(fname,"gvals","prey_equilib","y1_equilib","y2_equilib",...
    "prey_loop_min","prey_loop_max","y1_loop_min","y1_loop_max",...
    "y2_loop_min","y2_loop_max","period","stab");


Istab = find(stab == 1);
Inonstab = find(stab == 0);

    
fsz = 20;

color_prey = [0,0,1]; % blue
color_u1 =  [1,0,0]; % red
color_u2 = [1,0.8,0.2]; % dark yellow

figure(4);clf; hold on;

if isempty(Inonstab) == 0 

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
    xlabel('Juv. Pred. Damage Parameter, g','FontSize',fsz);
    ylabel('Population Size','Fontsize',fsz);
    legend(hand([1,5,9]))
    title(sprintf("Maturation Age = %.2f",tau),'FontSize',fsz);
    grid on

else
    plot(gvals(Istab),prey_equilib(Istab),'Linewidth',2,'color',color_prey,'DisplayName','prey');
    plot(gvals(Istab),y1_equilib(Istab),'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
    plot(gvals(Istab),y2_equilib(Istab),'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
    set(gca,'Fontsize',fsz)
    set(gcf,'color','white')
    legend('Location','northeast') 
    xlabel('Juv. Pred. Damage Parameter, g','FontSize',fsz);
    ylabel('Population Size','Fontsize',fsz);
    title(sprintf("Maturation Age = %.2f",tau),'FontSize',fsz);
    grid on
end

end
