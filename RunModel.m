function u = RunModel(tau,g)

% Plot graphs is graphix_flag == 1, else do not plot
graphix_flag = 0;

%%
T = 0.5e3; % max time

%% Parameters
[L,h_age,dt,age_vector,AgeNode,par_ifun,r,a,k,bb,s,rho,Mm,...
    z,birth_par,birth_exp_par,death_exp_par,mu_b_par] = assign_parameters(tau);

[I1,I2,mu_base,birth_rate_base] = assign_parameter_functions(par_ifun,...
    tau,L,mu_b_par,death_exp_par,birth_par,birth_exp_par);
% Save parameters
save("model_parameters.mat","L","tau","h_age","dt","par_ifun","r","a","k",...
    "bb","s","z","Mm","rho", "mu_b_par", "birth_par", "birth_exp_par", "death_exp_par");
dat = load("model_parameters.mat");

%%
t = 0:dt:T; % time vector

N_age = length(age_vector);
Nt = length(t);
u = zeros(N_age,1);
u1 = zeros(Nt,1);
u2 = zeros(Nt,1);
prey = zeros(Nt,1);


% Boundary conditions
prey(1) = 0.5;
% % Initial age distribution for the predator
initial_age_density = @(x) 0.1*heaviside(tau-x) + 0.05*heaviside(x-tau);
% Initial condition: the initial age distribution
u = initial_age_density(age_vector(:));

h05 = 0.5*h_age;

u1(1) = h05*(u(1)+u(AgeNode))+h_age*sum(u(2:AgeNode-1));
u2(1) = h05*(u(AgeNode)+u(end))+h_age*sum(u(AgeNode+1:end-1));

for i=1:Nt-1
    [prey(i+1),u,u1(i+1),u2(i+1)] = time_step_helper(prey(i),u,u1(i),u2(i),g,h_age,dt,age_vector,...
        AgeNode,r,a,k,bb,s,rho,z,Mm,I1,I2,mu_base,birth_rate_base);
end

prey_end = prey(end);

fname = sprintf('Data/Initial_guess_g%.2f.mat',g);
save(fname,'prey_end','u','g');


%% graphics
if graphix_flag == 1
    fsz = 20; % fontsize
    % population vs time
    figure(1);clf; hold on;
    plot(t,prey,'LineWidth',2,'displayname','Prey, x')
    plot(t,u1,'LineWidth',2,'displayname','Juvenile Predator, y_1')
    plot(t,u2,'LineWidth',2,'displayname','Adult Predator, y_2')
    legend('Location','northeast') 
    set(gca,'Fontsize',fsz)
    xlabel('time, t','Fontsize',fsz)
    ylabel('Population Size','Fontsize',fsz)
    title([' g = ', num2str(g) ', \tau* = ',num2str(tau) ] )

    % phase space
    I = floor(Nt*0.9)+1:Nt;

    figure(2);clf;
    plot3(prey(I),u1(I),u2(I),'Linewidth',3)
    grid on
    xlabel('Prey, x','Fontsize',fsz);
    ylabel('Juvenile predator, y_1','Fontsize',fsz)
    zlabel('Adult predator, y_2','Fontsize',fsz)
    set(gca,'Fontsize',fsz)
    title([' g = ', num2str(g) ', \tau* = ',num2str(tau) ] )
end
end