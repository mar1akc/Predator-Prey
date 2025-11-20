function PlotDiagram()
fsz = 24; % fontsize for graphics
lwt = 1; % linewidth
tau_text = "Maturation age, \tau*";
g_text = "Consumption rate, g";
%% Select the kind of diagram you want to plot
% flag = 0: image plot of three kinds of attractors on the (tau,g)-plane
% flag = 1: phase diagram of three types of attractors on the (tau,g)-plane
% flag = 2: bifurcation diagram at the given g
% flag = 3 for g = 0.52;
flag = 4; 
tau = 1.5; % value of tau for the bifurcation diagram

%% The format of data files
%     fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
%     fprintf('tau = %.4f, g = %.4f,max_abs_eval = %.4e\n',tau,g,max(abs(evals)));
%     save(fname,'tau','g','x','u','u1','u2','evals');

%     fname = sprintf("Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g);
%     save(fname,"prey_equilib","prey","u","g","evals");
% If file is missing, the program will skip it and keep running

tau_vals = 0.1:0.05:2.0;
glist = [0.12,0.36,0.38,0.44,0.50,0.52]; % data are available for this list
g = 0.50;

gvals = linspace(0,1,101);
dirname = ["/Users/mariacameron/Dropbox/AgeStructured_nu1_Doron",...
    "/Users/mariacameron/Dropbox/AgeStructured_nu1_Doron_work",...
    "/Users/mariacameron/Dropbox/Work/My_Programs/PredatorPrey/AgeStructured_smallstep_nu1"];
dirname_idx = zeros(length(tau_vals),1);
dirname_idx(1:18) = 2;
dirname_idx(19:28) = 1;
dirname_idx(29:end) = 3;

N_gvals = length(gvals);
N_tauvals = length(tau_vals);
%% plot the diagram on tau-g-plane: image
if flag == 0
    
    Bfun = zeros(N_gvals,N_tauvals);
    Itau = ones(N_tauvals,1);
    
    for jtau = 1 : N_tauvals
        tau = tau_vals(jtau);
        g = 0;
        fname = strcat(dirname(dirname_idx(jtau)),sprintf('/Data/equilibrium_tau%.4f_g%.4f.mat',tau,g));
        % fprintf("dirname = %s\n",dirname(dirname_idx(jtau)));
        % fprintf("fname = %s\n",fname);
        if isfile(fname)
            % fprintf("jtau = %d, tau = %d, dirname_idx = %d\n",jtau,tau,dirname_idx(jtau))
            max_abs_evals_eq = -ones(N_gvals,1);
            max_abs_evals_cyc = -ones(N_gvals,1);
            for j = 1 : N_gvals
                g = gvals(j);
                % equilibrium point
                fname = strcat(dirname(dirname_idx(jtau)),sprintf('/Data/equilibrium_tau%.4f_g%.4f.mat',tau,g));
                if isfile(fname)
                    data_eq = load(fname);
                    u1 = data_eq.u1;
                    u2 = data_eq.u2;
                    if min(u1,u2) > 1e-9
                        evals_eq = data_eq.evals;
                        max_abs_evals_eq(j) = max(abs(evals_eq));
                        % periodic solution
                        fname = strcat(dirname(dirname_idx(jtau)),sprintf("/Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g));
                        if isfile(fname)
                            data_cyc = load(fname);
                            evals_cyc = data_cyc.evals;
                            max_abs_evals_cyc(j) = max(abs(evals_cyc));
                        end
                        if max_abs_evals_eq(j)<=1 
                            Bfun(j,jtau) = 1;
                        else
                            if max_abs_evals_cyc(j) <=1
                                Bfun(j,jtau) = 2;
                            end
                        end
                    end
                end
            end
        else
            Itau(jtau) = 0;
        end

    end
    Iremove = find(Itau == 0);
    Bfun(:,Iremove) = [];
    tau_vals(Iremove) = [];
    
    figure(1); clf; hold on; 
    imagesc(tau_vals,gvals,Bfun);
    set(gca,'Fontsize',fsz);


    xlabel("Maturation Age, \tau^*",'FontSize',fsz);
    ylabel(g_text,'Fontsize',fsz);
end
%% plot the diagram on tau-g-plane: curves
if flag == 1
    tau_vals = (0.1:0.05:2.0)';
    gvals = linspace(0,1,101)';
    gstep = gvals(2)-gvals(1);

    N_gvals = length(gvals);
    N_tauvals = length(tau_vals);
    
    Ebdry = -ones(N_tauvals,1);
    Cbdry = -ones(N_tauvals,1);
    Itau = ones(N_tauvals,1);
    
    for jtau = 1 : N_tauvals
        tau = tau_vals(jtau);
        g = 0;
        fname = strcat(dirname(dirname_idx(jtau)),sprintf('/Data/equilibrium_tau%.4f_g%.4f.mat',tau,g));
        if isfile(fname)
            max_abs_evals_eq = -ones(N_gvals,1);
            max_abs_evals_cyc = -ones(N_gvals,1);
            Ig = ones(N_gvals,1);
            for j = 1 : N_gvals
                g = gvals(j);
                % equilibrium point
                fname = strcat(dirname(dirname_idx(jtau)),sprintf('/Data/equilibrium_tau%.4f_g%.4f.mat',tau,g));

                if isfile(fname)
                    data_eq = load(fname);
                    u1 = data_eq.u1;
                    u2 = data_eq.u2;
                    if min(u1,u2) > 1e-9
                        evals_eq = data_eq.evals;
                        max_abs_evals_eq(j) = max(abs(evals_eq));
                        % periodic solution
                        fname = strcat(dirname(dirname_idx(jtau)),sprintf("/Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g));
                        if isfile(fname)
                            data_cyc = load(fname);
                            evals_cyc = data_cyc.evals;
                            max_abs_evals_cyc(j) = max(abs(evals_cyc));
                        end
                        Ebdry(jtau) = gvals(j) + 0.5*gstep;
                        if max_abs_evals_eq(j)>1  && max_abs_evals_cyc(j)<1
                            Cbdry(jtau) = gvals(j) + 0.5*gstep;
                        end
                   end
                else
                    Ig(j) = 0;
                end
            end
        else
            Itau(jtau) = 0;
        end

    end
    Iremove = find(Itau == 0);
    tau_vals(Iremove) = [];
    Cbdry(Iremove) = [];
    Ebdry(Iremove) = [];

    inan = find(~isfinite(Ebdry));
    Ebdry(max(inan)) = 1;
    
    figure(2); clf; hold on; 
    plot(tau_vals,Ebdry,'linewidth',lwt);
    plot(tau_vals,Cbdry,'linewidth',lwt);
    grid on;
    axis([min(tau_vals),max(tau_vals),min(gvals),max(gvals)])

    
    % Shade the zones
    % First, triangulate zones
    % Zone 1, below Cbdry
    j1 = find(Cbdry < 0, 1, 'last' );
    v1 = [tau_vals(j1)-Cbdry(j1)*(tau_vals(j1+1)-tau_vals(j1))/(Cbdry(j1+1)-Cbdry(j1)),0];
    vv = [tau_vals((j1+1):end),Cbdry((j1+1):end)];
    v2 = [tau_vals(end),min(gvals)];
    verts = [v1;vv;v2];
    patch(verts(:,1),verts(:,2),'red');
    alpha(0.4);

    % Zone 3, above Ebdry
    j3 = find(Ebdry > 1,1,'last');
    vv_e = [tau_vals(j3:end),Ebdry(j3:end)];
    v5 = [max(tau_vals),max(gvals)];
    verts = [vv_e;v5];
    patch(verts(:,1),verts(:,2),'blue');
    alpha(0.2);

    % Zone 2, between Cbdry and Ebdry
    % j2 = find(abs(Cbdry(j1+1:end) - Ebdry(j1+1:end)) < gstep,1,'first');
    % j2 = j1+j2;
    % vv = [tau_vals((j1+1):j2),Cbdry((j1+1):j2)];
    % j3 = find(Ebdry > 1,1,'last');
    % vvv = [tau_vals(j2:-1:j3),Ebdry(j2:-1:j3)];
    v3 = [min(tau_vals),max(gvals)];
    v4 = [min(tau_vals),min(gvals)];
    verts = [v4;v1;vv;v2;flipud(vv_e);v3;v4];
    patch(verts(:,1),verts(:,2),'green');
    alpha(0.4);

    % Zone 3, above Ebdry
    vv = [tau_vals(j3:end),Ebdry(j3:end)];
    v5 = [max(tau_vals),max(gvals)];
    verts = [vv;v5];
    patch(verts(:,1),verts(:,2),'blue');
    alpha(0.2);

    set(gca,'Fontsize',fsz);
    xlabel("Maturation Age, \tau*",'FontSize',fsz);
    ylabel(g_text,'Fontsize',fsz);
    saveas(gcf,"Figures/Zones","epsc")

    % Save regions' boundaries for plot_zones.m
    fbdry_name = "ZonesPDEsmallstep.mat";
    save(fbdry_name,"tau_vals","Ebdry","Cbdry","j1","j3");


end



%% plot the bifurcation diagram for given tau
%     fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
%     fprintf('tau = %.4f, g = %.4f,max_abs_eval = %.4e\n',tau,g,max(abs(evals)));
%     save(fname,'tau','g','x','u','u1','u2','evals');

%     fname = sprintf("Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g);
%     save(fname,"prey_equilib","prey","u","g","evals");
[~,jtau] = min(abs(tau_vals - tau));
if flag == 2
    fname_par = strcat(dirname(dirname_idx(jtau)),"/model_parameters.mat");
    dat = load(fname_par);
    % save("model_parameters.mat","L","tau","h_age","dt","par_ifun","r","a",
    % "k","bb","s","z","Mm","rho", "mu_b_par", "birth_par", "birth_exp_par", "death_exp_par");
    par = [dat.L, dat.tau, dat.h_age, dat.dt, dat.par_ifun, dat.r, dat.a, dat.k,...
        dat.bb, dat.s, dat.z, dat.Mm, dat.rho,dat.mu_b_par,dat.birth_par,...
        dat.birth_exp_par,dat.death_exp_par];

    fname = sprintf("Data/FindCycleOutput_tau%.4f.mat",tau);
    data = load(fname);
    gvals = data.gvals; 
    u1_eq = data.u1_eq; 
    u2_eq = data.u2_eq; 
    prey_eq = data.prey_eq; 
    u1min = data.u1_loop_min;
    u1max = data.u1_loop_max;
    u2min = data.u2_loop_min;
    u2max = data.u2_loop_max;
    pmin = data.prey_loop_min;
    pmax = data.prey_loop_max;

    max_abs_evals_eq = -ones(length(gvals),1);
    Ig = ones(length(gvals),1);
    for j = 1 : length(gvals)
        g = gvals(j);
        % equilibrium point
        fname = strcat(dirname(dirname_idx(jtau)),sprintf('/Data/equilibrium_tau%.4f_g%.4f.mat',tau,g));
        if isfile(fname)
            data_eq = load(fname);
            u1 = data_eq.u1;
            u2 = data_eq.u2;
            if min(u1,u2) > 1e-9
                evals_eq = data_eq.evals;
                max_abs_evals_eq(j) = max(abs(evals_eq));
            else
                Ig(j) = 0;
            end
        else
            Ig(j) = 0;
        end
    end
I_eq_stable = Ig == 1 & max_abs_evals_eq <= 1 & max_abs_evals_eq >= 0;
I_eq_unstable = find(Ig == 0);
Iset = zeros(length(gvals),1);
Iset(I_eq_stable) = 1;
i = 1;
legend_count = 1;

color_prey = [0,0,1]; % blue
color_u1 =  [1,0,0]; % red
color_u2 = [1,0.8,0.2]; % dark yellow

figure(4);clf; hold on;

N_gvals = length(gvals);
while i < N_gvals
    stab_type = Iset(i);
    i_start = i;
    while i < N_gvals && Iset(i) == stab_type
        i = i+1;
    end
    if Iset(i) ~=stab_type
        i_finish = i - 1;
    else
        i_finish = N_gvals;
    end

    if stab_type == 1
        lstyle = '-';
    else
        lstype = '--';
    end
    fprintf("legend_count = %d, stab_type = %d,i_start = %d, ifinish = %d\n",legend_count,stab_type,i_start,i_finish);
    
    I = (i_start:i_finish);
    if stab_type == 1
        hand(legend_count) = plot(gvals(I),prey_eq(I),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(legend_count+1) = plot(gvals(I),u1_eq(I),'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
        hand(legend_count+2) = plot(gvals(I),u2_eq(I),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    else
        hand(legend_count) = plot(gvals(I),prey_eq(I),'--','Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(legend_count+1) = plot(gvals(I),u1_eq(I),'--','Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
        hand(legend_count+2) = plot(gvals(I),u2_eq(I),'--','Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    end 
    legend_count = legend_count + 3;
   

    if stab_type == 0
        if I(1) > 1
            I0 = I(1)-1;
            I = [I0,I];
            u1min(I0) = u1_eq(I0);
            u2min(I0) = u2_eq(I0);
            u1max(I0) = u1_eq(I0);
            u2max(I0) = u2_eq(I0);
            pmin(I0) = prey_eq(I0);
            pmax(I0) = prey_eq(I0);
        end
        if I(end) < N_tauvals
            I0 = I(end) + 1;
            I = [I,I0];
            u1min(I0) = u1_eq(I0);
            u2min(I0) = u2_eq(I0);
            u1max(I0) = u1_eq(I0);
            u2max(I0) = u2_eq(I0);
            pmin(I0) = prey_eq(I0);
            pmax(I0) = prey_eq(I0);
        end
        % hand(legend_count) = plot(gvals(I),pmin(I),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        % hand(legend_count+1) = plot(gvals(I),pmax(I),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        % 
        % hand(legend_count+2) = plot(gvals(I),u1min(I),'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
        % hand(legend_count+3) = plot(gvals(I),u1max(I),'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
        % 
        % hand(legend_count+4) = plot(gvals(I),u2min(I),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
        % hand(legend_count+5) = plot(gvals(I),u2max(I),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
         hand(legend_count) = plot(gvals,pmin,'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(legend_count+1) = plot(gvals,pmax,'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');

        hand(legend_count+2) = plot(gvals,u1min,'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
        hand(legend_count+3) = plot(gvals,u1max,'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');

        hand(legend_count+4) = plot(gvals,u2min,'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
        hand(legend_count+5) = plot(gvals,u2max,'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
       
        legend_count = legend_count + 6;
    end
end
if Iset(1) == 1
    Ilegend = [1,2,3];
else
    Ilegend = [4,6,8];
end
    set(gca,'Fontsize',fsz);
    xlabel(g_text,'FontSize',fsz);
    ylabel('Population Size','Fontsize',fsz);
    legend(hand(Ilegend))
    grid on;
    % xlim([min(gvals(Ig==1)),max(gvals(Ig==1))]);
    title(strcat("\tau* = ",sprintf("%.2f",tau)),'FontSize',fsz);
end


%% 
if flag == 3
    f_name = strcat(dirname,sprintf("FindCycleOutput_g%.4f.mat",g));
    % save(strcat(dirname,"FindCycleOutput.mat"),"g","tau_vals","u1_eq","u2_eq","prey_eq","u1_loop_min",...
    % "u1_loop_max","u2_loop_min","u2_loop_max","prey_loop_min","prey_loop_max",...
    % "period");
    data = load(f_name);
    g = data.g;
    tau_vals = data.tau_vals;
    N_tauvals = length(tau_vals);
    prey_eq = data.prey_eq;
    u1_eq = data.u1_eq;
    u2_eq = data.u2_eq;
    u1max = data.u1_loop_max;
    u1min = data.u1_loop_min;
    u2max = data.u2_loop_max;
    u2min = data.u2_loop_min;
    pmax = data.prey_loop_max;
    pmin = data.prey_loop_min;

    tol = 1e-2;
    ind = find(pmax-pmin <1e-2);

Iset = ones(N_tauvals,1);
Iset(ind) = 0;
i = 1;
legend_count = 1;

color_prey = [0,0,1]; % blue
color_u1 =  [1,0,0]; % red
color_u2 = [1,0.8,0.2]; % dark yellow

figure(4);clf; hold on;

while i < N_tauvals
    stab_type = Iset(i);
    i_start = i;
    while i < N_tauvals && Iset(i) == stab_type
        i = i+1;
    end
    if Iset(i) ~=stab_type
        i_finish = i - 1;
    else
        i_finish = N_tauvals;
    end

    if stab_type == 1
        lstyle = '-';
    else
        lstype = '--';
    end
    fprintf("legend_count = %d, stab_type = %d,i_start = %d, ifinish = %d\n",legend_count,stab_type,i_start,i_finish);
    
    I = (i_start:i_finish);
    if stab_type == 1
        hand(legend_count) = plot(tau_vals(I),prey_eq(I),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(legend_count+1) = plot(tau_vals(I),u1_eq(I),'Linewidth',lwt,'color',color_u1,'DisplayName','Juvenile predator');
        hand(legend_count+2) = plot(tau_vals(I),u2_eq(I),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    else
        hand(legend_count) = plot(tau_vals(I),prey_eq(I),'--','Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(legend_count+1) = plot(tau_vals(I),u1_eq(I),'--','Linewidth',lwt,'color',color_u1,'DisplayName','Juvenile predator');
        hand(legend_count+2) = plot(tau_vals(I),u2_eq(I),'--','Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    end 
    legend_count = legend_count + 3;
   

    if stab_type == 0
        if I(1) > 1
            I0 = I(1)-1;
            I = [I0,I];
            u1min(I0) = u1_eq(I0);
            u2min(I0) = u2_eq(I0);
            u1max(I0) = u1_eq(I0);
            u2max(I0) = u2_eq(I0);
            pmin(I0) = prey_eq(I0);
            pmax(I0) = prey_eq(I0);
        end
        if I(end) < N_tauvals
            I0 = I(end) + 1;
            I = [I,I0];
            u1min(I0) = u1_eq(I0);
            u2min(I0) = u2_eq(I0);
            u1max(I0) = u1_eq(I0);
            u2max(I0) = u2_eq(I0);
            pmin(I0) = prey_eq(I0);
            pmax(I0) = prey_eq(I0);
        end
        hand(legend_count) = plot(tau_vals(I),pmin(I),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(legend_count+1) = plot(tau_vals(I),pmax(I),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');

        hand(legend_count+2) = plot(tau_vals(I),u1min(I),'Linewidth',lwt,'color',color_u1,'DisplayName','Juvenile predator');
        hand(legend_count+3) = plot(tau_vals(I),u1max(I),'Linewidth',lwt,'color',color_u1,'DisplayName','Juvenile predator');

        hand(legend_count+4) = plot(tau_vals(I),u2min(I),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
        hand(legend_count+5) = plot(tau_vals(I),u2max(I),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
        
        legend_count = legend_count + 6;
    end
end
if Iset(1) == 1
    Ilegend = [1,2,3];
else
    Ilegend = [4,6,8];
end
    set(gca,'Fontsize',fsz);
    xlabel(tau_text,'FontSize',fsz);
    ylabel('Population Size','Fontsize',fsz);
    legend(hand(Ilegend))
    grid on;
    title(strcat("g = ",sprintf("%.2f",g)),'FontSize',fsz);
end

%%
if flag == 4
figure(4);clf; hold on;

    for jtau = 1 : N_tauvals
        tau = tau_vals(jtau);
    %% plot the bifurcation diagram for given tau
%     fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
%     fprintf('tau = %.4f, g = %.4f,max_abs_eval = %.4e\n',tau,g,max(abs(evals)));
%     save(fname,'tau','g','x','u','u1','u2','evals');

%     fname = sprintf("Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g);
%     save(fname,"prey_equilib","prey","u","g","evals");
[~,jtau] = min(abs(tau_vals - tau));
    fname_par = strcat(dirname(dirname_idx(jtau)),"/model_parameters.mat");
    dat = load(fname_par);
    % save("model_parameters.mat","L","tau","h_age","dt","par_ifun","r","a",
    % "k","bb","s","z","Mm","rho", "mu_b_par", "birth_par", "birth_exp_par", "death_exp_par");
    par = [dat.L, dat.tau, dat.h_age, dat.dt, dat.par_ifun, dat.r, dat.a, dat.k,...
        dat.bb, dat.s, dat.z, dat.Mm, dat.rho,dat.mu_b_par,dat.birth_par,...
        dat.birth_exp_par,dat.death_exp_par];

    fname = sprintf("Data_FindCycleOutput/FindCycleOutput_tau%.4f.mat",tau);
    data = load(fname);
    gvals = data.gvals; 
    u1_eq = data.u1_eq; 
    u2_eq = data.u2_eq; 
    prey_eq = data.prey_eq; 
    u1min = data.u1_loop_min;
    u1max = data.u1_loop_max;
    u2min = data.u2_loop_min;
    u2max = data.u2_loop_max;
    pmin = data.prey_loop_min;
    pmax = data.prey_loop_max;

N_gvals = length(gvals);
tvec = tau*ones(size(gvals));

    max_abs_evals_eq = -ones(length(gvals),1);
    Ig = ones(length(gvals),1);
    for j = 1 : length(gvals)
        g = gvals(j);
        % equilibrium point
        fname = strcat(dirname(dirname_idx(jtau)),sprintf('/Data/equilibrium_tau%.4f_g%.4f.mat',tau,g));
        if isfile(fname)
            data_eq = load(fname);
            u1 = data_eq.u1;
            u2 = data_eq.u2;
            if min(u1,u2) > 1e-9
                evals_eq = data_eq.evals;
                max_abs_evals_eq(j) = max(abs(evals_eq));
            else
                Ig(j) = 0;
            end
        else
            Ig(j) = 0;
        end
    end
I_eq_stable = Ig == 1 & max_abs_evals_eq <= 1 & max_abs_evals_eq >= 0;
I_eq_unstable = find(Ig == 0);
Iset = zeros(length(gvals),1);
Iset(I_eq_stable) = 1;
i = 1;
legend_count = 1;

color_prey = [0,0,1]; % blue
color_u1 =  [1,0,0]; % red
color_u2 = [0,0.5,0]; %[1,0.8,0.2]; % dark yellow


N_gvals = length(gvals);
while i < N_gvals
    stab_type = Iset(i);
    i_start = i;
    while i < N_gvals && Iset(i) == stab_type
        i = i+1;
    end
    if Iset(i) ~=stab_type
        i_finish = i - 1;
    else
        i_finish = N_gvals;
    end

    if stab_type == 1
        lstyle = '-';
    else
        lstype = '--';
    end
    fprintf("legend_count = %d, stab_type = %d,i_start = %d, ifinish = %d\n",legend_count,stab_type,i_start,i_finish);
    
    I = (i_start:i_finish);
    if stab_type == 1
        hand(legend_count) = plot3(gvals(I),tvec(I),prey_eq(I),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(legend_count+1) = plot3(gvals(I),tvec(I),u1_eq(I),'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
        hand(legend_count+2) = plot3(gvals(I),tvec(I),u2_eq(I),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    % else
    %     hand(legend_count) = plot3(tvec(I),gvals(I),prey_eq(I),'--','Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
    %     hand(legend_count+1) = plot3(tvec(I),gvals(I),u1_eq(I),'--','Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
    %     hand(legend_count+2) = plot3(tvec(I),gvals(I),u2_eq(I),'--','Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    end 
    legend_count = legend_count + 3;
   

    if stab_type == 0
        if I(1) > 1
            I0 = I(1)-1;
            I = [I0,I];
            u1min(I0) = u1_eq(I0);
            u2min(I0) = u2_eq(I0);
            u1max(I0) = u1_eq(I0);
            u2max(I0) = u2_eq(I0);
            pmin(I0) = prey_eq(I0);
            pmax(I0) = prey_eq(I0);
        end
        if I(end) < N_tauvals
            I0 = I(end) + 1;
            I = [I,I0];
            u1min(I0) = u1_eq(I0);
            u2min(I0) = u2_eq(I0);
            u1max(I0) = u1_eq(I0);
            u2max(I0) = u2_eq(I0);
            pmin(I0) = prey_eq(I0);
            pmax(I0) = prey_eq(I0);
        end
        % hand(legend_count) = plot(gvals(I),pmin(I),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        % hand(legend_count+1) = plot(gvals(I),pmax(I),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        % 
        % hand(legend_count+2) = plot(gvals(I),u1min(I),'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
        % hand(legend_count+3) = plot(gvals(I),u1max(I),'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
        % 
        % hand(legend_count+4) = plot(gvals(I),u2min(I),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
        % hand(legend_count+5) = plot(gvals(I),u2max(I),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
         hand(legend_count) = plot3(gvals,tvec,pmin,'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(legend_count+1) = plot3(gvals,tvec,pmax,'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');

        hand(legend_count+2) = plot3( gvals,tvec,u1min,'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
        hand(legend_count+3) = plot3(gvals,tvec,u1max,'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');

        hand(legend_count+4) = plot3(gvals,tvec,u2min,'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
        hand(legend_count+5) = plot3( gvals,tvec,u2max,'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
       
        legend_count = legend_count + 6;
    end
end
if Iset(1) == 1
    Ilegend = [1,2,3];
else
    Ilegend = [4,6,8];
end
    set(gca,'Fontsize',fsz);
    ylabel("\tau*",'FontSize',fsz);
    xlabel("g",'FontSize',fsz);
    zlabel('Population Size','Fontsize',fsz);
    legend(hand(Ilegend))
    grid on;
    % xlim([min(gvals(Ig==1)),max(gvals(Ig==1))]);
    % title(strcat("\tau* = ",sprintf("%.2f",tau)),'FontSize',fsz);
    view(3)
end



end
end
% %% graphics
% if graphix_flag == 1
% 
%     fsz = 20;
% 
%     color_prey = [0,0,1]; % blue
%     color_u1 =  [1,0,0]; % red
%     color_u2 = [1,0.8,0.2]; % dark yellow
%     figure(1);clf; hold on;
%     hand = zeros(1,9);
%     hand(1) = plot(tau_vals,prey_eq,'Linewidth',2,'color',color_prey,'DisplayName','prey');
%     hand(2) = plot(tau_vals,prey_loop_min,'Linewidth',2,'color',color_prey,'DisplayName','prey');
%     hand(3) = plot(tau_vals,prey_loop_max,'Linewidth',2,'color',color_prey,'DisplayName','prey');
%     hand(4) = plot(tau_vals,u1_eq,'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
%     hand(5) = plot(tau_vals,u1_loop_min,'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
%     hand(6) = plot(tau_vals,u1_loop_max,'Linewidth',2,'color',color_u1,'DisplayName','Juv. predator');
%     hand(7) = plot(tau_vals,u2_eq,'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
%     hand(8) = plot(tau_vals,u2_loop_min,'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
%     hand(9) = plot(tau_vals,u2_loop_max,'Linewidth',2,'color',color_u2,'DisplayName','Adult predator');
%     set(gca,'Fontsize',fsz);
%     xlabel('Maturation age, \tau^*','FontSize',fsz);
%     ylabel('Population size','Fontsize',fsz);
%     legend(hand(1:3:7))
% end
% 
% 
% 
% end


