function PlotDiagram()
fsz = 18; % fontsize for graphics
lwt = 3; % linewidth
g_text = "Consumption Rate of Juvenile Predator, g";
%% Select the kind of diagram you want to plot
% flag = 0: image plot of three kinds of attractors on the (tau,g)-plane
% flag = 1: phase diagram of three types of attractors on the (tau,g)-plane
% flag = 2: bifurcation diagram at the given tau
flag = 1; 
tau = 1; % value of tau for the bifurcation diagram

%% The format of data files
%     fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
%     fprintf('tau = %.4f, g = %.4f,max_abs_eval = %.4e\n',tau,g,max(abs(evals)));
%     save(fname,'tau','g','x','u','u1','u2','evals');

%     fname = sprintf("Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g);
%     save(fname,"prey_equilib","prey","u","g","evals");
% If file is missing, the program will skip it and keep running

tauvals = 0.1:0.1:2.6;
gvals = linspace(0,1,101);

N_gvals = length(gvals);
N_tauvals = length(tauvals);
%% plot the diagram on tau-g-plane: image
if flag == 0
    
    Bfun = zeros(N_gvals,N_tauvals);
    Itau = ones(N_tauvals,1);
    
    for jtau = 1 : N_tauvals
        tau = tauvals(jtau);
        g = 0;
        fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
        if isfile(fname)
            max_abs_evals_eq = -ones(N_gvals,1);
            max_abs_evals_cyc = -ones(N_gvals,1);
            for j = 1 : N_gvals
                g = gvals(j);
                % equilibrium point
                fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
                if isfile(fname)
                    data_eq = load(fname);
                    u1 = data_eq.u1;
                    u2 = data_eq.u2;
                    if min(u1,u2) > 1e-9
                        evals_eq = data_eq.evals;
                        max_abs_evals_eq(j) = max(abs(evals_eq));
                        % periodic solution
                        fname = sprintf("Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g);
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
    tauvals(Iremove) = [];
    
    figure(1); clf; hold on; 
    imagesc(tauvals,gvals,Bfun);
    set(gca,'Fontsize',fsz);


    xlabel("Maturation Age, \tau^*",'FontSize',fsz);
    ylabel(g_text,'Fontsize',fsz);
end
%% plot the diagram on tau-g-plane: curves
if flag == 1
    tauvals = (0.1:0.01:2.6)';
    gvals = linspace(0,1,101)';
    gstep = gvals(2)-gvals(1);

    N_gvals = length(gvals);
    N_tauvals = length(tauvals);
    
    Ebdry = -ones(N_tauvals,1);
    Cbdry = -ones(N_tauvals,1);
    Itau = ones(N_tauvals,1);
    
    for jtau = 1 : N_tauvals
        tau = tauvals(jtau);
        g = 0;
        fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
        if isfile(fname)
            max_abs_evals_eq = -ones(N_gvals,1);
            max_abs_evals_cyc = -ones(N_gvals,1);
            Ig = ones(N_gvals,1);
            for j = 1 : N_gvals
                g = gvals(j);
                % equilibrium point
                fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
                if isfile(fname)
                    data_eq = load(fname);
                    u1 = data_eq.u1;
                    u2 = data_eq.u2;
                    if min(u1,u2) > 1e-9
                        evals_eq = data_eq.evals;
                        max_abs_evals_eq(j) = max(abs(evals_eq));
                        % periodic solution
                        fname = sprintf("Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g);
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
%             I_eq_stable = find(max_abs_evals_eq <= 1 & max_abs_evals_eq >= 0);
%             I_eq_unstable = find(max_abs_evals_eq > 1 & max_abs_evals_cyc <= 1);
%             if min(Ig) == 1
%                 Ebdry(jtau) = NaN;
%             end
%             if ~isempty(I_eq_unstable)
%                 Istar = max(I_eq_unstable);
%                 Cbdry(jtau) = gvals(Istar) + 0.5*gstep;
%             end
        else
            Itau(jtau) = 0;
        end

    end
    Iremove = find(Itau == 0);
    tauvals(Iremove) = [];
    Cbdry(Iremove) = [];
    Ebdry(Iremove) = [];

    inan = find(~isfinite(Ebdry));
    Ebdry(max(inan)) = 1;
    
    figure(2); clf; hold on; 
    plot(tauvals,Ebdry,'linewidth',lwt);
    plot(tauvals,Cbdry,'linewidth',lwt);
    grid on;
    axis([min(tauvals),max(tauvals),min(gvals),max(gvals)])

    
    % Shade the zones
    % First, triangulate zones
    % Zone 1, below Cbdry
    j1 = find(Cbdry < 0, 1, 'last' );
    v1 = [tauvals(j1)-Cbdry(j1)*(tauvals(j1+1)-tauvals(j1))/(Cbdry(j1+1)-Cbdry(j1)),0];
    vv = [tauvals((j1+1):end),Cbdry((j1+1):end)];
    v2 = [tauvals(end),min(gvals)];
    verts = [v1;vv;v2];
    patch(verts(:,1),verts(:,2),'red');
    alpha(0.4);
    
    % Zone 2, between Cbdry and Ebdry
    j2 = find(abs(Cbdry(j1+1:end) - Ebdry(j1+1:end)) < gstep,1,'first');
    j2 = j1+j2;
    vv = [tauvals((j1+1):j2),Cbdry((j1+1):j2)];
    j3 = find(Ebdry > 1,1,'last');
    vvv = [tauvals(j2:-1:j3),Ebdry(j2:-1:j3)];
    v3 = [min(tauvals),max(gvals)];
    v4 = [min(tauvals),min(gvals)];
    verts = [v1;vv;vvv;v3;v4];
    patch(verts(:,1),verts(:,2),'green');
    alpha(0.4);

    % Zone 3, above Ebdry
    vv = [tauvals(j3:end),Ebdry(j3:end)];
    v5 = [max(tauvals),max(gvals)];
    verts = [vv;v5];
    patch(verts(:,1),verts(:,2),'blue');
    alpha(0.2);

    set(gca,'Fontsize',fsz);
    xlabel("Maturation Age, \tau*",'FontSize',fsz);
    ylabel(g_text,'Fontsize',fsz);
    saveas(gcf,"Figures/Zones","epsc")

    % Save regions' boundaries for plot_zones.m
    fbdry_name = "ZonesPDE.mat";
    save(fbdry_name,"tauvals","Ebdry","Cbdry","j2","j1");


end



%% plot the bifurcation diagram for given tau
%     fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
%     fprintf('tau = %.4f, g = %.4f,max_abs_eval = %.4e\n',tau,g,max(abs(evals)));
%     save(fname,'tau','g','x','u','u1','u2','evals');

%     fname = sprintf("Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g);
%     save(fname,"prey_equilib","prey","u","g","evals");

if flag == 2

    dat = load("model_parameters.mat");
    par = [dat.L, dat.tau, dat.h_age, dat.dt, dat.par_ifun, dat.r, dat.a, dat.k,...
    dat.bb, dat.s, dat.z, dat.Mm, dat.rho,dat.mu_b_par,dat.birth_par,...
    dat.birth_exp_par,dat.death_exp_par];
    %
%     tau = 0.2;
    par(2) = tau;
    %
    prey_eq = zeros(N_gvals,1);
    u1_eq = zeros(N_gvals,1);
    u2_eq = zeros(N_gvals,1);
    u1max = zeros(N_gvals,1);
    u1min = zeros(N_gvals,1);
    u2max = zeros(N_gvals,1);
    u2min = zeros(N_gvals,1);
    pmax = zeros(N_gvals,1);
    pmin = zeros(N_gvals,1);
    Icycle = zeros(N_gvals,1);
    max_abs_evals_eq = -ones(N_gvals,1);
    max_abs_evals_cyc = -ones(N_gvals,1);
    Ig = ones(N_gvals,1);
    for j = 1 : N_gvals
        g = gvals(j);
        % equilibrium point

        fname = sprintf('Data/equilibrium_tau%.4f_g%.4f.mat',tau,g);
        if isfile(fname)
            data_eq = load(fname);
            u1 = data_eq.u1;
            u2 = data_eq.u2;
            if min(u1,u2) > 1e-9
                evals_eq = data_eq.evals;
                max_abs_evals_eq(j) = max(abs(evals_eq));
                prey_eq(j) = data_eq.x;
                u1_eq(j) = data_eq.u1;
                u2_eq(j) = data_eq.u2;
                % periodic solution
                fname = sprintf("Data/LimitCycle_tau%.4f_g%.4f.mat",tau,g);
                if isfile(fname)
                    data_cyc = load(fname);
                    evals_cyc = data_cyc.evals;
                    max_abs_evals_cyc(j) = max(abs(evals_cyc));
                    Icycle(j) = 1;
                    u = data_cyc.u;
                    prey = data_cyc.prey;
                    [u1min(j),u1max(j),u2min(j),u2max(j),pmin(j),pmax(j),period] = ...
                            find_loop_min_max(prey_eq(j),prey,u,g,par);
                end
            else
                Ig(j) = 0;
            end
        else
            Ig(j) = 0;
        end
    end
    % g >= 0.47: equilibrium is stable
    % g <= 0.46: equilibrium is unstable
    % %% graphics
    % fsz = 20;
    % 
I_eq_stable = find(Ig == 1 & max_abs_evals_eq <= 1 & max_abs_evals_eq >= 0);
I_eq_unstable = find(Ig == 1 & max_abs_evals_eq > 1 & max_abs_evals_cyc <= 1);
if ~isempty(I_eq_unstable)
    Iend = I_eq_unstable(end);
    I = Iend+1;
    I_eq_unstable = [I_eq_unstable;I];
    u1min(I) = u1_eq(I);
    u2min(I) = u2_eq(I);
    u1max(I) = u1_eq(I);
    u2max(I) = u2_eq(I);
    pmin(I) = prey_eq(I);
    pmax(I) = prey_eq(I);
end

    color_prey = [0,0,1]; % blue
    color_u1 =  [1,0,0]; % red
    color_u2 = [1,0.8,0.2]; % dark yellow


    figure(4);clf; hold on;

    Ilegend = [];
    hand = zeros(1,9);
    if ~isempty(I_eq_stable)
        hand(1) = plot(gvals(I_eq_stable),prey_eq(I_eq_stable),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(5) = plot(gvals(I_eq_stable),u1_eq(I_eq_stable),'Linewidth',lwt,'color',color_u1,'DisplayName','Juvenile predator');
        hand(9) = plot(gvals(I_eq_stable),u2_eq(I_eq_stable),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
        Ilegend = [1,5,9];
    end
    if ~isempty(I_eq_unstable )
        hand(2) = plot(gvals(I_eq_unstable),prey_eq(I_eq_unstable),'--','Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(3) = plot(gvals(I_eq_unstable),pmin(I_eq_unstable),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
        hand(4) = plot(gvals(I_eq_unstable),pmax(I_eq_unstable),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
    
        hand(6) = plot(gvals(I_eq_unstable),u1_eq(I_eq_unstable),'--','Linewidth',lwt,'color',color_u1,'DisplayName','Juvenile predator');
        hand(7) = plot(gvals(I_eq_unstable),u1min(I_eq_unstable),'Linewidth',lwt,'color',color_u1,'DisplayName','Juvenile predator');
        hand(8) = plot(gvals(I_eq_unstable),u1max(I_eq_unstable),'Linewidth',lwt,'color',color_u1,'DisplayName','Juvenile predator');
    
        hand(10) = plot(gvals(I_eq_unstable),u2_eq(I_eq_unstable),'--','Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
        hand(11) = plot(gvals(I_eq_unstable),u2min(I_eq_unstable),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
        hand(12) = plot(gvals(I_eq_unstable),u2max(I_eq_unstable),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
        if isempty(Ilegend)
            Ilegend = [3,7,11];
        end
    end
    if ~isempty(I_eq_stable) && ~isempty(I_eq_unstable )
        set(gca,'Fontsize',fsz);
        xlabel(g_text,'FontSize',fsz);
        ylabel('Population Size','Fontsize',fsz);
        legend(hand(Ilegend))
    else
        if ~isempty(I_eq_stable) || ~isempty(I_eq_unstable )
            set(gca,'Fontsize',fsz);
            xlabel(g_text,'FontSize',fsz);
            ylabel('Population Size','Fontsize',fsz);
            legend();
        end
    end
    grid on;
    xlim([min(gvals(Ig==1)),max(gvals(Ig==1))]);
    title(strcat("\tau* = ",sprintf("%.2f",tau)),'FontSize',fsz);
end

end

