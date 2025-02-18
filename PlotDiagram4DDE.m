function PlotDiagram4DDE()
% flag = 0: image plot for dynamical behavior zones
% flag = 1: shaded zones of dynamical behavior
% flag = 2: bifurcation diagram for a particular tau
fsz = 18; % fontsize for graphics
lwt = 3; % linewidth
g_text = "Consumption Rate of Juvenile Predator, g";

tau_bif_diag = 1.2; % tau for the bifurcation diagram if flag = 2

fsz = 24; % fontsize for graphics
lwt = 3; % linewidth

par_ifun = 100;

flag = 2; 

% fname = sprintf("DataODE_DDE/DDE_BifurDiag_par_ifun%.1f_tau%.2f.mat",par_ifun,tau);

% save(fname,"gvals","prey_equilib","y1_equilib","y2_equilib",...
%     "prey_loop_min","prey_loop_max","y1_loop_min","y1_loop_max",...
%     "y2_loop_min","y2_loop_max","period","stab");

tauvals = (0.1:0.01:2.28)';
gvals = linspace(0,1,101)';

N_gvals = length(gvals);
N_tauvals = length(tauvals);
%% plot the diagram on tau-g-plane: image
if flag == 0
    
    Bfun = zeros(N_gvals,N_tauvals);
    Itau = ones(N_tauvals,1);
    for jtau = 1 : N_tauvals
        tau = tauvals(jtau);
        fname = sprintf("DataODE_DDE/DDE_BifurDiag_par_ifun%.1f_tau%.2f.mat",par_ifun,tau);
        if isfile(fname)
            data = load(fname);
            g = data.gvals;
            jg = length(g);
            stab = data.stab;
            Istab = find(stab == 1);
            Inonstab = find(stab == 0); 
            Istab = jg-length(Istab)+1:jg;
            Inonstab = 1:length(Inonstab);
            Bfun(Istab,jtau) = 1;
            Bfun(Inonstab,jtau) = 2;
        else
            Itau = 0;
        end
    end
    Iremove = find(Itau == 0);
    Bfun(:,Iremove) = [];
    tauvals(Iremove) = [];
    
    figure(1); clf; hold on; 
    imagesc(tauvals,gvals,Bfun);
    set(gca,'Fontsize',fsz);
    xlabel("Maturation Age, \tau^*",'FontSize',fsz);
    ylabel("Juv. Pred. Damage Coeff., g",'Fontsize',fsz);
end
%%
%% plot the diagram on tau-g-plane: curves and shaded zones
if flag == 1
    gstep = gvals(2) - gvals(1);
    Ebdry = ones(N_tauvals,1);
    Cbdry = zeros(N_tauvals,1);
    Xbdry = zeros(N_tauvals,1);
    for jtau = 1 : N_tauvals
        tau = tauvals(jtau);
        fname = sprintf("DataODE_DDE/DDE_BifurDiag_par_ifun%.1f_tau%.2f.mat",par_ifun,tau);
        if isfile(fname)
            data = load(fname);
            g = data.gvals;
            stab = data.stab;
            Istab = find(stab == 1);
            Inonstab = find(stab == 0);
            if max(Istab) == N_gvals
                Ebdry(jtau) = max(gvals);
            else
                Xbdry(jtau) = max(g) + 0.5*gstep;
                Ebdry(jtau) = max(g(Istab)) + 0.5*gstep;
                if ~isempty(Inonstab)
                    Cbdry(jtau) = max(g(Inonstab)) + 0.5*gstep;
                else
                end
            end
        end
    end
    
    figure(1); clf; hold on; 
    % Shade the zones
    % First, triangulate zones
    % Zone 1, below Cbdry : oscillations   
    j1 = find(Cbdry == 0, 1, 'last' );
    vv = [tauvals((j1):end),Cbdry((j1):end)];
    v1 = [tauvals(end),0];
    verts = [vv;v1];
    patch(verts(:,1),verts(:,2),'red');
    alpha(0.4);
    
    % Zone 2, between Cbdry and Ebdry: coexistence attractor
    vv2 = [tauvals,Ebdry];
    vv3 = [tauvals(j1:-1:1),zeros(size(tauvals(j1:-1:1)))];
    verts = [vv2;vv(end:-1:1,:);vv3];
    patch(verts(:,1),verts(:,2),'green');
    alpha(0.4);

    % Zone 3, extinction
    j2 = find(Ebdry == gvals(end),1,'last');
    vv4 = [tauvals(j2:end),Ebdry(j2:end)];
    v5 = [max(tauvals),max(gvals)];
    vv5 = [tauvals(end-1:-1:(j2+1)),max(gvals)*ones(size(tauvals(end-1:-1:(j2+1))))];
    verts = [vv4;v5;vv5];
    patch(verts(:,1),verts(:,2),'blue');
    alpha(0.2);


    plot(tauvals(j2:end),Ebdry(j2:end),'linewidth',lwt,'color',[0,0.5,0]);
    plot(tauvals(j1:end),Cbdry(j1:end),'linewidth',lwt,'color',[0.5,0,0]);

    fbdry_name = sprintf("ZonesDDE_par_ifun%d.mat",par_ifun);
    save(fbdry_name,"tauvals","Ebdry","Cbdry","j2","j1");

    
    set(gca,'Fontsize',fsz);
    xlabel("Maturation Age, \tau^*",'FontSize',fsz);
    ylabel(g_text,'Fontsize',fsz);
    saveas(gcf,"Zones","epsc")

end
%% plot bifurcation diagram
if flag == 2
    tau = tau_bif_diag;
    fname = sprintf("DataODE_DDE/DDE_BifurDiag_par_ifun%.1f_tau%.2f.mat",par_ifun,tau);
    data = load(fname);
    gvals = data.gvals;

    stab = data.stab;

    % save(fname,"gvals","prey_equilib","y1_equilib","y2_equilib",...
    %     "prey_loop_min","prey_loop_max","y1_loop_min","y1_loop_max",...
    %     "y2_loop_min","y2_loop_max","period","stab");
    prey_equilib = data.prey_equilib;
    prey_loop_min = data.prey_loop_min;
    prey_loop_max = data.prey_loop_max;
    y1_equilib = data.y1_equilib;
    y1_loop_min = data.y1_loop_min;
    y1_loop_max = data.y1_loop_max;
    y2_equilib = data.y2_equilib;
    y2_loop_min = data.y2_loop_min;
    y2_loop_max = data.y2_loop_max;

    Istab = find(stab == 1);
    Inonstab = find(stab == 0); 


    

color_prey = [0,0,1]; % blue
color_u1 =  [1,0,0]; % red
color_u2 = [1,0.8,0.2]; % dark yellow
if isempty(Inonstab) == 0 
    figure(1);clf; hold on;
    hand = zeros(1,9);
    hand(1) = plot(gvals(Istab),prey_equilib(Istab),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
    hand(2) = plot(gvals(Inonstab),prey_equilib(Inonstab),'--','Linewidth',lwt,'color',color_prey,'DisplayName','prey');
    hand(3) = plot(gvals,prey_loop_min,'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
    hand(4) = plot(gvals,prey_loop_max,'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');

    hand(5) = plot(gvals(Istab),y1_equilib(Istab),'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
    hand(6) = plot(gvals(Inonstab),y1_equilib(Inonstab),'--','Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
    hand(7) = plot(gvals,y1_loop_min,'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
    hand(8) = plot(gvals,y1_loop_max,'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');

    hand(9) = plot(gvals(Istab),y2_equilib(Istab),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    hand(10) = plot(gvals(Inonstab),y2_equilib(Inonstab),'--','Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    hand(11) = plot(gvals,y2_loop_min,'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    hand(12) = plot(gvals,y2_loop_max,'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    set(gca,'Fontsize',fsz);
    xlabel(g_text,'FontSize',fsz);%,'Interpreter','latex');
    ylabel('Population Size','Fontsize',fsz);%,'Interpreter','latex');
    legend(hand([1,5,9]));%,'Interpreter','latex')
    title(strcat('\tau* = ',sprintf(" %.1f",tau)));%,'Interpreter','latex')
    grid on

else
    figure(1);clf; hold on;
    plot(gvals(Istab),prey_equilib(Istab),'Linewidth',lwt,'color',color_prey,'DisplayName','Prey');
    plot(gvals(Istab),y1_equilib(Istab),'Linewidth',lwt,'color',color_u1,'DisplayName','Juv. predator');
    plot(gvals(Istab),y2_equilib(Istab),'Linewidth',lwt,'color',color_u2,'DisplayName','Adult predator');
    set(gca,'Fontsize',fsz)
    set(gcf,'color','white')
    legend('Location','northeast','Interpreter','latex') 
    xlabel(g_text,'FontSize',fsz);%,'Interpreter','latex');
    ylabel('Population Size','Fontsize',fsz);%,'Interpreter','latex');
    title(strcat('\tau* = ',sprintf(" %.1f",tau)));%,'Interpreter','latex')
    grid on
end

end
end