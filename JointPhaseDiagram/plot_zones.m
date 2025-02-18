function plot_zones()
data_ODE = load("ZonesODE_par_ifun100.mat");
data_DDE = load("ZonesDDE_par_ifun100.mat");
data_S = load("ZonesPDE_par_ifun100.mat");
data_100 = load("Zones_par_ifun100.mat");
data_10 = load("Zones_par_ifun10.mat");
data_1 = load("Zones_par_ifun1.mat");

g_text = "Consump. Rate of Juv. Pred., g";

%
figure()
hold on
grid
fsz = 18;
%
Code = data_ODE.Cbdry;
Eode = data_ODE.Ebdry;
j1 = data_ODE.j1;
j2 = data_ODE.j2;
tau_ode = data_ODE.tauvals;
%
Cdde = data_DDE.Cbdry;
Edde = data_DDE.Ebdry;
j1dde = data_DDE.j1;
j2dde = data_DDE.j2;
tau_dde = data_DDE.tauvals;
%
C1 = data_1.Cbdry;
E1 = data_1.Ebdry;
tau1 = data_1.tauvals;
%
C10 = data_10.Cbdry;
E10 = data_10.Ebdry;
tau10 = data_10.tauvals;
%
C100 = data_100.Cbdry;
E100 = data_100.Ebdry;
tau100 = data_100.tauvals;
%
C100S = data_S.Cbdry;
E100S = data_S.Ebdry;
tau100S = data_S.tauvals;

%
ode_red = [1,0,0];
dde_red = [0.5,0,0];
red1 = [0.7,0,0];
red10 = [1,0,0.5];
red100 = [1,0.5,0];
redS = [0.3,0,0];
%
ode_blue = [0,0,1];
dde_blue = [0,1,1];
blue1 = [0,0,0.5];
blue10 = [0.5,0,1];
blue100 = [0,0.4,0.4];
blueS = [0,0,0.3];
%
plot(tau_ode(j1:end),Code(j1:end),'Linewidth',2,'color',ode_red,'displayname','ODE: periodic coexistence attractor')
plot(tau_ode(j1:end),Eode(j1:end),'Linewidth',2,'color',ode_blue,'displayname','ODE: predator-free attractor')
plot(tau_dde,Cdde,'Linewidth',2,'color',dde_red,'displayname','DDE: periodic coexistence attractor')
plot(tau_dde,Edde,'Linewidth',2,'color',dde_blue,'displayname','DDE: predator-free attractor')
plot(tau1,C1,'Linewidth',2,'color',red1,'displayname','Periodic coexistence attractor, $\nu = 1$')
plot(tau1,E1,'Linewidth',2,'color',blue1,'displayname','predator-free attractor, $\nu = 1$')
plot(tau10,C10,'Linewidth',2,'color',red10,'displayname','Periodic coexistence attractor, $\nu = 10$')
plot(tau10,E10,'Linewidth',2,'color',blue10,'displayname','predator-free attractor, $\nu = 10$')
plot(tau100,C100,'Linewidth',2,'color',red100,'displayname','Periodic coexistence attractor, $\nu = 100$')
plot(tau100,E100,'Linewidth',2,'color',blue100,'displayname','predator-free attractor, $\nu = 100$')
plot(tau100S,C100S,'--','Linewidth',2,'color',redS,'displayname','Periodic coexistence attractor, $\nu = 100$')
plot(tau100S,E100S,'--','Linewidth',2,'color',blueS,'displayname','predator-free attractor, $\nu = 100$')
% lgd = legend('Location','northeast','Interpreter','latex')
lgd.FontSize = 12;
set(gca,'Fontsize',fsz)
xlabel("Maturation age, \tau*",'Fontsize',fsz);%,'Interpreter','latex')
ylabel(g_text,'Fontsize',fsz);%,'Interpreter','latex')
taumax = min([tau_ode(end),tau1(end),tau10(end),tau100(end)]);
axis([0.0,taumax,0.01,0.99])
end










