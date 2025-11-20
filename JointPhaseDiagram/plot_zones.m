function plot_zones()
data_ODE = load("ZonesODEnu100.mat");
data_DDE = load("ZonesDDEnu100.mat");
data_100 = load("ZonesPDEnu100.mat");
data_1 = load("ZonesPDEnu1.mat");
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
tau1 = data_1.tau_vals;
%
C100 = data_100.Cbdry;
E100 = data_100.Ebdry;
tau100 = data_100.tau_vals;
%
ode_red = [1,0,0];
dde_red = [0.5,0,0];
red1 = [0.7,0,0];
red100 = [1,0.5,0];
%
ode_blue = [0,0,1];
dde_blue = [0,1,1];
blue1 = [0,0,0.5];
blue10 = [0.5,0,1];
blue100 = [0,0.4,0.4];
blueS = [0,0,0.3];
bluesm = [0,76/255,153/255];
%
verts = [[tau_ode(j1:end),Code(j1:end)];[2,0]];
patch(verts(:,1),verts(:,2),'red');
alpha(0.1); 

verts = [[tau_dde(18:end),Cdde(18:end)];[2,0]];
patch(verts(:,1),verts(:,2),'red');
alpha(0.1); 

verts = [[tau1,C1];[2,0]];
patch(verts(:,1),verts(:,2),'red');
alpha(0.1); 

verts = [[tau100,C100];[2,0]];
patch(verts(:,1),verts(:,2),'red');
alpha(0.1); 

verts = [[tau_ode(18:end),Eode(18:end)];[2,1]];
patch(verts(:,1),verts(:,2),'blue');
alpha(0.1); 

verts = [[tau_dde(18:end),Edde(18:end)];[2,1]];
patch(verts(:,1),verts(:,2),'blue');
alpha(0.1); 

verts = [[tau1,E1];[2,1]];
patch(verts(:,1),verts(:,2),'blue');
alpha(0.1); 

verts = [[tau100(18:end),E100(18:end)];[2,1]];
patch(verts(:,1),verts(:,2),'blue');
alpha(0.1); 

lwt = 3;
plot(tau_ode(j1:end),Code(j1:end),'Linewidth',lwt,'color',ode_red,'displayname','ODE: periodic coexistence attractor')
plot(tau_ode(18:end),Eode(18:end),'Linewidth',lwt,'color',ode_blue,'displayname','ODE: predator-free attractor')
plot(tau_dde(18:end),Cdde(18:end),'Linewidth',lwt,'color',dde_red,'displayname','DDE: periodic coexistence attractor')
plot(tau_dde(18:end),Edde(18:end),'Linewidth',lwt,'color',dde_blue,'displayname','DDE: predator-free attractor')
plot(tau1,C1,'Linewidth',lwt,'color',red1,'displayname','Periodic coexistence attractor, $\nu = 1$')
plot(tau1,E1,'Linewidth',lwt,'color',blue1,'displayname','predator-free attractor, $\nu = 1$')
plot(tau100,C100,'Linewidth',lwt,'color',red100,'displayname','Periodic coexistence attractor, $\nu = 100$')
plot(tau100(18:end),E100(18:end),'Linewidth',lwt,'color',blue100,'displayname','predator-free attractor, $\nu = 100$')
set(gca,'Fontsize',fsz)
xlabel("Maturation age, \tau*",'Fontsize',fsz);
ylabel(g_text,'Fontsize',fsz);
taumax = min([tau_ode(end),tau1(end),tau100(end)]);
axis([0.1,taumax,0,1])
end










