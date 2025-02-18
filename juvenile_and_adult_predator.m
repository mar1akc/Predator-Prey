%%
function  [u1,u2] = juvenile_and_adult_predator(u,par)
tau = par(2); % maturation age
h_age = par(3); % step in age
AgeNode = floor(tau/h_age+1);
h05 = 0.5*h_age;
% juvenile predator at time 0
u1 = h05*(u(1)+u(AgeNode))+h_age*sum(u(2:AgeNode-1));
% adult predatir at time 0
u2 = h05*(u(AgeNode)+u(end))+h_age*sum(u(AgeNode+1:end-1));
end
