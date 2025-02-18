function [u1min,u1max,u2min,u2max,pmin,pmax,period] = ...
    find_loop_min_max(prey_equilib,prey,u,g,par)
% make one step to escape from the Poincare section
[u1,u2] = juvenile_and_adult_predator(u,par);
u1min = u1;
u1max = u1;
u2min = u2;
u2max = u2;
pmin = prey;
pmax = prey;
[prey,u,u1,u2] = PreyPredTimeStep(prey,u,g,par); 
u1min = min(u1,u1min);
u1max = max(u1,u1max);
u2min = min(u2,u2min);
u2max = max(u2,u2max);
pmin = min(prey,pmin);
pmax = max(prey,pmax);

% We needed to detect the crossing of 
% the Poincare section in the same direction
dt = par(4);
t = 0;
prey_sign = sign(prey - prey_equilib);
while prey_sign == sign(prey - prey_equilib)
    prey_prev = prey;
    u_prev = u;
   [prey,u,u1,u2] = PreyPredTimeStep(prey,u,g,par);
   t = t + dt;
   u1min = min(u1,u1min);
    u1max = max(u1,u1max);
    u2min = min(u2,u2min);
    u2max = max(u2,u2max);
    pmin = min(prey,pmin);
    pmax = max(prey,pmax);

end
while prey_sign ~= sign(prey - prey_equilib)
    prey_prev = prey;
    u_prev = u;
   [prey,u,~,~] = PreyPredTimeStep(prey,u,g,par);
   t = t + dt;
    u1min = min(u1,u1min);
    u1max = max(u1,u1max);
    u2min = min(u2,u2min);
    u2max = max(u2,u2max);
    pmin = min(prey,pmin);
    pmax = max(prey,pmax);
end   
% find the point of intersection of the trajectory with the hyperplane
aux = (prey_equilib - prey_prev)/(prey - prey_prev);
uloop = u_prev + (u - u_prev)*aux;
period = t - dt + dt*aux;
end
