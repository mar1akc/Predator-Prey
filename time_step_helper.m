function [prey_new,unew,u1,u2] = time_step_helper(prey,u,u1,u2,g,h_age,dt,age_vector,...
    AgeNode,r,a,k,bb,s,rho,z,Mm,I1,I2,mu_base,birth_rate_base)

h05 = 0.5*h_age;

u1(1) = h05*(u(1)+u(AgeNode))+h_age*sum(u(2:AgeNode-1));
u2(1) = h05*(u(AgeNode)+u(end))+h_age*sum(u(AgeNode+1:end-1));

prey_new = dt*prey*(r-a*prey+s*u1-bb*u2)+prey;

unew = u;

death_rate =  mu_base(age_vector(:)) + Mm*exp(-rho*prey) + g*prey*I1(age_vector(:));
birth_rate =  birth_rate_base(age_vector(:))*(1-exp(-z*prey)) + k*prey*I2(age_vector(:));

unew(1) = h05*(u(1)*birth_rate(1) + u(end)*birth_rate(end)) + ...
    h_age*sum(u(2:end-1).*birth_rate(2:end-1));
unew(2:end) = u(1:end-1).*(1-dt*death_rate(2:end));

u1 = h05*(unew(1)+unew(AgeNode))+h_age*sum(unew(2:AgeNode-1));
u2 = h05*(unew(AgeNode)+u(end))+h_age*sum(unew(AgeNode+1:end-1));
end