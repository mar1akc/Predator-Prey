function[prey,u,u1,u2] = PreyPredTimeStep(prey,u,g,par)

[L,tau,h_age,dt,age_vector,AgeNode,par_ifun,r,a,k,bb,s,rho,z,...
    Mm,mu_b_par,birth_par,birth_exp_par,death_exp_par] ...
    = unroll_parameters(par);


[I1,I2,mu_base,birth_rate_base] = assign_parameter_functions(par_ifun,...
    tau,L,mu_b_par,death_exp_par,birth_par,birth_exp_par);

h05 = 0.5*h_age;

u1(1) = h05*(u(1)+u(AgeNode))+h_age*sum(u(2:AgeNode-1));
u2(1) = h05*(u(AgeNode)+u(end))+h_age*sum(u(AgeNode+1:end-1));

[prey,u,u1,u2] = time_step_helper(prey,u,u1,u2,g,h_age,dt,age_vector,...
    AgeNode,r,a,k,bb,s,rho,z,Mm,I1,I2,mu_base,birth_rate_base);
end

