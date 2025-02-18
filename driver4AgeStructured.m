function driver4AgeStructured()
close all
g = 1;
for tau = 0.1:0.1:3.0
    fprintf("tau = %d\n",tau);
    u = RunModel(tau,g);
    umax = max(u);
    while umax < 1e-10 
        g = g - 0.01; 
        u = RunModel(tau,g);
        umax = max(u); 
        fprintf("g = %d, umax = %d\n",g,umax); 
    end
    Find_Equilib_driver(g);
    Find_periodic_solution_driver(g);
end
end