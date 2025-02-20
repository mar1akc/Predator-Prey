function driver4DDEmodel()

for tau = 0.1:0.1:3.0
    fprintf("tau = %d\n",tau);
    EquilibriumSolutionsDDE(tau);
end
end