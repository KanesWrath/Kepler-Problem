%% The True Error of the Newton - Ralphson Method for Kepler's Equation.

function [count, epsilon_t] = Keplers_Equation_Error(E_star, E_0, e, M)

E_n = E_0 - Kepler_Equation(M, E_0, e)/Kepler_Equation_Derivative(e, E_0);
i = 1;
epsilon_a = abs(E_n - E_0)/abs(E_n);
epsilon_t(i) = abs(E_n - E_star)/abs(E_star);
count(i) = i;

%while (epsilon_a > 10^-4)
for i = 1:5
    
    E_0 = E_n;
    E_n = E_0 - Kepler_Equation(M, E_0, e)/Kepler_Equation_Derivative(e, E_0);
    i = i + 1;
    epsilon_a = abs(E_n - E_0)/abs(E_n);
    epsilon_t(i) = abs(E_n - E_star)/abs(E_star);
    count(i) = i;
 
end

end
    
