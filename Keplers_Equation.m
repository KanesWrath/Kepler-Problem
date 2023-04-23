%% Newton - Ralphson Method of Kepler's Equation.

function [E_n] = Keplers_Equation(E_0, e, M)

E_n = E_0 - Kepler_Equation(M, E_0, e)/Kepler_Equation_Derivative(e, E_0);

i = 1;
count(i) = i;

for k = 1:30
    
    E_0 = E_n;
    E_n = E_0 - Kepler_Equation(M, E_n, e)/...
        Kepler_Equation_Derivative(e, E_n);
    i = i + 1;
    count(i) = i;
    
end

end
