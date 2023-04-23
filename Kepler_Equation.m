%% Kepler Equation.

function [KP] = Kepler_Equation(M, E_n, e)

KP = M - E_n + e*sin(E_n);

end
