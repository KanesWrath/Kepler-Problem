%% Derivative of Kepler Equation.

function [KP_prime] = Kepler_Equation_Derivative(e, E_n)

KP_prime = e*cos(E_n) - 1;

end

