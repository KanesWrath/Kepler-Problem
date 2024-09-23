import numpy as np
import matplotlib.pyplot as plt
import time

# Start Computation Timer.

start_time = time.time()

# TLE Data.

e = 0.0002714                       # Eccentricity.
M = 2.5                             # Mean Anomaly, radians.
n = 1.0027118 * (2 * np.pi / 86400) # Mean Motion, rad/s.

# Constants.

mu = 3.986004415e5                  # Gravitational Parameter for Earth, km^3/s^2.

# Initial Guess for Eccentric Anomaly, E_0 based on the Mean Anomaly M.

if -np.pi < M < 0 or M > np.pi:

    E_0 = M - e

else:

    E_0 = M - e

# Find Semimajor Axis.

a = (mu / (n ** 2)) ** (1 / 3)      # Semimajor Axis, km.

# Define Kepler's Equation.

def kepler_equation(M, E_n, e):

    KP = M - E_n + e * np.sin(E_n)

    return KP

# Derivative of the Kepler Equation.

def kepler_equation_derivative(E_n, e):

    KP_prime = e * np.cos(E_n) - 1

    return KP_prime

# Use Newton-Raphson Method to Solve for E_n.

def newton_raphson(M, e, E_0, tol = 1e-10, max_iter = 50):

    E_n = E_0

    for count in range(max_iter):

        KP = kepler_equation(M, E_n, e)
        KP_prime = kepler_equation_derivative(e, E_n)
        delta = KP / KP_prime
        E_n = E_n - delta

        if abs(delta) < tol:

            break

    return E_n, count

E_n, count = newton_raphson(M, e, E_0)

# Find True Anomaly.

nu = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E_n / 2)) * (180 / np.pi)

# Ensure nu is Positive.

if nu < 0:

    nu += 360

# # Calculate the Error during Iterations (for Plot).

# def kepler_equation_error(E_star, E_0, e, M):

#     epsilon_t = []

#     for _ in range(5):

#         E_star, _ = newton_raphson(M, e, E_0)
#         epsilon_t.append(abs(kepler_equation(M, E_star, e)))

#     return list(range(1, 6)), epsilon_t

# Function to Calculate the True Error of the Newton-Raphson Method for Kepler's Equation.

def keplers_equation_error(E_star, E_0, e, M):

    # First Newton-Raphson Step.

    E_n = E_0 - kepler_equation(M, E_0, e) / kepler_equation_derivative(e, E_0)
    i = 0
    epsilon_a = abs(E_n - E_0) / abs(E_n)
    epsilon_t = [abs(E_n - E_star) / abs(E_star)]
    count = [i + 1]

    # Perform 5 iterations.

    for i in range(5):

        E_0 = E_n
        E_n = E_0 - kepler_equation(M, E_0, e) / kepler_equation_derivative(e, E_0)
        epsilon_a = abs(E_n - E_0) / abs(E_n)
        epsilon_t.append(abs(E_n - E_star) / abs(E_star))
        count.append(i + 2)

    return count, epsilon_t

count, epsilon_t = keplers_equation_error(E_n, E_0, e, M)

# Plot Newton-Raphson Error.

plt.figure(1)
plt.plot(count, epsilon_t)
# plt.axis([0, 5, 0, 6e-12])
plt.grid(True)
plt.title('Newton-Raphson Error')
plt.xlabel('Iteration')
plt.ylabel('Error')

# Mean Anomaly vs. True Anomaly Plot.

E = np.linspace(-np.pi, np.pi, 60)
M_values = E + e * np.sin(E)
nu_values = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2)) * (180 / np.pi)

E_deg = E * (180 / np.pi)
M_deg = M_values * (180 / np.pi)

plt.figure(2)
plt.plot(E_deg, M_deg)
plt.grid(True)
plt.title('XX, Mean Anomaly vs. True Anomaly')
plt.legend(['e = '])
plt.xlabel('True Anomaly, ν (degrees)')
plt.ylabel('Mean Anomaly, M (degrees)')
plt.xticks(np.arange(-180, 181, 30), labels=[f'{x}°' for x in np.arange(-180, 181, 30)])
plt.yticks(np.arange(-180, 181, 30), labels=[f'{x}°' for x in np.arange(-180, 181, 30)])

# Mean Anomaly vs. Eccentric Anomaly Plot.

plt.figure(3)
plt.plot(E_deg, M_deg)
plt.grid(True)
plt.title('XX, Mean Anomaly vs. Eccentric Anomaly')
plt.legend(['e = '])
plt.xlabel('Eccentric Anomaly, E (degrees)')
plt.ylabel('Mean Anomaly, M (degrees)')

# Display the Plots.

plt.show()

# End Computation Time.

end_time = time.time()
elapsed_time = end_time - start_time
print(f'Calculation Time = {elapsed_time:.6f} [s]')