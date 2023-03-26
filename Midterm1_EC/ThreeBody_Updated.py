from matplotlib import pyplot as plt
import numpy as np

# Function to calculate state dot

def Three_B_Fun(state, G, M, N):

    X = state[0:N]
    Y = state[N:2*N]
    U = state[2*N:3*N]
    V = state[3*N:4*N]

    # Change in vel is acceleration
    Udot = np.zeros(N)
    Vdot = np.zeros(N)

    for i in range(N):
        for j in range(N):
            if not (i == j):
                r_ij = np.sqrt(np.power(X[i]-X[j], 2) + np.power(Y[i]-Y[j], 2))

                udot_j = -G * M[j] * (X[i]-X[j]) / np.power(r_ij, 3)
                vdot_j = -G * M[j] * (Y[i]-Y[j]) / np.power(r_ij, 3)

                Udot[i] = Udot[i] + udot_j
                Vdot[i] = Vdot[i] + vdot_j

    state_dot = np.append(U, [V, Udot, Vdot])
    return state_dot


# constants
Msun = 2e30 # kg
G = 6.67e-11
M = np.divide(np.array([1, 1, 1]), G)
AU = 1.5*np.power(10, 11)

# initializations
year_2_sec = 365*24*3600
N = 3
steps = 1001
end_time = 6.28
x = np.array([-0.97000436, 0., 0.97000436])  # x1 = -x3, x2 = 0
y = np.array([0.24208753, 0., -0.24208753])  # y1 = -y3, y2 = 0
vx = np.array([0.4662036850, -0.933240737, 0.4662036850])  # v1x = v3x
vy = np.array([0.4323657300, -0.86473146, 0.4323657300])  # v1y = v3y
X0 = np.append(x, [y, vx, vy])

tspan = np.linspace(0, end_time, steps)
del_t = tspan[1] - tspan[0]

Output = [X0]

for i in range(steps):

    X_current = Output[i]

    X_current_dot = Three_B_Fun(X_current, G, M, N)

    new_X = X_current + np.multiply(X_current_dot, del_t)

    Output.append(new_X)

plt.axes()

for i in range(N):
    plt.scatter(Output[:][i], Output[:][i+3])
plt.show()

print(Output)