import matplotlib.pyplot as plt
import numpy as np

#Differentiation
def funct_diff(x_n: np.array(float), u: np.float64):
    dx_n = np.zeros(degree_of_equation)
    dx_n[0] = x_n[1]
    dx_n[1] = -x_n[0] + (1 - x_n[0]**2 * x_n[1] + u)
    return dx_n

#Methods
#Euler method
def step_in_Euler(x, h, u):  #coordinates, delta, control_operation
    x_next = x + funct_diff(x, u) * h 
    return x_next

#Improved Euler method
def step_in_ImpEuler(x, h, u): 
    k1 = funct_diff(x, u)
    x_next = x + k1 * h
    k2 = funct_diff(x_next, u)
    x_next = x + (h * (k2 + k1) / 2)
    return x_next

#Runge-Kutta method
def step_in_RungKutt(x, h, u):
    k1 = funct_diff(x, u)
    k2 = funct_diff(x + h * (k1/2), u + h/2)
    k3 = funct_diff(x + h * (k2/2), u + h/2)
    k4 = funct_diff(x + h * k3, u + h)
    x_next = x + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return x_next

#Get parameters using methods above
#Euler parameteres
def get_Euler_params(x_n):
    Euler = [np.zeros(number_of_steps), np.zeros(number_of_steps)] #Coordinates for return
    for one_step in range(number_of_steps):
        x_n_new = step_in_Euler(x_n, h, u)
        x_n = x_n_new
        Euler[0][one_step] = x_n[0]
        Euler[1][one_step] = x_n[1]
    return Euler
    
#Improved Euler parameteres
def get_ImpEuler_params(x_n):
    ImpEuler = [np.zeros(number_of_steps), np.zeros(number_of_steps)] #Coordinates for return
    for one_step in range(number_of_steps):   
        x_n_new = step_in_ImpEuler(x_n, h, u)
        x_n = x_n_new
        ImpEuler[0][one_step] = x_n[0]
        ImpEuler[1][one_step] = x_n[1]
    return ImpEuler

#Runge-Kutta parameteres
def get_RungKutt_params(x_n):
    RungeKutt = [np.zeros(number_of_steps), np.zeros(number_of_steps)] #Coordinates for return
    for one_step in range(number_of_steps):   
        x_n_new = step_in_RungKutt(x_n, h, u)
        x_n = x_n_new
        RungeKutt[0][one_step] = x_n[0]
        RungeKutt[1][one_step] = x_n[1]
    return RungeKutt

# Plot
def plot_func(x, y, z, t):
    plt.subplot2grid ((3, 3),(2, 0), colspan= 2)
    plt.plot(x[0], x[1], label = 'Euler')
    plt.plot(y[0], y[1], label = 'Improved Euler')
    plt.plot(z[0], z[1], label = 'Runge Kutt')
    plt.title('x1, x2')
    plt.xlabel('x1')
    plt.ylabel('x2')
    plt.legend()
    plt.grid() 

    plt.subplot2grid ((3, 3),(0, 0), colspan= 2)
    plt.plot(t, x[0], label = 'Euler')
    plt.plot(t, y[0], label = 'Improved Euler')
    plt.plot(t, z[0], label = 'Runge Kutt')
    plt.title('x1, time')
    plt.ylabel('x1')
    plt.legend()
    plt.grid() 

    plt.subplot2grid ((3, 3),(1, 0), colspan= 2)
    plt.plot(t, x[1], label = 'Euler')
    plt.plot(t, y[1], label = 'Improved Euler')
    plt.plot(t, z[1], label = 'Runge Kutt')
    plt.title('x2, time')
    plt.ylabel('x2')
    plt.legend()
    plt.grid()

    plt.show()



degree_of_equation = 2
time_zero,  time_end,  h = 0.0,  10.0,  0.1  #Time scope(sec) and delta(h = 0.1 or 0.01 or 0.001)

number_of_steps = int((time_end - time_zero)/h)
time_axis = np.arange(time_zero, time_end, h) 

x_n = np.zeros(degree_of_equation) #Zero coordinates
u = 5.12345 #Control operation

Euler = get_Euler_params(x_n) #Gets x1 and x2
ImpEuler = get_ImpEuler_params(x_n)
RungKutt = get_RungKutt_params(x_n)

plot_func(Euler, ImpEuler, RungKutt, time_axis)
