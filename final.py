#! /usr/bin/env python

"""
File: final.py
Copyright (c) 2016 Michael Seaman
License: MIT

Description: Implementing the ODE Approximator Runge-Kutta 4 to approximate the three
coupled ODEs forming the Rossler Attractor:

x' = -y + -z
y' = x + .2y
z' = .2 + z(x - c)

Where we can fine tune the value c
"""

from unittest import TestCase
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class Rossler():
    """
    Solves and plots the rossler equations for different values of c and 
    initial conditions x0, y0, and z0

    x' = -y + -z
    y' = x + .2y
    z' = .2 + z(x - c)
    """
    def __init__(self, c, dt=0.001, T0 = 250, T=500, x0 = 0, y0 = 0, z0 = 0):
        """
        Takes float64 for c, the non-linear fine tuning, a float64 for dt which
        is the stepsize the run() call uses, T0 is an int from which plots against
        t begin from, and T is the amount of total time to pass.
        x0, y0, and z0 can also be varied.
        """
        self.c = c
        self.dt = float(dt)
        self.T0 = T0
        self.T = T
        self.t = np.arange(0, T+dt, dt, dtype = 'float64')
        self.x = np.zeros(len(self.t), dtype = 'float64') + x0
        self.y = np.zeros(len(self.t), dtype = 'float64') + y0
        self.z = np.zeros(len(self.t), dtype = 'float64') + z0

    def dxdt(self, x, y, z, t):
        return -1 * y + -1 * z

    def dydt(self, x, y, z, t):
        return x + .2 * y

    def dzdt(self, x, y, z, t):
        return .2 + z * (x - self.c)

    def run(self):
        """
        Runs a total of T / dt steps of the Runge Kutta 4 approximation 
        using the Rossler equations. Does not return output; instead, simply
        updates the contents of arrays self.x, self.y, and self.z
        """   
        dxdt = self.dxdt
        dydt = self.dydt
        dzdt = self.dzdt
        dt = self.dt
        x = self.x
        y = self.y
        z = self.z
        t = self.t

        for i in np.arange(0, len(t) - 1):
            k1_x = dt * dxdt(x[i], y[i], z[i], t[i])
            k1_y = dt * dydt(x[i], y[i], z[i], t[i])
            k1_z = dt * dzdt(x[i], y[i], z[i], t[i])

            k2_x = dt * dxdt(x[i] + .5 * k1_x, y[i] + .5 * k1_y, z[i] + .5 * k1_z, t[i] + .5 * dt)
            k2_y = dt * dydt(x[i] + .5 * k1_x, y[i] + .5 * k1_y, z[i] + .5 * k1_z, t[i] + .5 * dt)
            k2_z = dt * dzdt(x[i] + .5 * k1_x, y[i] + .5 * k1_y, z[i] + .5 * k1_z, t[i] + .5 * dt)

            k3_x = dt * dxdt(x[i] + .5 * k2_x, y[i] + .5 * k2_y, z[i] + .5 * k2_z, t[i] + .5 * dt)
            k3_y = dt * dydt(x[i] + .5 * k2_x, y[i] + .5 * k2_y, z[i] + .5 * k2_z, t[i] + .5 * dt)
            k3_z = dt * dzdt(x[i] + .5 * k2_x, y[i] + .5 * k2_y, z[i] + .5 * k2_z, t[i] + .5 * dt)

            k4_x = dt * dxdt(x[i] + k3_x, y[i] + k3_y, z[i] + k3_z, t[i+1])
            k4_y = dt * dydt(x[i] + k3_x, y[i] + k3_y, z[i] + k3_z, t[i+1])
            k4_z = dt * dzdt(x[i] + k3_x, y[i] + k3_y, z[i] + k3_z, t[i+1])

            x[i+1] = x[i] + (k1_x + 2*k2_x + 2*k3_x + k4_x) / 6
            y[i+1] = y[i] + (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6
            z[i+1] = z[i] + (k1_z + 2*k2_z + 2*k3_z + k4_z) / 6



    # def generate_parametric_graph(self):
    #     """
    #     Generates images ('.png's) of the plots of the Runge Kutta's
    #     solutions as well as the exact solutions.
    #     """

    #     approx_values = self.rk4_output
    #     u_list_approx = approx_values[1]
    #     v_list_approx = approx_values[2]

    #     fig, ax = plt.subplots(nrows = 1, ncols = 1)
    #     ax.plot(u_list_approx, v_list_approx, 'r')
    #     ax.set_xlabel("x(t)")
    #     ax.set_ylabel("y(t)")
    #     ax.set_title("Sombrero approximation (parametric) for n = " + str(self.n) + " and F = " + str(self.F))
    #     ax.grid(True)
    #     fig.savefig(self.__class__.__name__ + '_parametric'+ '_%3d' % (self.F) + '_%3d' % (self.x_0)  +  '_%3d' % (self.y_0) +  '_%3d.png' % (self.n))
    #     plt.close(fig)

    # def generate_x_graph(self):
    #     """
    #     Generates images ('.png's) of the plots of the Runge Kutta's
    #     solutions as well as the exact solutions.
    #     """

    #     approx_values = self.rk4_output
    #     x_list_approx = approx_values[1]
    #     t_list = approx_values[0]

    #     fig, ax = plt.subplots(nrows = 1, ncols = 1)
    #     ax.plot(t_list, x_list_approx, 'r')
    #     ax.set_xlabel("t")
    #     ax.set_ylabel("x(t)")
    #     ax.set_title("Sombrero approximation for x(t) for n = " + str(self.n) + " and F = " + str(self.F))
    #     ax.grid(True)
    #     fig.savefig(self.__class__.__name__ + '_x' + '_%3d' % (self.F) + '_%3d' % (self.x_0)  +  '_%3d' % (self.y_0) +  '_%3d.png' % (self.n))
    #     plt.close(fig)

    # def generate_poincare_graph(self, step, theta = 0 , xLim = (-1.5, 1.5), yLim = (-1.1, 1.1)):
    #     """
    #     Generates images ('.png's) of the plots of the Runge Kutta's
    #     solutions as well as the exact solutions.
    #     """
    #     approx_values = self.rk4_output
    #     u_list_approx = approx_values[1][int(theta*1000) :: step]
    #     v_list_approx = approx_values[2][int(theta*1000) :: step]

    #     #This method is more elegant, but rounding errors keep getting in the way
    #     #value_mask = np.where(approx_values[0] == step + theta)
    #     #u_list_approx = approx_values[1][value_mask]
    #     #v_list_approx = approx_values[2][value_mask]

    #     fig, ax = plt.subplots(nrows = 1, ncols = 1)
    #     ax.plot(u_list_approx[0::3], v_list_approx[0::3], 'r.')
    #     ax.plot(u_list_approx[1::3], v_list_approx[1::3], 'b.')
    #     ax.plot(u_list_approx[2::3], v_list_approx[2::3], 'g.')
    #     ax.set_xlabel("x")
    #     ax.set_ylabel("x'")
    #     ax.set_ylim(yLim)
    #     ax.set_xlim(xLim)
    #     ax.set_title("Sombrero approximation for theta = ." + str(theta))
    #     ax.grid(True)
    #     fig.savefig('frame%02d.png' % (100*theta))
    #     plt.close(fig)


class test_Rossler(TestCase):
    def test_linear(self, dt = .5, T = 500, x_prime = 2, y_prime = -1, z_prime = 1.5):
        """
        Runs the RK approximation implemented in the Rossler's run function, but
        with specifying different differential equations such that dx, dy, and dz / dt
        all are constant, and produce a linear result. RK 4 approximates linear functions
        perfectly by definition
        """
        test_R = Rossler(5, dt = dt, T = T)

        apt0 = len(test_R.t) == T // dt + 1 and test_R.t[-1] == T
        msg = "t was not properly initilaized"
        assert apt0, msg

        def constant_dxdt(x, y, z, t):
            return x_prime

        def constant_dydt(x, y, z, t):
            return y_prime

        def constant_dzdt(x, y, z, t):
            return z_prime

        test_R.dxdt = constant_dxdt
        test_R.dydt = constant_dydt
        test_R.dzdt = constant_dzdt

        test_vx = dt * x_prime * np.arange(0, T // dt + 1)
        test_vy = dt * y_prime * np.arange(0, T // dt + 1)
        test_vz = dt * z_prime * np.arange(0, T // dt + 1)

        test_R.run()


        aptx = np.array_equal(test_R.x, test_vx)
        msgx = "x was not correctly approximated"
        assert aptx, msgx

        apty = np.array_equal(test_R.y, test_vy)
        msgy = "y was not correctly approximated"
        assert apty, msgy

        aptz = np.array_equal(test_R.z, test_vz)
        msgz = "z was not correctly approximated"
        assert aptz, msgz
