#!/usr/bin/env python3

import argparse
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

class Problem:
    def __init__(self, N, L, T):
        self.N = N
        self.dx = 1/N
        
        self.T = T

        self.M = int((T * N**2) / L)
        self.dt = T / self.M

        self.u = np.array([[0] * (N + 1)] * (self.M + 1), dtype=float)
        self.gabarito = np.array([[0] * (self.N + 1)] * (self.M + 1), dtype=float) # array M x N 

    def initial_condition(self):
        pass

    def frontier_condition(self):
        pass

    def heat_source(self, x, t):
        pass

    def solve(self):
        self.initial_condition()
        self.frontier_condition()

        for k in range(1, self.M + 1):
            for i in range(1, self.N):
                self.u[k][i] = self.u[k - 1][i] 
                self.u[k][i] += self.dt * ((self.u[k - 1][i - 1] - 2*self.u[k - 1][i] + self.u[k - 1][i + 1])/(self.dx ** 2) + self.heat_source(i*self.dx, (k - 1)*self.dt))

    def error(self):
        '''
        Erro local de truncamento
        '''

        tau = np.array([[0] * (self.N + 1)] * (self.M + 1), dtype=float)

        for k in range(1, self.M + 1):
            for i in range(1, self.N):
                tau[k-1][i] = (self.gabarito[k][i] - self.gabarito[k-1][i]) / self.dt
                tau[k-1][i] -= (self.gabarito[k-1][i-1] - 2*self.gabarito[k-1][i] + self.gabarito[k-1][i+1]) / (self.dx**2)
                tau[k-1][i] -= self.heat_source(i*self.dx, (k-1)*self.dt)


        self.e = np.array([[0] * (self.N + 1)] * (self.M + 1), dtype=float)

        for k in range(1, self.M + 1):
            for i in range(1, self.N):
                self.e[k][i] = self.e[k-1][i]
                self.e[k][i] += self.dt * ((self.e[k-1][i-1] - 2*self.e[k-1][i] + self.e[k-1][i+1]) / (self.dx**2) + tau[k-1][i])


    def plot_solution(self):
        x_axis = np.arange(0, 1, self.dx)

        step = int(self.M/10)
        y_axis = self.u[::step]

        plot = go.Figure()
        for i in range(len(y_axis)):
            plot.add_trace(go.Scatter(
                x = x_axis,
                y = y_axis[i],
                name = "tk = {:.1f}".format(i*0.1),
                line=dict(width=2)
            ))
        plot.update_layout(xaxis_title='Comprimento',
                            yaxis_title='Temperatura',
                            title="Solução")

        plot.show()

    def plot_gabarito(self):
        x_axis = np.arange(0, 1, self.dx)

        step = int(self.M/10)
        y_axis = self.gabarito[::step]

        plot = go.Figure()
        for i in range(len(y_axis)):
            plot.add_trace(go.Scatter(
                x = x_axis,
                y = y_axis[i],
                line=dict(width=2),
                name = "tk = {:.1f}".format(i*0.1),
            ))

        plot.update_layout(xaxis_title='Comprimento',
                            yaxis_title='Temperatura',
                            title="Gabarito")

        plot.show()

    def plot_heatmap(self):
        x_axis = np.arange(0, self.T, self.dt)

        y_axis = np.arange(0, 1, self.dx)
        
        fig1 = go.Figure(data=go.Heatmap(
                        x = x_axis,
                        y = y_axis,
                        z = self.u,
                        type = 'heatmap',
                        colorscale = 'Blues'))
        fig1.show()



class Item_A1(Problem):
    def __init__(self, N, L, T):
        super().__init__(N, L, T)

    def initial_condition(self):
        self.u[0] = [0] * len(self.u[0])

    def frontier_condition(self):
        for i in self.u:
            i[0] = 0
            i[-1] = 0

    def heat_source(self, x, t):
        return 10*(x**2) * (x - 1) - (60*x*t) + 20*t

    def funcao_gabarito(self):
        for k in range(0, self.M + 1):
                for i in range(0, self.N + 1):
                    self.gabarito[k][i] = 10 * (k*self.dt) * ((i*self.dx)**2) * ((i*self.dx) - 1)

class Item_A2(Problem):
    def __init__(self, N, L, T):
        super().__init__(N, L, T)

    def initial_condition(self):
        for i in range(0, self.N + 1):
            self.u[0][i] = (i*self.dx)**2 * (1 - (i*self.dx))**2

    def frontier_condition(self):
        for i in self.u:
            i[0] = 0
            i[-1] = 0

    def heat_source(self, x, t):
        return 10 * np.cos(10*t) * (x**2) * (1-x)**2 - (1 + np.sin(10* t)) * (12*(x**2) - 12*x + 2)

    def funcao_gabarito(self):
        for k in range(0, self.M + 1):
                for i in range(0, self.N + 1):
                    self.gabarito[k][i] = (1 + np.sin(10*k*self.dt)) * ((i*self.dx)**2) * (1 - (i*self.dx))**2

class Item_B(Problem):
    def __init__(self, N, L, T):
        super().__init__(N, L, T)

    def initial_condition(self):
        for i in range(0, self.N + 1):
            self.u[0][i] = np.exp(-i*self.dx)

    def frontier_condition(self):
        for k in range(1, self.M + 1):
            self.u[k][0] = np.exp(k*self.dt)
            self.u[k][-1] = np.exp((k*self.dt) - 1) * np.cos(5*k*self.dt)

    def heat_source(self, x, t):
        return np.exp(t-x) * ((25*(t**2)) * np.cos(5*t*x) - 5*(x + 2*t) * np.sin(5*t*x))

    def funcao_gabarito(self):
        for k in range(0, self.M + 1):
            for i in range(0, self.N + 1):
                self.gabarito[k][i] = np.exp((k*self.dt - i*self.dx)) * np.cos(5 * k*self.dt * i*self.dx)


class Item_C(Problem):
    def __init__(self, N, L, T):
        super().__init__(N, L, T)

    def initial_condition(self):
        for x in self.u[0]:
            self.x = 0

    def frontier_condition(self):
        for i in self.u:
            self.u[0] = 0
            self.u[-1] = 0

    def heat_source(self, x, t):
        if ((0.25 - self.dx/2) <= x <= (0.25 + self.dx/2)):
            return 10000 * (1 - 2*t**2) * (1/self.dx)
        else:
            return 0


# Args
N = 100
l = 0.25
T = 1

# a1 = Item_A1(N, l, T)
# a1.solve()
# a1.plot_solution()
# a1.plot_heatmap()

# a2 = Item_A2(N, l, T)
# a2.solve()
# a2.plot_solution()
# # a2.plot_heatmap()

# b = Item_B(N, l, T)
# b.solve()
# b.plot_solution()

# b.funcao_gabarito()
# b.plot_gabarito()

c = Item_C(N, l, T)
c.solve()
c.plot_solution()



