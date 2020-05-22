#!/usr/bin/env python3

import time
import numpy as np
import plotly.graph_objects as go

from .tarefas import Problem


def finite_difference_method(Item: Problem):
    Item.initial_condition()
    Item.frontier_condition()

    for k in range(1, Item.M + 1):
        for i in range(1, Item.N):
            Item.u[k][i] = Item.u[k - 1][i] 
            Item.u[k][i] += Item.dt * ((Item.u[k - 1][i - 1] - 2*Item.u[k - 1][i] + Item.u[k - 1][i + 1])/(Item.dx ** 2) + Item.heat_source(i*Item.dx, (k - 1)*Item.dt))

    return Item.u

def trucation_error(Item: Problem):
    tau = np.array([[0] * (Item.N + 1)] * (Item.M + 1), dtype=float)
    
    for k in range(1, Item.M + 1):
        for i in range(1, Item.N):
            tau[k-1][i] = (Item.gabarito[k][i] - Item.gabarito[k-1][i]) / Item.dt
            tau[k-1][i] -= (Item.gabarito[k-1][i-1] - 2*Item.gabarito[k-1][i] + Item.gabarito[k-1][i+1]) / (Item.dx**2)
            tau[k-1][i] -= Item.heat_source(i*Item.dx, (k-1)*Item.dt)

    e = np.array([[0] * (Item.N + 1)] * (Item.M + 1), dtype=float)

    for k in range(1, Item.M + 1):
        for i in range(1, Item.N):
            e[k][i] = e[k-1][i]
            e[k][i] += Item.dt * ((e[k-1][i-1] - 2*e[k-1][i] + e[k-1][i+1]) / (Item.dx**2) + tau[k-1][i])

    return e

def plot(u, dx, M, title="Solução"):
    x_axis = np.arange(0, 1, dx)

    step = int(M/10)
    y_axis = u[::step]
    
    plot = go.Figure()
    for i in range(len(y_axis)):
        plot.add_trace(go.Scatter(
            x = x_axis,
            y = y_axis[i],
            line = dict(width=2),
            name = "Tk = {:.1f}".format(i*0.1)
        ))
    plot.update_layout(xaxis_title='Comprimento',
                        yaxis_title='Temperatura',
                        title=title,
                        paper_bgcolor="white",
                        plot_bgcolor="white",
                        )

    plot.layout['yaxis'].update(dict(showgrid=True, gridcolor='#e6e6e6'))

    plot.show()


def create_A_matrix(l, size):
    """
    Creates a representation of A with two arrays

    @param l    Lambda value
    @param size Size of one dimension of A
    """
    A = np.array([1 + 2*l] * size, dtype=float)
    B = np.array([-l] * size, dtype=float)

    B[0] = 0

    return A, B

def LDL_decomposition(A, B):
    L = np.array([0] * len(A), dtype=float)
    D = np.array([0] * len(A), dtype=float)

    D[0] = A[0]
    L[0] = 0

    for i in range(1, len(A)):
        L[i] = B[i] / D[i - 1]
        D[i] = A[i] - (L[i]**2 * D[i - 1])
    
    return L, D


def implicit_euler_method(Item: Problem, L: np.array, D: np.array):
    '''
    Solution is stored on Problem.u matrix
    '''
    Item.initial_condition()
    Item.frontier_condition()

    y = np.array([0] * (Item.N - 1), dtype=float)

    def calc_b(k):
        b = np.array([0] * (Item.N - 1), dtype=float)
        
        b[0] = Item.u[k-1][1] + Item.dt * Item.heat_source(1 * Item.dx, (k) * Item.dt) + Item.Lambda * Item.u[k][0]
        
        for i in range(1, Item.N-2):
            b[i] = Item.u[k-1][i+1] + Item.dt * Item.heat_source((i+1) * Item.dx, (k) * Item.dt)
        
        b[Item.N - 2] = Item.u[k-1][Item.N-1] + Item.dt * Item.heat_source((Item.N-1) * Item.dx, (k) * Item.dt) + Item.Lambda * Item.u[k][-1]
        
        return b

    for k in range(1, Item.M + 1):
        b = calc_b(k)

        y[0] = b[0] / D[0]
        for i in range(1, Item.N-1):
            y[i] = ((b[i] - L[i] * (y[i-1] * D[i-1])) / D[i])

        Item.u[k][Item.N-1] = y[-1]
        for i in range(Item.N-2, 0, -1):
            Item.u[k][i] = y[i-1] - (L[i] * Item.u[k][i+1])



def crank_nicolson_method(Item: Problem, L, D):
    Item.initial_condition()
    Item.frontier_condition()

    y = np.array([0] * (Item.N - 1), dtype=float)

    def calc_b(k):
        b = np.array([0] * (Item.N - 1), dtype=float)
        
        b[0] = (1-Item.Lambda)*Item.u[k-1][1] \
            + (Item.Lambda/2)*(Item.u[k-1][2] + Item.u[k-1][0]) \
            + (Item.dt/2) * (Item.heat_source(1*Item.dx, (k)*Item.dt) + Item.heat_source(1*Item.dx, (k-1)*Item.dt)) \
            + (Item.Lambda/2) * Item.u[k][0]
        
        for i in range(1, Item.N-2):
            b[i] = (1-Item.Lambda)*Item.u[k-1][i+1] \
                + (Item.Lambda/2)*(Item.u[k-1][i+2] + Item.u[k-1][i]) \
                + (Item.dt/2) * (Item.heat_source((i+1)*Item.dx, (k)*Item.dt) + Item.heat_source((i+1)*Item.dx, (k-1)*Item.dt))
        
        b[Item.N-2] = (1-Item.Lambda)*Item.u[k-1][Item.N-1] \
                + (Item.Lambda/2)*(Item.u[k-1][Item.N-2] + Item.u[k-1][Item.N]) \
                + (Item.dt/2) * (Item.heat_source((Item.N-1)*Item.dx, (k)*Item.dt) + Item.heat_source((Item.N-1)*Item.dx, (k-1)*Item.dt)) \
                + (Item.Lambda/2) * Item.u[k][-1]
        
        return b

    for k in range(1, Item.M + 1):
        b = calc_b(k)

        y[0] = b[0] / D[0]
        for i in range(1, Item.N-1):
            y[i] = ((b[i] - L[i] * (y[i-1] * D[i-1])) / D[i])

        Item.u[k][Item.N-1] = y[-1]
        for i in range(Item.N-2, 0, -1):
            Item.u[k][i] = y[i-1] - (L[i] * Item.u[k][i+1])





'''
b[i] = u[i] + dt*f(xi,tk+1)

z[i] = (b[i] - l[i] * z[i-1]) 

y[1] = b[1]
y[i] = [(b[i] - l[i] * (y[i-1] * d[i-1]))]  / d[i]  i=2...(N-1)

u[i] = y[i] - ( l[i+1]*u[i+1])  i=(N-1)...1


u[0] = g1 
u[N] = g2


'''

def error(Item: Problem):
    tau = np.array([[0] * (Item.N + 1)] * (Item.M + 1), dtype=float)
    
    for k in range(1, Item.M + 1):
        for i in range(1, Item.N):
            tau[k-1][i] = (Item.gabarito[k][i] - Item.gabarito[k-1][i]) / Item.dt
            tau[k-1][i] -= (Item.gabarito[k][i-1] - 2*Item.gabarito[k][i] + Item.gabarito[k][i+1]) / (Item.dx**2)
            tau[k-1][i] -= Item.heat_source(i*Item.dx, (k-1)*Item.dt)

    e = np.array([[0] * (Item.N + 1)] * (Item.M + 1), dtype=float)
    # !
    for k in range(1, Item.M + 1):
        for i in range(1, Item.N):
            e[k][i] = e[k-1][i]
            e[k][i] += Item.dt * ((e[k-1][i-1] - 2*e[k-1][i] + e[k-1][i+1]) / (Item.dx**2) + tau[k-1][i])

    return e