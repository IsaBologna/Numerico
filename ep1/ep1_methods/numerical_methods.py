#!/usr/bin/env python3

import time
import os
import numpy as np
import plotly.graph_objects as go

from .tarefas import Problem


def finite_difference_method(Item: Problem):
    '''Finite difference method implementation

    Args:
    -----
        Item: 
    
    Returns:
    -------
        U solution matrix
    '''
    print("Resolvendo o problema pelo método de diferenças finitas...")

    # Inicialização da matriz com condições iniciais de fronteira do problema
    Item.initial_condition()
    Item.frontier_condition()

    begin_time = time.time()
    # Implementação da equação (11)
    for k in range(1, Item.M + 1):
        for i in range(1, Item.N):
            Item.u[k][i] = Item.u[k - 1][i] 
            Item.u[k][i] += Item.dt * ((Item.u[k - 1][i - 1] - 2*Item.u[k - 1][i] + Item.u[k - 1][i + 1])/(Item.dx ** 2) + Item.heat_source(i*Item.dx, (k - 1)*Item.dt))

    elapsed_time = time.time() - begin_time
    print("Tempo para a solucao: {:.4f} segundos".format(elapsed_time))

    return Item.u

def finite_difference_method_error(Item: Problem, time):
    '''Calculates the error between exact and approximate solution

    Args:
    -----
        Item:
        time: Exact time to get the error value

    Returns:
    --------
        Max absolute error value at that time
    '''

    # Erro local de truncamento
    tau = np.array([[0] * (Item.N + 1)] * (Item.M + 1), dtype=float)
    for k in range(1, Item.M + 1):
        for i in range(1, Item.N):
            tau[k-1][i] = (Item.gabarito[k][i] - Item.gabarito[k-1][i]) / Item.dt
            tau[k-1][i] -= (Item.gabarito[k-1][i-1] - 2*Item.gabarito[k-1][i] + Item.gabarito[k-1][i+1]) / (Item.dx**2)
            tau[k-1][i] -= Item.heat_source(i*Item.dx, (k-1)*Item.dt)

    # Erro entre as soluções
    e = np.array([[0] * (Item.N + 1)] * (Item.M + 1), dtype=float)
    for k in range(1, Item.M + 1):
        for i in range(1, Item.N):
            e[k][i] = e[k-1][i]
            e[k][i] += Item.dt * ((e[k-1][i-1] - 2*e[k-1][i] + e[k-1][i+1]) / (Item.dx**2) + tau[k-1][i])

    return max(abs(e[time]))

def plot_solution(Item: Problem):
    '''Plots solution U every t=0.1s

    Args:
    -----
        Item: 
    '''

    # Valores do eixo x
    x_axis = np.arange(0, 1, Item.dx)

    # Valores do eixo y: temperatura a cada 0.1s
    step = int(Item.M/10)
    y_axis = Item.u[::step]
    
    plot_title = "Solução item {}, N = {}, M = {}, Lambda = {}".format(Item.item_name, Item.N, Item.M, Item.Lambda)
    # Criação do gráfico
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
                        title=plot_title,
                        paper_bgcolor="white",
                        plot_bgcolor="white",
                        )

    # Plot grid
    plot.layout['yaxis'].update(dict(showgrid=True, gridcolor='#e6e6e6'))

    plot.show()

    # Salva gráfico .html
    if not os.path.exists("images"):
        os.mkdir("images")

    file_name = "Item_{}_{}_{}.html".format(Item.item_name, Item.N, str(Item.Lambda).replace('.', ''))
    plot.write_html("images/" + file_name)
    print("Gráfico {} salvo na pasta images/".format(file_name))


def create_A_matrix(l, size):
    '''Creates a representation of A with two arrays

    Args:
    -----
        l:    Lambda value
        size: Size of one dimension of A

    Returns:
    --------
        Matrix A represented as two arrays (A, B)
    '''

    A = np.array([1 + 2*l] * size, dtype=float)
    B = np.array([-l] * size, dtype=float)

    B[0] = 0  # valor l1 = 0

    return A, B

def LDL_decomposition(A, B):
    '''LDLT decomposition of a matrix represented by two arrays

    Args:
    -----
        A:  diagonal of the original matrix
        B:  subdiagonal of the original matrix

    Returns:
    -------
        LDL decomposition represented as two arrays (L, D)
    '''

    L = np.array([0] * len(A), dtype=float)
    D = np.array([0] * len(A), dtype=float)

    D[0] = A[0]
    L[0] = 0

    # Calculo dos valores 'Di' e 'Li' de cada matriz
    for i in range(1, len(A)):
        L[i] = B[i] / D[i - 1]
        D[i] = A[i] - (L[i]**2 * D[i - 1])
    
    return L, D


def implicit_euler_method(Item: Problem, L: np.array, D: np.array):
    '''Implicit Euler method implementation

    Args:
    -----
        Item:
        L:  array from the LDL decomposition
        D:  array from the LDL decomposition
    '''
    print("Resolvendo o problema pelo método de Euler implicito...")

    # Inicialização da matriz com condições iniciais de fronteira do problema
    Item.initial_condition()
    Item.frontier_condition()

    y = np.array([0] * (Item.N - 1), dtype=float)

    def calc_b(k):
        '''
        Função que calcula a matriz coluna do lado direito do sistema em função de k
        '''
        b = np.array([0] * (Item.N - 1), dtype=float)
        
        # Primeiro valor da matriz, considerando g1(tk)
        b[0] = Item.u[k-1][1] + Item.dt * Item.heat_source(1 * Item.dx, (k) * Item.dt) + Item.Lambda * Item.u[k][0]
        
        for i in range(1, Item.N-2):
            b[i] = Item.u[k-1][i+1] + Item.dt * Item.heat_source((i+1) * Item.dx, (k) * Item.dt)
        
        # Ultimo valor da matriz, considerando g2(tk)
        b[Item.N - 2] = Item.u[k-1][Item.N-1] + Item.dt * Item.heat_source((Item.N-1) * Item.dx, (k) * Item.dt) + Item.Lambda * Item.u[k][-1]
        
        return b

    begin_time = time.time()
    # Loop principal da solução do sistema "[L][D][Lt] [x] = [b]" para cada tk, k=1...M
    for k in range(1, Item.M + 1):
        b = calc_b(k) # calcula o vetor correspondente ao lado direito do sistema

        # Primeiro, resolvemos [L][D] [y] = [b]
        y[0] = b[0] / D[0]
        for i in range(1, Item.N-1):
            y[i] = ((b[i] - L[i] * (y[i-1] * D[i-1])) / D[i])

        # Por ultimo, [Lt] [x] = [y] para achar a solução x
        Item.u[k][Item.N-1] = y[-1]
        for i in range(Item.N-2, 0, -1):
            Item.u[k][i] = y[i-1] - (L[i] * Item.u[k][i+1])

    elapsed_time = time.time() - begin_time
    print("Tempo para a solucao: {:.4f} segundos".format(elapsed_time))


def crank_nicolson_method(Item: Problem, L, D):
    '''Crank-Nicolson method implementation

    Args:
    -----
        Item:
        L:  array from the LDL decomposition
        D:  array from the LDL decomposition
    '''
    print("Resolvendo o problema pelo método de Crank-Nicolson ...")

    Item.initial_condition()
    Item.frontier_condition()

    y = np.array([0] * (Item.N - 1), dtype=float)

    def calc_b(k):
        '''
        Função que calcula a matriz coluna do lado direito do sistema em função de k
        '''
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

    begin_time = time.time()

    # Loop principal da solução do sistema "[L][D][Lt] [x] = [b]" para cada tk, k=1...M
    for k in range(1, Item.M + 1):
        b = calc_b(k)

        # Primeiro, resolvemos [L][D] [y] = [b]
        y[0] = b[0] / D[0]
        for i in range(1, Item.N-1):
            y[i] = ((b[i] - L[i] * (y[i-1] * D[i-1])) / D[i])

        # Por ultimo, [Lt] [x] = [y] para achar a solução x
        Item.u[k][Item.N-1] = y[-1]
        for i in range(Item.N-2, 0, -1):
            Item.u[k][i] = y[i-1] - (L[i] * Item.u[k][i+1])


    elapsed_time = time.time() - begin_time
    print("Tempo para a solucao: {:.4f} segundos".format(elapsed_time))



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
    
    for k in range(1, Item.M + 1):
        for i in range(1, Item.N):
            e[k][i] = Item.gabarito[k][i] - Item.u[k][i] + tau[k][i]

    return max(abs(e[time]))