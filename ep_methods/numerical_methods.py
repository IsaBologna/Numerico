#!/usr/bin/env python3

import os
import math
import numpy as np
import plotly.graph_objects as go

from .tarefas import Problem

# ****** METODO CRANK-NICOLSON ******

def crank_nicolson_method(Item: Problem, s):
    '''
    Crank-Nicolson method implementation
    ----
    Arg:
        s:  índice do vetor pk para a conta do heat_source
    '''
    # print("Resolvendo o problema pelo método de Crank-Nicolson ...")

    Item.initial_condition()
    Item.frontier_condition()

    y = np.array([0] * (Item.N - 1), dtype=float)

    def LDL_decomposition():
        '''
        LDLT decomposition of a matrix represented by two arrays
        '''

        #* Creates a representation of A with two arrays (A, B)
        A = np.array([1 + 2 * (Item.Lambda/2)] * Item.N, dtype=float)
        B = np.array([-(Item.Lambda/2)] * Item.N, dtype=float)
        B[0] = 0  # valor l1 = 0

        L = np.array([0] * len(A), dtype=float)
        D = np.array([0] * len(A), dtype=float)

        D[0] = A[0]
        L[0] = 0

        # Calculo dos valores 'Di' e 'Li' de cada matriz
        for i in range(1, len(A)):
            L[i] = B[i] / D[i - 1]
            D[i] = A[i] - (L[i]**2 * D[i - 1])
        
        return L, D

    L, D = LDL_decomposition()

    def calc_b(k):
        '''
        Função que calcula a matriz coluna do lado direito do sistema em função de k
        '''
        b = np.array([0] * (Item.N - 1), dtype=float)
        
        b[0] = (1-Item.Lambda)*Item.u[k-1][1] \
            + (Item.Lambda/2)*(Item.u[k-1][2] + Item.u[k-1][0]) \
            + (Item.dt/2) * (Item.heat_source(1*Item.dx, (k)*Item.dt, s) + Item.heat_source(1*Item.dx, (k-1)*Item.dt, s)) \
            + (Item.Lambda/2) * Item.u[k][0]
        
        for i in range(1, Item.N-2):
            b[i] = (1-Item.Lambda)*Item.u[k-1][i+1] \
                + (Item.Lambda/2)*(Item.u[k-1][i+2] + Item.u[k-1][i]) \
                + (Item.dt/2) * (Item.heat_source((i+1)*Item.dx, (k)*Item.dt, s) + Item.heat_source((i+1)*Item.dx, (k-1)*Item.dt, s))
        
        b[Item.N-2] = (1-Item.Lambda)*Item.u[k-1][Item.N-1] \
                + (Item.Lambda/2)*(Item.u[k-1][Item.N-2] + Item.u[k-1][Item.N]) \
                + (Item.dt/2) * (Item.heat_source((Item.N-1)*Item.dx, (k)*Item.dt, s) + Item.heat_source((Item.N-1)*Item.dx, (k-1)*Item.dt, s)) \
                + (Item.Lambda/2) * Item.u[k][-1]
        
        return b

    #* Loop principal da solução do sistema "[L][D][Lt] [x] = [b]" para cada tk, k=1...M
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



#**** Cálculo dos Vetores uk
def matrix_uk(Item:Problem):
    ''' 
    Cria a matriz uk(T,xi) 
    ----
    Aproximações da solução usando o método de Crank-Nicolson para cada fonte de calor em pk
    '''
    for s in range(0, Item.nf):
        crank_nicolson_method(Item, s)
        Item.uk[s] = Item.u[-1]
        # print(Item.uk)
        Item.u = np.zeros((Item.N+1,Item.M+1))  # zerar elementos u

#* Resolução Sistema Normal
def solve_normal_system(Item:Problem): 
    '''
    Resolução do sistema normal para encontrar os coeficientes ak
    ----

    '''
    
    #* Matrizes para o sistema da forma Ax = B
    A = np.zeros((Item.nf, Item.nf))
    B = np.zeros(Item.nf)

    #* Calculo produto interno <u,v> = Σ(𝑖=1,𝑁−1) u(xi)v(xi) 
    for i in range(0, Item.nf):
        B[i] = np.vdot(Item.gabarito, Item.uk[i])
        for j in range(0,Item.nf): 
            A[i][j] = np.vdot(Item.uk[i], Item.uk[j])
    # print(A)
    # print(B)

    #* Decomposição LDLt
    D = np.zeros((Item.nf, Item.nf))
    L = np.zeros((Item.nf, Item.nf))

    for i in range(0, Item.nf):
        for j in range(0, i+1):
            if i == j:
                L[i][i] = 1.0  # Diagonal principal de L = 1

                ld = 0.0
                for k in range(0, j):
                    ld += L[j][k]**2 * D[k][k]
                D[j][j] = A[j][j] - ld

            else:
                ldl = 0.0
                for k in range(0, j):
                    ldl += L[i][k] * D[k][k] * L[j][k]
                L[i][j] = (A[i][j] - ldl) / D[j][j]
    # print("Matriz L:")
    # print(L)
    # print(" ")
    # print("Matriz D:")
    # print(D)
    # print(" ")

    #* Resolver o sistema
    ''' 
    Loop principal de resolução do sistema L.D.Lt * ak = b
    ---
    Eq:
        Lt.a = y -> a = y.L
        L.D.y = b
    Solve:
        yi = (bi - Σ(k=1,i-1) Lik * (yk * dkk) ) / dii (I)     , i = 1...n
        ai = yi - Σ(k=i+1,n) lki . ak       (II)    , i = n...1
    '''
    # (I):
    y = np.zeros(Item.nf)
    for i in range(0, Item.nf):
        sum_y = 0.0
        for k in range(0, i):
            sum_y += (L[i][k] * (y[k] * D[k][k]))
        y[i] = ( B[i] - sum_y)/ D[i][i]
    
    # (II):
    for i in range(Item.nf-1, -1, -1): 
        sum_x = 0.0
        for k in range(i+1, Item.nf): # só entra em i<nf
            sum_x += L[k][i] * Item.a[k]
        Item.a[i] = (y[i] - sum_x)

#***** Erro Quadrático *****
def quatratic_error(Item: Problem):
    '''
    Calculo do Erro Quadratico
    '''
    sum_mmq = 0
    for i in range(1, Item.N - 1):
        calc_solution = 0.0
        for k in range(0, Item.nf):
            calc_solution += Item.a[k] * Item.uk[k][i]

        sum_mmq += (Item.gabarito[i]- calc_solution) ** 2

    return math.sqrt(Item.dx * sum_mmq)


def plot_solution(Item: Problem):
    '''
    Grafico contendo a solução calculada e a proveniente do arquivo, em T = 1
    '''

    # Valores do eixo X
    x_axis = np.arange(0, 1, Item.dx)

    plot_title = "Solução item {}, com N = {}".format(Item.item_name, Item.N)

    # Criação do gráfico
    plot = go.Figure()
    
    #* Dados de uT provenientes do arquivo fornecido
    plot.add_trace(go.Scatter(
        x = x_axis,
        y = Item.gabarito,
        line = dict(width=6),
        name = "Arquivo"
    ))

    #* Solucao calculada com os coeficientes das intesidades recuperadas
    _solucao_calculada = np.zeros(Item.N+1)
    for k in range(0, Item.nf):
        _solucao_calculada += Item.a[k] * Item.uk[k]  # uT = Σ(ak * uk(T))

    plot.add_trace(go.Scatter(
        x = x_axis,
        y = _solucao_calculada,
        line = dict(width=2),
        name = "Solução Recuperada"
    ))

    plot.update_layout(
        xaxis_title = 'Comprimento',
        yaxis_title = 'Temperatura',
        title = plot_title,
        paper_bgcolor = "white",
        plot_bgcolor = "white"
    )

    # Plot grid
    plot.layout['yaxis'].update(dict(showgrid=True, gridcolor='#e6e6e6'))

    plot.show()

     # Salva gráfico .html
    if not os.path.exists("images"):
        os.mkdir("images")

    file_name = "Item_{}_{}.html".format(Item.item_name, Item.N)
    plot.write_html("images/" + file_name)
    print("Gráfico {} salvo na pasta images/".format(file_name))
