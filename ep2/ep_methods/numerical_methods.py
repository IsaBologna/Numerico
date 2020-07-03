import time
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

    #? tirar da função
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
    # print("Tempo para a solucao: {:.4f} segundos".format(elapsed_time))



'''
b[i] = u[i] + dt*f(xi,tk+1)

z[i] = (b[i] - l[i] * z[i-1]) 

y[1] = b[1]
y[i] = [(b[i] - l[i] * (y[i-1] * d[i-1]))]  / d[i]  i=2...(N-1)

u[i] = y[i] - ( l[i+1]*u[i+1])  i=(N-1)...1


u[0] = g1 
u[N] = g2


'''
#**** Cálculo dos Vetores uk
def matrix_uk(Item:Problem):
    ''' Create the matrix uk(T,xi) 
    '''
    for s in range(0,Item.nf):
        crank_nicolson_method(Item, s)
        Item.uk[s] = Item.u[Item.T]
        Item.u = np.zeros((Item.N+1,Item.M+1)) #? zerar elementos u

#* Resolução Sistema Normal
def solve_normal_system(Item:Problem): 
#todo
    '''
    '''
    
    A = np.zeros((Item.nf, Item.nf))
    B = np.zeros(Item.nf)

    Item.exact_solution() #calcula uT a partir de uk

    #* Calculo produto interno <u,v> = Σ(𝑖=1,𝑁−1) u(xi)v(xi)

    # Matriz A 
    for i in range(0,Item.nf):
        for j in range(0,Item.nf): #? dava pra fazer melhor acho
            A[i][j] = np.vdot(Item.uk[i], Item.uk[j])
    # print(A)

    # vector B
    for i in range(0,Item.nf):
        B[i] = np.vdot(Item.gabarito[i], Item.uk[i])
    # print(B)
    

    #* Decomposição LDLt
    D = np.zeros((Item.nf, Item.nf))
    L = np.zeros((Item.nf, Item.nf))

    ld = 0
    ldl = 0
    for i in range(0,Item.nf): # inverter i e j nas equações (prometo que faz sentido... ou n)
        for j in range(0,Item.nf): #??

            for k in range(0,i-1): # somatorias
                ld += L[i][k]**2 * D[k][k]
                ldl += L[j][k] * D[k][k] * L[i][k]

            D[i][i] = A[i][i] - ld 
            L[j][i] = (A[j][i] - ldl) / D[i][i]
    
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
    for i in range(0,Item.nf):
        sum_y = 0.0
        for k in range(0,i-1):
            sum_y += (L[i][k] * (y[k] * D[k][k]))
        y[i] = ( B[i] - sum_y)/ D[i][i]
    
    # (II):
    for i in range(Item.nf-1,-1,-1): #wtf 
        sum_x = 0.0
        for k in range(i,Item.nf): # só entra em i>nf
            sum_x += L[k][i] * Item.a[i]
        Item.a[i] = (y[i] - sum_x)

#***** Erro Quadrático *****
def quatratic_error(Item: Problem): #? Por que ta aqui e não dentro da classe?
    sum_mmq = 0.0
    for k in range(0, Item.nf):
        for i in range(0, Item.N-1 ):        
            # print(Item.uk[k][i])
            sum_mmq += ( Item.gabarito[k][i] - (Item.a[k] * Item.uk[k][i]) )**2
    return math.sqrt(Item.dx * sum_mmq)
#todo pensar nos plots