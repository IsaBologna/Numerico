import time
import os
import numpy as np
import plotly.graph_objects as go

from .tarefas import Problem

# ! precisa modificar -> matriz A diferente
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

# * METODO CRANK-NICOLSON 

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

def create_normal_system(uk): 
    '''Creates a representation of A with two arrays

    Args:
    -----
        uk:     vectors calculated by Nicolson method, with pk points
        uT:     given vector ?

    Returns:
    --------
        Matrix A:   <uk,uk>
        Vector B:   <uT,uk> size k
    '''

    A = np.array([[0] * nf] * nf, dtype=float)
    B = np.array([0] * nf, dtype=float)

    # Calculo produto interno <u,v> = Σ(𝑖=1,𝑁−1) u(xi)v(xi)
    for i in range(1,nf+1): # de 1 a nf
        for j in range(1,nf+1):
            

    return A, B
