#!/usr/bin/env python3

import numpy as np
import math
from abc import abstractmethod

class Problem:
    '''Base class for multiple itens definitions

    Attributes
    ----------
        N: defines dx space discretization
        M: defines dt time discretization
        Lambda: 
        T: integration time period [0, T]
        u: solution matrix
        gabarito: exact solution matrix
    ''' 


    def __init__(self, N, T, p:np.ndarray): #? tirar lambda e M pois N = M sempre. Incluir p
        '''Inits Problem with chosen N, M, Lambda and T
        '''
        self.N = N
        self.M = N
        self.T = T

        # Calulo de dx e dt
        self.dx = 1 / N
        self.dt = T / N
        
        self.Lambda = (self.dt)/(self.dx**2)

        #Forças Pontuais
        self.p = p
        self.nf = len(p)

        # Inicializacao matrizes M x N
        self.u = np.array([[0] * (self.N + 1)] * (self.M + 1), dtype=float) # auxiliar para o método de Nicolson
        
        self.uk = np.array([[0] * (self.N + 1)] * (self.nf), dtype=float) # vetores uk(T,xi) 
        self.gabarito = np.array([0] * (self.N + 1), dtype=float)  # calculada a partir de self.uk
        
        #todo vetor de coeficientes ak (intensidade das forças temporais)
        self.a = np.array([0] * (self.nf), dtype=float)

    def r(self,t):
        '''
        Função da variação temporal das forçantes r(t) = 10(1 + cos(5t))
        '''
        return 10*(1 + math.cos(5*t))


    def initial_condition(self): 
        '''
        Applies initial problem condition on solution matrix [u]
        '''
        self.u[0] = [0] * len(self.u[0])


    def frontier_condition(self): 
        '''
        Applies frontier condition on solution matrix [u]
        '''
        for i in self.u:
            i[0] = 0
            i[-1] = 0


    def heat_source(self, x, t, k):
        '''Heat source (pseudo) function f(x,t) = r(t) * gh(x)
        
        Args:
        -----
            x:   Postition
            t:   Time
            k:   Índice do vetor de pontos pk para o cálculo de gh(x)
        ''' 
        
        if ((self.p[k] - self.dx/2) <= x <= (self.p[k] + self.dx/2)):
            return self.r(t) * (1/self.dx) # gh(x) = 1/h 
        else: 
            return 0

    @abstractmethod
    def exact_solution(self):
        raise NotImplementedError

    

# Em todos os testes utilizaremos T = 1 e 

class Teste_A(Problem): #? rename
    def __init__(self, N, T, p:np.ndarray):
        super().__init__(N, T, p)
        self.item_name = 'A'

    # r(self,t)
 
    def exact_solution(self):
        # gabarito = a1 * u1 + a2 * u2 ...
        self.gabarito = 7.0 * self.uk[0]


class Teste_B(Problem): #? rename
    def __init__(self, N, T, p:np.ndarray):
        super().__init__(N, T, p)
        self.item_name = 'B'

    # r(self,t)
 
    def exact_solution(self):
        # gabarito = a1 * u1 + a2 * u2 ...
        # 2.3u1(T, xi) + 3.7u2(T, xi) + 0.3u3(T, xi) + 4.2u4(T, xi)
        self.gabarito = 2.3 * self.uk[0] + 3.7 * self.uk[1] + 0.3 * self.uk[2] + 4.2 * self.uk[3]

class Teste_C(Problem): #? rename
    def __init__(self, N, T, p:np.ndarray):
        super().__init__(N, T, p)
        self.item_name = 'C'

    # r(self,t)
 
    def exact_solution(self, uT:np.ndarray):
        # self.gabarito = np.zeros(self.N) # muda pra um vetor 1xN
        coef=int(2048/self.N)
        for i in range(0,self.N):
            self.gabarito[i] = uT[i*coef]

class Teste_D(Problem): #? rename
    def __init__(self, N, T, p:np.ndarray):
        super().__init__(N, T, p)
        self.item_name = 'D'

    # r(self,t)
 
    def exact_solution(self, uT:np.ndarray):
        # self.gabarito = np.zeros(self.N) # muda pra um vetor 1xN
        coef=int(2048/self.N)
        for i in range(0,self.N):
            self.gabarito[i] = uT[i*coef]
