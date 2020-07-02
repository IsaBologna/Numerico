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
        self.M = N #!
        self.T = T

        # Calulo de dx e dt
        self.dx = 1 / N
        self.dt = T / N
        
        self.Lambda = (self.dt)/(self.dx**2)

        #Forças Pontuais
        self.p = p
        self.nf = len(p)

        # Inicializacao matrizes M x N
        self.u = np.array([[0] * (self.N + 1)] * (self.M + 1), dtype=float)
        self.gabarito = np.array([[0] * (self.N + 1)] * (self.M + 1), dtype=float)  # calculada a partir de self.u
        
        #todo vetor de coeficientes ak (intensidade das forças temporais)
        self.a = np.array([0] * (self.nf+1), dtype=float)

    def r(self,t):
        '''
        função da variação temporal das forçantes r(t) = 10(1 + cos(5t))
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


    @abstractmethod
    def heat_source(self, x, t, p):
        '''Heat source function f(x,t) = r(t) * gh(x)
        
        Args:
        -----
            x:   Postition
            t:   Time
        -----
            pk: forçantes pontuais
            (arg? ou pk em init e k como arg?)
        '''
        #? loop para todas as fontes pk? (arg?) 
        # p = 0.25
        if ((p - self.dx/2) <= x <= (p + self.dx/2)):
            return self.r(t) * (1/self.dx) # gh(x) = 1/h 
        else:
            return 0


    @abstractmethod
    def exact_solution(self):
        raise NotImplementedError

    

# Em todos os testes utilizaremos T = 1 e 

class Teste(Problem): #? rename
    def __init__(self, N, T):
        super().__init__(N, T)
        self.item_name = 'A'

    # r(self,t)
 
    def exact_solution(self):
        # gabarito = a1 * u1 + a2 * u2 ... 
        return


   