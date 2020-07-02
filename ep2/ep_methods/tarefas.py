#!/usr/bin/env python3

import numpy as np
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


    def __init__(self, N, T): #? tirar lambda e M pois N = M sempre
        '''Inits Problem with chosen N, M, Lambda and T
        '''
        self.N = N
        # self.M = M
        # self.Lambda = Lambda
        self.T = T

        # Calulo de dx e dt
        self.dx = 1 / N
        self.dt = T / N

        # Inicializacao matrizes N x N
        self.u = np.array([[0] * (self.N + 1)] * (self.N + 1), dtype=float)
        # self.gabarito = np.array([[0] * (self.N + 1)] * (self.N + 1), dtype=float)  #? Não tem mais solução gabarito

    @abstractmethod
    def initial_condition(self): 
        '''
        Applies initial problem condition on solution matrix [u]
        '''
        raise NotImplementedError


    @abstractmethod
    def frontier_condition(self): 
        '''
        Applies frontier condition on solution matrix [u]
        '''
        raise NotImplementedError


    @abstractmethod
    def heat_source(self, x, t):
        '''Heat source function f(x,t)
        
        Args:
        -----
            x:   Postition
            t:   Time

            pk: forçantes pontuais
            nf: size of pk ?
            r: function r(t)
        '''
        raise NotImplementedError


    @abstractmethod
    def exact_solution(self):
        raise NotImplementedError



# Em todos os testes utilizaremos T = 1 e r(t) = 10(1 + cos(5t))

class Teste_A(Problem):
    def __init__(self, N, T):
        super().__init__(N, T)
        self.item_name = 'A'

    def initial_condition(self):
        self.u[0] = [0] * len(self.u[0])
    
    def frontier_condition(self):
        for i in self.u:
            i[0] = 0
            i[-1] = 0

    def heat_source(self, x, t):
        # array p, function r(t)
        return 
    
    def exact_solution(self):

        return
