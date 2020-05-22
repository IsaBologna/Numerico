#!/usr/bin/env python3

import numpy as np
from abc import abstractmethod

class Problem:
    def __init__(self, N, M, Lambda, T):
        self.N = N
        self.M = M
        self.Lambda = Lambda
        self.T = T

        self.dx = 1 / N
        self.dt = T / M

        self.u = np.array([[0] * (self.N + 1)] * (self.M + 1), dtype=float)
        self.gabarito = np.array([[0] * (self.N + 1)] * (self.M + 1), dtype=float) # array M x N 

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
        '''
        Heat source function f(x,t)
        
        @param x    Postition
        @param t    Time
        '''
        raise NotImplementedError


    @abstractmethod
    def exact_solution(self):
        raise NotImplementedError




class Item_A1(Problem):
    def __init__(self, N, M, Lambda, T):
        super().__init__(N, M, Lambda, T)

    def initial_condition(self):
        self.u[0] = [0] * len(self.u[0])
    
    def frontier_condition(self):
        for i in self.u:
            i[0] = 0
            i[-1] = 0

    def heat_source(self, x, t):
        return 10*(x**2) * (x - 1) - (60*x*t) + 20*t

    def exact_solution(self):
        for k in range(0, self.M + 1):
                for i in range(0, self.N + 1):
                    self.gabarito[k][i] = 10 * (k*self.dt) * ((i*self.dx)**2) * ((i*self.dx) - 1)

class Item_A2(Problem):
    def __init__(self, N, M, Lambda, T):
        super().__init__(N, M, Lambda, T)

    def initial_condition(self):
        for i in range(0, self.N + 1):
            self.u[0][i] = (i*self.dx)**2 * (1 - (i*self.dx))**2

    def frontier_condition(self):
        for i in self.u:
            i[0] = 0
            i[-1] = 0

    def heat_source(self, x, t):
        return 10 * np.cos(10*t) * (x**2) * (1-x)**2 - (1 + np.sin(10* t)) * (12*(x**2) - 12*x + 2)

    def exact_solution(self):
        for k in range(0, self.M + 1):
                for i in range(0, self.N + 1):
                    self.gabarito[k][i] = (1 + np.sin(10*k*self.dt)) * ((i*self.dx)**2) * (1 - (i*self.dx))**2


class Item_B(Problem):
    def __init__(self, N, M, Lambda, T):
        super().__init__(N, M, Lambda, T)

    def initial_condition(self):
        for i in range(0, self.N + 1):
            self.u[0][i] = np.exp(-i*self.dx)

    def frontier_condition(self):
        for k in range(1, self.M + 1):
            self.u[k][0] = np.exp(k*self.dt)
            self.u[k][-1] = np.exp((k*self.dt) - 1) * np.cos(5*k*self.dt)

    def heat_source(self, x, t):
        return np.exp(t-x) * ((25*(t**2)) * np.cos(5*t*x) - 5*(x + 2*t) * np.sin(5*t*x))

    def exact_solution(self):
        for k in range(0, self.M + 1):
            for i in range(0, self.N + 1):
                self.gabarito[k][i] = np.exp((k*self.dt - i*self.dx)) * np.cos(5 * k*self.dt * i*self.dx)


class Item_C(Problem):
    def __init__(self, N, M, Lambda, T):
        super().__init__(N, M, Lambda, T)

    def initial_condition(self):
        self.u[0] = [0] * len(self.u[0])
        

    def frontier_condition(self):
        for i in self.u:
            i[0] = 0
            i[-1] = 0

    def heat_source(self, x, t):
        if ((0.25 - self.dx/2) <= x <= (0.25 + self.dx/2)):
            return 10000 * (1 - 2*t**2) * (1/self.dx)
        else:
            return 0
