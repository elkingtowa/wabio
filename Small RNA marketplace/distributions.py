import random
import matplotlib
import matplotlib.pyplot as plt

N  = 3 # Default size of population
mu = 100. # Default mean of population's wealth

def sample(distribution, N=N, mu=mu):
    "Sample from the distribution N times, then normalize results to have mean mu."
    return normalize([distribution() for _ in range(N)], mu * N)

def constant(mu=mu):          return mu
def uniform(mu=mu, width=mu): return random.uniform(mu-width/2, mu+width/2)
def gauss(mu=mu, sigma=mu/3): return random.gauss(mu, sigma) 
def beta(alpha=2, beta=3):    return random.betavariate(alpha, beta)
def pareto(alpha=4):          return random.paretovariate(alpha)
    
def normalize(numbers, total):
    "Scale the numbers so that they add up to total."
    factor = total / float(sum(numbers))
    return [x * factor for x in numbers]