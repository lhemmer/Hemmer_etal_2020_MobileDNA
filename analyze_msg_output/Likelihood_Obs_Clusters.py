#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: lucashemmer
"""

# libraries and functions to import
import csv
import pandas as pd
from math import factorial
from math import log, exp
import numpy as np
from scipy.optimize import minimize
from scipy.stats.distributions import chi2

####load data or just copy from below   
#train = pd.read_csv("/Users/lucashemmer/Downloads/Virilis_cluster_table.csv")

#dys_non = train["dys.nondys"].tolist()
dys_non = ['dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'dys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys', 'nondys']
#chrom_x = train["x.total"].tolist()
chrom_x = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#chrom_3 = train["3.total"].tolist()
chrom_3 = [0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 28, 0, 0, 15, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#prog_tot = train["total.prog"].tolist()
prog_tot = [2, 10, 2, 3, 16, 5, 49, 3, 3, 4, 2, 10, 14, 13, 24, 2, 14, 2, 1, 2, 32, 5, 1, 1, 1, 1, 1, 5, 11, 4, 1, 10, 6, 3, 13, 10, 2, 9, 6, 1, 1, 1, 1, 3, 10, 2, 1, 2, 3, 32, 2, 6, 30, 2, 11, 1, 4, 1, 1, 35, 1, 1, 2, 11, 3, 4, 40, 4, 8, 10, 4, 2, 10, 13, 8, 29, 19, 14, 9, 12, 20, 12, 17, 20, 26, 11, 9, 11, 9, 8, 10, 9, 9]


#### binomial equations

## prob cluster not observed
def prob_cluster_n_obs(beta,n):
    if n < 4:
        return float(1)
    else:
        p_no_obs_clust = ((factorial(n)/(factorial(0)*factorial(n)))*(beta**0)*((1-beta)**n)) + \
        ((factorial(n)/(factorial(1)*factorial(n-1)))*(beta**1)*((1-beta)**(n-1))) + \
        ((factorial(n)/(factorial(2)*factorial(n-2)))*(beta**2)*((1-beta)**(n-2))) + \
        ((factorial(n)/(factorial(3)*factorial(n-3)))*(beta**3)*((1-beta)**(n-3))) 
    return p_no_obs_clust

## prob cluster observed
def prob_obs_cluster(beta,n,k):
    prob_cluster = ((factorial(n)/(factorial(k)*factorial(n-k)))*(beta**k)*((1-beta)**(n-k)))
    return prob_cluster


#### parts of likelhood equation

## likelihood of no cluster           
def ln_no_cluster(alpha,beta,n):
    #ll = log(1-alpha) + log(alpha*prob_cluster_n_obs(beta, n))
    ll = log((1-alpha) + alpha*prob_cluster_n_obs(beta, n))
    return ll 

#likelihood of cluster
def ln_cluster(alpha,beta,n,k):
    ll = log(alpha*prob_obs_cluster(beta, n, k))
    #ll = log(alpha) + log(prob_obs_cluster(beta, n, k))
    return ll    

#if x and 3 have no clusters
def ln_neither_cluster(alpha,beta,n):
    ll = 2*ln_no_cluster(alpha,beta,n)
    return ll

#if one has a cluster, define k in likelihood function as chrom_x[i] or chrom_3[i]
def ln_one_cluster(alpha,beta,n,k):
    ll = ln_no_cluster(alpha,beta,n) + ln_cluster(alpha,beta,n,k)
    return ll

## likelihood calc
def ln_likelihood_calc(alpha,beta,n,cX,c3):
    if cX == 0 and c3 == 0:
        LL = ln_neither_cluster(alpha,beta,n)
    elif cX > 0 and c3 == 0:
        k = cX
        LL = ln_one_cluster(alpha,beta,n,k) 
    elif cX == 0 and c3 > 0:
        k = c3
        LL = ln_one_cluster(alpha,beta,n,k)
    return -LL


#### likelihood based on number of parameters

## two paramters, alpha and beta
def ln_likelihood_2(params):
    assert len(params) == 2
    LL = []
    alpha = params[0]
    beta = params[1]
    for i in range(len(prog_tot)):
        n = prog_tot[i]
        cX = chrom_x[i]
        c3 = chrom_3[i]
        if alpha <= 0.0:
            alpha = 0.00000000001
        elif alpha >= 1.0:
            alpha = 0.99999999999
        if beta <= 0.0:
            beta = 0.00000000001
        elif beta >= 1.0:
            beta = 0.99999999999
        LL.append( ln_likelihood_calc(alpha,beta,n,cX,c3) )
    return sum(LL)

## three paramters, alpha_dys, alpha_non, beta
def ln_likelihood_3a(params):
    assert len(params) == 3
    LL = [] 
    beta = params[2]
    for i in range(len(prog_tot)):
        n = prog_tot[i]
        cX = chrom_x[i]
        c3 = chrom_3[i]
        if dys_non[i] == 'dys':
            alpha = params[0]
            #LL.append( ln_likelihood_calc(alpha,beta,n,cX,c3) )
        elif dys_non[i] == 'nondys':
            alpha = params[1]
        if alpha <= 0.0:
            alpha = 0.00000000001
        elif alpha >= 1.0:
            alpha = 0.99999999999
        if beta <= 0.0:
            beta = 0.00000000001
        elif beta >= 1.0:
            beta = 0.99999999999
        #print(alpha)
        LL.append( ln_likelihood_calc(alpha,beta,n,cX,c3) )    
    return sum(LL)
  
# three parameters, alpha, beta_dys, beta_non    
def ln_likelihood_3b(params):
    assert len(params) == 3
    LL = [] 
    alpha = params[0]
    for i in range(len(prog_tot)):
        n = prog_tot[i]
        cX = chrom_x[i]
        c3 = chrom_3[i]
        if dys_non[i] == 'dys':
            beta = params[1]
            #LL.append( ln_likelihood_calc(alpha,beta,n,cX,c3) )
        elif dys_non[i] == 'nondys':
            beta = params[2]
        if alpha <= 0.0:
            alpha = 0.00000000001
        elif alpha >= 1.0:
            alpha = 0.99999999999
        if beta <= 0.0:
            beta = 0.00000000001
        elif beta >= 1.0:
            beta = 0.99999999999
        LL.append( ln_likelihood_calc(alpha,beta,n,cX,c3) )    
    return sum(LL)

# four parameters, alpha_dys, alpha_non, beta_dys, beta_non  
def ln_likelihood_4(params):
    assert len(params) == 4
    LL = [] 
    for i in range(len(prog_tot)):
        n = prog_tot[i]
        cX = chrom_x[i]
        c3 = chrom_3[i]
        if dys_non[i] == 'dys':
            alpha = params[0]
            beta = params[2]
        elif dys_non[i] == 'nondys':
            alpha = params[1]
            beta = params[3]
            #alphas.append(alpha)
            #betas.append(beta)
        if alpha <= 0.0:
            alpha = 0.00000000001
        elif alpha >= 1.0:
            alpha = 0.99999999999
        if beta <= 0.0:
            beta = 0.00000000001
        elif beta >= 1.0:
            beta = 0.99999999999
        LL.append( ln_likelihood_calc(alpha,beta,n,cX,c3) )
        #print(alpha)
    return sum(LL)


# Make a list of initial parameter guesses    
# alpha is probability of cluster, beta is frequency of cluster if cluster is observed
params2 = [0.01, 0.01] # alpha, beta
params3a = [0.01, 0.01, 0.01] #alpha_dys, alpha_non, beta
params3b = [0.01, 0.01, 0.01] #alpha, beta_dys, beta_non 
params4 = [0.01, 0.01, 0.01, 0.01] #alpha_dys, alpha_non, beta_dys, beta_non 


#### Run the minimizer, fun = -log likelihood, x: array values are alpha and beta
L2 = minimize(ln_likelihood_2, params2)#, method='nelder-mead')
#print(L2)
print("In the two parameter model,")
print("Alpha =",round(L2.x[0],3))
print("Beta =",round(L2.x[1],3))

#### Run the minimizer, fun = -log likelihood, x: array values are alpha_dys, alpha_non, beta
L3a = minimize(ln_likelihood_3a, params3a)#, method='nelder-mead')
#print(L3a)
print("In the three parameter model with two alpha parameter estimates,")
print("Alpha in dysgenic flies =",round(L3a.x[0],3))
print("Alpha in non-dysgenic flies =",round(L3a.x[1],3))
print("Beta =",round(L3a.x[2],3),"\n")

#### Run the minimizer, fun = -log likelihood, x: array values are alpha, beta_dys, beta_non
L3b = minimize(ln_likelihood_3b, params3b)#, method='nelder-mead')
#print(L3b)
print("In the three parameter model with two beta parameter estimates,")
print("Alpha  =",round(L3b.x[0],3))
print("Beta in dysgenic flies =",round(L3b.x[1],3))
print("Beta in non-dysgenic flies =",round(L3b.x[2],3),"\n")

#### Run the minimizer, fun = -log likelihood, x: array values are alpha_dys, alpha_non, beta_dys, beta_non
L4 = minimize(ln_likelihood_4, params4) #method='nelder-mead')
#print(L4)
print("In the four parameter model with separate alpha and beta parameter estimates,")
print("Alpha in dysgenic flies =",round(L4.x[0],3))
print("Alpha in non-dysgenic flies =",round(L4.x[1],3))
print("Beta in dysgenic flies =",round(L4.x[1],3))
print("Beta in non-dysgenic flies =",round(L4.x[2],3),"\n")


def likelihood_ratio(llmin, llmax):
    return(2*(llmax-llmin))
    
def return_statement(p):
    if p < 0.05:
        return(print("We favor the more model with additional estimated parameters\n"))
    else:
        return(print("We do not favor the the model with additional estimated parameters\n"))
    
#### Likelihood ratio test 1, single alpha/beta vs two alphas
LR = likelihood_ratio(L3a.fun,L2.fun)
p = chi2.sf(LR, 1) # L2 has 1 DoF more than L1
print("The p-value when comparing the 3 parameter / 2 alpha model to the two parameter model is")
print('p: %.20f' % p)
return_statement(p)

#### Likelihood ratio test 2, single alpha/beta vs two betas
LR = likelihood_ratio(L3b.fun,L2.fun)
p = chi2.sf(LR, 1) # L2 has 1 DoF more than L1
print("The p-value when comparing the 3 parameter / 2 beta model to the two parameter model is")
print('p: %.20f' % p)
return_statement(p)

#### Likelihood ratio test 2, one alpha / two betas vs two alphas / two betas
LR = likelihood_ratio(L4.fun,L3b.fun)
p = chi2.sf(LR, 1) # L2 has 1 DoF more than L1
print("The p-value when comparing the 4 parameter model to the 3 parameter / 2 beta model is")
print('p: %.20f' % p)
return_statement(p)




