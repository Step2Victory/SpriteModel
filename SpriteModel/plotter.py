# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 00:53:38 2022

@author: Step2Victory
"""
import matplotlib.pyplot as plt
import numpy as np




def plot1d(file):
    f = open(file, 'r')
    text = f.readlines()
    for i in range(len(text)):
        y = list(map(float, text[i].split()))
        x = np.linspace(0, 1, len(y))
        plt.plot(x,y)
        plt.xlabel('r')
        plt.ylabel('$H_\phi$')
        plt.legend('numerical')
        
    
    

x = np.linspace(0, 1, 1000)
y = 8 * x
plt.plot(x,y)
plt.legend('true')

plot1d('1d.txt')
        