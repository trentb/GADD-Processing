# -*- coding: utf-8 -*-
"""
Created on Sat Dec 27 21:10:43 2014

@author: Trent
"""

filename = dat #variable
freq = 1000 #Hz
t = 1.000 #microns
d = 600 #microns
minpts = 10
maxpts = 10
maxi = 1

f=dat["swp"]["FREQUENCY_set"].index(freq)

M=[filename["swp"]["AC_set"],filename["meas"]["CD0"]["real"]()]