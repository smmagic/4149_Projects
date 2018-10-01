# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 12:42:55 2018

@author: Alex Snouffer
"""

import numpy

waterFName = 'WaterPropertyTables.csv'
airFName = 'AirPropertyTables.csv'

wtemp, wcond, wcp, whfg, wpr, wvisdyn, wvg, wvl = numpy.loadtxt(waterFName, delimiter = ',', skiprows = 2, unpack = True)
atemp, acond, acp, apr, adens, avisc = numpy.loadtxt(airFName, delimiter = ',', skiprows = 2, unpack = True)

def interpolation(value, x, y):
    interpValue = numpy.interp(value, x, y)
    return interpValue

## Variables for Problem
temp_Water = 17+273.15
water_Var = [wtemp, wcond, wcp, whfg, wpr, wvisdyn, wvg, wvl]
air_Var = [atemp, acond, acp, apr, adens, avisc]

water_interp = []

for n in range(0,7):
    i = interpolation(temp_Water, water_Var[0], water_Var[n+1])
    water_interp.append(i)
    print(i, "\n")


    
    
        
