# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 22:01:28 2018

@author: Alex Snouffer
"""

import sympy
import cantera
from matplotlib import pyplot
import numpy

def h_OutCompressor(n_c, hOI, hI):
    """
    """
    (n_comp, hOutIsentropic, hIn, 
     hOutActual) = sympy.symbols('n_comp hOutIsentropic hIn hOutActual')
    
    eEff = sympy.Eq(n_comp, (hOutIsentropic - hIn)/(hOutActual - hIn))
    e_h_Out = sympy.solve(eEff, hOutActual)
    h_Out = e_h_Out[0].evalf(subs = {n_comp: n_c, hOutIsentropic: hOI, 
                   hIn: hI})
    return h_Out

def h_OutTurbine(n_t, hOI, hI):
    """
    """
    (n_turbine, hOutIsentropic, hIn, 
     hOutActual) = sympy.symbols('n_turbine hOutIsentropic hIn hOutActual')
    
    eEff = sympy.Eq(n_turbine, (hIn - hOutActual)/(hIn - hOutIsentropic))
    e_h_Out = sympy.solve(eEff, hOutActual)
    h_Out = e_h_Out[0].evalf(subs = {n_turbine: n_t, hOutIsentropic: hOI, 
                   hIn: hI})
    return h_Out

def atm2Pa(P):
    P = P * 101325
    return P

def rev_irrev(hin, hout, sin, sout, Tin, Qin, Qout, mdotratio):
    To = 300
    irrev = mdotratio * (Qout + To * (sout - sin))
    rev = mdotratio * ((hin - hout) - To * (sin - sout) + Qin * (1 - (To/Tin)))
    rev_irrev_ans = [rev, irrev]
    return rev_irrev_ans
    
#   Brayton Cycle Paramters
n_comp_Brayton = 0.8
n_turbine_Brayton = 0.85
n_heatExchanger = 0.86
n_pump_Rankine = 0.9
n_turbine_Rankine = 0.9

W_NET = []
massRatio = []
Q_In = []
Q_Out = []
n_thermal = []

p5 = []
p6 = []
p7 = []
p8 = []
p9 = []
p1 = []
p2 = []
p3 = []
p4 = []

t5 = []
t6 = []
t7 = []
t8 = []
t9 = []
t1 = []
t2 = []
t3 = []
t4 = []

s5 = []
s6 = []
s7 = []
s8 = []
s9 = []
s1 = []
s2 = []
s3 = []
s4 = []

h5 = []
h6 = []
h7 = []
h8 = []
h9 = []
h1 = []
h2 = []
h3 = []
h4 = []

x3 = []

Wrev12 = []
Wrev23 = []
Wrev34 = []
Wrev41 = []
Wrev56 = []
Wrev67 = []
Wrev78 = []
Wrev89 = []

W12a = []
W23a = []
W34a = []
W41a = []
W56a = []
W67a = []
W78a = []
W89a = []

ri12 = []
ri23 = []
ri34 = []
ri41 = []
ri56 = []
ri67 = []
ri78 = []
ri89 = []

for pr in range(30, 201, 1):
    #Initial State of Brayton
    air5 = cantera.Solution('air.cti')
    p5.append(atm2Pa(1))
    t5.append(300)
    air5.TP = t5[-1], p5[-1]
    
    s5.append(air5.entropy_mass)
    h5.append(air5.enthalpy_mass)
    
    #State after Compressor
    air6 = cantera.Solution('air.cti')
    p6.append(atm2Pa(pr/10))
    air6.SP = s5[-1], p6[-1]
    t6.append(air6.T)
    
    h6_IS = air6.enthalpy_mass
    h6.append(h_OutCompressor(n_comp_Brayton, h6_IS, h5[-1]))
    air6.HP = h6[-1], p6[-1]
    s6.append(air6.entropy_mass)
    
    #State before Turbine (After Combustion)
    air7 = cantera.Solution('air.cti')
    p7.append(p6[-1])
    t7.append(1400)
    air7.TP = t7[-1], p7[-1]
    s7.append(air7.entropy_mass)

    h7.append(air7.enthalpy_mass)
    
    #State after Turbine
    air8 = cantera.Solution('air.cti')
    p8.append(p5[-1])
    air8.SP = s7[-1], p8[-1]
    h8_IS = air8.enthalpy_mass
    h8.append(h_OutTurbine(n_turbine_Brayton, h8_IS, h7[-1]))
    air8.HP = h8[-1], p8[-1]
    t8.append(air8.T)
    s8.append(air8.entropy_mass)
    
    #State After HRSG
    air9 = cantera.Solution('air.cti')
    p9.append(atm2Pa(1))
    t9.append(450)
    air9.TP = t9[-1], p9[-1]
    h9.append(air9.enthalpy_mass)
    s9.append(air9.entropy_mass)
    
    #State after the Condensor
    water1 = cantera.Water()
    p1.append(5 * 10**3)
    water1.PX = p1[-1], 0.0
    
    h1.append(water1.enthalpy_mass)
    s1.append(water1.entropy_mass)
    t1.append(water1.T)
    
    #State after Pump
    water2 = cantera.Water()
    p2.append(7 * 10**6)
    water2.SP = s1[-1], p2[-1]
    h2_IS = water2.enthalpy_mass
    h2.append(h_OutCompressor(n_pump_Rankine, h2_IS, h1[-1]))
    water2.HP = h2[-1], p2[-1]
    t2.append(water2.T)
    s2.append(water2.entropy_mass)
    
    #State after HRSG
    water3Perf = cantera.Water()
    p3.append(p2[-1])
    water3Perf.TP = t8[-1], p3[-1]
    h3perfect = water3Perf.enthalpy_mass
    h3.append(n_heatExchanger * (h3perfect - h2[-1]) + h2[-1])
    water3 = cantera.Water()
    water3.HP = h3[-1], p3[-1]
    
    s3.append(water3.entropy_mass)
    t3.append(water3.T)
    x3.append(water3.X)
    
    #State After Turbine Rankine
    water4 = cantera.Water()
    p4.append(p1[-1])
    water4.SP = s3[-1], p4[-1]
    h4_IS = water4.enthalpy_mass
    h4.append(h_OutTurbine(n_turbine_Rankine, h4_IS, h3[-1]))
    water4.HP = h4[-1], p4[-1]
    t4.append(water4.T)
    s4.append(water4.entropy_mass)
    
    W_BraytonComp = h5[-1] - h6[-1]
    W_BraytonTurbine = h7[-1] - h8[-1]
     
    massRatio.append((h8[-1] - h9[-1]) / (h3[-1] - h2[-1]))
    
    W_RankinePump = massRatio[-1] * (h1[-1] - h2[-1])
    W_RankineTurbine = massRatio[-1] * (h3[-1] - h4[-1])
    
    W_Top = W_BraytonTurbine + W_BraytonComp 
    W_Bottom = W_RankineTurbine + W_RankinePump
    W_NET.append(W_Top + W_Bottom)
    
    Q_In.append(h7[-1] - h6[-1])
    Q_Out.append(h4[-1] - h1[-1])
    
    n_thermal.append(W_NET[-1]/Q_In[-1])
    
    ri12.append(rev_irrev(h1[-1], h2[-1], s1[-1], s2[-1], t1[-1], 0, 0, massRatio[-1]))
    ri23.append(rev_irrev(h2[-1], h3[-1], s2[-1], s3[-1], t2[-1], (h3[-1] - h2[-1]), 0, massRatio[-1]))
    ri34.append(rev_irrev(h3[-1], h4[-1], s3[-1], s4[-1], t3[-1], 0, 0, massRatio[-1]))
    ri41.append(rev_irrev(h4[-1], h1[-1], s4[-1], s1[-1], t4[-1], 0, -Q_Out[-1], massRatio[-1]))
    
    ri56.append(rev_irrev(h5[-1], h6[-1], s5[-1], s6[-1], t5[-1], 0, 0, 1))
    ri67.append(rev_irrev(h6[-1], h7[-1], s6[-1], s7[-1], t6[-1], Q_In[-1], 0, 1))
    ri78.append(rev_irrev(h7[-1], h8[-1], s7[-1], s8[-1], t7[-1], 0, 0, 1))
    ri89.append(rev_irrev(h8[-1], h9[-1], s8[-1], s9[-1], t8[-1], 0, (h8[-1] - h9[-1]), 1))
    
    
pr = numpy.linspace(3.0, 20.0, 171)
pyplot.plot(pr, n_thermal)
bestPR = n_thermal.index(max(n_thermal))
print('The max effeciency is: ', max(n_thermal))
print('The most effecient Pressure Ratio is: ', (bestPR + 30)/10.0)


