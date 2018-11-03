# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 22:01:28 2018

@author: Alex Snouffer
"""

import sympy
import cantera

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
    
    eEff = sympy.Eq(n_turbine, (hIn - hOutIsentropic)/(hIn - hOutActual))
    e_h_Out = sympy.solve(eEff, hOutActual)
    h_Out = e_h_Out[0].evalf(subs = {n_turbine: n_t, hOutIsentropic: hOI, 
                   hIn: hI})
    return h_Out

def atm2Pa(P):
    P = P * 101325
    return P

#   Brayton Cycle Paramters
n_comp_Brayton = 0.8
n_turbine_Brayton = 0.85

for pr in range(3,20):
    airCompIn_Brayton = cantera.Oxygen()
    airPresCompIn_Brayton = atm2Pa(1)
    airTempCompIn_Brayton = 300
    airCompIn_Brayton.TP = airTempCompIn_Brayton, airPresCompIn_Brayton
    
    airEntropyComp_Brayton = airCompIn_Brayton.entropy_mass
    airEnthalpyCompIn_Brayton = airCompIn_Brayton.enthalpy_mass
    
    airCompOut_Brayton = cantera.Oxygen()
    airPresCompOut_Brayton = atm2Pa(1*pr)
    airCompOut_Brayton.SP = airEntropyComp_Brayton, airPresCompOut_Brayton
    
    airIsEnthalpyCompOut_Brayton = airCompOut_Brayton.enthalpy_mass
    
    airEnthalpyCompOut_Brayton = h_OutCompressor(n_comp_Brayton, 
                                                 airIsEnthalpyCompOut_Brayton,
                                                 airEnthalpyCompIn_Brayton)

    print(airEnthalpyCompOut_Brayton)
    
    airTurbineIn_Brayton = cantera.Oxygen()
    airPresTurbineIn_Brayton = airPresCompOut_Brayton
    airTempTurbineIn_Brayton = 1400
    airTurbineIn_Brayton.TP = airTempTurbineIn_Brayton,airPresTurbineIn_Brayton
    
    airEntropyTurbine_Brayton = airTurbineIn_Brayton.entropy_mass
    airEnthalpyTurbineIn_Brayton = airTurbineIn_Brayton.enthalpy_mass
    
    airTurbineOut_Brayton = cantera.Oxygen()
    airPresTurbineOut_Brayton = airPresCompIn_Brayton
    airTurbineOut_Brayton.SP = airEntropyTurbine_Brayton, airPresTurbineOut_Brayton
    
    airIsEnthalpyTurbineOut_Brayton = airTurbineOut_Brayton.enthalpy_mass
    
    airEnthalpyTurbineOut_Brayton = h_OutTurbine(n_turbine_Brayton, 
                                                 airIsEnthalpyTurbineOut_Brayton,
                                                 airEnthalpyTurbineIn_Brayton)
    
    print(airEnthalpyTurbineOut_Brayton,'\n')
    