# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 13:08:56 2018

@author: Alex Snouffer
"""

import cantera
import sympy

(m_turbine, m_pump, m_ratio, W, h_inTurbine, h_outTurbine,
 h_inPump, h_outPump, h_isentropicTurbine, h_isentropicPump,
 n_Turbine, n_Pump) = sympy.symbols('m_turbine m_pump m_ratio W h_inTurbine h_outTurbine h_inPump h_outPump h_isentropicTurbine h_isentropicPump n_Turbine n_Pump')

eBalanceTurbine = sympy.Eq(W, m_turbine * (h_inTurbine - h_isentropicTurbine))
print('\n',eBalanceTurbine,'\n')

eBalancePump = sympy.Eq(W, m_pump * (h_isentropicPump - h_inPump))
print(eBalancePump,'\n')

eCombined = sympy.Eq(-eBalanceTurbine.lhs + eBalancePump.lhs,
                     -eBalanceTurbine.rhs + eBalancePump.rhs)
print(eCombined,'\n')

mass_Ratio = sympy.Eq(m_ratio,
                      (eBalanceTurbine.rhs / m_turbine)/
                      (eBalancePump.rhs / m_pump))
print(mass_Ratio,'\n')

inSteamValue_Turbine = cantera.Water()
P_inTurbine = 10 * 10**5
T_inTurbine = 240 + 273
inSteamValue_Turbine.TP = T_inTurbine, P_inTurbine
hTi = inSteamValue_Turbine.enthalpy_mass

outWaterValue_Turbine = cantera.Water()
P_OutTurbine = 1.5 * 10**5
S_OutTurbine = inSteamValue_Turbine.entropy_mass
outWaterValue_Turbine.SP = S_OutTurbine, P_OutTurbine
hTo = outWaterValue_Turbine.enthalpy_mass

inWaterValue_Pump = cantera.Water()
T_inPump = 15 + 273
P_inPump = 1 * 10**5
inWaterValue_Pump.TP = T_inPump, P_inPump
hPi = inWaterValue_Pump.enthalpy_mass

outWaterValue_Pump = cantera.Water()
P_outPump = 60 * 10**5
S_outPump = inWaterValue_Pump.entropy_mass
outWaterValue_Pump.SP = S_outPump, P_outPump
hPo = outWaterValue_Pump.enthalpy_mass

massRatioSol = sympy.solve(mass_Ratio, m_ratio)
idealRatio = massRatioSol[0].evalf(subs={h_inTurbine: hTi,
                         h_isentropicTurbine: hTo,
                         h_inPump: hPi, h_isentropicPump: hPo})
print('The ideal Ratio is:', idealRatio,'\n')

eTurbineEff = sympy.Eq((h_inTurbine - h_outTurbine)/(h_inTurbine -
                       h_isentropicTurbine), n_Turbine)
print(eTurbineEff,'\n')

ePumpEff = sympy.Eq((h_isentropicPump - h_inPump)/(h_outPump -
                       h_inPump), n_Pump)
print(ePumpEff,'\n')

h_OutActTurbine = sympy.solve(eTurbineEff, h_outTurbine)
print(h_OutActTurbine,'\n')

h_OutActPump = sympy.solve(ePumpEff, h_outPump)
print(h_OutActPump,'\n')

eActualRatio = mass_Ratio.subs([(h_isentropicTurbine, h_OutActTurbine[0]), 
                                (h_isentropicPump, h_OutActPump[0])])
print(eActualRatio.simplify(),'\n')

solActualRatio = sympy.solve(eActualRatio, m_ratio)
print(solActualRatio,'\n')

ActualRatio = solActualRatio[0].evalf(subs = {n_Pump: 0.75, n_Turbine: 0.85,
                            h_inTurbine: hTi, h_isentropicTurbine: hTo,
                            h_inPump: hPi, h_isentropicPump: hPo})
print('The Actual Effeciency of the System is:', ActualRatio)

IdealRatio = solActualRatio[0].evalf(subs = {n_Pump: 1.0, n_Turbine: 1.0,
                            h_inTurbine: hTi, h_isentropicTurbine: hTo, 
                            h_inPump: hPi, h_isentropicPump: hPo})
print('The Ideal Effeciency of the System is: ',IdealRatio)

percentageRatio = ((IdealRatio - ActualRatio) / IdealRatio) * 100
print('The Percent Decrease from the Ideal to the Actual system is: ',
      percentageRatio)
