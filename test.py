import numpy as np
from psychrolib import *

# Establecer sistema SI
SetUnitSystem(SI)

# Constantes
P_atm = 101325  # Pa
Tdb_fixed = 180  # °C (fija para el test)

# Ticks de Humidity Ratio en g/kg
tick_positions = np.linspace(0, 70, 8)
for W_g in tick_positions:
    try:
        W = W_g / 1000  # a kg/kg
        Psat = GetSatVapPres(Tdb_fixed)
        Pv = (W * P_atm) / (0.622 + W)
        RH = Pv / Psat
        if RH > 1.0:
            RH = 1.0  # limitar a 100%
        Tdp_val = GetTDewPointFromRelHum(Tdb_fixed, RH)
    except Exception as e:
