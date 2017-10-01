from psychroutils import *
from psychrothermo import *
from psychrochart import *
from psychroprocess import *
from pyomo.environ import Block

def heat(element='heat', **fixedVars):
    b = Block()
    b.element = element
    b.type1 = "stream"
    b.type2 = "heat"
    b.H    = myVar(0   , (-big,big), "kJ")      # extensive quantities
    b.A    = myVar(0   , (   0,  0), "kg")
    b.V    = myVar(0   , (   0,  0), "kg")
    b.W    = myVar(0   , (   0,  0), "kg")
    for (var, val) in fixedVars.items(): setattr(b, var+"_fixed", myCon(getattr(b, var) == val))
    return b

def water(element='water', **fixedVars):
    b = Block()
    b.element = element
    b.type1 = "stream"
    b.type2 = "water"
    b.T = myVar(300, (173, 523), "K")
    b.t = myVar( 25, (-100, 250), "°C")
    b.H    = myVar(0   , (-big,big), "kJ")      # extensive quantities
    b.A    = myVar(0   , (   0,  0), "kg")
    b.V    = myVar(0   , (   0,  0), "kg")
    b.W    = myVar(0   , (-big,big), "kg")
    b.T_convert = myCon(b.T == Tref + b.t)
    b.H_def = myCon(b.H == (H2O.l.H0(b.T) - H2O.l.H0(Tref))/H2O.M * b.W)
    for (var, val) in fixedVars.items(): setattr(b, var+"_fixed", myCon(getattr(b, var) == val))
    return b

def vapor(element='vapor', **fixedVars):
    b = Block()
    b.element = element
    b.type1 = "stream"
    b.type2 = "vapor"
    b.T = myVar(300, (173, 523), "K")
    b.t = myVar( 25, (-100, 250), "°C")
    b.H    = myVar(0   , (-big,big), "kJ")      # extensive quantities
    b.A    = myVar(0   , (   0,  0), "kg")
    b.V    = myVar(0   , (-big,big), "kg")
    b.W    = myVar(0   , (   0,  0), "kg")
    b.T_convert = myCon(b.T == Tref + b.t)
    b.H_def = myCon(b.H == (H2O.g.H0(b.T) - H2O.l.H0(Tref))/H2O.M * b.V)
    for (var, val) in fixedVars.items(): setattr(b, var+"_fixed", myCon(getattr(b, var) == val))
    return b

def humid_air(element, patm=101325, **fixedVars):
    ref = {  'A': 1.0,
             'H': 45.54956808796108,
             'T': 298.15,
             'Td': 283.6258531107488,
             'Tw': 289.41173683833046,
             'V': 0.007912721306430524,
             'W': 0,
             'cp': 29.602790962509385,
             'element': 'ha',
             'h': 45.54956808796108,
             'm3': 0.8587547216582929,
             'psat': 3169.7468549523596,
             'pvap': 1267.898741980944,
             'rh': 0.4,
             't': 25.0,
             'td': 10.475853110748833,
             'tw': 16.261736838330496,
             'v': 0.007912721306430524,
             'vs': 0.8587547216582929}
    b = Block()
    b.patm = patm
    b.element = element
    b.type1 = "stream"
    b.type2 = "humid air"

    b.t    = myVar(ref["t"]   , (-100, 250), "°C")          # intensive quantities
    b.tw   = myVar(ref["tw"]  , (-100, 250), "°C")
    b.td   = myVar(ref["td"]  , (-100, 250), "°C")
    b.T    = myVar(ref["T"]   , ( 173, 523), "K" )
    b.Tw   = myVar(ref["Tw"]  , ( 173, 523), "K" )
    b.Td   = myVar(ref["Td"]  , ( 173, 523), "K" )
    b.psat = myVar(ref["psat"], (   0, big), "Pa")
    b.pvap = myVar(ref["pvap"], (   0, big), "Pa")
    b.rh   = myVar(ref["rh"]  , (   0,   1), "-")
    b.v    = myVar(ref["v"]   , (   0, big), "kg/kga")
    b.h    = myVar(ref["h"]   , (-big, big), "kJ/kga")
    b.vs   = myVar(ref["vs"]  , (   0, big), "m³/kga")
    b.m3   = myVar(ref["m3"]  , (   0, big), "m³/kga")
    b.cp   = myVar(ref["cp"]  , (   0, big), "J/mola/K")
    b.H    = myVar(ref["H"]   , (-big, big), "kJ")          # extensive quantities
    b.A    = myVar(ref["A"]   , (   0, big), "kg")
    b.V    = myVar(ref["V"]   , (   0, big), "kg")
    b.W    = myVar(ref["W"]   , (   0,   0), "kg")

    bA1 = 1
    b.T_convert  = myCon(b.T == Tref + b.t)
    b.Tw_convert = myCon(b.Tw == Tref + b.tw)
    b.Td_convert = myCon(b.Td == Tref + b.td)
    b.Tw_def     = myCon( H2O.psat(b.Tw) / patm - b.v / H2O.M / (b.v / H2O.M + bA1 / air.M) == b.cp / H2O.L(Tref) * (b.T - b.Tw))
    b.Td_def     = myCon( b.v / H2O.M / (bA1 / air.M + b.v / H2O.M) * patm == H2O.psat(b.Td) )
    b.psat_def   = myCon( b.psat == H2O.psat(b.T) )
    b.pvap_def   = myCon( b.pvap == b.v / H2O.M / (bA1 / air.M + b.v / H2O.M) * patm )
    b.rh_def     = myCon( b.psat * b.rh == b.v / H2O.M / (bA1 / air.M + b.v / H2O.M) * patm )
    b.h_def      = myCon( bA1 * b.h == bA1 / air.M * (air.H0(b.T) - air.H0(Tref)) + b.v / H2O.M * (H2O.g.H0(b.T) - H2O.l.H0(Tref)) )
    b.vs_def     = myCon( b.vs == (b.v / H2O.M + bA1 / air.M) * Rgas * b.T / patm / bA1 * 1000)
    b.m3_def     = myCon( b.m3 == b.vs * b.A)
    b.cp_def     = myCon( bA1 * b.cp / air.M == bA1 / air.M * air.cp0(b.T) + b.v / H2O.M * H2O.g.cp0(b.T))
    b.V_def      = myCon(b.V == b.v * b.A)
    b.H_def      = myCon(b.H == b.h * b.A)
    for (var, val) in fixedVars.items(): setattr(b, var+"_fixed", myCon(getattr(b, var) == val))

    return b
