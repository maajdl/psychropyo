from psychrothermo import *
from psychrometry import *
from psychroutils import *
from pyomo.environ import Block

def mix(processname, ha, hb, hc):
    b = Block()
    b.processname = processname
    #b.ha = ha
    #b.hb = hb
    #b.hc = hc
    b.airbal   = myCon(ha.A + hb.A == hc.A)
    b.waterbal = myCon(ha.V + ha.W + hb.V + hb.W == hc.V + hc.W)
    b.heatbal  = myCon(ha.H + hb.H == hc.H)
    return b

def negStream(bb):
    b = Block()
    #b.bb = bb
    b.stream = "-" + bb.stream
    b.A    = myVar(0   , (-big,big), "kg")
    b.V    = myVar(0   , (-big,big), "kg")
    b.W    = myVar(0   , (-big,big), "kg")
    b.H    = myVar(0   , (-big,big), "kJ")
    b.A_neg = myCon(b.A == -bb.A)
    b.V_neg = myCon(b.V == -bb.V)
    b.W_neg = myCon(b.W == -bb.W)
    b.H_neg = myCon(b.H == -bb.H)
    return b
Block.__neg__ = negStream

def mix2(processname, ha, hb, hc):
    b = Block()
    b.processname = processname
    #b.ha = ha
    #b.hb = hb
    #b.hc = hc
    b.airbal   = myCon(ha.A + hb.A + hc.A == 0)
    b.waterbal = myCon(ha.V + ha.W + hb.V + hb.W + hc.V + hc.W == 0)
    b.heatbal  = myCon(ha.H + hb.H + hc.H == 0)
    return b