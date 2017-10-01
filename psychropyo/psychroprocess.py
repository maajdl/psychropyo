from psychrothermo import *
from psychrometry import *
from psychroutils import *
from pyomo.environ import Block

def mix(element, ha, hb, hc):
    b = Block()
    b.element = element
    b.type1 = "process"
    b.type2 = "mix"
    b.airbal   = myCon(ha.A + hb.A == hc.A)
    #b.waterbal = myCon(ha.W + hb.W == hc.W)
    #b.vaporbal = myCon(ha.V + hb.V == hc.V)
    b.vaporbal = myCon(ha.V + hb.V + ha.W + hb.W == hc.V + hc.W)
    b.heatbal  = myCon(ha.H + hb.H == hc.H)
    return b

def split(element, ha, hb, hc, xb):
    b = Block()
    b.element = element
    b.type1 = "process"
    b.type2 = "split"

    b.airb     = myCon(hb.A == ha.A * xb)
    b.waterb   = myCon(hb.W == ha.W * xb)
    b.vaporb   = myCon(hb.V == ha.V * xb)
    b.heatb    = myCon(hb.H == ha.H * xb)
    b.m3b      = myCon(hb.m3 == ha.m3 * xb)

    xc = 1-xb
    b.airc     = myCon(hc.A == ha.A * xc)
    b.waterc   = myCon(hc.W == ha.W * xc)
    b.vaporc   = myCon(hc.V == ha.V * xc)
    b.heatc    = myCon(hc.H == ha.H * xc)
    b.m3c      = myCon(hc.m3 == ha.m3 * xc)

    return b

def exchange(element, ha, ht, hc, hA, tout):
    b = Block()
    b.element = element
    b.type1 = "process"
    b.type2 = "exchange"
    b.exch = myCon(ht.H == hA * (tout - hc.t))
    b.mix = mix(element+'mix', ha, ht, hc)
    return b


def flowsheet(element, *elements):
    b = Block()
    b.element = element
    b.type1 = "process"
    b.type2 = "flowsheet"
    for e in elements:
        setattr(b, element + "_" + e.element, e)
    return b
