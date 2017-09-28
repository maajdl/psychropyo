from psychrothermo import *
from psychrometry import *
from psychroutils import *
from pyomo.environ import Block

def mix(element, ha, hb, hc):
    b = Block()
    b.element = element
    b.airbal   = myCon(ha.A + hb.A == hc.A)
    b.waterbal = myCon(ha.V + ha.W + hb.V + hb.W == hc.V + hc.W)
    b.heatbal  = myCon(ha.H + hb.H == hc.H)
    return b

def split(element, ha, hb, hc, xb):
    b = Block()
    b.element = element

    b.airb     = myCon(hb.A == ha.A * xb)
    b.waterb   = myCon(hb.W == ha.W * xb)
    b.vaporb   = myCon(hb.V == ha.V * xb)
    b.heatb    = myCon(hb.H == ha.H * xb)

    xc = 1-xb
    b.airc     = myCon(hc.A == ha.A * xc)
    b.waterc   = myCon(hc.W == ha.W * xc)
    b.vaporc   = myCon(hc.V == ha.V * xc)
    b.heatc    = myCon(hc.H == ha.H * xc)

    return b

def flowsheet(element, *elements):
    b = Block()
    b.element = element
    for element in elements:
        setattr(b, element.element, element)
    return b
