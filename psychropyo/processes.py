from psychrothermo import air, H2O, Rgas, Tref, pref
from psychrometry import *
from pyomo.environ import Var, Constraint, Block, ConcreteModel, SolverFactory, sqrt, log, exp


def ipoptsolve(block):
    model       = ConcreteModel()
    model.block = block
    opt = SolverFactory('ipopt')
    block.result = opt.solve(model)
    if block.result.Solver[0].Status.key != 'ok': print(block.fixedVars)
    return block

Block.solve = ipoptsolve

def mix(ha, hb, hc):
    b = Block()
    b.ha = ha
    b.hb = hb
    b.hc = hc
    b.airbal   = myCon(ha.A + hb.A == hc.A)
    b.waterbal = myCon(ha.A * ha.Wv + hb.A * hb.Wv == hc.A * hc.Wv)
    b.heatbal  = myCon(ha.A * ha.hm + hb.A * hb.hm == hc.A * hc.hm)
    return b