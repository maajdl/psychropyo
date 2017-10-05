from pyomo.environ import Var, Constraint, Block, ConcreteModel, SolverFactory, Objective, VarList, ConstraintList
import pandas as pd

big = 1E12

def myVar(initialize, bounds, units):
    v = Var(bounds=bounds, initialize=initialize)
    v.units = units
    return v

def myCon(expr):
    return Constraint(expr=expr)


def ipoptsolve(block):
    model       = ConcreteModel()
    model.block = block
    opt = SolverFactory('ipopt')
    block.result = opt.solve(model)
    status = block.result.Solver[0].Status.key
    if status != 'ok': print("Solver returned : ", status)
    return block
Block.solve = ipoptsolve

def valuelist(block, props):
    return [getattr(block,p).value for p in props]
Block.valuelist = valuelist

def streamDic(ha):
    dic = {v.local_name: v.value for v in ha.component_objects(Var, descend_into=False)}
    dic.update({'psyname':ha.psyname})
    return dic
Block.streamDic = streamDic

def streamTable(lha, fmt='{:3.4g}'):
    pd.options.display.float_format = fmt.format
    df = [streamDic(ha) for ha in lha if ha.psytype1=="stream"]
    di = [ha.psyname for ha in lha if ha.psytype1=="stream"]
    df = pd.DataFrame(df, index=di)
    df = df.drop('psyname', 1)
    df = df.sort_index()
    return df

def constraintTable(ha, fmt='{:3.4g}'):
    pd.options.display.float_format = fmt.format
    df = pd.DataFrame({"constraint": {c.local_name: c.body() for c in ha.component_objects(Constraint)}})
    return df
Block.constraintTable = constraintTable

from IPython.display import display as Idisplay
from IPython.display import HTML

def Ihtml(str):
    Idisplay(HTML(str))
