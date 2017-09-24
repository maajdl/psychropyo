from psychrothermo import air, H2O, Rgas, Tref, pref
from pyomo.environ import Var, Constraint, Block, ConcreteModel, SolverFactory, sqrt, log, exp
import pandas as pd
import holoviews as hv
import numpy as np
from IPython.display import display as Idisplay
import pickle
import os

def myVar(initialize, bounds, units):
    v = Var(bounds=bounds, initialize=initialize)
    v.units = units
    return v

def myCon(expr):
    return Constraint(expr=expr)

def humid_air(patm=101325, statename="", **fixedVars):
    ref = {  'A'   : 1.0,
             'Gv'  : -135.70142061190091,
             'T'   : 298.15,
             'Td'  : 283.62585311806976,
             'Tw'  : 289.4117368385543,
             'Wv'  : 0.007912721381495895,
             'cpha': 29.602790966545705,
             'hm'  : 45.54956827926809,
             'psat': 3169.7468549523596,
             'pvap': 1267.898741980944,
             'rh'  : 0.4,
             't'   : 25.0,
             'td'  : 10.475853118069793,
             'tw'  : 16.261736838554324,
             'vs'  : 0.8587547217602346 }
    big = 1E12
    b = Block()
    b.patm = patm
    b.statename = statename
    b.fixedVars = fixedVars

    b.t    = myVar(ref["t"]   , (-100, 250), "°C")
    b.tw   = myVar(ref["tw"]  , (-100, 250), "°C")
    b.td   = myVar(ref["td"]  , (-100, 250), "°C")
    b.T    = myVar(ref["T"]   , ( 173, 523), "K" )
    b.Tw   = myVar(ref["Tw"]  , ( 173, 523), "K" )
    b.Td   = myVar(ref["Td"]  , ( 173, 523), "K" )
    b.psat = myVar(ref["psat"], (   0, big), "Pa")
    b.pvap = myVar(ref["pvap"], (   0, big), "Pa")
    b.Gv   = myVar(ref["Gv"]  , (-big, big), "J/mol")
    b.A    = myVar(ref["A"]   , (   0, big), "kg")
    b.Wv   = myVar(ref["Wv"]  , (   0, big), "kg/kga")
    b.rh   = myVar(ref["rh"]  , (   0,   1), "-")
    b.hm   = myVar(ref["hm"]  , (-big, big), "kJ/kga")
    b.vs   = myVar(ref["vs"]  , (   0, big), "m³/kga")
    b.cpha = myVar(ref["cpha"], (   0, big), "J/mola/K")

    b.T_convert  = myCon(b.T == Tref + b.t)
    b.Tw_convert = myCon(b.Tw == Tref + b.tw)
    b.Td_convert = myCon(b.Td == Tref + b.td)

    bA1 = 1

    b.rh_def   = myCon( b.psat * b.rh == b.Wv / H2O.M / (bA1 / air.M + b.Wv / H2O.M) * patm )
    b.psat_def = myCon( b.psat == H2O.psat(b.T) )
    b.pvap_def = myCon( b.pvap == b.Wv / H2O.M / (bA1 / air.M + b.Wv / H2O.M) * patm )
    b.Gv_def   = myCon( b.Gv == b.Wv / H2O.M * (H2O.l.G0(b.T) + Rgas * b.T * log(b.pvap / b.psat)) )
    b.Td_def   = myCon( b.Wv / H2O.M / (bA1 / air.M + b.Wv / H2O.M) * patm == H2O.psat(b.Td) )
    b.hm_def   = myCon( bA1 * b.hm ==
                       bA1 / air.M * (air.H0(b.T) - air.H0(Tref)) + b.Wv / H2O.M * (H2O.g.H0(b.T) - H2O.l.H0(Tref)) )
    b.cpha_def = myCon( bA1 * b.cpha / air.M ==
                       bA1 / air.M * air.cp0(b.T) + b.Wv / H2O.M * H2O.g.cp0(b.T))
    b.Tw_def   = myCon( H2O.psat(b.Tw) / patm - b.Wv / H2O.M / (b.Wv / H2O.M + bA1 / air.M) ==
                       b.cpha / H2O.L(Tref) * (b.T - b.Tw))
    b.vs_def   = myCon( b.vs == (b.Wv / H2O.M + bA1 / air.M) * Rgas * b.T / patm / bA1 * 1000)

    for (var, val) in fixedVars.items():
        setattr(b, var+"_fixed", myCon(getattr(b, var) == val))

    return b

def ipoptsolve(block):
    model       = ConcreteModel()
    model.block = block
    opt = SolverFactory('ipopt')
    block.result = opt.solve(model)
    if block.result.Solver[0].Status.key != 'ok': print(block.fixedVars)
    return block

Block.solve = ipoptsolve

def valuelist(block, props):
    return [getattr(block,p).value for p in props]

Block.valuelist = valuelist

def state(ha):
    dic = {v.local_name: v.value for v in ha.component_objects(Var)}
    dic.update({'statename':ha.statename})
    return dic
Block.state = state

def statetable(lha, fmt='{:3.4g}'):
    pd.options.display.float_format = fmt.format
    df = [state(ha) for ha in lha]
    df = pd.DataFrame(df)
    df.index = df.statename
    df = df.drop('statename', 1)
    return df

def constraintTable(ha, fmt='{:3.4g}'):
    pd.options.display.float_format = fmt.format
    df = pd.DataFrame({"constraint": {c.local_name: c.body() for c in ha.component_objects(Constraint)}})
    return df
Block.constraintTable = constraintTable

def display(df, fmt='{:3.4g}'):
    pd.options.display.float_format = fmt.format
    Idisplay(df)
pd.DataFrame.display = display

class psy_curves:
    def __init__(self, patm=101325):
        self.patm = patm
        self.pickleFile = os.path.join(os.path.dirname(__file__), 'psy_curves_' + str(patm) + '.p')

    def read(self):
        try:
            self = pickle.load(open(self.pickleFile, "rb"))
            return self
        except:
            return self.write()

    def write(self):
        t_range  = np.arange(    0,   52,    2)
        rh_range = [0.001]+np.arange( 0.1, 1.1, 0.1).tolist()
        hm_range = np.arange(   10,  150,   10)
        tw_range = np.arange(    0,   40,    5)
        vs_range = np.arange( 0.78, 1.02, 0.02)
        patm = self.patm
        self.rh = {rh:[humid_air(patm=patm, A=1, rh=rh, t =t ).solve().valuelist(['t','Wv']) for t  in t_range ] for rh in rh_range}
        self.hm = {hm:[humid_air(patm=patm, A=1, hm=hm, rh=rh).solve().valuelist(['t','Wv']) for rh in rh_range] for hm in hm_range}
        self.tw = {tw:[humid_air(patm=patm, A=1, tw=tw, rh=rh).solve().valuelist(['t','Wv']) for rh in rh_range] for tw in tw_range}
        self.vs = {vs:[humid_air(patm=patm, A=1, vs=vs, rh=rh).solve().valuelist(['t','Wv']) for rh in rh_range] for vs in vs_range}
        pickle.dump(self, open(self.pickleFile, "wb"))
        return self

    def chart(self, points=None, paths=None, rh_alpha=0.5, hm_alpha=0.5, tw_alpha=0.0, vs_alpha=0.0):
        if points is None: points = []
        if paths is None: paths = []
        rh = hv.Overlay( [hv.Curve(v) for (k,v) in self.rh.items()] )
        hm = hv.Overlay( [hv.Curve(v) for (k,v) in self.hm.items()] )
        tw = hv.Overlay( [hv.Curve(v) for (k,v) in self.tw.items()] )
        vs = hv.Overlay( [hv.Curve(v) for (k,v) in self.vs.items()] )

        txyn = [(point.t.value, point.Wv.value, point.statename)   for point in points]
        po = hv.Points(txyn) * hv.Overlay([ hv.Text(x+1,y+0.001,n) for (x,y,n) in txyn ])

        pa = hv.Path([[(point.t.value, point.Wv.value) for point in path] for path in paths] )

        gr = hv.Overlay(  [hv.Curve([( 0, y),(50,   y)]) for y in np.arange(0, 0.042, 0.002)]+
                            [hv.Curve([( x, 0),( x,0.04)]) for x in np.arange(0,    52,     2)] )

        g = rh({'Curve' : {'style': {'color':'red'   , 'alpha':rh_alpha, 'hover_line_alpha':1}}}) *\
            hm({'Curve' : {'style': {'color':'blue'  , 'alpha':hm_alpha, 'hover_line_alpha':1}}}) *\
            tw({'Curve' : {'style': {'color':'green' , 'alpha':tw_alpha}}}) *\
            vs({'Curve' : {'style': {'color':'grey'  , 'alpha':vs_alpha}}}) *\
            po({'Points': {'style': {'color':'black' , 'size':5 }}})        *\
            pa({'Path'  : {'style': {'color':'black'}}})                    *\
            gr({'Curve' : {'style': {'color':'black' , 'alpha':0.1}}})

        g = g.redim(x ={'range': (0, 50   ),'label': 'Dry bulb temperature  [°C]'})
        g = g.redim(y ={'range': (0, 0.030),'label': 'Humidity ratio     [kg/kg]'})

        g = g({'Curve':{'style': {'line_width':0}}})
        g = g({'Curve':{'plot':{'tools':['hover'], 'toolbar':'above', 'yaxis':'right'}}})
        g = g({'Overlay':{'plot': {'width': 675, 'height':556, 'title_format': 'Psychrometric Chart (' + str(self.patm) + ' Pa)'}} })
        g = g({'Text':{'style': {'text_font_size':'8pt'}}})

        return g
