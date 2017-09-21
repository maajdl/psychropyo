from psychrothermo import air, H2O, Rgas, Tref, pref
from pyomo.environ import Var, Constraint, Block, ConcreteModel, SolverFactory, sqrt, log, exp
import pandas as pd
import holoviews as hv
import numpy as np
from IPython.display import display as Idisplay
import pickle
import os

def myVar(initialize, bounds, units, fmt):
    v = Var(bounds=bounds, initialize=initialize)
    v.units = units
    v.fmt = fmt
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
    
    b.t    = myVar(ref["t"]   , (-100, 250), "°C",       "{:8.3f}")
    b.tw   = myVar(ref["tw"]  , (-100, 250), "°C",       "{:8.3f}")
    b.td   = myVar(ref["td"]  , (-100, 250), "°C",       "{:8.3f}")
    b.T    = myVar(ref["T"]   , ( 173, 523), "K" ,       "{:8.3f}")
    b.Tw   = myVar(ref["Tw"]  , ( 173, 523), "K" ,       "{:8.3f}")
    b.Td   = myVar(ref["Td"]  , ( 173, 523), "K" ,       "{:8.3f}")
    b.psat = myVar(ref["psat"], (   0, big), "Pa",       "{:8.1f}")
    b.pvap = myVar(ref["pvap"], (   0, big), "Pa",       "{:8.1f}")
    b.Gv   = myVar(ref["Gv"]  , (-big, big), "J/mol",    "{:8.1f}")
    b.A    = myVar(ref["A"]   , (   0, big), "kg",       "{:8.3f}")
    b.Wv   = myVar(ref["Wv"]  , (   0, big), "kg",       "{:8.4f}")
    b.rh   = myVar(ref["rh"]  , (   0,   1), "-",        "{:8.3f}")
    b.hm   = myVar(ref["hm"]  , (-big, big), "kJ/kga",   "{:8.1f}")
    b.vs   = myVar(ref["vs"]  , (   0, big), "m³/kga",   "{:8.3f}")
    b.cpha = myVar(ref["cpha"], (   0, big), "J/mola/K", "{:8.3f}")

    b.T_convert  = myCon(b.T == Tref + b.t)
    b.Tw_convert = myCon(b.Tw == Tref + b.tw)
    b.Td_convert = myCon(b.Td == Tref + b.td)

    b.rh_def   = myCon( b.psat * b.rh == b.Wv / H2O.M / (b.A / air.M + b.Wv / H2O.M) * patm )
    b.psat_def = myCon( b.psat == H2O.psat(b.T) )
    b.pvap_def = myCon( b.pvap == b.Wv / H2O.M / (b.A / air.M + b.Wv / H2O.M) * patm )
    b.Gv_def   = myCon( b.Gv == b.Wv / H2O.M * (H2O.l.G0(b.T) + Rgas * b.T * log(b.pvap / b.psat)) )
    b.Td_def   = myCon( b.Wv / H2O.M / (b.A / air.M + b.Wv / H2O.M) * patm == H2O.psat(b.Td) )
    b.hm_def   = myCon( b.A * b.hm ==
                       b.A / air.M * (air.H0(b.T) - air.H0(Tref)) + b.Wv / H2O.M * (H2O.g.H0(b.T) - H2O.l.H0(Tref)) )
    b.cpha_def = myCon( b.A * b.cpha / air.M == 
                       b.A / air.M * air.cp0(b.T) + b.Wv / H2O.M * H2O.g.cp0(b.T))
    b.Tw_def   = myCon( H2O.psat(b.Tw) / patm - b.Wv / H2O.M / (b.Wv / H2O.M + b.A / air.M) ==
                       b.cpha / H2O.L(Tref) * (b.T - b.Tw))
    b.vs_def   = myCon( b.vs == (b.Wv / H2O.M + b.A / air.M) * Rgas * b.T / patm / b.A * 1000)

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
        
    def read(self):
        try:
            self = pickle.load(open(self.pickleFile, "rb"))
            return self
        except:
            return self.write()
    
    def chart(self, lha=[]):
        curveStyle   = {'Curve':{'style': {'line_width':0, 'hover_line_alpha':1}}}
        curvePlot    = {'Curve':{'plot':{'width': 675, 'height':556, 'tools':['hover'], 'toolbar':'above', 'yaxis':'right'}}}
        overlayPlot  = {'Overlay':{'plot': {'title_format': 'Psychrometric Chart'}} }
        gridStyle      = {'Curve':{'style': {'line_width':0, 'color':'black', 'alpha':0.2}}}
        
        process = [(ha.t.value, ha.Wv.value)  for ha in lha]
        p = hv.Overlay( [hv.Curve(v, kdims=['T'], vdims=['Wv']) for (k,v) in self.rh.items()] )
        q = hv.Overlay( [hv.Curve(v, kdims=['T'], vdims=['Wv']) for (k,v) in self.hm.items()] )
        r = hv.Overlay( [hv.Curve(v, kdims=['T'], vdims=['Wv']) for (k,v) in self.tw.items()] )
        s = hv.Overlay( [hv.Curve(v, kdims=['T'], vdims=['Wv']) for (k,v) in self.vs.items()] )
        t =              hv.Curve(process, kdims=['T'], vdims=['Wv'])
        u =              hv.Points(process, kdims=['T', 'Wv'])

        grid = hv.Overlay(
                [hv.Curve([( 0, y),(50,   y)], kdims=['T'], vdims=['Wv']) for y in np.arange(0, 0.042, 0.002)]+
                [hv.Curve([( x, 0),( x,0.04)], kdims=['T'], vdims=['Wv']) for x in np.arange(0,    52,     2)] )

        def psychart(rh,hm,tw,vs):
            pcolor  = {'Curve': {'style': {'color':'red'   , 'alpha':rh}}}
            qcolor  = {'Curve': {'style': {'color':'blue'  , 'alpha':hm}}}
            rcolor  = {'Curve': {'style': {'color':'green' , 'alpha':tw}}}
            scolor  = {'Curve': {'style': {'color':'grey'  , 'alpha':vs}}}
            tcolor  = {'Curve': {'style': {'color':'black' , 'alpha':1 }}}
            ucolor  = {'Points':{'style': {'size':3 , 'alpha':1 }}}
            
            g = p(pcolor) * q(qcolor) * r(rcolor) * s(scolor) * t(tcolor) * u(ucolor) * grid(gridStyle)
            g = g.redim(T ={'range': (0, 50   ),'label': 'Dry bulb temperature  [°C]'})
            g = g.redim(Wv={'range': (0, 0.030),'label': 'Humidity ratio  [kg/kg]'})
            return g
        
        dmap = hv.DynamicMap(psychart, kdims=['rh','hm','tw','vs'], )
        slider = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        g = dmap.redim.values(rh=slider,hm=slider, tw=slider, vs=slider)
        chart = g(curveStyle)(curvePlot)(overlayPlot)

        return chart
