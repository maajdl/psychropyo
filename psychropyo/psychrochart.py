from psychropyo.psychrostream import *
import numpy as np
import pickle
import os

import holoviews as hv
hv.extension('bokeh')

class PsychrometricChart:
    def __init__(self, patm=101325):
        self.patm = patm
        self.pickleFile = os.path.join(os.path.dirname(__file__), 'PsychrometricCharts', 'PsychrometricChart_' + str(patm) + '.p')

    def read(self):
        try:
            self = pickle.load(open(self.pickleFile, "rb"))
        except:
            print("Preparing chart file: " + self.pickleFile)
            self = self.write(open(self.pickleFile, "wb"))
        return self

    def write(self, file):
        t_range  = np.arange(    0,   52,    2)
        rh_range = [0.001]+np.arange( 0.1, 1.1, 0.1).tolist()
        hs_range = np.arange(   10,  150,   10)
        tw_range = np.arange(    0,   40,    5)
        vs_range = np.arange( 0.78, 1.02, 0.02)
        patm = self.patm
        self.rh = {rh:[humid_air('', patm=patm, A=1, rh=rh, t =t).solve().valuelist(['t', 'v']) for t in t_range] for rh in rh_range}
        self.hs = {hs:[humid_air('', patm=patm, A=1, h =hs, rh=rh).solve().valuelist(['t', 'v']) for rh in rh_range] for hs in hs_range}
        self.tw = {tw:[humid_air('', patm=patm, A=1, tw=tw, rh=rh).solve().valuelist(['t', 'v']) for rh in rh_range] for tw in tw_range}
        self.vs = {vs:[humid_air('', patm=patm, A=1, vs=vs, rh=rh).solve().valuelist(['t', 'v']) for rh in rh_range] for vs in vs_range}
        pickle.dump(self, file)
        return self

    def chart(self, paths=None, rh_alpha=0.5, hs_alpha=0.5, tw_alpha=0.0, vs_alpha=0.0, width=675, height=556, xmax=50, ymax=0.030):
        if paths is None: paths = []
        points = list(set([po for po in [po for pa in paths for po in pa]]))
        rh = hv.Overlay( [hv.Curve(v) for (k,v) in self.rh.items()] )
        hs = hv.Overlay( [hv.Curve(v) for (k,v) in self.hs.items()] )
        tw = hv.Overlay( [hv.Curve(v) for (k,v) in self.tw.items()] )
        vs = hv.Overlay( [hv.Curve(v) for (k,v) in self.vs.items()] )

        txyn = [(point.t.value, point.v.value, point.psyname)   for point in points]
        po = hv.Points(txyn) * hv.Overlay([ hv.Text(x+1,y+0.001,n) for (x,y,n) in txyn ])

        pa = hv.Path([[(point.t.value, point.v.value) for point in path] for path in paths] )

        gr = hv.Overlay(  [hv.Curve([( 0, y),(50,   y)]) for y in np.arange(0, 0.042, 0.002)]+
                            [hv.Curve([( x, 0),( x,0.04)]) for x in np.arange(0,    52,     2)] )

        g = rh({'Curve' : {'style': {'color':'red'   , 'alpha':rh_alpha, 'hover_line_alpha':1}}}) *\
            hs({'Curve' : {'style': {'color':'blue'  , 'alpha':hs_alpha, 'hover_line_alpha':1}}}) *\
            tw({'Curve' : {'style': {'color':'green' , 'alpha':tw_alpha}}}) *\
            vs({'Curve' : {'style': {'color':'grey'  , 'alpha':vs_alpha}}}) *\
            po({'Points': {'style': {'color':'black' , 'size':5 }}})        *\
            pa({'Path'  : {'style': {'color':'black'}}})                    *\
            gr({'Curve' : {'style': {'color':'black' , 'alpha':0.1}}})

        g = g.redim(x ={'range': (0, xmax),'label': 'Dry bulb temperature  [°C]'})
        g = g.redim(y ={'range': (0, ymax),'label': 'Humidity ratio     [kg/kg]'})

        g = g({'Curve':{'style': {'line_width':0}}})
        g = g({'Curve':{'plot':{'tools':['hover'], 'toolbar':'above', 'yaxis':'right'}}})
        g = g({'Overlay':{'plot': {'width': width, 'height':height, 'title_format': 'Psychrometric Chart (' + str(self.patm) + ' Pa)'}} })
        g = g({'Text':{'style': {'text_font_size':'8pt'}}})

        return g
