from psychropyo.psychroutils import *
from pyomo.environ import Block

def mix(psyname, ha, hb, hc):
    b = Block()
    b.streams = [ha,hb,hc]
    b.psyname = psyname
    b.psytype1 = "process"
    b.psytype2 = "mix"
    b.airbal   = myCon(ha.A + hb.A == hc.A)
    b.waterbal = myCon(ha.V + hb.V + ha.W + hb.W == hc.V + hc.W)
    b.heatbal  = myCon(ha.H + hb.H == hc.H)
    return b

def condense(psyname, ha, hb, hc):
    b = Block()
    b.streams = [ha,hb,hc]
    b.psyname = psyname
    b.psytype1 = "process"
    b.psytype2 = "condense"
    b.mix = mix(psyname + "_mix", ha, hb, hc)
    b.maxrh    = myCon(hc.rh <= 1)
    b.thermEqu = myCon(hb.t == hc.t)
    b.reactEqu = myCon(hb.W * (hc.rh - 1) == 0)     # ;>)  equilibrium if (no water) or (saturation)
    return b

def split(psyname, ha, hb, hc, **fixedVars):
    b = Block()
    b.streams = [ha,hb,hc]
    b.psyname = psyname
    b.psytype1 = "process"
    b.psytype2 = "split"
    b.xb = myVar(0   , (0,1), ".")
    b.xc = myVar(0   , (0,1), ".")
    b.xbxc     = myCon(b.xb + b.xc == 1)

    b.airb     = myCon(hb.A  == ha.A  * b.xb)
    b.waterb   = myCon(hb.W  == ha.W  * b.xb)
    b.vaporb   = myCon(hb.V  == ha.V  * b.xb)
    b.heatb    = myCon(hb.H  == ha.H  * b.xb)

    b.airc     = myCon(hc.A  == ha.A  * b.xc)
    b.waterc   = myCon(hc.W  == ha.W  * b.xc)
    b.vaporc   = myCon(hc.V  == ha.V  * b.xc)
    b.heatc    = myCon(hc.H  == ha.H  * b.xc)

    for (var, val) in fixedVars.items(): setattr(b, var + "_fixed", myCon(getattr(b, var) == val))
    return b

def walllosses(psyname, ha, hb, hc, hA, tout):
    b = Block()
    b.streams = [ha,hb,hc]
    b.psyname = psyname
    b.psytype1 = "process"
    b.psytype2 = "walllosses"
    b.exch = myCon(hb.H == hA * (tout - hc.t))
    b.mix = mix(psyname+'mix', ha, hb, hc)
    return b

def countercurrent(psyname, ha1, ha2, hb1, hb2, eff=1):
    # ok for ideal eff=1, simplistic otherwise
    b = Block()
    b.streams = [ha1, ha2, hb1, hb2]
    b.psyname = psyname
    b.psytype1 = "process"
    b.psytype2 = "countercurrent"
    b.haA   = myCon(ha1.A == ha2.A)
    b.haW   = myCon(ha1.W == ha2.W)
    b.haV   = myCon(ha1.V == ha2.V)
    b.hbA   = myCon(hb1.A == hb2.A)
    b.hbW   = myCon(hb1.W == hb2.W)
    b.hbV   = myCon(hb1.V == hb2.V)
    #b.exch  = myCon((ha1.t - hb2.t)*(ha2.t - hb1.t) == 0)
    b.exch  = myCon((ha2.t - ha1.t - eff*(hb1.t - ha1.t))*(hb2.t - hb1.t - eff*(ha1.t - hb1.t)) == 0)
    b.hbal  = myCon(ha1.H + hb1.H == ha2.H + hb2.H)
    return b

def flowsheetsection(psyname, *devices):
    b = Block()
    b.psyname = psyname
    b.psytype1 = "process"
    b.psytype2 = "flowsheetsection"

    for device in devices: setattr(b, device.psyname, device)
    streams = [stream for device in devices for stream in device.streams ]

    externalstreams = set([s for s in streams if streams.count(s) == 1])
    internalstreams = set([s for s in streams if streams.count(s) != 1])
    for stream in internalstreams: setattr(b, stream.psyname, stream)
    b.streams = externalstreams
    b.internalstreams = internalstreams
    return b

def flowsheet(psyname, *devices):
    b = Block()
    b.streams = set([])
    b.psyname = psyname
    b.psytype1 = "process"
    b.psytype2 = "flowsheet"

    for device in devices: setattr(b, device.psyname, device)
    streams = set([stream for device in devices for stream in device.streams ])
    for stream in streams: setattr(b, stream.psyname, stream)
    b.streams = streams
    return b

