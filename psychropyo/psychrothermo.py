from pyomo.environ import sqrt, log, exp

Rgas = 8.3144598     # J/mol/K
pref = 101325        # Pa
Tref = 273.15        # K

class air:
    '''
    Thermodynamic data of water from US Bureau of Mines
        for standard pressure (101325 Pa)'''
    M = 28.8503972
    def cp0(T):
        '''Specific heat, J/mol/K  '''
        return 4.184 * (9.9044 + 126531./(T*T) - 77.24/sqrt(T) + 0.00042677*T - 9.257100000000002e-8*T*T)
    def H0(T):
        '''Enthalpy, J/mol  '''
        return 4.184 * (0 + 9.9044*T + 0.000213385*T*T -126531.0/T -154.48*sqrt(T) -3.0857e-08*T*T*T + 121.0)
    def S0(T):
        '''Entropy, J/mol/K  '''
        return 4.184*(-16.28923809866967-126531.0/(2*T*T)+154.48/sqrt(T)+2*0.000213385*T+(3*-3.0857e-08*T*T)/2+9.9044*log(T))
    def G0(T):
        '''Gibbs free energy, J/mol  '''
        return air.H0(T) - T * air.S0(T)

class H2O:
    '''
    Thermodynamic data of water
        H, S, G0, cp from US Bureau of Mines
        psat from IAPWS'''
    M = 18.01528
    class g:
        '''Vapor phase for standard pressure (101325 Pa)'''
        def cp0(T):
            '''Specific heat, J/mol/K  '''
            return 4.184 * (4.7301998 + 15523./(T*T) + 29.5200005/sqrt(T) + 0.0049613*T - 7.72539e-7*T*T)
        def H0(T):
            '''Enthalpy, J/mol  '''
            return 4.184 * (-57795.0 + 4.7301998*T + 0.00248065*T*T -15523.0/T + 59.040001*sqrt(T) -2.57513e-07*T*T*T -2591.0)
        def S0(T):
            '''Entropy, J/mol/K  '''
            return 4.184*(20.21690179199231-15523.0/(2*T*T)-59.040001/sqrt(T)+2*0.00248065*T+(3*-2.57513e-07*T*T)/2+4.7301998*log(T))
        def G0(T):
            '''Gibbs free energy, J/mol  '''
            return H2O.g.H0(T) - T * H2O.g.S0(T)
    class l:
        '''Liquid phase for standard pressure (101325 Pa)'''
        def cp0(T):
            '''Specific heat, J/mol/K  '''
            return 4.184 * (558.898 + 3.85798e6/(T*T) - 7248.38/sqrt(T) - 0.6867*T + 0.00045269099999999997*T*T)
        def H0(T):
            '''Enthalpy, J/mol  '''
            return 4.184 * (-68315.0 + 558.898*T -0.34335*T*T -3857980.0/T -14496.76*sqrt(T) + 0.000150897*T*T*T + 123142.0)
        def S0(T):
            '''Entropy, J/mol/K  '''
            return 4.184*(-3800.901885684364-3857980.0/(2*T*T)+14496.76/sqrt(T)+2*-0.34335*T+(3*0.000150897*T*T)/2+558.898*log(T))
        def G0(T):
            '''Gibbs free energy, J/mol  '''
            return H2O.l.H0(T) - T * H2O.l.S0(T)
    def L(T):
        '''Latent heat of vaporization'''
        return H2O.g.H0(T)-H2O.l.H0(T)
    def psat(T):
        '''
        Saturation pressure
            IAPWS R7-97(2012) equation 30
            http://www.iapws.org/relguide/IF97-Rev.pdf
            https://github.com/jjgomera/iapws/blob/master/iapws/iapws97.py
            https://en.wikipedia.org/wiki/Vapour_pressure_of_water'''
        ni = [1167.0521452767,-724213.16703206,-17.073846940092,12020.82470247,-3232555.0322333,14.91510861353,-4823.2657361591,405113.40542057,-0.23855557567849, 650.17534844798]
        v = T + ni[8] / (T - ni[9])
        A = 1 * v ** 2 + ni[0] * v + ni[1]
        B = ni[2] * v ** 2 + ni[3] * v + ni[4]
        C = ni[5] * v ** 2 + ni[6] * v + ni[7]
        return 1000000*(2 * C / (-B + (B*B - 4 * A * C) ** 0.5)) ** 4
    def psatIdeal(T):
        '''Saturation pressure from ideal gas assumption'''
        return pref * exp((H2O.l.G0(T)-H2O.g.G0(T))/Rgas/T)
    def psatBuck(T):
        '''Saturation pressure from Buck equation'''
        return 611.21 * exp( (19.84282 - T/234.5) * ((T - 273.15)/(-16.01 + T)) )