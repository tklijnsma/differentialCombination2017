#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys, itertools

import ROOT

from math import sqrt, log, exp, pi, asin
from array import array


########################################
# Main
########################################

def main():
    ReadCouplingFile()


    # All widths precalculated, see:
    # https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/74x-root6/data/lhc-hxswg/couplings/Gamma_Hgammagamma.txt

    # Contents of said file:
    # M_H      G_gaga    G_tt/G_gaga  G_bb/G_gaga  G_WW/G_gaga  G_tb/G_gaga  G_tW/G_gaga  G_bW/G_gaga G_ll/G_gaga  G_tl/G_gaga  G_bl/G_gaga  G_lW/G_gaga
    # 123.  0.86609E-05  0.72602E-01  0.20885E-04   1.5934     -0.18391E-02 -0.68026      0.86351E-02 0.25148E-04 -0.19851E-02  0.45838E-04  0.93221E-02
    # 124.  0.89615E-05  0.72048E-01  0.20152E-04   1.5911     -0.18029E-02 -0.67717      0.84914E-02 0.24292E-04 -0.19461E-02  0.44253E-04  0.91672E-02
    # 125.  0.92729E-05  0.71482E-01  0.19444E-04   1.5887     -0.17672E-02 -0.67400      0.83495E-02 0.23464E-04 -0.19077E-02  0.42720E-04  0.90143E-02
    # 126.  0.95959E-05  0.70904E-01  0.18760E-04   1.5863     -0.17319E-02 -0.67074      0.82093E-02 0.22663E-04 -0.18696E-02  0.41239E-04  0.88634E-02
    # 127.  0.99309E-05  0.70312E-01  0.18099E-04   1.5838     -0.16970E-02 -0.66741      0.80708E-02 0.21888E-04 -0.18321E-02  0.39806E-04  0.87144E-02



    print '\n\n' + '='*70
    print 'EXCLUDING THE CHARM'

    manualBRcalculator = ManualBRcalculator(
        # includeCharm = True
        includeCharm = False
        )

    couplingCombinations = GetCombinations( [ 't', 'b', 'W' ]  )

    for mH in [ 123., 124., 125., 126., 127. ]:
        print '\n' + '-'*50 + '\nSetting mH to {0}'.format(mH)
        manualBRcalculator.SetMH( mH )

        for coupling1, coupling2 in couplingCombinations:

            if coupling1 == coupling2:
                A = getattr( manualBRcalculator, 'A'+coupling1 )()
                term = norm( A*A ) / manualBRcalculator.norm_squared_SM
            else:
                A1 = getattr( manualBRcalculator, 'A'+coupling1 )()
                A2 = getattr( manualBRcalculator, 'A'+coupling2 )()

                term = 2.*( A1 * A2 ).real / manualBRcalculator.norm_squared_SM

            termFromFile = GetWidthFromFile( coupling1+coupling2, mH )
            print '{0}{1}: calc: {2:+.5f} , in file: {3:+.5f}'.format( coupling1, coupling2, term, termFromFile )







    print '\n\n' + '='*70
    print 'INCLUDING THE CHARM'

    manualBRcalculator = ManualBRcalculator(
        includeCharm = True
        # includeCharm = False
        )

    couplingCombinations = GetCombinations( [ 't', 'b', 'W', 'c' ]  )

    for mH in [ 123., 124., 125., 126., 127. ]:
        print '\n' + '-'*50 + '\nSetting mH to {0}'.format(mH)
        manualBRcalculator.SetMH( mH )

        for coupling1, coupling2 in couplingCombinations:

            if coupling1 == coupling2:
                A = getattr( manualBRcalculator, 'A'+coupling1 )()
                term = norm( A*A ) / manualBRcalculator.norm_squared_SM
            else:
                A1 = getattr( manualBRcalculator, 'A'+coupling1 )()
                A2 = getattr( manualBRcalculator, 'A'+coupling2 )()

                term = 2.*( A1 * A2 ).real / manualBRcalculator.norm_squared_SM

            termFromFile = GetWidthFromFile( coupling1+coupling2, mH )
            print '{0}{1}: calc: {2:+.5f} , in file: {3:+.5f}'.format( coupling1, coupling2, term, termFromFile )




    sys.exit()


    print '\n\n' + '='*70
    print 'INCLUDING THE CHARM'

    for mH in [ 123., 124., 125., 126., 127. ]:

        print '\n' + '-'*50 + '\nSetting mH to {0}'.format(mH)
        manualBRcalculator.SetMH( mH )

        print 'total width = {0}'.format( manualBRcalculator.norm_squared_SM )
        print 'tt width = {0}'.format( manualBRcalculator.evaluateWidth('t') )
        print 'bb width = {0}'.format( manualBRcalculator.evaluateWidth('b') )
        print 'WW width = {0}'.format( manualBRcalculator.evaluateWidth('w') )

        if manualBRcalculator.includeCharm:
            print 'cc width = {0}'.format( manualBRcalculator.evaluateWidth('c') )

        print ''
        print 'tW interference width = {0}'.format( 2.*( manualBRcalculator.At() * manualBRcalculator.AW() ).real / manualBRcalculator.norm_squared_SM )
        print 'bW interference width = {0}'.format( 2.*( manualBRcalculator.Ab() * manualBRcalculator.AW() ).real / manualBRcalculator.norm_squared_SM )
        print 'cW interference width = {0}'.format( 2.*( manualBRcalculator.Ac() * manualBRcalculator.AW() ).real / manualBRcalculator.norm_squared_SM )


        print '\n( 2 Re( AW * Ab ) ) / ( 2 Re( AW * Ac ) ) = {0}'.format(
            ( 2.*( manualBRcalculator.Ab() * manualBRcalculator.AW() ).real )
            /
            ( 2.*( manualBRcalculator.Ac() * manualBRcalculator.AW() ).real )
            )


    # # ======================================
    # # Old evaluation method
   
    # mH        = ROOT.RooRealVar( 'mH', 'mH', 1.0  )
    # kappa_t   = ROOT.RooRealVar( 'kappa_t', 'kappa_t', 1.0  )
    # kappa_W   = ROOT.RooRealVar( 'kappa_W', 'kappa_W', 1.0  )
    # mb        = ROOT.RooRealVar( 'mb', 'mb', 1.0  )
    # kappa_b   = ROOT.RooRealVar( 'kappa_b', 'kappa_b', 1.0  )

    # Scaling_hgg_old = ROOT.RooScaleHGamGamLOSM(
    #     'Scaling_hgg_old',
    #     'Scaling_hgg_old',
    #     mH, kappa_t, kappa_W, mb, kappa_b
    #     )

    # print Scaling_hgg_old.getVal()




class ManualBRcalculator():
    def __init__( self, includeCharm=False ):
        
        self.includeCharm = includeCharm

        # self.mt = 172.5
        # self.mb = 5
        # self.mc = 1.25

        self.mW = 80.4

        self.mH = 125.


        self.kappa_t = 1.
        self.kappa_b = 1.
        self.kappa_c = 1.
        self.kappa_W = 1.



        # ======================================
        # Higgs mass as RooRealVar for the spline evaluation

        self.mH_roo = ROOT.RooRealVar( 'mH', 'mH', self.mH  )
        ROOT.SetOwnership( self.mH_roo, False )


        # ======================================
        # Running mb

        self.mbSpline = textToSpline(
            name     = 'mb',
            filename = os.path.join(
                    os.environ['CMSSW_BASE'],
                    'src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg',
                    # 'running_constants.txt'
                    'running_fermion_masses.txt'
                    ),
            xv       = self.mH_roo,
            ycol     = 2
            )
        ROOT.SetOwnership( self.mbSpline, False )


        # ======================================
        # Running mc

        self.mcSpline = textToSpline(
            name     = 'mc',
            filename = os.path.join(
                    os.environ['CMSSW_BASE'],
                    'src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg',
                    # 'running_constants.txt'
                    'running_fermion_masses.txt'
                    ),
            xv       = self.mH_roo,
            ycol     = 1
            )
        ROOT.SetOwnership( self.mcSpline, False )


        # ======================================
        # Running mt

        self.mtSpline = textToSpline(
            name     = 'mt',
            filename = os.path.join(
                    os.environ['CMSSW_BASE'],
                    'src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg',
                    # 'running_constants.txt'
                    'running_fermion_masses.txt'
                    ),
            xv       = self.mH_roo,
            ycol     = 3
            )
        ROOT.SetOwnership( self.mtSpline, False )



        self.mbRunning = True
        self.mtRunning = True
        self.mcRunning = True

        self.norm_squared_SM = self.evaluate( SM_evaluation = True )




    # ======================================
    # Quark masses as functions


    def SetMH( self, mH ):
        self.mH = mH
        self.mH_roo.setVal( self.mH )

        kappa_t_SAVED = self.kappa_t
        kappa_b_SAVED = self.kappa_b
        kappa_c_SAVED = self.kappa_c
        kappa_W_SAVED = self.kappa_W

        self.kappa_t = 1.
        self.kappa_b = 1.
        self.kappa_c = 1.
        self.kappa_W = 1.

        self.norm_squared_SM = self.evaluate( SM_evaluation = True )

        self.kappa_t = kappa_t_SAVED
        self.kappa_b = kappa_b_SAVED
        self.kappa_c = kappa_c_SAVED
        self.kappa_W = kappa_W_SAVED



    def mb( self ):
        if self.mbRunning:
            self.mH_roo.setVal( self.mH )
            return self.mbSpline.getVal()
        else:
            return 4.18

    def mt( self ):
        if self.mtRunning:
            self.mH_roo.setVal( self.mH )
            return self.mtSpline.getVal()
        else:       
            return 172.5

    def mc( self ):
        if self.mcRunning:
            self.mH_roo.setVal( self.mH )
            return self.mcSpline.getVal()
        else:
            return 1.25



    def evaluate( self, SM_evaluation = False ):

        Nc = 3.

        ret = self.At() + self.Ab() + self.AW()

        if self.includeCharm:
            ret += self.Ac()

        ret = norm_squared( ret )

        if SM_evaluation:
            return ret
        else:
            return ret / self.norm_squared_SM


    def evaluateWidth( self, particle = 't' ):
        particle = particle.lower()
        if not particle in [ 't', 'b', 'w', 'c' ]:
            raise RuntimeError( 'Particle \'{0}\' is not implemented'.format(particle) )

        Nc = 3.
        if particle == 't':
            return norm_squared( self.At() ) / self.norm_squared_SM
        if particle == 'b':
            return norm_squared( self.Ab() ) / self.norm_squared_SM
        if particle == 'w':
            return norm_squared( self.AW() ) / self.norm_squared_SM

        if particle == 'c':
            if not self.includeCharm:
                print 'Warning: charm was not included in SM width calculation but width was requested'
            return norm_squared( self.Ac() ) / self.norm_squared_SM


    def At( self ):
        return self.kappa_t * 3 * (2./3.)**2 * F_half( self.mt(), self.mH )

    def Ab( self ):
        return self.kappa_b * 3 * (-1./3.)**2 * F_half( self.mb(), self.mH )

    def AW( self ):
        return self.kappa_W * F_one( self.mW, self.mH )

    def Ac( self ):
        return self.kappa_c * 3 * (2./3.)**2 * F_half( self.mc(), self.mH )




def F_one( m, mH = 125. ):
    tau = 4* (m**2) / (mH**2)
    return 2. + 3.*tau + 3.*tau*(2.-tau)*f(tau)

def F_half( m, mH = 125. ):
    tau = 4.* (m**2) / (mH**2)
    return -2.*tau*( 1. + (1.-tau)*f(tau) )


def f( tau ):

    if tau >= 1.:
        ret = ( asin( sqrt( 1./tau ) ) )**2
        return ret

    else:

        temp = ( 1. + sqrt(1.-tau) ) / ( 1. - sqrt(1.-tau) )

        ret = log( temp ) + 1j * pi
        ret = -0.25 * ret**2

        return ret


def norm_squared( c ):
    return c.real**2 + c.imag**2

def norm( c ):
    return sqrt( norm_squared(c) )


def textToSpline(
    name,
    filename,
    xv,
    ycol     = 1,
    xcol     = 0,
    algo     = "CSPLINE"
    ):

    x = []
    y = []

    with open( filename, 'r' ) as fP:
        lines = [ l for l in fP.readlines() ]

    for line in lines[1:]:
        if len(line.strip()) == 0: continue
        cols = line.split();
        x.append(float(cols[xcol]))
        y.append(float(cols[ycol]))

    spline = ROOT.RooSpline1D(
        name,
        "file %s, x=%d, y=%d" % (filename,xcol,ycol),
        xv,
        len(x),
        array('d', x),
        array('d', y),
        algo
        )

    return spline



def GetCombinations(
        couplings = [ 't', 'b', 'W', 'l' ]
        ):

    couplingCombinations = []
    couplingCombinations.extend( [ [ coupling, coupling ] for coupling in couplings ] )
    couplingCombinations.extend( [ list(couplingTuple) for couplingTuple in itertools.combinations( couplings, 2 ) ] )

    return couplingCombinations



# ======================================
# 


GetWidthFromFile = 0

def ReadCouplingFile():

    filename = os.path.join(
        os.environ['CMSSW_BASE'],
        'src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/couplings/Gamma_Hgammagamma.txt',
        )

    with open( filename, 'r' ) as fP:
        lines = fP.readlines()
    lines = [ l.split() for l in lines ]

    headers = lines.pop(0)
    couplings = [ 't', 'b', 'W', 'l' ]
    couplingCombinations = []
    couplingCombinations.extend( [ [ coupling, coupling ] for coupling in couplings ] )
    couplingCombinations.extend( [ list(couplingTuple) for couplingTuple in itertools.combinations( couplings, 2 ) ] )

    columnDict = {}
    for coupling1, coupling2 in couplingCombinations:
        for iHeader, header in enumerate(headers):
            if coupling1+coupling2 in header or coupling2+coupling1 in header:
                columnDict[ coupling1+coupling2 ] = iHeader
                columnDict[ coupling2+coupling1 ] = iHeader
                break
        else:
            raise RuntimeError(
                'Could not determine a column for {0}{1}'.format( coupling1, coupling2 )
                )


    mHlist = [ float(l[0]) for l in lines ]


    def reader( channel, mH ):

        if not channel in columnDict:
            # raise RuntimeError(
            #     'Channel \'{0}\' not in columnDict'.format(channel)
            #     )
            return 0.0
        iCol = columnDict[channel]


        if not mH in mHlist:
            raise RuntimeError(
                'Values in file not known for mH = {0}'.format(mH)
                )
        iRow = mHlist.index( mH )

        return float(lines[iRow][iCol])

    global GetWidthFromFile
    GetWidthFromFile = reader















########################################
# End of Main
########################################
if __name__ == "__main__":
    main()