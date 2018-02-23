from HiggsAnalysis.CombinedLimit.PhysicsModel import *

from array import array
from time import sleep

from MethodHandler import MethodHandler
methodhandler = MethodHandler([
    'physicsModels.CouplingModelMethods.YieldParameters',
    'physicsModels.CouplingModelMethods.Parametrization',
    'physicsModels.CouplingModelMethods.CouplingDependentBRs',
    'physicsModels.CouplingModelMethods.TheoryUncertainties',
    'physicsModels.CouplingModelMethods.OtherStudies',
    ])

class CouplingModel( PhysicsModel ):
    ''' Model used to unfold differential distributions '''

    class CouplingModelError(Exception):
        pass

    def __init__(self):
        PhysicsModel.__init__(self)
        self.mHRange=[]
        self.mass = 125.
        self.verbose = 1

        self.Persistence = []
        self.theories = []
        self.SMbinBoundaries = []
        self.SMXS = []

        self.isOnlyHgg = False
        self.isOnlyHZZ = False

        self.includeLinearTerms   = False
        self.splitggH             = False
        self.MakeLumiScalable     = False
        self.FitBR                = False
        self.ProfileTotalXS       = False
        self.FitOnlyNormalization = False
        self.FitRatioOfBRs        = False
        self.DoBRUncertainties    = False
        self.ReweightedXS         = None

        self.theoryUncertaintiesPassed = False
        self.correlationMatrixPassed   = False
        self.covarianceMatrixPassed    = False
        self.skipOverflowBinTheoryUncertainty = False

        self.manualExpBinBoundaries = []
        self.skipBins = []
        self.skipOverflowBin = False

        self.SM_HIGG_DECAYS = [ "hww", "hzz", "hgg", "htt", "hbb", 'hzg', 'hmm', 'hcc', 'hgluglu' ]
        self.SM_HIGG_DECAYS.append( 'hss' )


    def setPhysicsOptions( self, physOptions ):
        self.chapter( 'Starting model.setPhysicsOptions()' )

        for physOption in physOptions:
            optionName, optionValue = physOption.split('=',1)
            if False: pass
            elif optionName == 'splitggH':
                self.splitggH = eval(optionValue)
            elif optionName == 'isOnlyHZZ':
                self.isOnlyHZZ = eval(optionValue)
            elif optionName == 'isOnlyHgg':
                self.isOnlyHgg = eval(optionValue)
            elif optionName == 'binBoundaries':
                self.manualExpBinBoundaries = [ float(i) for i in optionValue.split(',') ]
            elif optionName == 'verbose':
                self.verbose = int(optionValue)
            elif optionName == 'linearTerms':
                self.includeLinearTerms = eval(optionValue)
            elif optionName == 'higgsMassRange':
                self.mHRange = [ float(v) for v in optionValue.split(',') ]

            elif optionName == 'theoryUncertainties':
                self.theoryUncertainties = self.readErrorFile( optionValue )
                self.theoryUncertaintiesPassed = True
            elif optionName == 'correlationMatrix':
                self.correlationMatrix = self.readCorrelationMatrixFile( optionValue )
                self.correlationMatrixPassed = True
            elif optionName == 'covarianceMatrix':
                self.covarianceMatrix = self.readCorrelationMatrixFile( optionValue )
                self.covarianceMatrixPassed = True

            elif optionName == 'lumiScale':
                self.MakeLumiScalable = eval(optionValue)
            elif optionName == 'ReweightCrossSections':
                self.ReweightedXS = [ float(v) for v in optionValue.split(',') ]
            elif optionName == 'FitBR':
                self.FitBR = eval(optionValue)
            elif optionName == 'ProfileTotalXS':
                self.ProfileTotalXS = eval(optionValue)
            elif optionName == 'FitOnlyNormalization':
                self.FitOnlyNormalization = eval(optionValue)
            elif optionName == 'FitRatioOfBRs':
                self.FitRatioOfBRs = eval(optionValue)
            elif optionName == 'DoBRUncertainties':
                self.DoBRUncertainties = eval(optionValue)

            elif optionName == 'theory':
                # Syntax: --PO theory=[ct=1,cg=1,file=some/path/.../] , brackets optional
                # Delete enclosing brackets if present
                if optionValue.startswith('['): optionValue = optionValue[1:]
                if optionValue.endswith(']'): optionValue = optionValue[:-1]
                theoryDict = {                
                    'couplings'      : {},
                    'binBoundaries'  : [],
                    'crosssection'   : [],
                    'ratios'         : [],
                    }
                for component in optionValue.split(','):
                    couplingName, couplingValue = component.split('=',1)
                    # Only exception is file, which is the file to read the shape from
                    if couplingName == 'file':
                        theoryFile = couplingValue
                        theoryDict['binBoundaries'], theoryDict['ratios'], theoryDict['crosssection'] = self.readTheoryFile(theoryFile)
                        continue
                    theoryDict['couplings'][couplingName] = float(couplingValue)
                self.theories.append( theoryDict )
            elif optionName == 'SM':
                # Delete enclosing brackets if present
                if optionValue.startswith('['): optionValue = optionValue[1:]
                if optionValue.endswith(']'): optionValue = optionValue[:-1]
                self.SMDict = {                
                    'couplings'      : {},
                    'binBoundaries'  : [],
                    'crosssection'   : [],
                    }
                for component in optionValue.split(','):
                    couplingName, couplingValue = component.split('=',1)
                    # Only exception is file, which is the file to read the shape from
                    if couplingName == 'file':
                        self.SMDict['binBoundaries'], dummy, self.SMDict['crosssection'] = self.readTheoryFile( couplingValue )
                        continue
                    self.SMDict['couplings'][couplingName] = float(couplingValue)
                # SM cross section extended with underflow and overflow bin
                # self.SMXS = [ self.SMDict['crosssection'][0] ] + self.SMDict['crosssection'] + [ self.SMDict['crosssection'][-1] ]
                self.SMXS = [ 0.0 ] + self.SMDict['crosssection']

            else:
                print 'Unknown option \'{0}\'; skipping'.format(optionName)
                continue

        # Minor checks before continuing
        if self.FitBR and self.FitRatioOfBRs:
            raise self.CouplingModelError( 'Options FitBR and FitRatioOfBRs are exclusive!' )
        if not(self.ReweightedXS is None):
            if not( len(self.ReweightedXS) == len(self.manualExpBinBoundaries)-1 ):
                raise self.CouplingModelError('Dim mismatch: len(ReweightCrossSections)={0}, but len(binBoundaries)-1={1}'.format(len(self.ReweightedXS), len(self.manualExpBinBoundaries)-1 ))


    def doParametersOfInterest(self):
        self.chapter('Starting model.doParametersOfInterest()')
        self.setMH()

        if self.FitOnlyNormalization:
            self.makeParametrizationsFromTheory()
            self.MakeTotalXSExpressions()
            self.modelBuilder.doVar('bkg_modifier[1.0]')
            self.modelBuilder.out.var('bkg_modifier').setConstant(True)
            self.modelBuilder.doVar('hgg_xH_modifier[1.0]')
            self.modelBuilder.out.var('hgg_xH_modifier').setConstant(True)
            self.modelBuilder.doVar('hzz_xH_modifier[1.0]')
            self.modelBuilder.out.var('hzz_xH_modifier').setConstant(True)
            self.modelBuilder.out.defineSet( 'POI', ','.join(self.couplings) )
        else:
            self.figureOutBinning()
            self.makeParametrizationsFromTheory()
            if self.FitBR:
                self.MakeWidthExpressions()
            if self.ProfileTotalXS:
                self.MakeTotalXSExpressions()
            if self.FitRatioOfBRs:
                self.MakeRatioOfBRsExpressions()
            self.defineYieldParameters()
            self.addTheoryUncertaintyNuisances()

        self.chapter( 'Starting model.getYieldScale()' )

    def figureOutBinning( self ):
        bins    = self.DC.bins
        signals = self.DC.signals
        signals = [ s for s in signals if not 'OutsideAcceptance' in s ]
        productionModes = list(set([ s.split('_')[0] for s in signals ]))
        if self.splitggH: signals = [ s for s in signals if 'ggH' in s ]
        signals.sort( key = lambda n: float(
            n.split('_')[2].replace('p','.').replace('m','-').replace('GT','').replace('GE','').replace('LT','').replace('LE','') ) )
        productionMode, observableName = signals[0].split('_')[:2]

        if len(self.manualExpBinBoundaries) == 0:
            raise self.CouplingModelError(
                'Specify \'--PO binBoundaries=0,15,...\' on the command line.'
                'Automatic boundary determination no longer supported.'
                )
        expBinBoundaries = self.manualExpBinBoundaries
        nExpBins = len(expBinBoundaries)-1
        if self.verbose: print 'Taking manually specified bin boundaries:', expBinBoundaries

        self.allProcessBinBoundaries = []
        for signal in signals:
            if 'OutsideAcceptance' in signal: continue
            for bound in signal.split('_')[2:]:
                bound = bound.replace('p','.').replace('m','-').replace('GT','').replace('GE','').replace('LT','').replace('LE','')
                self.allProcessBinBoundaries.append( float(bound) )

        # Attach to class so they are accesible in other methods
        self.expBinBoundaries = expBinBoundaries
        self.nExpBins         = nExpBins
        self.signals          = signals
        self.productionMode   = productionMode
        self.observableName   = observableName

        expBinBoundarySet = []
        for i, expBinBoundary in enumerate(self.expBinBoundaries):
            if self.skipOverflowBin and expBinBoundary == expBinBoundaries[-1]: continue
            self.modelBuilder.doVar( 'expBinBound{0}[{1}]'.format( i, expBinBoundary ) )
            expBinBoundarySet.append( 'expBinBound{0}'.format(i) )
        self.modelBuilder.out.defineSet( 'expBinBoundaries', ','.join(expBinBoundarySet) )


    def setMH( self ):
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                print 'MH will be assumed to be', self.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.mass)
        else:
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
            else:
                print 'MH (not there before) will be assumed to be', self.mass
                self.modelBuilder.doVar("MH[%g]" % self.mass)

    def chapter( self, text, indent=0 ):
        if self.verbose:
            print '\n{tabs}{line}\n{tabs}{text}\n'.format(
                tabs = '    ' * indent,
                line = '-' * 70,
                text = text
                )

    def BinIsHgg( self, bin ):
        # If the input is hardcoded to be only hzz or hgg:
        if self.isOnlyHgg:
            return True
        if self.isOnlyHZZ:
            return False
        
        # This is not very robust coding
        if bin.startswith('hgg_'):
            # if self.verbose: print '    Bin \'{0}\' is classified as a hgg bin!'.format( bin )
            return True
        elif bin.startswith('hzz_'):
            # if self.verbose: print '    Bin \'{0}\' is classified as a hzz bin!'.format( bin )
            return False
        else:
            raise self.CouplingModelError( 'Bin \'{0}\' can not be classified as either a hgg or hzz bin'.format(bin) )

#____________________________________________________________________
# Finalize
methodhandler.make_methods(CouplingModel)
couplingModel = CouplingModel()