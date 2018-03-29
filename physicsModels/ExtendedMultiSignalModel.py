import os
from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder

class ExtendedMultiSignalModel(MultiSignalModel):
    """docstring for ExtendedMultiSignalModel"""
    def __init__(self):
        MultiSignalModel.__init__(self)
        self.exprs = []
        self.get_splines = False
        
    def setPhysicsOptions(self, physOptions):
        for po in physOptions:
            if po.startswith("get_splines"):
                self.get_splines = True
            if po.startswith("factory="):
                self.exprs.append(po.replace('factory=',''))
        MultiSignalModel.setPhysicsOptions(self, physOptions)

    def make_splines(self):
        print 'Making splines'
        # Load spline for HZZ BR (function of MH)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        datadir = os.environ['CMSSW_BASE'] + '/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg'
        self.SMH.textToSpline( 'BR_hzz', os.path.join( datadir, 'sm/br/BR4.txt' ), ycol=11 );
        self.spline_SMBR_hzz = self.modelBuilder.out.function('BR_hzz')

        # Seems to work:
        # self.modelBuilder.out.var('MH').setVal(125.09)
        # spline_SMBR_hzz.Print()

        # Load spling for Hgg BR
        self.spline_SMBR_hgg = self.modelBuilder.out.function('fbr_13TeV')

    def doParametersOfInterest(self):  
        if self.get_splines:
            self.make_splines()
        for expr in self.exprs:
            expr = expr.replace('#spline_hgg', self.spline_SMBR_hgg.GetName())
            expr = expr.replace('#spline_hzz', self.spline_SMBR_hzz.GetName())
            print 'Processing factory {0}'.format(expr)
            self.modelBuilder.factory_(expr)
        MultiSignalModel.doParametersOfInterest(self)

extendedMultiSignalModel = ExtendedMultiSignalModel()
