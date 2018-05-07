import os, re, sys
from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder


class Variable():
    def __init__(self, factory_str):
        self.factory_str = factory_str
        match_without_ranges = re.search(r'(.*)\[([\d\.\-]+)\]', factory_str)
        if match_without_ranges:
            match = match_without_ranges
            self.name = match.group(1)
            self.val = match.group(2)
            self.x_min = None
            self.x_max = None
        else:
            match_with_ranges = re.search(r'(.*)\[([\d\.\-]+),([\d\.\-]+),([\d\.\-]+)\]', factory_str)
            if match_with_ranges:
                match = match_with_ranges
                self.name = match.group(1)
                self.val = match.group(2)
                self.x_min = match.group(3)
                self.x_max = match.group(4)

class Expression():
    def __init__(self, factory_str):
        self.factory_str = factory_str
        match = re.search(r'\w+::(.*?)\(', factory_str)
        if not match:
            print 'ERROR: could not determine name from expr str \'{0}\''.format(factory_str)
            self.name = 'None'
        else:
            self.name = match.group(1)
        

class ExtendedMultiSignalModel(PhysicsModel):

    def __init__(self):
        self.mHRange = []
        self.poiMap  = []
        self.get_splines = False

        self.pois = []
        self.variables = []
        self.exprs = []

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:

            if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"

            if po.startswith("get_splines"):
                self.get_splines = True

            if po.startswith("verbose"):
                self.verbose = True

            if po.startswith("map="):
                (maplist,poi) = po.replace("map=","").split(":",1)
                maps = maplist.split(",")
                poiname = re.sub("\[.*","", poi)
                if self.verbose:  print "Mapping ",poiname," to ",maps," patterns"
                self.poiMap.append((poiname, maps))

            if po.startswith("expr="):
                self.exprs.append(
                    Expression(po.replace('expr=',''))
                    )

            if po.startswith('poi=') or po.startswith('variable='):
                variable = Variable(po.replace('poi=', '').replace('variable=', ''))
                self.variables.append(variable)
                if po.startswith('poi='):
                    self.pois.append(variable)
 
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

    def process_mass(self):
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                print 'MH will be assumed to be', self.options.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
        else:
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
            else:
                print 'MH (not there before) will be assumed to be', self.options.mass
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)

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
        """Create POI and other parameters, and define the POI set."""

        if self.get_splines:
            self.make_splines()

        # poiNames = []
        # # first do all non-factory statements, so all params are defined
        # for pn,pf in self.pois.items():
        #     poiNames.append(pn)
        #     self.modelBuilder.doVar(pf)
        # # then do all factory statements (so vars are already defined)
        # for pf in self.factories:
        #     self.modelBuilder.factory_(pf)

        for variable in self.variables:
            print '\nCreating variable {0} (factory str: {1})'.format(variable.name, variable.factory_str)
            self.modelBuilder.doVar(variable.factory_str)
            print 'Test evaluation:'
            self.modelBuilder.out.var(variable.name).Print()

        for expr in self.exprs:
            factory_str = expr.factory_str
            factory_str = factory_str.replace('#spline_hgg', self.spline_SMBR_hgg.GetName())
            factory_str = factory_str.replace('#spline_hzz', self.spline_SMBR_hzz.GetName())
            print '\nCreating expr {0} (factory str: {1})'.format(expr.name, expr.factory_str)
            self.modelBuilder.factory_(factory_str)
            print 'Test evaluation:'
            self.modelBuilder.out.function(expr.name).Print()

        poi_names = ['MH']
        for poi in self.pois:
            poi_names.append(poi.name)
        print 'The following variables are considered POIs:', ', '.join(poi_names)
        self.modelBuilder.doSet("POI",",".join(poi_names))

        self.process_mass()


    def getYieldScale(self,bin,process):
        string = "%s/%s" % (bin,process)
        poi = 1
        for p, maps in self.poiMap:
            for map_str in maps:
                if re.match(map_str, string): poi = p
        print "Will scale ", string, " by ", poi
        if poi in ["1","0"]: return int(poi)
        return poi;


# class ExtendedMultiSignalModel(MultiSignalModel):
#     """docstring for ExtendedMultiSignalModel"""
#     def __init__(self):
#         MultiSignalModel.__init__(self)
#         self.variables = []
#         self.exprs = []
#         self.extended_pois = []
#         self.get_splines = False

#     def setPhysicsOptions(self, physOptions):
#         for po in physOptions:
#             if po.startswith("get_splines"):
#                 self.get_splines = True
#             if po.startswith("factory="):
#                 self.exprs.append(po.replace('factory=',''))
#             if po.startswith('poi=') or po.startswith('variable='):
#                 variable = Variable(po.replace('poi=', ''))
#                 self.variables.append(variable)
#                 if po.startswith('poi='):
#                     self.extended_pois.append(variable)
#         MultiSignalModel.setPhysicsOptions(self, physOptions)

#     def make_splines(self):
#         print 'Making splines'
#         # Load spline for HZZ BR (function of MH)
#         self.SMH = SMHiggsBuilder(self.modelBuilder)
#         datadir = os.environ['CMSSW_BASE'] + '/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg'
#         self.SMH.textToSpline( 'BR_hzz', os.path.join( datadir, 'sm/br/BR4.txt' ), ycol=11 );
#         self.spline_SMBR_hzz = self.modelBuilder.out.function('BR_hzz')

#         # Seems to work:
#         # self.modelBuilder.out.var('MH').setVal(125.09)
#         # spline_SMBR_hzz.Print()

#         # Load spling for Hgg BR
#         self.spline_SMBR_hgg = self.modelBuilder.out.function('fbr_13TeV')

#     def doParametersOfInterest(self):  
#         if self.get_splines:
#             self.make_splines()

#         for variable in self.variables:
#             print 'Processing variable {0} (factory_str: {1})'.format(variable.name, variable.factory_str)
#             self.modelBuilder.doVar(variable.factory_str)

#         for expr in self.exprs:
#             expr = expr.replace('#spline_hgg', self.spline_SMBR_hgg.GetName())
#             expr = expr.replace('#spline_hzz', self.spline_SMBR_hzz.GetName())
#             print 'Processing factory {0}'.format(expr)
#             self.modelBuilder.factory_(expr)

#         MultiSignalModel.doParametersOfInterest(self)

#         # Add own variables to POI set
#         poi_set = self.modelBuilder.out.set('POI')
#         for variable in self.extended_pois:
#             roovariable = self.modelBuilder.out.var(variable.name)
#             print 'Adding poi {0} to the POI set'.format(roovariable.GetName())
#             poi_set.add(roovariable)

#         print 'Checking r_hgg and r_hzz'
#         print self.modelBuilder.out.var('r_hgg')
#         print self.modelBuilder.out.var('r_hzz')
#         sys.exit()

#         print 'Leaving doParametersOfInterest'

#     def getYieldScale(self, bin, process):
#         print 'Calling self.getYieldScale(bin={0}, process={1})'.format(bin, process)
#         MultiSignalModel.getYieldScale(self, bin, process)


extendedMultiSignalModel = ExtendedMultiSignalModel()
