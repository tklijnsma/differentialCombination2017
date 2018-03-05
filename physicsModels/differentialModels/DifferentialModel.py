
from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder

from DifferentialModelUtils import *


class DifferentialModel(MultiSignalModel):
    """docstring for DifferentialModel"""
    def __init__(self):
        MultiSignalModel.__init__(self)
        self.exprs = []
        self.binning = None
        self.verbose = True
        self.make_underflow = False

    def get_processes(self):
        processes = [ s for s in self.DC.signals if not('OutsideAcceptance' in s) and not('xH' in s) ]
        processes.sort(key=left_bound_of_process)
        return processes

    def setPhysicsOptions(self, physOptions):
        # Run expressions before maps
        MultiSignalModel.setPhysicsOptions(self, physOptions)
        for i_po, po in enumerate(physOptions):
            if po.startswith('expr='):
                self.exprs.append(po.replace('expr=',''))
            if po.startswith('binning='):
                self.binning = [ float(i) for i in po.replace('binning=','').split(',') ]
            if po.startswith('make_underflow'):
                self.make_underflow = True

    def doParametersOfInterest(self):
        for expr in self.exprs:
            self.modelBuilder.factory_(expr)

        self.yield_parameters = make_yield_parameters(self.get_processes(), self.binning, self.make_underflow)
        for y_name in self.yield_parameters:
            expr = y_name + '[1.0,-1.0,4.0]'
            print 'factory: {0}'.format(expr)
            self.modelBuilder.doVar(expr)
        self.modelBuilder.out.defineSet('yieldParameters', ','.join(self.yield_parameters))
        self.modelBuilder.out.defineSet('POI', ','.join(self.yield_parameters))
        # self.modelBuilder.doVar('one[1.0]')
        # MultiSignalModel.doParametersOfInterest(self)
        print 'Leaving doParametersOfInterest'

    def getYieldScale(self, bin, process):
        if process in self.get_processes():
            if self.binning is None:
                y = 'r_' + process
            else:
                i_bin = find_i_bin(process, self.binning)
                if i_bin == -1:
                    y = 1
                else:
                    y = self.yield_parameters[i_bin]
        else:
            y = 1
        if self.verbose: print 'Scaling proc:{0}/bin:{1} with {2}'.format(process, bin, y)
        return y


class LumiScaleDifferentialModel(DifferentialModel):
    """docstring for LumiScaleDifferentialModel"""
    def __init__(self):
        DifferentialModel.__init__(self)
        self.map_to_new_y = {}

    def doParametersOfInterest(self):
        print '[in LumiScaleDifferentialModel:doParametersOfInterest]'
        DifferentialModel.doParametersOfInterest(self)
        self.modelBuilder.doVar('lumiScale[8.356546]')
        for y in self.yield_parameters:
            new_y = 'lumiScale_times_' + y
            expr = 'prod::{0}(lumiScale,{1})'.format(new_y, y)
            print 'factory: {0}'.format(expr)
            self.modelBuilder.factory_(expr)
            self.map_to_new_y[y] = new_y

    def getYieldScale(self, bin, process):
        y = DifferentialModel.getYieldScale(self, bin, process, verbose=False)
        if y == 1:
            y = 'lumiScale'
        elif y in self.yield_parameters:
            y = self.map_to_new_y[y]
        else:
            raise ValueError(
                'Do not know what to do with process={0}, '
                'bin={1}, y={2}'
                .format(process, bin, y)
                )
        print 'Scaling proc:{0}/bin:{1} with {2}'.format(process, bin, y)
        return y

