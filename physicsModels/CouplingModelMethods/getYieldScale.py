from physicsModels.MethodHandler import flag_as_method
import sys, re

class BinProcInterpreter(object):
    """docstring for BinProcInterpreter"""
    def __init__(self, model, bin, proc):
        super(BinProcInterpreter, self).__init__()
        self.model = model
        self.bin = bin
        self.proc = proc
        # self.debug = True
        self.debug = False
        
        self.is_hzz = False
        self.is_hgg = False
        self.is_hbb = False

        self.is_smH = False
        self.is_ggH = False
        self.is_xH  = False
        self.is_OOA = False

        self.is_signal = False
        self.is_bkg = False

        self.left = None
        self.right = None

        self.analyze()


    def analyze(self):
        # Reconstructed
        if self.model.isOnlyHZZ or self.bin.startswith('hzz_'):
            self.is_hzz = True
        if self.model.isOnlyHgg or self.bin.startswith('hgg_'):
            self.is_hgg = True
        if getattr(self.model, 'isOnlyHbb', False) or self.bin.startswith('hbb_'):
            self.is_hbb = True

        # Generator
        if self.proc.startswith('smH'):
            self.is_smH = True
        elif self.proc.startswith('ggH'):
            self.is_ggH = True
        elif self.proc.startswith('xH'):
            self.is_xH = True

        if self.model.DC.isSignal[self.proc]:
            # This considers OOA signal, but that's okay
            self.is_signal = True
        else:
            self.is_bkg = True

        if 'OutsideAcceptance' in self.proc:
            self.is_OOA = True
        elif 'INC' in self.proc:
            # flag?
            pass
        elif self.is_smH or self.is_ggH or self.is_xH:
            self.get_boundaries()


    def decay_channel(self):
        if self.is_hzz:
            dc = 'hzz'
        elif self.is_hgg:
            dc = 'hgg'
        elif self.is_hbb:
            dc = 'hbb'
        else:
            raise RuntimeError(
                'Tried to get the decay channel from {0}/{1}, but no decay channel was determined.'
                .format(self.bin, self.proc)
                )
        return dc

    def get_boundaries(self):
        match_regular_bin  = re.search( r'([\dmp]+)_([\dmp]+)', self.proc )
        if match_regular_bin:
            self.left  = self.str_to_num(match_regular_bin.group(1))
            self.right = self.str_to_num(match_regular_bin.group(2))
        else:
            match_overflow_bin = re.search( r'[GTLE]+([\dmp]+)', self.proc )
            if match_overflow_bin:
                self.left  = self.str_to_num(match_overflow_bin.group(1))
            else:
                print 'Process {0} has no clearly defined range; get_boundaries failed'.format(self.proc)
                return
        if self.debug:
            print 'get_boundaries called with self.proc = {0}; left = {1}, right = {2}'.format(self.proc, self.left, self.right)

    def str_to_num(self, s):
        num = s.replace('p','.').replace('m','-')
        return float(num)



@flag_as_method
def getYieldScale( self, bin, process ):
    if self.verbose: sys.stdout.write('proc = {0:26} | bin = {1:26} | '.format(process[:26], bin[:26]))

    # Run a general analyzer over the bin/process combination
    binproc = BinProcInterpreter(self, bin, process)

    r = 1 # Return / yield parameter

    # For fitting only the total XS
    if self.FitOnlyNormalization:
        if binproc.is_smH or binproc.is_ggH:
            if self.BRs_kappa_dependent:
                r = 'r_{0}_smH_INC'.format(binproc.decay_channel())
            else:
                r = 'r_smH_INC'
        else:
            r = 1

        if self.verbose: sys.stdout.write('r = {0}\n'.format(r))
        return r

    # Retrieve the right YieldParameterContainer
    yieldParameterContainer = (
        self.yieldParameters_per_decay_channel[binproc.decay_channel()] \
        if self.distinguish_between_decay_channels() else self.yieldParameters
        )

    if binproc.is_bkg:
        r = yieldParameterContainer.bkg_yieldParameter.name
    elif binproc.is_OOA:
        r = yieldParameterContainer.OutsideAcceptance_yieldParameter.name
    elif binproc.is_ggH or binproc.is_smH or binproc.is_xH:
        r = yieldParameterContainer.get_match_for_binproc(binproc)

    if self.verbose: sys.stdout.write('r = {0}\n'.format(r))
    return r

    # if self.FitOnlyNormalization:
    #     if ( self.splitggH and 'ggH' in process ) or ( 'smH' in process ):
    #         return 'totalXSmodifier'
    #     else:
    #         return 1

    # # Retrieve the right YieldParameterContainer
    # if self.distinguish_between_decay_channels():
    #     match = re.match( r'(h[a-zA-Z]+)_', bin )
    #     if not match:
    #         raise RuntimeError( 'Cannot determine decay channel for bin {0}'.format(bin) )
    #     decayChannel = match.group(1)
    #     yieldParameterContainer = self.yieldParameters_per_decay_channel[decayChannel]
    # else:
    #     yieldParameterContainer = self.yieldParameters


    # # bkg
    # if not self.DC.isSignal[process]:
    #     yieldParameter = yieldParameterContainer.bkg_yieldParameter.name

    # # OutsideAcceptance
    # elif 'OutsideAcceptance' in process:
    #     yieldParameter = yieldParameterContainer.OutsideAcceptance_yieldParameter.name

    # # signal or xH
    # else:
    #     if self.splitggH and 'xH' in process:
    #         yieldParameter = yieldParameterContainer.find_corresponding_xH_yieldParameter(process)
    #     elif self.FitOnlyNormalization:
    #         yieldParameter = 'globalTotalXSmodifier'
    #     elif ( self.splitggH and 'ggH' in process ) or ( 'smH' in process ):
    #         yieldParameter = yieldParameterContainer.find_corresponding_ggH_yieldParameter(process)
    #     else:
    #         raise RuntimeError('Failure for process \'{0}\': Production process is not \'xH\', \'ggH\' or \'smH\''.format(process))

    # if self.verbose:
    #     print '    --> Scaling with \'{0}\''.format( yieldParameter )
        
    #     # print '          test print:'
    #     # try:
    #     #     self.modelBuilder.out.var(yieldParameter).Print()
    #     # except ReferenceError:
    #     #     try:
    #     #         self.modelBuilder.out.function(yieldParameter).Print()
    #     #     except ReferenceError:
    #     #         print '          yieldParameter \'{0}\' does not seem to be in the ws!!'.format(yieldParameter)

    # return yieldParameter
