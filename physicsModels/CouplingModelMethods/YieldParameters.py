from physicsModels.MethodHandler import flag_as_method
import sys, re
import physicsModels.RooFactoryInterface as RooFactoryInterface

class Container:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class YieldParameterContainer(Container):
    """docstring for YieldParameterContainer"""
    instances = []
    bkg_yieldParameter = None
    debug = True

    @staticmethod
    def all_ggH_yieldParameters():
        ret = []
        for instance in YieldParameterContainer.instances:
            ret.extend(instance.ggH_yieldParameters)
        return ret

    @staticmethod
    def all_xH_yieldParameters():
        ret = []
        for instance in YieldParameterContainer.instances:
            ret.extend(instance.xH_yieldParameters)
        return ret

    @staticmethod
    def all_OutsideAcceptance_yieldParameters():
        ret = []
        for instance in YieldParameterContainer.instances:
            ret.append(instance.OutsideAcceptance_yieldParameter)
        return ret

    def __init__(self, decayChannel=None):
        self.instances.append(self) # Bit dangerous since the objects are not garbage collected now
        self.decayChannel = decayChannel
        self.ggH_yieldParameters = []
        self.xH_yieldParameters = []
        self.OutsideAcceptance_yieldParameter = None

    def find_corresponding_ggH_yieldParameter(self, process):
        left, right = self.get_range_from_process(process)
        return self.search_yieldParameter_in_left_right( left, right, self.ggH_yieldParameters )

    def find_corresponding_xH_yieldParameter(self, process):
        left, right = self.get_range_from_process(process)
        return self.search_yieldParameter_in_left_right( left, right, self.xH_yieldParameters )

    def get_range_from_process(self, process):
        matchRegularBin  = re.search( r'([\dmp]+)_([\dmp]+)', process )
        if matchRegularBin:
            left  = self.str_to_num(matchRegularBin.group(1))
            right = self.str_to_num(matchRegularBin.group(2))
        else:
            matchOverflowBin = re.search( r'[GTLE]+([\dmp]+)', process )
            if matchOverflowBin:
                left  = self.str_to_num(matchOverflowBin.group(1))
                right = None
            else:
                raise RuntimeError( 'Process {0} has no clearly defined range'.format(process) )
        if self.debug:
            print 'get_range_from_process called with process = {0}; left = {1}, right = {2}'.format(process, left, right)
        return ( left, right )

    def search_yieldParameter_in_left_right(self, left, right, some_list):
        if right is None:
            yp = some_list[-1]
            # Changed second comp. to "<"
            # Example: binning [ ..., 350, 600 ], input_yp r_ggH_GT600
            # Then right=None, left=600, and "left(600)<=yp.right(600) == True"
            # To prevent the last True, use strictly smaller than
            if yp.left <= left and left < yp.right:
                return yp.name
            else:
                # raise RuntimeError((
                #     'Process {0} should be scaled by the overflow yieldParameter, '
                #     'but the yieldParameter has range {1} to {2}, '
                #     'where as the process left bound is {3}'
                #         ).format(process, yp.left, yp.right, left)
                #     )
                # Treat as bkg (considers signal processes that are not included in the range)
                return self.bkg_yieldParameter.name

        if self.debug:
            print 'search_yieldParameter_in_left_right called with left={0}, right = {1}'.format(left, right)
            for e in some_list:
                print '   ',e.name

        for yp in some_list:
            if left >= yp.left and right <= yp.right:
                return yp.name
        else:
            if self.debug: print '    Could not find a match in this list; treating as bkg'
            return self.bkg_yieldParameter.name

    def str_to_num(self, s):
        num = s.replace('p','.').replace('m','-')
        return float(num)


@flag_as_method
def defineYieldParameters(self):
    self.chapter( 'Starting model.defineYieldParameters()' )

    self.SMXSInsideExperimentalBins = []

    # Define a constant 1.0
    self.modelBuilder.doVar('one[1]')
    self.modelBuilder.out.var('one').setConstant(True)


    # ======================================
    # Create the parametrizations for the exp binning (by integrating over theory bins)

    expRooParametrizations = []
    for iExpBin in xrange(self.nExpBins):
        expBoundLeft  = self.expBinBoundaries[iExpBin]
        expBoundRight = self.expBinBoundaries[iExpBin+1]
        if self.verbose:
            print '\n' + '- '*30
            print 'Processing bin {0}, from {1} to {2}'.format(iExpBin, expBoundLeft, expBoundRight)
        expRooParametrization = self.make_parametrization_for_experimental_bin(expBoundLeft, expBoundRight)
        print
        self.commit_parseable_to_ws(expRooParametrization)
        expRooParametrizations.append(expRooParametrization)
    self.modelBuilder.out.defineSet( 'parametrizations_exp', ','.join([ p.name for p in expRooParametrizations ]) )


    # ======================================
    # The actual yield parameters differ in the case of bkg, xH, ggH, OutsideAcceptance, hgg/hzz etc.
    # Define yield parameters for all cases

    if self.verbose:
        print '\n' + '- '*30
        print 'Defining final yield parameters and adding modifiers\n'

    # Open up RooProducts for all yield parameters; need more parameters in the case of distinction between decay channels
    if self.distinguish_between_decay_channels():
        decayChannels = self.get_decay_channels()
        self.yieldParameters_per_decay_channel = {}
        for decayChannel in decayChannels:
            c = YieldParameterContainer()
            for rooParametrization in expRooParametrizations:
                ggH_yieldParameter = RooFactoryInterface.RooProduct('r_{0}_{1}'.format(decayChannel, rooParametrization.binStr))
                ggH_yieldParameter.add_variable(rooParametrization.name)
                ggH_yieldParameter.left = rooParametrization.left
                ggH_yieldParameter.right = rooParametrization.right
                c.ggH_yieldParameters.append(ggH_yieldParameter)

                xH_yieldParameter = RooFactoryInterface.RooProduct('r_{0}_{1}'.format(decayChannel, rooParametrization.binStr.replace('ggH','xH')))
                xH_yieldParameter.left = rooParametrization.left
                xH_yieldParameter.right = rooParametrization.right
                c.xH_yieldParameters.append(xH_yieldParameter)

            c.OutsideAcceptance_yieldParameter = RooFactoryInterface.RooProduct('r_{0}_OutsideAcceptance'.format(decayChannel))
            self.yieldParameters_per_decay_channel[decayChannel] = c
    else:
        c = YieldParameterContainer()
        for rooParametrization in expRooParametrizations:
            ggH_yieldParameter = RooFactoryInterface.RooProduct('r_{0}'.format(rooParametrization.binStr))
            ggH_yieldParameter.add_variable(rooParametrization.name)
            ggH_yieldParameter.left = rooParametrization.left
            ggH_yieldParameter.right = rooParametrization.right
            c.ggH_yieldParameters.append(ggH_yieldParameter)

            xH_yieldParameter = RooFactoryInterface.RooProduct('r_{0}'.format(rooParametrization.binStr.replace('ggH','xH')))
            xH_yieldParameter.left = rooParametrization.left
            xH_yieldParameter.right = rooParametrization.right
            c.xH_yieldParameters.append(xH_yieldParameter)

        c.OutsideAcceptance_yieldParameter = RooFactoryInterface.RooProduct('r_OutsideAcceptance')
        self.yieldParameters = c

    # bkg modifier is the same for all
    # one_yieldParameter = RooFactoryInterface.RooProduct('one')
    bkg_yieldParameter = RooFactoryInterface.RooProduct('r_bkg')
    YieldParameterContainer.bkg_yieldParameter = bkg_yieldParameter
    # YieldParameterContainer.one_yieldParameter = one_yieldParameter


    # ======================================
    # Add modifiers to these yield parameters

    if self.MakeLumiScalable:
        self.modelBuilder.doVar('lumiScale[8.356546]')
        # Luminosity modifier works on all parameters
        bkg_yieldParameter.add_variable('lumiScale')
        # one_yieldParameter.add_variable('lumiScale')
        for ggH_yieldParameter in YieldParameterContainer.all_ggH_yieldParameters():
            ggH_yieldParameter.add_variable('lumiScale')
        for xH_yieldParameter in YieldParameterContainer.all_xH_yieldParameters():
            xH_yieldParameter.add_variable('lumiScale')
        for OutsideAcceptance_yieldParameter in YieldParameterContainer.all_OutsideAcceptance_yieldParameters():
            OutsideAcceptance_yieldParameter.add_variable('lumiScale')

    if not(self.ReweightedXS is None):
        # Reweighting only ggH now
        for expectedXS, rooParametrization in zip( self.ReweightedXS, expRooParametrizations ):
            reweightor = RooFactoryInterface.RooFormulaVar( rooParametrization.name.replace('parametrization','reweightor') )
            # reweightor.formula = '{0}/{1}'.format(expectedXS, rooParametrization.SMXS)
            reweightor.formula = '{0}/{1}'.format(rooParametrization.SMXS, expectedXS)
            rooParametrization.reweightor = reweightor
            self.commit_parseable_to_ws(reweightor)
        for ggH_yieldParameter in YieldParameterContainer.all_ggH_yieldParameters():
            parametrization_name = [ v for v in ggH_yieldParameter.variables if 'parametrization' in v ][0]
            reweightor_name = parametrization_name.replace('parametrization', 'reweightor')
            ggH_yieldParameter.add_variable(reweightor_name)

    if self.FitBR:
        # Add a previously defined scale parameter to the ggH yieldParameters
        scaleParameter = { 'hgg' : 'scalingBR_hggModifier', 'hzz' : 'scalingBR_hzzModifier', 'hbb' : 'scalingBR_hbbModifier' }
        for decayChannel in decayChannels:
            yieldParameterContainer = self.yieldParameters_per_decay_channel[decayChannel]
            for ggH_yieldParameter in yieldParameterContainer.ggH_yieldParameters:
                ggH_yieldParameter.add_variable(scaleParameter[decayChannel])
        for xH_yieldParameter in YieldParameterContainer.all_xH_yieldParameters():
            xH_yieldParameter.add_variable('scalingBR_xHModifier')

    if self.ProfileTotalXS:
        for ggH_yieldParameter in YieldParameterContainer.all_ggH_yieldParameters():
            ggH_yieldParameter.add_variable('totalXSmodifier')


    # ======================================
    # Import in ws

    if self.verbose:
        print '\n' + '- '*30
        print 'Importing in ws and checking\n'

    # self.commit_parseable_to_ws(one_yieldParameter)
    self.commit_parseable_to_ws(bkg_yieldParameter)
    for ggH_yieldParameter in YieldParameterContainer.all_ggH_yieldParameters():
        self.commit_parseable_to_ws(ggH_yieldParameter)
    for xH_yieldParameter in YieldParameterContainer.all_xH_yieldParameters():
        self.commit_parseable_to_ws(xH_yieldParameter)
    for OutsideAcceptance_yieldParameter in YieldParameterContainer.all_OutsideAcceptance_yieldParameters():
        self.commit_parseable_to_ws(OutsideAcceptance_yieldParameter)

    self.modelBuilder.out.defineSet( 'all_ggH_yieldParameters', ','.join([ p.name for p in YieldParameterContainer.all_ggH_yieldParameters() ]) )
    self.modelBuilder.out.defineSet( 'all_xH_yieldParameters', ','.join([ p.name for p in YieldParameterContainer.all_xH_yieldParameters() ]) )
    self.modelBuilder.out.defineSet( 'all_OutsideAcceptance_yieldParameters', ','.join([ p.name for p in YieldParameterContainer.all_OutsideAcceptance_yieldParameters() ]) )

    SMXSset = []
    for iExpBin in xrange(self.nExpBins):
        SMXS          = self.SMXSInsideExperimentalBins[iExpBin]
        expBoundLeft  = self.expBinBoundaries[iExpBin]
        expBoundRight = self.expBinBoundaries[iExpBin+1]
        name = 'SMXS_{0}_{1}'.format(expBoundLeft, expBoundRight)
        self.modelBuilder.doVar('{0}[{1}]'.format(name, SMXS))
        SMXSset.append(name)
    self.modelBuilder.out.defineSet('SMXS', ','.join([p for p in SMXSset]))



@flag_as_method
def getYieldScale( self, bin, process ):

    if self.verbose:
        print 'Getting scale for process = {0:21}, bin = {1}'.format( process, bin )

    # Retrieve the right YieldParameterContainer
    if self.distinguish_between_decay_channels():
        match = re.match( r'(h[a-zA-Z]+)_', bin )
        if not match:
            raise RuntimeError( 'Cannot determine decay channel for bin {0}'.format(bin) )
        decayChannel = match.group(1)
        yieldParameterContainer = self.yieldParameters_per_decay_channel[decayChannel]
    else:
        yieldParameterContainer = self.yieldParameters


    # bkg
    if not self.DC.isSignal[process]:
        yieldParameter = yieldParameterContainer.bkg_yieldParameter.name

    # OutsideAcceptance
    elif 'OutsideAcceptance' in process:
        yieldParameter = yieldParameterContainer.OutsideAcceptance_yieldParameter.name

    # signal or xH
    else:
        if self.splitggH and 'xH' in process:
            yieldParameter = yieldParameterContainer.find_corresponding_xH_yieldParameter(process)
        elif self.FitOnlyNormalization:
            yieldParameter = 'globalTotalXSmodifier'
        elif ( self.splitggH and 'ggH' in process ) or ( 'smH' in process ):
            yieldParameter = yieldParameterContainer.find_corresponding_ggH_yieldParameter(process)
        else:
            raise RuntimeError('Failure for process \'{0}\': Production process is not \'xH\', \'ggH\' or \'smH\''.format(process))

    if self.verbose:
        print '    --> Scaling with \'{0}\''.format( yieldParameter )
        
        # print '          test print:'
        # try:
        #     self.modelBuilder.out.var(yieldParameter).Print()
        # except ReferenceError:
        #     try:
        #         self.modelBuilder.out.function(yieldParameter).Print()
        #     except ReferenceError:
        #         print '          yieldParameter \'{0}\' does not seem to be in the ws!!'.format(yieldParameter)

    return yieldParameter


@flag_as_method
def make_parametrization_for_experimental_bin(self, leftBound, rightBound):
    binStr = self.get_binStr(leftBound, rightBound)

    # Determine which parametrizations enter into the experimental bin
    theoryIndices, theoryBinBoundaries = find_contained_theory_bins(leftBound, rightBound, self.theoryBinBoundaries)
    nTheoryBins = len(theoryIndices)
    theoryBinWidths = [ theoryBinBoundaries[i+1] - theoryBinBoundaries[i] for i in xrange(nTheoryBins) ]

    # Calculate SMXS in the bin
    SMXS = 0.0
    for width, theoryIndex in zip( theoryBinWidths, theoryIndices ):
        SMXS += width * self.SMXS[theoryIndex]
    self.SMXSInsideExperimentalBins.append(SMXS)

    # Make debug printout
    if self.verbose:
        print '\nCalculated SMXS:'
        print '  {0:7} | {1:15} | {2:11} | {3:11} | {4:15}'.format(
            'index', 'SMXS/GeV in bin', 'left', 'right', 'SMXS in bin')
        for i in xrange(nTheoryBins):
            print '  {0:7d} | {1:15.5f} | {2:11.5f} | {3:11.5f} | {4:15.5f}'.format(
                theoryIndices[i],
                self.SMXS[theoryIndices[i]],
                theoryBinBoundaries[i], theoryBinBoundaries[i+1],
                theoryBinWidths[i] * self.SMXS[theoryIndices[i]]
                )
        print ' '*57 + '---------------- +'
        print ' '*57 + str(SMXS)

    # Calculate combined parametrization for the experimental bin
    weights = [ theoryBinWidths[i] * self.SMXS[theoryIndices[i]] / SMXS for i in xrange(nTheoryBins) ]
    average_coefficients = RooFactoryInterface.get_average_coefficients(
        weights,
        [ self.rooParametrizations[i] for i in theoryIndices ]
        )

    # Create average rooParametrization and import
    rooParametrization = RooFactoryInterface.RooParametrization('parametrization_{0}'.format(binStr), verbose=True)
    rooParametrization.variables.extend(self.couplings)
    for coefficient, couplingList in zip( average_coefficients, self.couplingCombinations ):
        rooParametrization.add_term( coefficient, couplingList )

    if self.verbose:
        print 'The unmodified polynomial for bin {0} is:'.format(binStr)
        print rooParametrization.name,'=',rooParametrization.get_formula()

    # Do later
    # self.modelBuilder.factory_( rooParametrization.parse() )
    # if self.verbose:
    #     print '\nTest evaluation of {0}:'.format(rooParametrization.name)
    #     for coupling in self.couplings:
    #         self.modelBuilder.out.var(coupling).Print()
    #     self.modelBuilder.out.function(rooParametrization.name).Print()
    #     print ''

    # Attach for easy access
    rooParametrization.binStr = binStr
    rooParametrization.SMXS   = SMXS
    rooParametrization.left   = leftBound
    rooParametrization.right  = rightBound

    return rooParametrization

@flag_as_method
def commit_parseable_to_ws(self, parameter):
    _v = parameter.verbose
    parameter.verbose = True
    self.modelBuilder.factory_(parameter.parse())
    parameter.verbose = _v
    if self.verbose:
        print 'Test evaluation of {0}:'.format(parameter.name)
        self.modelBuilder.out.function(parameter.name).Print()
        print

@flag_as_method
def get_binStr(self, leftBound, rightBound):
    # Determine if bin is the overflow bin
    isOverflowBin = False
    if leftBound == self.allProcessBinBoundaries[-1]:
        isOverflowBin = True

    # Get a name string for the bin
    if isOverflowBin:
        binStr = '{0}_{1}_GT{2}'.format(
            self.productionMode, self.observableName,
            leftBound if not leftBound.is_integer() else int(leftBound)
            )
    else:
        binStr =  '{0}_{1}_{2}_{3}'.format(
            self.productionMode, self.observableName,
            leftBound  if not leftBound.is_integer() else int(leftBound),
            rightBound if not rightBound.is_integer() else int(rightBound)
            )

    return binStr

@flag_as_method
def distinguish_between_decay_channels(self):
    if self.FitBR:
        return True
    return False

@flag_as_method
def get_decay_channels(self):
    decayChannels = []
    for b in self.DC.bins:
        match = re.match( r'(h[a-zA-Z]+)_', b )
        if match:
            decayChannels.append(match.group(1))
    decayChannels = list(set(decayChannels))

    return decayChannels


def find_contained_theory_bins(left, right, theory_bin_boundaries):
    """Returns indices and boundaries of a bin-boundary list between left and right"""

    contained_indices    = []
    contained_boundaries = []

    # Add left extrapolation
    # The parametrization for this extrapolation is currently set to 0.0, but can be set
    # to e.g. the same as the first defined theory bin
    if left < theory_bin_boundaries[0]:
        contained_indices.append(0)
        contained_boundaries.append(left)

    for i_theory_bin in xrange(len(theory_bin_boundaries)-1):

        # Left-most bin
        if (
            left >= theory_bin_boundaries[i_theory_bin] and
            left < theory_bin_boundaries[i_theory_bin+1]
            ):
            contained_boundaries.append( left )
            contained_indices.append( i_theory_bin+1 )

        # Bins in between; only stricly 'between' (and not 'on') the left and right bounds
        if (
            left < theory_bin_boundaries[i_theory_bin] and
            theory_bin_boundaries[i_theory_bin] < right
            ):
            contained_boundaries.append( theory_bin_boundaries[i_theory_bin] )
            contained_indices.append( i_theory_bin+1 )

        # Right-most bin
        if (
            right > theory_bin_boundaries[i_theory_bin] and
            right <= theory_bin_boundaries[i_theory_bin+1]
            ):
            contained_boundaries.append( right )
            break

    if right > theory_bin_boundaries[-1]:
        contained_boundaries.append(theory_bin_boundaries[-1])

    return contained_indices, contained_boundaries

