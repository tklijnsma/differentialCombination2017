# from CouplingModel import CouplingModel, add_method

import sys, re
import RooFactoryInterface

#____________________________________________________________________
class Container:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

#____________________________________________________________________
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
            if not( yp.left <= left and left <= yp.right ):
                # raise RuntimeError((
                #     'Process {0} should be scaled by the overflow yieldParameter, '
                #     'but the yieldParameter has range {1} to {2}, '
                #     'where as the process left bound is {3}'
                #         ).format(process, yp.left, yp.right, left)
                #     )
                return 1
            else:
                return yp.name

        if self.debug:
            print 'search_yieldParameter_in_left_right called with left={0}, right = {1}'.format(left, right)
            for e in some_list:
                print '   ',e.name

        for yp in some_list:
            if left >= yp.left and right <= yp.right:
                return yp.name
        else:
            if self.debug: print '    Could not find a match in this list'
            return 1

    def str_to_num(self, s):
        num = s.replace('p','.').replace('m','-')
        return float(num)


#____________________________________________________________________
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
    bkg_yieldParameter = RooFactoryInterface.RooProduct('r_bkg')
    YieldParameterContainer.bkg_yieldParameter = bkg_yieldParameter


    # ======================================
    # Add modifiers to these yield parameters

    if self.MakeLumiScalable:
        # Luminosity modifier works on all parameters
        bkg_yieldParameter.add_variable('lumiScale')
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
            reweightor.formula = '{0}/{1}'.format(expectedXS, rooParametrization.SMXS)
            rooParametrization.reweightor = reweightor
            self.commit_parseable_to_ws(reweightor)

        for ggH_yieldParameter in YieldParameterContainer.all_ggH_yieldParameters():
            parametrization_name = [ v for v in ggH_yieldParameter.variables if 'parametrization' in v ][0]
            reweightor_name = parametrization_name.replace('parametrization', 'reweightor')
            ggH_yieldParameter.add_variable(reweightor_name)

    if self.FitBR:
        # Add a previously defined scale parameter to the ggH yieldParameters
        scaleParameter = { 'hgg' : 'scalingBR_hggModifier', 'hzz' : 'scalingBR_hzzModifier' }
        for decayChannel in decayChannels:
            yieldParameterContainer = self.yieldParameters_per_decay_channel[decayChannel]
            for ggH_yieldParameter in yieldParameterContainer.ggH_yieldParameters:
                ggH_yieldParameter.add_variable(scaleParameter[decayChannel])

        for xH_yieldParameter in YieldParameterContainer.all_xH_yieldParameters():
            xH_yieldParameter.add_variable('scalingBR_xHModifier')


    # ======================================
    # Import in ws

    if self.verbose:
        print '\n' + '- '*30
        print 'Importing in ws and checking\n'

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


#____________________________________________________________________
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
            raise RuntimeError( 'Failure for process \'{0}\': Production process is not \'xH\', \'ggH\' or \'smH\''.format(process) )

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


#____________________________________________________________________
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


#____________________________________________________________________
def commit_parseable_to_ws(self, parameter):
    _v = parameter.verbose
    parameter.verbose = True
    self.modelBuilder.factory_(parameter.parse())
    parameter.verbose = _v
    if self.verbose:
        print 'Test evaluation of {0}:'.format(parameter.name)
        self.modelBuilder.out.function(parameter.name).Print()
        print


#____________________________________________________________________
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


#____________________________________________________________________
def distinguish_between_decay_channels(self):
    if self.FitBR:
        return True
    return False


#____________________________________________________________________
def get_decay_channels(self):
    
    decayChannels = []
    for b in self.DC.bins:
        match = re.match( r'(h[a-zA-Z]+)_', b )
        if match:
            decayChannels.append(match.group(1))
    decayChannels = list(set(decayChannels))

    return decayChannels


#____________________________________________________________________
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


########################################
# Old code
########################################

#____________________________________________________________________
def defineYieldParameters_old( self ):
    self.chapter( 'Starting model.defineYieldParameters()' )

    # Define some variables in local for convenience
    nExpBins            = self.nExpBins
    expBinBoundaries    = self.expBinBoundaries
    nTheoryBins         = self.nTheoryBins
    theoryBinBoundaries = self.theoryBinBoundaries
    parametrizations    = self.parametrizations


    # ======================================
    # Extra variables for further studies

    self.hgg_yieldParameterNames = []
    self.hzz_yieldParameterNames = []

    if self.MakeLumiScalable:
        self.modelBuilder.doVar( 'lumiScale[1.0]' )
        self.modelBuilder.out.var('lumiScale').setConstant()


    # ======================================
    # Start of loop over experimental bins

    SMXSInsideExperimentalBins = []
    self.yieldParameterNames = []
    for iExpBin in xrange(nExpBins):

        expBoundLeft  = expBinBoundaries[iExpBin]
        expBoundRight = expBinBoundaries[iExpBin+1]

        # Determine if bin is the overflow bin
        isOverflowBin = False
        if expBoundLeft == self.allProcessBinBoundaries[-1]:
            isOverflowBin = True

        # if iExpBin == nExpBins-1:
        #     expBinStr     =  '{0}_{1}_GT{2}'.format(
        #         self.productionMode, self.observableName,
        #         expBoundLeft if not expBoundLeft.is_integer() else int(expBoundLeft)
        #         )
        if isOverflowBin:
            expBinStr     = '{0}_{1}_GT{2}'.format(
                self.productionMode, self.observableName,
                expBoundLeft if not expBoundLeft.is_integer() else int(expBoundLeft)
                )
        else:
            expBinStr     =  '{0}_{1}_{2}_{3}'.format(
                self.productionMode, self.observableName,
                expBoundLeft  if not expBoundLeft.is_integer() else int(expBoundLeft),
                expBoundRight if not expBoundRight.is_integer() else int(expBoundRight)
                )


        if self.verbose:
            print '\n' + '='*80
            print 'Processing bin {0}: {1}'.format( iExpBin, expBinStr )

        theoryBinBoundariesInsideExperimentalBin = []
        parametrizationIndices                   = []

        # Extrapolation on left of theory spectrum
        if expBoundLeft < theoryBinBoundaries[0]:
            theoryBinBoundariesInsideExperimentalBin.append( expBoundLeft )
            parametrizationIndices.append( 0 )


        # ======================================
        # Determine which parametrizations enter into the experimental bin

        for iTheoryBin in xrange(nTheoryBins):

            # Left-most bin
            if (
                expBoundLeft >= theoryBinBoundaries[iTheoryBin] and
                expBoundLeft < theoryBinBoundaries[iTheoryBin+1]
                ):
                theoryBinBoundariesInsideExperimentalBin.append( expBoundLeft )
                parametrizationIndices.append( iTheoryBin+1 )

            # Bins in between; only stricly 'between' (and not 'on') the left and right bounds
            if (
                expBoundLeft < theoryBinBoundaries[iTheoryBin] and
                theoryBinBoundaries[iTheoryBin] < expBoundRight
                ):
                theoryBinBoundariesInsideExperimentalBin.append( theoryBinBoundaries[iTheoryBin] )
                parametrizationIndices.append( iTheoryBin+1 )

            # Right-most bin
            if (
                expBoundRight > theoryBinBoundaries[iTheoryBin] and
                expBoundRight <= theoryBinBoundaries[iTheoryBin+1]
                ):
                theoryBinBoundariesInsideExperimentalBin.append( expBoundRight )
                break

        # Extrapolation on right of theory spectrum
        # Sep28: Stop doing this
        #        Use only the parametrizations defined by theorists
        if expBoundRight > theoryBinBoundaries[-1]:
            # theoryBinBoundariesInsideExperimentalBin.append( expBoundRight )
            # parametrizationIndices.append( nTheoryBins+1 )
            theoryBinBoundariesInsideExperimentalBin.append( theoryBinBoundaries[-1] )


        nTheoryBinsInsideExperimentalBin = len(theoryBinBoundariesInsideExperimentalBin)-1
        theoryBinWidthsInsideExperimentalBin = [
            theoryBinBoundariesInsideExperimentalBin[i+1] - theoryBinBoundariesInsideExperimentalBin[i] for i in xrange(nTheoryBinsInsideExperimentalBin)
            ]


        if self.verbose:
            # print '\nProcessing experimental bin {0}'.format(expBinStr)
            print 'Found the following theory bin boundaries inside experimental bin {0} to {1}:'.format(
                expBoundLeft, expBoundRight )
            print '  ', theoryBinBoundariesInsideExperimentalBin
            print '   The following parametrizations will be used to evaluate the cross section:'
            print '  ', parametrizationIndices


        if self.verbose:
            print '\n' + '-'*60
            print 'Calculating cross section'

        # Calculate total cross section (*not* /GeV) in experimental bin
        SMXSInsideExperimentalBin = 0.
        for iParametrization, binWidth in zip( parametrizationIndices, theoryBinWidthsInsideExperimentalBin ):
            SMXSInsideExperimentalBin += self.SMXS[iParametrization] * binWidth

            if self.verbose > 1:
                print '\n    Adding contribution of theory bin {0} to SM cross section in exp bin {1}'.format( iParametrization, expBinStr )
                print '    contribution = {2:+8.3f} ( self.SMXS[iParametrization] = {0:+8.3f} * binWidth = {1:+8.3f} )'.format(
                    self.SMXS[iParametrization] ,
                    binWidth ,
                    self.SMXS[iParametrization] * binWidth ,
                    )

        if self.verbose:
            print '\nTotal cross section in {1} is {0}'.format( SMXSInsideExperimentalBin, expBinStr )


        if len(theoryBinWidthsInsideExperimentalBin) == 0:
            print '  Did not find any theoretical bin boundaries inside this experimental bin'
            print '   (i.e. can not do a parametrization here)'
            print '   Yield parameter will be 1.0'
            self.modelBuilder.doVar( 'r_{0}[1.0]'.format( expBinStr ))
            continue

        SMXSInsideExperimentalBins.append( SMXSInsideExperimentalBin )


        if self.verbose:
            print '\n' + '-'*60
            print 'Creating expression for yield parameter'

        componentWeights = []
        for iComponent in xrange(self.nComponents):
            parametrizationWeights = []
            for iParametrization, binWidth in zip( parametrizationIndices, theoryBinWidthsInsideExperimentalBin ):
                weight = (
                    ( binWidth * self.SMXS[iParametrization] )
                    /
                    SMXSInsideExperimentalBin
                    )
                parametrizationWeights.append(weight)

                if self.verbose > 1:
                    print '\n  Weight for component {0}, parametrization {1}'.format( iComponent, iParametrization )
                    print '    weight = {0:+8.3f} ( ( binWidth[{1:+8.3f}] * SMXS/GeV[{2:+8.3f}] ) / SMXS_expBin[{3:+8.3f}] )'.format(
                        weight, binWidth, self.SMXS[iParametrization], SMXSInsideExperimentalBin
                        )

            componentWeights.append( parametrizationWeights )


        # ======================================
        # Calculate the weighted average of components

        if self.verbose:
            print '\nCalculating weighted average components'

        averageComponents = []
        for iComponent in xrange(self.nComponents):

            if self.verbose > 1:
                print '\n  Component {0} ('.format(iComponent), self.couplingCombinations[iComponent], '):'

            averageComponent = 0.
            parametrizationWeights = componentWeights[iComponent]
            for iParametrization, weight in zip( parametrizationIndices, parametrizationWeights ):
                averageComponent += weight * parametrizations[iParametrization][iComponent]

                if self.verbose > 1:
                    print '    Parametrization {0:3}: product = {1:+8.3f} ( weight[{2:+8.3f}] * parameterValue[{3:+8.3f}]'.format(
                        iParametrization, 
                        weight * parametrizations[iParametrization][iComponent],
                        weight,
                        parametrizations[iParametrization][iComponent]
                        )

            if self.verbose > 1:
                print '  average Component {0} = {1}'.format( iComponent, averageComponent )

            averageComponents.append( averageComponent )


        if self.verbose:
            print '\nThe determined polynomial for bin {0} is:\n'.format(expBinStr)

            line = 'r_{0} = '.format(expBinStr)
            for iComponent in xrange(self.nComponents):

                line += '{0}*{1} + '.format(
                    averageComponents[iComponent],
                    '*'.join(self.couplingCombinations[iComponent])
                    )
            print line [:-2]


        # ======================================
        # Importing into WS is somewhat delicate

        argumentIndices = { coupling : '@{0}'.format(iCoupling) for iCoupling, coupling in enumerate(self.couplings) }

        averageComponentNames = []
        productNames = []
        yieldParameterFormula = []
        for iAverageComponent, averageComponent in enumerate(averageComponents):

            couplingList = self.couplingCombinations[iAverageComponent]
            nCouplingsForThisParameter = len(couplingList)

            averageComponentName = 'averageComponent_{signal}_{whichCouplings}'.format(
                signal = expBinStr,
                whichCouplings = ''.join(couplingList) if nCouplingsForThisParameter > 0 else 'CONSTANT',
                )

            self.modelBuilder.doVar(
                '{0}[{1}]'.format(
                    averageComponentName,
                    averageComponent,
                    )
                )
            averageComponentNames.append( averageComponentName )


            # First add a "@number" entry in the argumentIndices dict, then use that in the formula
            argumentIndices[averageComponentName] = '@{0}'.format( len(argumentIndices) )

            yieldParameterFormulaComponent = argumentIndices[averageComponentName]
            for coupling in couplingList:
                yieldParameterFormulaComponent += '*{0}'.format( argumentIndices[coupling] )
            yieldParameterFormula.append( yieldParameterFormulaComponent )

        # Compile expression string for final yieldParameter
        yieldParameterExpression = 'expr::r_{signal}( "({formulaString})", {commaSeparatedParameters} )'.format(
            signal                   = expBinStr,
            formulaString            = '+'.join(yieldParameterFormula),
            commaSeparatedParameters = ','.join( self.couplings + averageComponentNames )
            )

        if self.verbose:
            print '\nFinal yield parameter expression: {0}'.format(yieldParameterExpression)
            print '  Overview:'
            for key, value in argumentIndices.iteritems():
                print '    {0:4} = {1}'.format( value, key )

        self.modelBuilder.factory_( yieldParameterExpression )

        if self.verbose:
            print '\nTest evaluation of yieldParameter:'
            for coupling in self.couplings:
                self.modelBuilder.out.var(coupling).Print()
            self.modelBuilder.out.function( 'r_{0}'.format(expBinStr) ).Print()
            print ''

        self.yieldParameterNames.append( 'r_{0}'.format(expBinStr) )


        # ======================================
        # Add other modifiers

        hgg_modifiers = [ 'r_{0}'.format(expBinStr) ]
        hzz_modifiers = [ 'r_{0}'.format(expBinStr) ]

        if self.MakeLumiScalable:
            hgg_modifiers.append( 'lumiScale' )
            hzz_modifiers.append( 'lumiScale' )

        if self.FitBR:
            hgg_modifiers.append( 'hggBRmodifier' )
            hzz_modifiers.append( 'hzzBRmodifier' )

        if self.ProfileTotalXS:
            hgg_modifiers.append( 'totalXSmodifier' )
            hzz_modifiers.append( 'totalXSmodifier' )

        if self.FitRatioOfBRs:
            hgg_modifiers.append( 'hgg_ratioBRmodifier' )
            hzz_modifiers.append( 'hzz_ratioBRmodifier' )

        self.modelBuilder.factory_( 'prod::r_hgg_{0}({1})'.format( expBinStr, ','.join(hgg_modifiers) ) )
        self.modelBuilder.factory_( 'prod::r_hzz_{0}({1})'.format( expBinStr, ','.join(hzz_modifiers) ) )

        self.hgg_yieldParameterNames.append( 'r_hgg_{0}'.format(expBinStr) )
        self.hzz_yieldParameterNames.append( 'r_hzz_{0}'.format(expBinStr) )


    # ======================================
    # Other yield parameters: xH, OutsideAcceptance, and background

    self.modelBuilder.doVar( 'one[1]' )
    self.modelBuilder.out.var('one').setConstant(True)

    # xH
    hgg_xH_modifier = [ 'one' ]
    hzz_xH_modifier = [ 'one' ]

    if self.MakeLumiScalable:
        hgg_xH_modifier.append( 'lumiScale' )
        hzz_xH_modifier.append( 'lumiScale' )
    if self.FitBR:
        hgg_xH_modifier.extend([ 'xH_modifier', 'hggBRmodifier' ])
        hzz_xH_modifier.extend([ 'xH_modifier', 'hzzBRmodifier' ])

    self.modelBuilder.factory_( 'prod::hgg_xH_modifier({0})'.format( ','.join(hgg_xH_modifier) ) )
    self.modelBuilder.factory_( 'prod::hzz_xH_modifier({0})'.format( ','.join(hzz_xH_modifier) ) )

    # OutsideAcceptance
    hgg_OutsideAcceptance_modifier = [ 'one' ]
    hzz_OutsideAcceptance_modifier = [ 'one' ]

    if self.FitBR:
        hgg_OutsideAcceptance_modifier.append( 'hggBRmodifier' )
        hzz_OutsideAcceptance_modifier.append( 'hzzBRmodifier' )
    if self.MakeLumiScalable:
        hgg_OutsideAcceptance_modifier.append( 'lumiScale' )
        hzz_OutsideAcceptance_modifier.append( 'lumiScale' )
    
    self.modelBuilder.factory_( 'prod::hgg_OutsideAcceptance_modifier({0})'.format( ','.join(hgg_OutsideAcceptance_modifier) ) )
    self.modelBuilder.factory_( 'prod::hzz_OutsideAcceptance_modifier({0})'.format( ','.join(hzz_OutsideAcceptance_modifier) ) )

    # bkg
    bkg_modifier = [ 'one' ]
    if self.MakeLumiScalable:
        bkg_modifier.append( 'lumiScale' )
    self.modelBuilder.factory_( 'prod::bkg_modifier({0})'.format( ','.join(bkg_modifier) ) )


    # ======================================
    # Save some sets to ws

    self.modelBuilder.out.defineSet( 'POI', ','.join(self.couplings) )
    self.modelBuilder.out.defineSet( 'couplings', ','.join(self.couplings) )
    self.modelBuilder.out.defineSet( 'yieldParameters', ','.join(self.yieldParameterNames) )
    self.modelBuilder.out.defineSet( 'hgg_yieldParameters', ','.join(self.hgg_yieldParameterNames) )
    self.modelBuilder.out.defineSet( 'hzz_yieldParameters', ','.join(self.hzz_yieldParameterNames) )

    self.modelBuilder.out.defineSet( 'modifiers', ','.join([
        'hgg_xH_modifier',
        'hzz_xH_modifier',
        'hgg_OutsideAcceptance_modifier',
        'hzz_OutsideAcceptance_modifier',
        'bkg_modifier',
        ]))

    self.SMXSInsideExperimentalBins = SMXSInsideExperimentalBins

    if self.FitRatioOfBRs:
        self.modelBuilder.out.set('POI').add( self.modelBuilder.out.var('hgg_ratioBRmodifier') )
        self.modelBuilder.out.set('POI').add( self.modelBuilder.out.var('ratio_BR_hgg_hzz') )


    if self.verbose:

        print '\n\n*** Done defining all yield parameters; control printout:'

        allVarsToPrint = self.hgg_yieldParameterNames + self.hzz_yieldParameterNames + [
            'hgg_xH_modifier',
            'hzz_xH_modifier',
            'hgg_OutsideAcceptance_modifier',
            'hzz_OutsideAcceptance_modifier',
            'bkg_modifier',
            ]

        for name in allVarsToPrint:
            # print ''
            self.modelBuilder.out.function(name).Print()

        print '\n'



#____________________________________________________________________
def getYieldScale_old( self, bin, process ):

    if self.verbose:
        print 'Getting scale for process = {0:21}, bin = {1}'.format( process, bin )

    # 'hgg_xH_modifier',
    # 'hzz_xH_modifier',
    # 'hgg_OutsideAcceptance_modifier',
    # 'hzz_OutsideAcceptance_modifier',
    # 'bkg_modifier',

    if not self.DC.isSignal[process]:
        yieldParameter = 'bkg_modifier'

    elif 'OutsideAcceptance' in process:
        yieldParameter = 'hgg_OutsideAcceptance_modifier' if self.BinIsHgg(bin) else 'hzz_OutsideAcceptance_modifier'

        if self.FitOnlyNormalization:
            yieldParameter = 'globalTotalXSmodifier'

    else:
        if self.splitggH and 'xH' in process:
            yieldParameter = 'hgg_xH_modifier' if self.BinIsHgg(bin) else 'hzz_xH_modifier'

        elif self.FitOnlyNormalization:
            yieldParameter = 'globalTotalXSmodifier'

        elif ( self.splitggH and 'ggH' in process ) or ( 'smH' in process ):
            # If it's a ggH process

            isOverflowBin = False
            isInsideUserBinBoundaries = False

            matchRegularBin  = re.search( r'([\dmp]+)_([\dmp]+)', process )
            matchOverflowBin = re.search( r'([GTLE]+)([\dmp]+)', process )

            if matchRegularBin:
                leftBound  = float(matchRegularBin.group(1).replace('m','-').replace('p','.'))
                rightBound = float(matchRegularBin.group(2).replace('m','-').replace('p','.'))
                for iBin in xrange(self.nExpBins):
                    if leftBound >= self.expBinBoundaries[iBin] and rightBound <= self.expBinBoundaries[iBin+1]:
                        yieldParameter = self.hgg_yieldParameterNames[iBin] if self.BinIsHgg(bin) else self.hzz_yieldParameterNames[iBin]
                        break
                else:
                    yieldParameter = 'hgg_OutsideAcceptance_modifier' if self.BinIsHgg(bin) else 'hzz_OutsideAcceptance_modifier'

            elif matchOverflowBin:
                # bound = float(matchOverflowBin.group(1).replace('m','-').replace('p','.'))
                lastYieldParameter = self.hgg_yieldParameterNames[-1] if self.BinIsHgg(bin) else self.hzz_yieldParameterNames[-1]
                if matchOverflowBin.group() in lastYieldParameter:
                    yieldParameter = lastYieldParameter
                else:
                    yieldParameter = 'hgg_OutsideAcceptance_modifier' if self.BinIsHgg(bin) else 'hzz_OutsideAcceptance_modifier'

            else:
                raise CouplingModelError( 'Failure for process \'{0}\': Could not extract any bin boundary information'.format(process) )

        else:
            raise CouplingModelError( 'Failure for process \'{0}\': Production process is not \'xH\', \'ggH\' or \'smH\''.format(process) )


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
