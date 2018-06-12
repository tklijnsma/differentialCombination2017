from physicsModels.MethodHandler import flag_as_method
import physicsModels.RooFactoryInterface as RooFactoryInterface
import numpy, itertools, sys

@flag_as_method
def makeParametrizationsFromTheory(self):
    self.ParamB = ParametrizationBuilder(self)
    self.ParamB.check_input_dimensions()

    if self.radial_coord_for_ctcg:
        self.ParamB.import_radial_ktkg()
    else:
        self.ParamB.import_couplings()

    self.ParamB.get_coupling_combinations()

    self.parametrizations = self.ParamB.get_parametrization_coefficients()
    self.ParamB.print_parametrizations(self.parametrizations)
    self.rooParametrizations = self.ParamB.get_rooParametrizations(self.parametrizations)
    self.ParamB.import_rooParametrizations(self.rooParametrizations)

    # Better not to use these, but rather the values in ParamB directly. Keep around for legacy
    self.couplingCombinations = self.ParamB.coupling_combinations
    self.nComponents          = self.ParamB.n_coefficients
    self.nCouplings           = self.ParamB.n_couplings
    self.couplings            = self.ParamB.couplings
    self.nTheoryBins          = self.ParamB.n_bins
    self.theoryBinBoundaries  = self.ParamB.bin_boundaries

    # Some sets for easy access after t2ws
    self.ParamB.import_bin_boundaries_as_set()
    self.modelBuilder.out.defineSet( 'parametrizations', ','.join([ p.name for p in self.rooParametrizations ]) )
    
    if self.radial_coord_for_ctcg:
        self.modelBuilder.out.defineSet('POI', 'r,theta')
    else:
        self.modelBuilder.out.defineSet('POI', ','.join(self.couplings))

    self.modelBuilder.out.defineSet( 'couplings', ','.join(self.couplings) )

@flag_as_method
def MakeTotalXSExpressions(self):
    self.chapter('Starting model.MakeTotalXSExpressions()')
    self.ParamB.import_expr_for_totalXS()
    self.ParamB.import_expr_for_totalXS_SM()

    if self.scale_with_mu_totalXS:
        self.ParamB.import_totalXS_modifier()
    elif self.ProfileTotalXS:
        self.ParamB.import_floating_totalXS_modifier()
    elif self.FitOnlyNormalization:
        self.ParamB.import_totalXS_modifier()

        if self.distinguish_between_decay_channels():
            # Will need a yield parameter specific per decay channel
            for decay_channel in self.get_decay_channels():
                yieldparameter = RooFactoryInterface.RooProduct('r_{0}_smH_INC'.format(decay_channel))
                yieldparameter.add_variable('totalXSmodifier')

                # Theory uncertainty
                if not(self.SMXS_of_input_ws is None) and len(self.SMXS_of_input_ws) == 1:
                    self.modelBuilder.factory_('expr::reweightor_smH_INC("@0/{0}", totalXS_SM, )'.format(self.SMXS_of_input_ws[0]))
                    self.modelBuilder.out.function('reweightor_smH_INC').Print()
                    yieldparameter.add_variable('reweightor_smH_INC')

                # BR modifier
                if self.BRs_kappa_dependent:
                    yieldparameter.add_variable(
                        { 'hgg' : 'brmodifier_hgg', 'hzz' : 'brmodifier_hzz', 'hbb' : 'brmodifier_hbb' }[decay_channel]
                        )

                print
                self.commit_parseable_to_ws(yieldparameter)

        else:
            yieldparameter = RooFactoryInterface.RooProduct('r_smH_INC')
            yieldparameter.add_variable('totalXSmodifier')
            if not(self.SMXS_of_input_ws is None) and len(self.SMXS_of_input_ws) == 1:
                self.modelBuilder.factory_('expr::reweightor_smH_INC("@0/{0}", totalXS_SM, )'.format(self.SMXS_of_input_ws[0]))
                self.modelBuilder.out.function('reweightor_smH_INC').Print()
                yieldparameter.add_variable('reweightor_smH_INC')
            print
            self.commit_parseable_to_ws(yieldparameter)


    if not(self.inc_xs_uncertainty is None):
        inc_xs_uncertainty = self.inc_xs_uncertainty
        if inc_xs_uncertainty < 1.: inc_xs_uncertainty += 1. # Uncertainty should be 1 + Delta

        errDict = {}
        for binName in self.DC.bins:
            errDict[binName] = {}
            for processName in self.DC.processes:
                errDict[binName][processName] = 0.
                if self.DC.isSignal[processName] and not 'OutsideAcceptance' in processName:
                    errDict[binName][processName] = inc_xs_uncertainty
               
        systematicName = 'theoryUncertainty_INC'
        self.DC.systs.append(
            ( systematicName, False, 'lnN', [], errDict )
            )
        if self.verbose:
            print 'Added nuisance \'{0}\''.format( systematicName )



class Theory(object):
    """docstring for Theory"""
    def __init__(self, theory_dict):
        super(Theory, self).__init__()
        self.theory_dict = theory_dict
        self.couplings = theory_dict['couplings']
        self.n_couplings = len(self.couplings)
        self.bin_boundaries = theory_dict['binBoundaries']
        self.n_bins = len(self.bin_boundaries)-1
        self.crosssection = theory_dict['crosssection']
        self.n_crosssection = len(self.crosssection)
        if 'ratios' in theory_dict:
            self.ratios = theory_dict['ratios']
            self.n_ratios = len(self.ratios)

    def repr(self):
        s = (
            '{0}:'
            '\n    couplings:            {1}'
            '\n    bin_boundaries ({2}): {3}'
            '\n    crosssection   ({4}): {5}'
            .format(
                self,
                self.couplings,
                self.n_bins, self.bin_boundaries,
                self.n_crosssection, self.crosssection
                )
            )
        return s


class ParametrizationBuilder(object):
    """docstring for ParametrizationBuilder"""
    def __init__(self, model):
        super(ParametrizationBuilder, self).__init__()
        self.model = model
        self.modelBuilder = model.modelBuilder # Convenience
        self.SM = Theory(model.SMDict)
        self.theories = [ Theory(theory_dict) for theory_dict in self.model.theories ]
        self.check_input_dimensions()

    def check_input_dimensions(self):
        first_theory = self.theories[0]
        for theory in self.theories:
            if (
                first_theory.n_couplings == theory.n_couplings
                and first_theory.n_bins == theory.n_bins
                and first_theory.n_crosssection == theory.n_crosssection
                and first_theory.bin_boundaries == theory.bin_boundaries
                and list(sorted(first_theory.couplings.keys())) == list(sorted(theory.couplings.keys()))
                ):
                continue
            else:
                raise ValueError(
                    'At least the following two theories have mismatching dimensions, bin boundaries or couplings:\n'
                    '{0}\nand:\n{1}'
                    .format(first_theory.repr(), theory.repr())
                    )

        self.couplings = first_theory.couplings.keys()
        self.bin_boundaries = first_theory.bin_boundaries
        self.n_couplings = first_theory.n_couplings
        self.n_bins = first_theory.n_bins
        self.n_crosssection = first_theory.n_crosssection

    def get_coupling_combinations(self):
        # Get list of squared terms and unique combinations (works as expected also for 1 couplings cases)
        coupling_combinations = []
        coupling_combinations.extend( [ [ coupling, coupling ] for coupling in self.couplings ] )
        coupling_combinations.extend( [ list(couplingTuple) for couplingTuple in itertools.combinations( self.couplings, 2 ) ] )
        if self.model.includeLinearTerms:
            # Include also singular coupling terms and a constant
            coupling_combinations.extend( [ [coupling] for coupling in self.couplings ] )
            coupling_combinations.append( [] )
        self.coupling_combinations = coupling_combinations
        self.n_coefficients = len(self.coupling_combinations)

        if len(self.theories) > self.n_coefficients:
            print 'Limiting number of theories to {0} (was {1})'.format(self.n_coefficients, len(self.theories))
            self.theories = self.theories[:self.n_coefficients]
        elif len(self.theories) < self.n_coefficients:
            raise ValueError(
                'Need at least {0} theories for {0} coefficients, but found {1}'
                .format(self.n_coefficients, len(self.theories))
                )

    def import_couplings(self):
        # First make sure there are RooRealVars for all the couplings in the workspace
        for coupling in self.couplings:
            self.modelBuilder.doVar(
                '{coupling}[{default},{down},{up}]'.format(
                    coupling = coupling,
                    default  = self.model.SMDict['couplings'][coupling],
                    down     = self.model.SMDict['couplings'][coupling] - 1000.0,
                    up       = self.model.SMDict['couplings'][coupling] + 1000.0,
                    )
                )

    def import_radial_ktkg(self):
        self.modelBuilder.doVar(
            'theta[0.0,{0},{1}]'
            .format(-0.75*3.14159265359, 0.25*3.14159265359)
            )
        self.modelBuilder.doVar(
            'r[1.0,0.0,5.0]'
            )
        self.modelBuilder.factory_(
            'expr::ct( "@0*TMath::Cos(@1)", r, theta )'
            )
        self.modelBuilder.factory_(
            'expr::cg( "@0*TMath::Sin(@1)", r, theta )'
            )

        print '\nImported couplings ct and cg as dependent on radial variables theta and r:'
        self.modelBuilder.out.var('theta').Print()
        self.modelBuilder.out.var('r').Print()
        self.modelBuilder.out.function('ct').Print()
        self.modelBuilder.out.function('cg').Print()


    def get_coupling_matrix(self):
        # Calculate the coupling matrix
        coupling_mat = []
        for theory in self.theories:
            row = []
            for couplings in self.coupling_combinations:
                product = 1.0
                for coupling in couplings:
                    product *= theory.couplings[coupling]
                row.append(product)
            coupling_mat.append(row)
        return coupling_mat

    def get_inverse_coupling_matrix(self):
        coupling_matrix = numpy.array(self.get_coupling_matrix())
        return numpy.linalg.inv(coupling_matrix)

    def get_parametrization_coefficients(self):
        inverse_coupling_matrix = self.get_inverse_coupling_matrix()
        # Parametrization for each bin
        parametrizations = []

        # Underflow (left extrapolation); always zero
        parametrizations.append( [ 0. for i in xrange(self.n_coefficients) ] )
        
        for i in xrange(self.n_bins):
            # Column vector of ratios per theory in bin i
            ratios = numpy.array([ [theory.ratios[i]] for theory in self.theories ])
            # inverse_coupling_matrix dot ratios, flattened to python list
            parametrization = list(itertools.chain.from_iterable(inverse_coupling_matrix.dot(ratios)))
            parametrizations.append(parametrization)
        self.parametrizations = parametrizations
        return parametrizations

    def print_parametrizations(self, parametrizations):
        if self.model.verbose:
            print '\nParametrizations per theory bin:'
            for i, parametrization in enumerate(parametrizations):
                print '{0:4}: {1}'.format( i, parametrization )
            print

    def get_rooParametrizations(self, parametrizations):
        rooParametrizations = []
        for iParametrization, parametrization in enumerate(parametrizations):
            rooParametrization = RooFactoryInterface.RooParametrization('parametrization{0}'.format(iParametrization), verbose=True)
            rooParametrization.variables.extend(self.couplings)
            for coefficient, couplings in zip(parametrization, self.coupling_combinations):
                rooParametrization.add_term(coefficient, couplings)
            rooParametrizations.append(rooParametrization)
        self.rooParametrizations = rooParametrizations
        return rooParametrizations

    def import_rooParametrizations(self, rooParametrizations):
        for rooParametrization in rooParametrizations:        
            self.modelBuilder.factory_(rooParametrization.parse())

    def import_bin_boundaries_as_set(self):
        names_in_set = []
        for i, bin_bound in enumerate(self.bin_boundaries):
            self.modelBuilder.doVar('theoryBinBound{0}[{1}]'.format(i, bin_bound))
            names_in_set.append( 'theoryBinBound{0}'.format(i) )
        self.modelBuilder.out.defineSet('theoryBinBoundaries', ','.join(names_in_set))

    def import_expr_for_totalXS(self):
        inc_xs_per_bin = []

        # Remember the added underflow parametrization
        underflow_bin_width = self.bin_boundaries[0] - 0.0
        self.modelBuilder.factory_(
            'expr::IncXS_underflow( "{0}*{1}*@0", parametrization0 )'
            .format(underflow_bin_width, 0.0) # or a better guess than zero
            )
        inc_xs_per_bin.append('IncXS_underflow')

        for i in xrange(self.n_bins):
            binWidth = self.bin_boundaries[i+1] - self.bin_boundaries[i]
            self.modelBuilder.factory_(
                'expr::IncXS{0}( "{1}*{2}*@0", parametrization{0} )'
                .format(i+1, binWidth, self.SM.crosssection[i])
                )
            inc_xs_per_bin.append( 'IncXS{0}'.format(i+1) )
        self.modelBuilder.factory_('sum::totalXS( {0} )'.format(','.join(inc_xs_per_bin)))
        self.modelBuilder.out.function('totalXS').Print()

    def import_expr_for_totalXS_SM(self):
        # Make sure couplings and MH are at SM
        self.modelBuilder.out.var('MH').setVal(125.)
        for coupling in self.couplings:
            self.modelBuilder.out.var(coupling).setVal(self.SM.couplings[coupling])

        # Get the inclusive SMXS
        self.modelBuilder.doVar( 'totalXS_SM[{0}]'.format(self.modelBuilder.out.function('totalXS').getVal()))
        self.modelBuilder.out.var('totalXS_SM').setConstant(True)
        self.modelBuilder.out.var('totalXS_SM').Print()

    def import_floating_totalXS_modifier(self):
        # Useful for letting the totalXS be profiled in the fit
        self.modelBuilder.doVar( 'r_totalXS[1.0,0.0,3.0]' )
        self.modelBuilder.out.var('r_totalXS').Print()
        # self.modelBuilder.factory_( 'expr::totalXSmodifier( "@0*@1/@2", r_totalXS, totalXS_SM, totalXS )' )
        self.modelBuilder.factory_( 'expr::totalXSmodifier( "@0*@1/@2", r_totalXS, totalXS, totalXS_SM )' )
        self.modelBuilder.out.function('totalXSmodifier').Print()

    def import_totalXS_modifier(self):
        # Useful for fitting only the normalization, and dropping shape information
        # self.modelBuilder.factory_( 'expr::totalXSmodifier( "@0/@1", totalXS_SM, totalXS )' )
        self.modelBuilder.factory_( 'expr::totalXSmodifier( "@0/@1", totalXS, totalXS_SM )' )
        self.modelBuilder.out.function('totalXSmodifier').Print()


    # Not yet implemented - keep old code for now.
    def make_parametrization_for_coarser_bin(self, left, right):
        binStr = self.model.get_binStr(left, right)

        # Determine which parametrizations enter into the experimental bin
        theory_indices, theory_bin_boundaries = find_contained_theory_bins(left, right, self.bin_boundaries)
        n_theory_bins = len(theory_indices)
        theory_bin_widths = [ theory_bin_boundaries[i+1] - theory_bin_boundaries[i] for i in xrange(n_theory_bins) ]

        # Calculate SMXS in the bin
        SMXS = 0.0
        for width, theory_index in zip(theory_bin_widths, theory_indices):
            SMXS += width * self.model.SMXS[theory_index]
        # self.model.SMXSInsideExperimentalBins.append(SMXS) # Should go somewhere else

        # Make debug printout
        if self.model.verbose:
            print '\nCalculated SMXS:'
            print '  {0:7} | {1:15} | {2:11} | {3:11} | {4:15}'.format(
                'index', 'SMXS/GeV in bin', 'left', 'right', 'SMXS in bin')
            for i in xrange(n_theory_bins):
                print '  {0:7d} | {1:15.5f} | {2:11.5f} | {3:11.5f} | {4:15.5f}'.format(
                    theory_indices[i],
                    self.SM.crosssection[theory_indices[i]],
                    theory_bin_boundaries[i], theory_bin_boundaries[i+1],
                    theory_bin_widths[i] * self.SM.crosssection[theory_indices[i]]
                    )
            print ' '*57 + '---------------- +'
            print ' '*57 + str(SMXS)

        # Calculate combined parametrization for the experimental bin
        weights = [ theory_bin_widths[i] * self.SM.crosssection[theory_indices[i]] / SMXS for i in xrange(n_theory_bins) ]
        average_coefficients = RooFactoryInterface.get_average_coefficients(
            weights,
            [ self.rooParametrizations[i] for i in theory_indices ]
            )

        # Create average rooParametrization and import
        rooParametrization = RooFactoryInterface.RooParametrization('parametrization_{0}'.format(binStr), verbose=True)
        rooParametrization.variables.extend(self.couplings)
        for coefficient, couplingList in zip( average_coefficients, self.coupling_combinations ):
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
        rooParametrization.left   = left
        rooParametrization.right  = right

        return rooParametrization


# NOT IN USE YET - this should replace code currently in YieldParameters.py
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
