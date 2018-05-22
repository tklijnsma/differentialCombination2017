import differentials
import core
from core import AttrDict
import ROOT
import logging
import numpy

class WSParametrization(object):
    """docstring for WSParametrization"""
    def __init__(self, ws_file):
        logging.debug('Initializing parametrization with file {0}'.format(ws_file))
        self.ws_file = ws_file
        with core.openroot(ws_file) as ws_fp:
            self.w = ws_fp.Get('w')
        self.old_style = False
        self.yield_parameters = []
        self.smxs = []

    def get_yield_parameters(self):
        yp_list = self.get_yield_parameter_arglist()
        for i in xrange(yp_list.getSize()):
            yp = yp_list[i]
            yp_name = yp.GetName()
            logging.debug('Found yield parameter: '.format(yp))

            if not(self.old_style) and yp_name.startswith('r_'):
                new_yp_name = yp_name.replace('r_', 'parametrization_')
                logging.info('Taking {0} instead of {1}'.format(new_yp_name, yp_name))
                new_yp = self.w.function(new_yp_name)
                if new_yp == None:
                    logging.error('Variable {0} does not exist in ws; taking {1} instead'.format(new_yp_name, yp_name))
                else:
                    yp = new_yp
            self.yield_parameters.append(yp)

    def get_yield_parameter_arglist(self):
        # if self.set_exists('all_ggH_yieldParameters'):
            # logging.debug('Found set called all_ggH_yieldParameters')
            # argset = self.w.set('all_ggH_yieldParameters')
        if self.set_exists('parametrizations_exp'):
            logging.debug('Found set called parametrizations_exp')
            argset = self.w.set('parametrizations_exp')
        elif self.set_exists('yieldParameters'):
            logging.debug('Found set called yieldParameters')
            self.old_style = True
            argset = self.w.set('yieldParameters')
        else:
            raise RuntimeError(
                'Sets \'{0}\' and \'{1}\' do not exist in {2}'
                .format('all_ggH_yieldParameters', 'yieldParameters', self.ws_file)
                )
            # raise RuntimeError(
            #     'Set \'{0}\' does not exist in {1}'
            #     .format(set_name, self.ws_file)
            #     )
        ROOT.SetOwnership(argset, False)
        arglist = ROOT.RooArgList(argset)
        ROOT.SetOwnership(arglist, False)
        return arglist

    def set_exists(self, set_name):
        logging.debug('Checking if set {0} exists'.format(set_name))
        s = self.w.set(set_name)
        logging.debug('  Raw repr of w.set: {0}'.format(s))
        if s == None:
            return False
        else:
            return True

    def set(self, name, value):
        logging.debug('Setting {0} to {1}'.format(name, value))
        roovar = self.w.var(name)
        if roovar == None:
            raise RuntimeError(
                'Variable \'{0}\' does not exist in {1}'
                .format(name, self.ws_file)
                )
        roovar.setVal(value)

    def set_smxs(self, smxs):
        self.smxs = smxs

    def get_mus_exp(self, **kwargs):
        if len(self.yield_parameters) == 0:
            self.get_yield_parameters()

        for name, value in kwargs.iteritems():
            self.set(name, value)
        mus = [ yp.getVal() for yp in self.yield_parameters ]
        return mus

    def get_xs_exp(self, **kwargs):
        if len(self.smxs)==0:
            raise RuntimeError('Need to set smxs to a list of xs first')
        mus = self.get_mus_exp(**kwargs)
        if not len(mus) == len(self.smxs):
            raise ValueError(
                'Found {0} yield parameters, and {1} smxs'
                .format(len(mus), len(self.smxs))
                )
        return [ mu * xs for mu, xs in zip(mus, self.smxs) ]




class Parametrization(object):
    """docstring for Parametrization"""
    def __init__(self):
        super(Parametrization, self).__init__()
        self.do_linear_terms = True
        self.variations = []
        self.parametrizations = []

        self.parametrize_by_matrix_inversion = False

        self.c1_name = 'c1'
        self.c2_name = 'c2'

        self.c1_SM = 1.0
        self.c2_SM = 1.0

        self.rebinner = None
        self.SM_set = False


    def evaluate(self, c1, c2):
        return [ p(c1, c2) if abs(p(c1, c2))>1e-12 else 0.0 for p in self.parametrizations ]

    def get_xs_exp(self, c1, c2):
        if self.rebinner is None:
            raise RuntimeError(
                'First define a rebinner so that the theory spectrum may be mapped to the experimental one'
                )
        if not self.SM_set: self.set_SM()
        xs_theory = self.evaluate(c1, c2)
        xs_exp = self.rebinner.rebin_values(xs_theory)
        return xs_exp

    def get_xs_exp_integrated_per_bin(self, c1, c2):
        xs_exp = self.get_xs_exp(c1, c2)
        xs_exp_integrated_per_bin = [ xs*bw for xs, bw in zip(xs_exp, self.bin_widths) ]
        return xs_exp_integrated_per_bin

    def get_mus_exp(self, c1, c2):
        xs_exp = self.get_xs_exp(c1, c2)
        mu_exp = [ xs / xs_SM for xs, xs_SM in zip(xs_exp, self.xs_exp_SM) ]
        return mu_exp

    def make_rebinner(self, theory_binning, exp_binning):
        self.rebinner = differentials.integral.Rebinner(
            bin_boundaries_old = theory_binning,
            bin_boundaries_new = exp_binning
            )
        self.bin_widths = [ r-l for l, r in zip(exp_binning[:-1], exp_binning[1:]) ]

    def set_SM(self):
        self.xs_theory_SM = self.evaluate(self.c1_SM, self.c2_SM)
        if not(self.rebinner is None):
            self.xs_exp_SM = self.rebinner.rebin_values(self.xs_theory_SM)
        self.SM_set = True


    def add_variation(self, c1, c2, crosssections):
        logging.debug('Entering variation for c1={0}, c2={1}'.format(c1, c2))
        self.variations.append(
            AttrDict(c1=c1, c2=c2, xs=crosssections)
            )

    def from_theory_dicts(self, theories):
        for theory in theories:
            c1 = theory[self.c1_name]
            c2 = theory[self.c2_name]
            crosssections = theory.crosssection
            self.add_variation(c1, c2, crosssections)
        self.parametrize()

    def parametrize(self):
        if self.parametrize_by_matrix_inversion:
            if self.do_linear_terms:
                n_needed = 6
            else:
                n_needed = 3
            n = len(self.variations)
            if n < n_needed:
                raise RuntimeError('Need {0} variations, but only {1} were given'.format(n_needed, n))
            elif n > n_needed:
                self.variations = self.variations[:n_needed]
                logging.warning(
                    'Taking only first {0} variations: {1}'
                    .format(
                        n_needed,
                        ', '.join([ '({0}={1},{2}={3})'.format(self.c1_name, v.c1, self.c2_name, v.c2) for v in self.variations ])
                        )
                    )

        self.n_bins = len(self.variations[0].xs)
        c1s = [ v.c1 for v in self.variations ]
        c2s = [ v.c2 for v in self.variations ]

        if self.parametrize_by_matrix_inversion:
            inv_coupling_matrix = self.get_inv_coupling_matrix(c1s, c2s)

        self.parametrizations = []
        for i_bin in xrange(self.n_bins):
            xs = [ v.xs[i_bin] for v in self.variations ]
            if self.parametrize_by_matrix_inversion:
                parametrization = self.get_paramatrization_by_matrix_inversion(inv_coupling_matrix, xs)
            else:
                parametrization = self.get_fitted_parametrization(c1s, c2s, xs)
            self.parametrizations.append(parametrization)


    def get_paramatrization_by_matrix_inversion(self, inv_coupling_matrix, xs):
        xs_column = numpy.array([ [x] for x in xs ])
        coefficients_column = inv_coupling_matrix.dot(xs_column)
        coefficients = [ coefficients_column[i][0] for i in xrange(len(xs)) ]
        return self.get_parabola_for_given_coefficients(coefficients)

    def get_inv_coupling_matrix(self, c1s, c2s):
        coupling_matrix = []
        
        if self.do_linear_terms:
            def get_row(c1, c2):
                # Order should be identical as in get_parabola_for_given_coefficients
                return [ c1**2, c2**2, c1*c2, c1, c2, 1. ]
        else:
            def get_row(c1, c2):
                return [ c1**2, c2**2, c1*c2 ]

        for c1, c2 in zip(c1s, c2s):
            coupling_matrix.append(get_row(c1, c2))

        coupling_matrix = numpy.array(coupling_matrix)
        inv_coupling_matrix = numpy.linalg.inv(coupling_matrix)
        return inv_coupling_matrix

    def get_parabola_for_given_coefficients(self, coefficients):
        if self.do_linear_terms:
            A, B, C, D, E, F = coefficients
            def parabola(c1, c2):
                return A*c1**2 + B*c2**2 + C*c1*c2 + D*c1 + E*c2 + F
        else:
            A, B, C = coefficients
            def parabola(c1, c2):
                return A*c1**2 + B*c2**2 + C*c1*c2
        return parabola


    def get_parabola(self):
        """
        Construct a simple parabolic function.
        Can take numpy arrays.
        """
        if self.do_linear_terms:
            def parabola(c1, c2, A, B, C, D, E, F):
                return A*c1**2 + B*c2**2 + C*c1*c2 + D*c1 + E*c2 + F
        else:
            def parabola(c1, c2, A, B, C):
                return A*c1**2 + B*c2**2 + C*c1*c2
        return parabola

    def get_parabola_for_scipy(self):
        """
        Returns a parabola that can be evaluated for a
        list of couplings, and returns an equal length list
        """
        parabola = self.get_parabola()
        def parabola_for_scipy(coupling_tuple, *coefficients):
            c1s, c2s = coupling_tuple
            return parabola(c1s, c2s, *coefficients)
        return parabola_for_scipy


    def get_fitted_parametrization(self, couplings1, couplings2, crosssections):
        logging.info('Trying to import scipy.optimize.curve_fit')
        from scipy.optimize import curve_fit

        coupling_tuple = ( numpy.array(couplings1), numpy.array(couplings2) )

        xs_tuple = numpy.array(crosssections)

        p0 = numpy.array([ 1.0 for i in xrange(6 if self.do_linear_terms else 3) ])

        fitted_coefficients, cov_mat = curve_fit(
            self.get_parabola_for_scipy(),
            coupling_tuple,
            xs_tuple,
            p0
            )

        fitted_parabola_base = self.get_parabola()
        def fitted_parabola(c1, c2):
            return fitted_parabola_base(c1, c2, *fitted_coefficients)
        return fitted_parabola



#     def ParametrizeByFitting_Scipy(
#             self,
#             containers,
#             couplingsToParametrize = [ 'kappab', 'kappac' ],
#             includeLinearTerms = True,
#             ):

#         print 'Trying to import scipy.optimize.curve_fit'
#         from scipy.optimize import curve_fit

#         self.parametrizedByFitting = True

#         self.couplings = couplingsToParametrize
#         self.nCouplings = len(self.couplings)

#         couplingCombinations = []
#         couplingCombinations += [ list(couplingTuple) for couplingTuple in itertools.combinations( self.couplings, 2 ) ]
#         couplingCombinations += [ [ coupling, coupling ] for coupling in self.couplings ]
#         if includeLinearTerms:
#             couplingCombinations += [ [coupling] for coupling in self.couplings ] + [[]]
#         self.couplingCombinations = couplingCombinations

#         nComponents = len(couplingCombinations)
#         self.nComponents = nComponents

#         if len(containers) < nComponents:
#             Commands.ThrowError( 'Need at least as much input ({0}) as desired components ({1})'.format( len(containers), nComponents ) )
#             sys.exit()

#         nBins = len(containers[0].binBoundaries) - 1
#         self.nBins = nBins
#         self.binBoundaries = containers[0].binBoundaries

#         self.nContainers = len(containers)

#         # ======================================
#         # Fitting


#         # Construct function to be minimized

#         couplingNameToIndex = { coupling : iCoupling for iCoupling, coupling in enumerate(self.couplings) }
#         couplingIndexToName = { iCoupling : coupling for iCoupling, coupling in enumerate(self.couplings) }

#         def parabolicParametrization( C, *coefficients ):
#             paramValue = 0.
#             for couplingList, coefficient in zip( self.couplingCombinations, coefficients ):
#                 component = coefficient
#                 for coupling in couplingList:
#                     component *= C[ couplingNameToIndex[coupling] ]
#                 paramValue += component
#             return paramValue
#         self.parabolicParametrizationFunction = parabolicParametrization


#         self.coefficientsPerBin = []

#         for iBin in xrange(self.nBins):
#             # if iBin == 0: continue

#             print '\n' + '='*50 + '\nFitting bin {0}'.format(iBin)

#             # Construct input for least_sq: ( c1, c2 ) per input container
#             couplingTuple = [ [] for coupling in self.couplings ]
#             for iCoupling, coupling in enumerate(self.couplings):
#                 for container in containers:
#                     couplingTuple[iCoupling].append( getattr(container, coupling) )
#             couplingTuple = tuple(couplingTuple)

#             # Construct input for least_sq: xs[iBin] per input container
#             xsTuple = tuple([ container.crosssection[iBin] for container in containers ])

#             # initial guess for coefficients
#             p0 = tuple([ 1. for i in xrange(self.nComponents) ])


#             print '\n    Overview of given input to minimizer:'
#             print '    couplingTuple :', couplingTuple
#             print '    xsTuple       :', xsTuple
#             print '    p0            :', p0

#             print '\n    Fitting...'
#             fittedCoefficients, covMat = curve_fit(
#                 parabolicParametrization,
#                 couplingTuple,
#                 xsTuple,
#                 p0
#                 )
#             fittedCoefficients = list(fittedCoefficients)
#             print '\n    fittedCoefficients:', fittedCoefficients


#             print '\n    Test evaluations'

#             for i in xrange(len(xsTuple)):

#                 couplings = tuple([ couplingTuple[iCoupling][i] for iCoupling in xrange(len(self.couplings)) ])
#                 xs        = xsTuple[i]

#                 xsParametrization = parabolicParametrization( couplings, *fittedCoefficients )

#                 line = []
#                 for iCoupling, couplingValue in enumerate(couplings):
#                     couplingName = couplingIndexToName[iCoupling]
#                     line.append( '{0:7} = {1:+5.1f}'.format( couplingName, couplingValue ) )

#                 line.append( 'xs_file = {0:10.4f}'.format(xs) )
#                 line.append( 'xs_param = {0:10.4f}'.format(xsParametrization) )

#                 print '    ' + ' | '.join(line)

#             self.coefficientsPerBin.append( fittedCoefficients )




#     def Evaluate( self, **kwargs ):

#         self.verbose = kwargs.get( 'verbose', False )
#         for key, value in kwargs.iteritems():
#             if self.verbose:
#                 if self.THROW_EVALUATION_WARNINGS:print '[in Parametrization] Setting {0} to {1}'.format( key, value )
#             setattr( self, key, value )


#         if self.parametrizedByMatrixInversion:

#             xsFromParametrization = []
#             for iBin in xrange(self.nBins):
#                 components = self.componentsPerBin[iBin]

#                 res = 0.
#                 for couplingCombination, component in zip( self.couplingCombinations, components ):
#                     product = 1.0
#                     for coupling in couplingCombination:
#                         product *= getattr( self, coupling )
#                     res += product * component

#                 xsFromParametrization.append( res )

#             self.THROW_EVALUATION_WARNINGS = False
#             return xsFromParametrization


#         elif self.ParametrizeByFitting:

#             if self.fitWithScipy:

#                 xs = []
#                 for coefficients in self.coefficientsPerBin:
#                     C = tuple([ getattr( self, coupling ) for coupling in self.couplings ])
#                     xs.append(
#                         self.parabolicParametrizationFunction( C, *coefficients )
#                         )
#                 self.THROW_EVALUATION_WARNINGS = False
#                 return xs

#             else:

#                 # Set RooRealVars to the right couplings
#                 for couplingRooRealVar in self.couplingRooRealVars:
#                     couplingRooRealVar.setVal(
#                         getattr( self, couplingRooRealVar.GetName() )
#                         )

#                 xs = []
#                 for fitEval in self.fitEvals:
#                     xs.append( fitEval.getVal() )

#                 self.THROW_EVALUATION_WARNINGS = False
#                 return xs

