import differentials
import core
from core import AttrDict
import ROOT
import logging
import numpy
import itertools


class Parabola(object):
    c3 = 1.
    def __init__(self, coefficients):
        super(Parabola, self).__init__()
        self.coefficients = coefficients
        if len(coefficients) == 6:
            self.A, self.B, self.C, self.D, self.E, self.F = coefficients
        else:
            self.A, self.B, self.C = coefficients
            self.D, self.E, self.F = (0., 0., 0.)

    def __call__(self, c1, c2):
        return self.A*c1**2 + self.B*c2**2 + self.C*c1*c2 + self.D*c1*Parabola.c3 + self.E*c2*Parabola.c3 + self.F*Parabola.c3*Parabola.c3


class ParabolaNDim(object):
    """docstring for ParabolaNDim"""
    def __init__(self, coefficients, coupling_combinations):
        super(ParabolaNDim, self).__init__()
        self.coefficients = coefficients
        self.coupling_combinations = coupling_combinations

    def __call__(self, **coupling_vals):
        r = 0.0
        for coeff, comb in zip(self.coefficients, self.coupling_combinations):
            r += coeff * prod([coupling_vals[cname] for cname in comb])
        return r


class Variation(object):
    """docstring for Variation"""
    print_n_xs = 3

    def __init__(self, crosssections, **coupling_vals):
        super(Variation, self).__init__()
        self.coupling_vals = coupling_vals
        self.xs = crosssections

    def __getitem__(self, name):
        return self.coupling_vals[name]

    def __repr__(self):
        names = self.coupling_vals.keys()
        names.sort()
        coupling_str = ' '.join(
            ['{0}={1:+.4f}'.format(name, self.coupling_vals[name]) for name in names]
            )
        xs_str = ', '.join(['{0:+.3f}'.format(xs) for xs in self.xs[:self.print_n_xs]])
        if len(self.xs) > self.print_n_xs:
            xs_str += ', ...'
        return '{0} (xs = {1})'.format(coupling_str, xs_str)


def prod(list):
    r = 1
    for x in list:
        r *= x
    return r


class WSParametrization(object):
    """docstring for WSParametrization"""
    def __init__(self, ws_file):
        logging.debug('Initializing parametrization with file {0}'.format(ws_file))
        self.ws_file = ws_file
        with core.openroot(ws_file) as ws_fp:
            self.w = ws_fp.Get('w')

        self.old_style = False

        self.try_to_include_brmodifiers = False

        self.smxs = []
        self.parametrizations = []
        self.ggH_yield_parameters = [] # Taking into account also brmodifiers for example

        self.smxs_xH = []
        self.xH_yield_parameters = []

        self.decay_channel = None
        self.is_dc_specific = False



    def arglist_to_pylist(self, arglist):
        pylist = []
        for i in xrange(arglist.getSize()):
            pylist.append(arglist[i])
        return pylist

    def prod_list(self, l1, l2):
        if not len(l1) == len(l2):
            raise ValueError(
                'Inconsistent list length for prod_list: Found {0} in l1, and {1} in l2'
                .format(len(l1), len(l2))
                )
        return [ e1 * e2 for e1, e2 in zip(l1, l2) ]


    def set_exists(self, set_name):
        logging.debug('Checking if set {0} exists'.format(set_name))
        s = self.w.set(set_name)
        logging.debug('  Raw repr of w.set: {0}'.format(s))
        if s == None:
            return False
        else:
            return True

    def get_argset_as_pylist(self, set_name):
        if not self.set_exists(set_name):
            raise ValueError('Workspace has no set {0}'.format(set_name))
        argset = self.w.set(set_name)
        ROOT.SetOwnership(argset, False)
        arglist = ROOT.RooArgList(argset)
        ROOT.SetOwnership(arglist, False)
        return self.arglist_to_pylist(arglist)


    def set(self, name, value):
        logging.debug('Setting {0} to {1}'.format(name, value))
        roovar = self.w.var(name)
        if roovar == None:
            raise RuntimeError(
                'Variable \'{0}\' does not exist in {1}'
                .format(name, self.ws_file)
                )
        roovar.setVal(value)

    def set_kwargs(self, kwargs):
        for name, value in kwargs.iteritems():
            self.set(name, value)

    def set_smxs(self, smxs):
        self.smxs = smxs

    def get_smxs_from_ws(self, set_name='SMXS'):
        xs_pars = self.get_argset_as_pylist(set_name)
        self.smxs = [ p.getVal() for p in xs_pars ]


    #____________________________________________________________________
    # For the parametrization only (ignoring other modifiers e.g. BR)

    def get_parametrizations_from_ggH_yield_parameters(self, ggH_yield_parameters):
        parametrizations = []
        for ggH_yp in ggH_yield_parameters:
            param_name = ggH_yp.GetName().replace('r_', 'parametrization_')
            parametrizations.append(self.w.function(param_name))
        return parametrizations

    def get_parametrizations(self):
        self.parametrizations = []

        if self.try_to_include_brmodifiers:
            if self.set_exists('all_ggH_yieldParameters'):
                # Get the parametrizations from the full ggH modifier (simply replace 'r_' by 'parametrization_' in the name)
                logging.debug('Found set called all_ggH_yieldParameters')
                ggH_yield_parameters = self.get_argset_as_pylist('all_ggH_yieldParameters')
                self.parametrizations = ggH_yield_parameters
            else:
                raise RuntimeError(
                    'Tried to include brmodifiers but something went wrong'
                    )
            return

        if self.set_exists('parametrizations_exp'):
            # Newest way, parametrizations binned for exp are stored directly
            logging.debug('Found set called parametrizations_exp')
            self.parametrizations = self.get_argset_as_pylist('parametrizations_exp')

        elif self.set_exists('all_ggH_yieldParameters'):
            # Get the parametrizations from the full ggH modifier (simply replace 'r_' by 'parametrization_' in the name)
            logging.debug('Found set called all_ggH_yieldParameters')
            ggH_yield_parameters = self.get_argset_as_pylist('all_ggH_yieldParameters')
            self.parametrizations = self.get_parametrizations_from_ggH_yield_parameters(ggH_yield_parameters)

        elif self.set_exists('yieldParameters'):
            logging.debug('Found set called yieldParameters')
            ggH_yield_parameters = self.get_argset_as_pylist('yieldParameters')
            self.parametrizations = self.get_parametrizations_from_ggH_yield_parameters(ggH_yield_parameters)

        else:
            raise RuntimeError(
                'Could not find any appropriately named set that could contain parametrizations in {0}'
                .format(self.ws_file)
                )

    def get_mus_exp(self, **kwargs):
        if len(self.parametrizations) == 0:
            self.get_parametrizations()
        self.set_kwargs(kwargs)
        return [ param.getVal() for param in self.parametrizations ]

    def get_xs_exp(self, **kwargs):
        if len(self.smxs)==0:
            raise RuntimeError('Need to set smxs to a list of xs first')
        mus = self.get_mus_exp(**kwargs)
        return self.prod_list(mus, self.smxs)

    def get_shape_exp(self, **kwargs):
        xss = self.get_xs_exp(**kwargs)
        shape = [ xs/sum(xss) for xs in xss ]
        return shape

    #____________________________________________________________________
    # Full ggH yield parameter; similar to the parametrization, but also
    # including br modifications etc.

    def select_decay_channel(self, dc):
        self.decay_channel = dc
        self.is_dc_specific = True

    def filter_for_decay_channel(self, list_of_pars):
        return [ p for p in list_of_pars if self.decay_channel in p.GetName() ]


    def get_mus_exp_ggH(self, **kwargs):
        if len(self.ggH_yield_parameters)==0:
            self.ggH_yield_parameters = self.get_argset_as_pylist('all_ggH_yieldParameters')
            if self.is_dc_specific:
                self.ggH_yield_parameters = self.filter_for_decay_channel(self.ggH_yield_parameters)
        self.set_kwargs(kwargs)
        return [ yp.getVal() for yp in self.ggH_yield_parameters ]

    def get_xs_exp_ggH(self, **kwargs):
        if len(self.smxs)==0:
            raise RuntimeError('Need to set smxs to a list of xs first')
        mus = self.get_mus_exp_ggH(**kwargs)
        return self.prod_list(mus, self.smxs)

    #____________________________________________________________________
    # xH specific

    def set_smxs_xH(self, smxs):
        self.smxs_xH = smxs

    def get_mus_exp_xH(self, **kwargs):
        if len(self.xH_yield_parameters) == 0:
            self.xH_yield_parameters = self.get_argset_as_pylist('all_xH_yieldParameters')
            if self.is_dc_specific:
                self.xH_yield_parameters = self.filter_for_decay_channel(self.xH_yield_parameters)
        self.set_kwargs(kwargs)
        return [ yp.getVal() for yp in self.xH_yield_parameters ]

    def get_xs_exp_xH(self, **kwargs):
        if len(self.smxs_xH)==0:
            raise RuntimeError('Need to set smxs_xH to a list of xs first')
        mus = self.get_mus_exp_xH(**kwargs)
        return self.prod_list(mus, self.smxs_xH)



class ParametrizationMultiDim(object):
    """docstring for ParametrizationMultiDim"""
    def __init__(self, n_couplings):
        super(ParametrizationMultiDim, self).__init__()
        self.n_couplings = n_couplings

        self.do_linear_terms = True
        self.variations = []
        self.parametrizations = []
        self.parametrize_by_matrix_inversion = False
        
        self.rebinner = None
        self.SM_set = False
        self.binning_set = False

    def evaluate(self, **coupling_vals):
        return [ p(**coupling_vals) if abs(p(**coupling_vals))>1e-12 else 0.0 for p in self.parametrizations ]

    def set_binning(self, binning):
        self.binning = binning
        self.bin_widths = [r-l for l, r in zip(self.binning[:-1], self.binning[1:])]
        self.binning_set = True

    def set_SM(self, **sm_coupling_vals):
        if not(self.SM_set):
            self.xs_theory_SM = self.evaluate(**sm_coupling_vals)
            if not(self.rebinner is None):
                self.xs_exp_SM = self.rebinner.rebin_values(self.xs_theory_SM)
            self.SM_set = True

    def incl_xs(self, **coupling_vals):
        xs_per_GeV = self.evaluate(**coupling_vals)
        return sum([width*xs for xs, width in zip(xs_per_GeV, self.bin_widths)])

    def make_rebinner(self, theory_binning, exp_binning):
        self.rebinner = differentials.integral.Rebinner(
            bin_boundaries_old = theory_binning,
            bin_boundaries_new = exp_binning
            )
        self.bin_widths = [ r-l for l, r in zip(exp_binning[:-1], exp_binning[1:]) ]

    def set_coupling_names(self, *names):
        if len(names) != self.n_couplings:
            raise ValueError(
                'Pass a list of names with length n_couplings={0}'.format(self.n_couplings)
                )
        self.coupling_names = names

    def from_theory_dicts(self, theories):
        for theory in theories:
            self.add_variation(
                theory.crosssection,
                **{ cname : theory[cname] for cname in self.coupling_names }
                )
        self.parametrize()

    def add_variation(self, crosssections, **coupling_vals):
        variation = Variation(crosssections, **coupling_vals)
        logging.debug('Entering variation for {0}'.format(variation))
        self.variations.append(variation)

    def get_coupling_combinations(self):
        # Get list of squared terms and unique combinations (works as expected also for 1 couplings cases)
        self.coupling_combinations = (
            # Squared terms
            [[ coupling, coupling ] for coupling in self.coupling_names]
            +
            # Interference terms
            [list(comb) for comb in itertools.combinations(self.coupling_names, 2)]
            )
        if self.do_linear_terms:
            # Include also singular coupling terms and a constant
            self.coupling_combinations.extend([ [coupling] for coupling in self.coupling_names ])
            self.coupling_combinations.append([])
        self.n_coefficients = len(self.coupling_combinations)

    def get_coupling_matrix(self):
        mat = []
        for variation in self.variations:
            mat.append(self.get_row(variation))
        return mat

    def get_row(self, variation):
        row = []
        for coupling_combination in self.coupling_combinations:
            row.append(
                prod([ variation[c] for c in coupling_combination ])
                )
        return row

    def get_inv_coupling_matrix(self):
        return numpy.linalg.inv(
            numpy.array(self.get_coupling_matrix())
            )

    def get_paramatrization_by_matrix_inversion(self, inv_coupling_matrix, xs):
        xs_column = numpy.array([ [x] for x in xs ])
        coefficients_column = inv_coupling_matrix.dot(xs_column)
        coefficients = [ coefficients_column[i][0] for i in xrange(len(xs)) ]
        return self.get_parabola_for_given_coefficients(coefficients)

    def get_parabola_for_given_coefficients(self, coefficients):
        return ParabolaNDim(coefficients, self.coupling_combinations)

    def cut_variations(self):
        if not self.parametrize_by_matrix_inversion:
            raise RuntimeError(
                'parametrize_by_matrix_inversion is set to False;'
                ' cutting variations makes no sense.'
                )
        n = len(self.variations)
        if n < self.n_coefficients:
            raise RuntimeError(
                'Need {0} variations, but only {1} were given'.format(self.n_coefficients, n)
                )
        elif n > self.n_coefficients:
            self.variations = self.variations[:self.n_coefficients]
            logging.warning('Taking only first {0} variations:'.format(self.n_coefficients))
            for variation in self.variations:
                logging.warning(variation)


    def parametrize(self):
        if not(self.parametrize_by_matrix_inversion):
            raise NotImplementedError(
                'self.parametrize_by_matrix_inversion=False, this is not implemented yet'
                )

        self.parametrizations = []
        self.n_bins = len(self.variations[0].xs)
        if self.parametrize_by_matrix_inversion:
            self.do_parametrize_by_matrix_inversion()

    def do_parametrize_by_matrix_inversion(self):
        self.get_coupling_combinations()
        self.cut_variations()
        inv_coup_mat = self.get_inv_coupling_matrix()
        for i_bin in xrange(self.n_bins):
            xs = [ v.xs[i_bin] for v in self.variations ]
            self.parametrizations.append(
                self.get_paramatrization_by_matrix_inversion(inv_coup_mat, xs)
                )

    def do_parametrize_by_fitting(self):
        raise NotImplementedError('TODO')
        for i_bin in xrange(self.n_bins):
            xs = [ v.xs[i_bin] for v in self.variations ]
            self.parametrizations.append(self.get_fitted_parametrization(c1s, c2s, xs))




class Parametrization2Dim(ParametrizationMultiDim):
    """docstring for Parametrization2Dim"""
    def __init__(self):
        super(Parametrization2Dim, self).__init__(2)
        self.c1_name = 'c1'
        self.c2_name = 'c2'
        self.c1_SM = 1.0
        self.c2_SM = 1.0


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

    # def get_parabola_for_given_coefficients(self, coefficients):
    #     # if self.do_linear_terms:
    #     #     A, B, C, D, E, F = coefficients
    #     #     def parabola(c1, c2):
    #     #         return A*c1**2 + B*c2**2 + C*c1*c2 + D*c1 + E*c2 + F
    #     # else:
    #     #     A, B, C = coefficients
    #     #     def parabola(c1, c2):
    #     #         return A*c1**2 + B*c2**2 + C*c1*c2
    #     # setattr(parabola, 'coefficients', coefficients)
    #     return Parabola(coefficients)

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

