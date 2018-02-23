import uuid
from copy import deepcopy

class RooFactoryInterface(object):
    """docstring for RooFactoryInterface"""
    def __init__(self, name=None, variables=None, verbose=False ):
        if name is None:
            self.name = uuid.uuid4()
        else:
            self.name = name

        self.verbose = verbose

        if variables is None or len(variables)==0:
            self.variables = []
        else:
            self.variables = variables

    def add_variable( self, variable_name ):
        self.variables.append(variable_name)

    def make_variable_dict( self ):
        self.variable_dict = {}
        for variable in self.variables:
            self.variable_dict[variable] = len(self.variable_dict)
    
    def copy(self):
        newcopy = deepcopy(self)
        return newcopy

    def plus( self, other ):
        self = plus(self, other, keepName=True)

    def times( self, other ):
        self = times(self, other, keepName=True)


class RooFormulaVar(RooFactoryInterface):
    """docstring for RooFormulaVar"""

    def __init__(self, *args, **kwargs ):
        super(RooFormulaVar, self).__init__( *args, **kwargs )
        self.formula = ''

    def is_empty(self):
        if len(self.formula) == 0:
            return True
        return False

    def get_formula(self):
        return self.formula

    def parse(self):
        self.make_variable_dict()
        formula = self.get_formula()

        # Replace names in formula with numbers
        for variable in self.variables:
            formula = formula.replace( '{' + variable + '}', '@' + str(self.variable_dict[variable]) )

        expr = 'expr::{0}("{1}", {2})'.format(
            self.name,
            formula,
            ','.join(self.variables)
            )

        if self.verbose:
            print 'Parsed expression for \'{0}\':'.format(self.name)
            print expr
        return expr



class RooProduct(RooFactoryInterface):
    """docstring"""

    def __init__(self, *args, **kwargs ):
        super(RooProduct, self).__init__( *args, **kwargs )

    def is_empty(self):
        if len(self.variables) == 0:
            return True
        return False

    def get_formula(self):
        formula = '{' + '}*{'.join(self.variables) + '}'
        return formula

    def parse(self):
        if self.is_empty():
            self.variables.append('one')
        expr = 'prod::{0}({1})'.format( self.name, ','.join(self.variables) )
        if self.verbose:
            print 'Parsed product for \'{0}\':'.format(self.name)
            print expr
        return expr


class RooAddition(RooFactoryInterface):
    """docstring"""

    def __init__(self, *args, **kwargs ):
        super(RooAddition, self).__init__( *args, **kwargs )

    def is_empty(self):
        if len(self.variables) == 0:
            return True
        return False

    def get_formula(self):
        formula = '{' + '}+{'.join(self.variables) + '}'
        return formula

    def parse(self):
        expr = 'sum::{0}({1})'.format( self.name, ','.join(self.variables) )
        if verbose:
            print 'Parsed addition for \'{0}\':'.format(self.name)
            print expr
        return expr




def plus( A, B, keepName=False ):
    if keepName:
        name = A.name
    else:
        A.name + '_plus_' + B.name

    if isinstance( A, RooAddition ) and isinstance( B, RooAddition ):
        C = RooAddition(name)
        C.variables = A.variables + B.variables
    else:
        C = RooFormulaVar(name)
        C.variables = A.variables + B.variables

        if A.is_empty() and B.is_empty():
            C.formula = ''
        elif A.is_empty():
            C.formula = B.get_formula()
        elif B.is_empty():
            C.formula = A.get_formula()
        else:
            C.formula = '({0})+({1})'.format( A.get_formula(), B.get_formula() )


def times( A, B, keepName=False ):
    if keepName:
        name = A.name
    else:
        A.name + '_plus_' + B.name

    if isinstance( A, RooProduct ) and isinstance( B, RooProduct ):
        C = RooProduct( A.name + '_times_' + B.name )
        C.variables = A.variables + B.variables
    else:
        C = RooFormulaVar( A.name + '_times_' + B.name )
        C.variables = A.variables + B.variables

        if A.is_empty() and B.is_empty():
            C.formula = ''
        elif A.is_empty():
            C.formula = B.get_formula()
        elif B.is_empty():
            C.formula = A.get_formula()
        else:
            C.formula = '({0})*({1})'.format( A.get_formula(), B.get_formula() )


########################################
# Specific for the polynomial parametrization
########################################

class RooParametrization(RooFormulaVar):
    """docstring for RooParametrization"""
    def __init__(self, *args, **kwargs ):
        super(RooParametrization, self).__init__( *args, **kwargs )
        self.terms = []

    def add_term(self, coefficient, couplings):
        self.terms.append( RooParametrizationTerm(coefficient, couplings) )

    def get_formula(self):
        return '+'.join([ t.get_formula() for t in self.terms ])


class RooParametrizationTerm(object):
    """docstring for RooParametrizationTerm"""
    def __init__(self, coefficient, couplings):
        self.coefficient = coefficient
        self.couplings   = couplings
        
    def get_formula(self):
        return '*'.join([ str(self.coefficient) ] + [ '{'+c+'}' for c in self.couplings ] )


def get_average_coefficients( weights, rooParametrizations, verbose=True ):
    
    def p(text):
        if verbose: print text

    n_terms = len(rooParametrizations[0].terms)

    average_coefficients = []
    p('\nCalculating average coefficients out of {0} parametrizations'.format(len(rooParametrizations)))
    for i_term in xrange(n_terms):
        p('  Coefficient {0}'.format(i_term))
        p('    {0:18} | {1:15} | {2:15} | {3:15}'.format('param. name', 'coefficient', 'weight', 'coeff*weight'))
        average_coefficient = 0.0
        for weight, rooParametrization in zip(weights, rooParametrizations):
            coefficient = rooParametrization.terms[i_term].coefficient
            average_coefficient += weight * coefficient
            p('    {0:18} | {1:15.4f} | {2:15.4f} | {3:15.4f}'.format(rooParametrization.name, coefficient, weight, weight*coefficient ))
        p( ' '*60 + '-'*16 + ' +' )
        p( ' '*61 + '{0:15.4f}'.format(average_coefficient) )
        average_coefficients.append(average_coefficient)

    return average_coefficients










