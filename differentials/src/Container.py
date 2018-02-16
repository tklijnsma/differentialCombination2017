class Container(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def list_attributes( self, onlyVariables=False ):
        if onlyVariables:
            return [ a for a in dir(self) if not a.startswith('__') and not a == 'ListAttributes' and not a == 'PrintAttributes' and not callable(a) ]
        else:
            return [ a for a in dir(self) if not a.startswith('__') and not a == 'ListAttributes' and not a == 'PrintAttributes' ]

    def print_attributes( self ):
        variables = self.list_attributes( onlyVariables=True )
        for variable in variables:
            print '{0:20} : {1}'.format( variable, getattr( self, variable ) )