class Container(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def ListAttributes( self ):
        return [ a for a in dir(self) if not a.startswith('__') and not a == 'ListAttributes' ]