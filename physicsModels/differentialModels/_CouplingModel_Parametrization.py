def makeParametrizationsFromTheory( self ):

    ########################################
    # Handling theory
    ########################################

    # Make sure list lengths do not differ between theories
    #   (Could also check whether binBoundaries matches exactly)
    for key in [ 'couplings', 'binBoundaries', 'crosssection' ]:
        lengthsOfLists = list(set([ len(theory[key]) for theory in self.theories ]))
        if len(lengthsOfLists) > 1:
            print 'WARNING: Found different list lengths for theory attribute \'{0}\': {1}'.format( key, lengthsOfLists )

    nCouplings                = len(self.theories[0]['couplings'])
    couplings                 = self.theories[0]['couplings'].keys()
    nTheoryBins               = len(self.theories[0]['binBoundaries']) - 1
    theoryBinBoundaries       = self.theories[0]['binBoundaries']

    # Get list of squared terms and unique combinations (works as expected also for 1 couplings cases)
    couplingCombinations = []
    couplingCombinations.extend( [ [ coupling, coupling ] for coupling in couplings ] )
    couplingCombinations.extend( [ list(couplingTuple) for couplingTuple in itertools.combinations( couplings, 2 ) ] )
    if self.includeLinearTerms:
        # Include also singular coupling terms and a constant
        couplingCombinations.extend( [ [coupling] for coupling in couplings ] )
        couplingCombinations.append( [] )


    # Number of theories necessary to perform a parametrization
    nComponents = sum(range(nCouplings+1))
    if self.includeLinearTerms:
        nComponents += len(couplings) + 1

    if not len(couplingCombinations) == nComponents:
        raise CouplingModelError( 'ERROR: amount of coupling combinations should be equal to number of expected parameters' )
        sys.exit()

    if len(self.theories) > nComponents:
        print '[FIXME!] Need to choose {0} theories with DIFFERENT COUPLINGS, otherwise matrix is singular!'.format( nComponents )
        print '\n{0} theories supplied but only {1} needed for parametrization; taking only the first {1} theories:'.format(
            len(self.theories), nComponents )
        self.theories = self.theories[:nComponents]
        for theory in self.theories:
            print '    ' + ', '.join([ '{0:10} = {1:10}'.format( cName, cValue ) for cName, cValue in theory['couplings'].iteritems() ])
    elif len(self.theories) < nComponents:
        raise CouplingModelError( 'ERROR: cannot parametrize because number of supplied theories is too small (need {0} theories, found {1})'.format(
            nComponents, len(self.theories) ) )


    # ======================================
    # Find the parametrizations for all nTheoryBins, plus 1 underflow and 1 overflow parametrization

    # First make sure there are RooRealVars for all the couplings in the workspace
    for coupling in couplings:
        self.modelBuilder.doVar(
            '{coupling}[{default},{down},{up}]'.format(
                coupling = coupling,
                default  = self.SMDict['couplings'][coupling],
                down     = self.SMDict['couplings'][coupling] - 1000.0,
                up       = self.SMDict['couplings'][coupling] + 1000.0,
                )
            )

    # Calculate the coupling matrix
    couplingMat = []
    for theory in self.theories:
        if self.includeLinearTerms:
            row = []
            for couplingList in couplingCombinations:
                product = 1.0
                for coupling in couplingList:
                    product *= theory['couplings'][coupling]
                row.append(product)
            couplingMat.append(row)
        else:
            couplingMat.append(
                [ theory['couplings'][c1]*theory['couplings'][c2] for c1, c2 in couplingCombinations ]
                )

    couplingMat = numpy.array( couplingMat )
    if self.verbose:
        print '\nSquared coupling terms:'
        print '  ', couplingCombinations
        print 'Found the following coupling matrix:'
        print couplingMat
    couplingMatInv = numpy.linalg.inv( couplingMat )


    # Create the actual parametrization for each theory bin
    parametrizations = []

    # Underflow (left extrapolation)
    parametrizations.append( [ 0. for i in xrange(nComponents) ] )

    # Add parametrization for each bin
    for iTheoryBin in xrange(nTheoryBins):
        ratios = numpy.array([ [theory['ratios'][iTheoryBin]] for theory in self.theories ])
        parametrization = list(itertools.chain.from_iterable( couplingMatInv.dot( ratios ) ))
        parametrizations.append( parametrization )

    # Use the SAME parametrization for the underflow bin as the first actually defined bin
    # Otherwise the zero contribution pulls down the integral
    #   Sep28: Actually just put 0.0. The underflow bin is considered to have no contribution
    # parametrizations[0] = deepcopy( parametrizations[1] )

    # Overflow (right extrapolation)
    #   Sep28: Stop adding the overflow parametrization.
    #          Just use whatever parametrizations there are in the overflow bin, or return 1.0
    # parametrizations.append( parametrizations[-1][:] )

    if self.verbose:
        print '\nParametrizations per theory bin:'
        for i, parametrization in enumerate(parametrizations):
            print '{0:4}: {1}'.format( i, parametrization )
        print

    # Create a "@number" string per coupling
    argCoupling = { coupling : '@{0}'.format(iCoupling) for iCoupling, coupling in enumerate(couplings) }

    parametrizationNames = []
    for iParametrization, parametrization in enumerate( parametrizations ):

        # Also import it to the workspace
        parametrizationName = 'parametrization{0}'.format( iParametrization )
        parametrizationNames.append( parametrizationName )

        parametrizationString = ''
        for parameter, couplingList in zip( parametrization, couplingCombinations ):
            if len(couplingList) == 2:
                coupling1, coupling2 = couplingList
                parametrizationString += '+{0}*{1}*{2}'.format(
                    parameter, argCoupling[coupling1], argCoupling[coupling2]
                    )
            elif len(couplingList) == 1:
                parametrizationString += '+{0}*{1}'.format(
                    parameter, argCoupling[couplingList[0]]
                    )
            elif len(couplingList) == 0:
                parametrizationString += '+{0}'.format( parameter )

        # Strip off the '+' at the beginning
        parametrizationString = parametrizationString[1:]

        parametrizationExpression = (
            'expr::{name}('
            '"({string})", '
            '{arguments} )'
            ).format(
                name      = parametrizationName,
                string    = parametrizationString,
                arguments = ','.join(couplings),
                )

        if self.verbose > 1:
            print 'Importing parametrization {0}:'.format( iParametrization )
            print parametrizationExpression

        self.modelBuilder.factory_( parametrizationExpression )

    self.modelBuilder.out.defineSet( 'parametrizations', ','.join(parametrizationNames) )

    # Attach to class

    self.parametrizations     = parametrizations

    self.couplingCombinations = couplingCombinations
    self.nComponents          = nComponents

    self.nCouplings           = nCouplings
    self.couplings            = couplings
    self.nTheoryBins          = nTheoryBins
    self.theoryBinBoundaries  = theoryBinBoundaries


    theoryBinBoundarySet = []
    for i, theoryBinBoundary in enumerate(self.theoryBinBoundaries):
        self.modelBuilder.doVar( 'theoryBinBound{0}[{1}]'.format( i, theoryBinBoundary ) )
        theoryBinBoundarySet.append( 'theoryBinBound{0}'.format(i) )
    self.modelBuilder.out.defineSet( 'theoryBinBoundaries', ','.join(theoryBinBoundarySet) )

