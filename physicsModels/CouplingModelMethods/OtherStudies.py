from physicsModels.MethodHandler import flag_as_method

@flag_as_method
def MakeRatioOfBRsExpressions( self ):
    self.chapter( 'Starting model.MakeRatioOfBRsExpressions()' )

    # Define the two floating pars
    self.modelBuilder.doVar( 'hgg_ratioBRmodifier[1.0,-2.0,4.0]' )
    self.modelBuilder.doVar( 'ratio_BR_hgg_hzz[0.086,0.0,0.5]' )

    # Load spline for HZZ BR (function of MH)
    self.SMH = SMHiggsBuilder(self.modelBuilder)
    datadir = os.environ['CMSSW_BASE'] + '/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg'
    self.SMH.textToSpline( 'BR_hzz', os.path.join( datadir, 'sm/br/BR4.txt' ), ycol=11 );
    spline_SMBR_hzz = self.modelBuilder.out.function('BR_hzz')

    # Seems to work:
    # self.modelBuilder.out.var('MH').setVal(125.09)
    # spline_SMBR_hzz.Print()

    # Load spling for Hgg BR
    spline_SMBR_hgg = self.modelBuilder.out.function('fbr_13TeV')


    # ======================================
    # Create 'hzz_BRmodifier', as a function of 'hgg_BRmodifier' and 'ratio_BR_hgg_hzz'

    hzz_modifier_expr = 'expr::{name}("{formula}",{commaSeparatedParameters})'.format(
        name                     = 'hzz_ratioBRmodifier',
        formula                  = '@0*(@1/@2)/@3',
        commaSeparatedParameters = ','.join([
            'hgg_ratioBRmodifier',
            spline_SMBR_hgg.GetName(),
            spline_SMBR_hzz.GetName(),
            'ratio_BR_hgg_hzz'
            ])
        )
    print 'Processing expr:'
    print '    ',hzz_modifier_expr
    self.modelBuilder.factory_( hzz_modifier_expr )


# @flag_as_method
def MakeTotalXSExpressions_old( self ):
    self.chapter( 'Starting model.MakeTotalXSExpressions()' )

    self.modelBuilder.doVar( 'r_totalXS[1.0,0.0,3.0]' )

    thingsToSum = []
    for iTheoryBin in xrange( self.nTheoryBins ):
        binWidth = self.theoryBinBoundaries[iTheoryBin+1] - self.theoryBinBoundaries[iTheoryBin]
        self.modelBuilder.factory_(
            'expr::IncXS{0}( "{1}*{2}*@0", parametrization{0} )'.format(
                iTheoryBin, binWidth, self.SMXS[iTheoryBin]
                )
            )
        thingsToSum.append( 'IncXS{0}'.format(iTheoryBin) )
    self.modelBuilder.factory_('sum::totalXS( {0} )'.format(','.join(thingsToSum)))

    # Make sure couplings are at SM
    for coupling in self.couplings:
        self.modelBuilder.out.var(coupling).setVal(self.SMDict['couplings'][coupling])

    # Set SM spectrum integral to the SM xs
    self.modelBuilder.doVar( 'totalXS_SM[{0}]'.format(self.modelBuilder.out.function('totalXS').getVal()))
    self.modelBuilder.out.var('totalXS_SM').setConstant(True)

    # Build the modifier
    self.modelBuilder.factory_( 'expr::totalXSmodifier( "@0*@1/@2", r_totalXS, totalXS_SM, totalXS )' )

    if self.FitOnlyNormalization:
        # Make also modifier if only fitting normalization
        self.modelBuilder.out.var('r_totalXS').setVal(1.0)
        self.modelBuilder.out.var('r_totalXS').setConstant(True)
        self.modelBuilder.factory_( 'expr::globalTotalXSmodifier( "@0/@1", totalXS, totalXS_SM )' )
        if self.verbose:
            self.modelBuilder.out.function('globalTotalXSmodifier').Print()

