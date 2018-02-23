from physicsModels.MethodHandler import flag_as_method
import physicsModels.RooFactoryInterface as RooFactoryInterface

from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder

@flag_as_method
def GetCouplingOrDefineNew( self, couplingName, possibyExistingCouplings=None, defaultValue=1. ):
    # Look for possibly existing couplings that also define the desired coupling

    if possibyExistingCouplings is None:
        possibyExistingCouplings = [ couplingName ]
    if not couplingName in possibyExistingCouplings:
        possibyExistingCouplings.append( couplingName )

    # ======================================
    # Get list of all components in workspace

    componentArgset   = self.modelBuilder.out.components()
    nComponents       = componentArgset.getSize()
    componentIterator = self.modelBuilder.out.componentIterator()
    allComponentsList = []
    for i in xrange(nComponents):
        element = componentIterator.Next()
        allComponentsList.append( element.GetName() )

    # ======================================
    # Check for possibly existing couplings

    for possibyExistingCoupling in possibyExistingCouplings:
        if possibyExistingCoupling in allComponentsList:
            if self.verbose:
                print 'Found pre-existing coupling \'{0}\', which will be used for \'{1}\''.format( possibyExistingCoupling, couplingName )
            usedCouplingName = possibyExistingCoupling
            break
    else:
        if self.verbose:
            print 'Creating new variable \'{0}\''.format( couplingName )
        self.modelBuilder.doVar( '{0}[{1}]'.format( couplingName, defaultValue ) )
        self.modelBuilder.out.var(couplingName).setConstant(True)
        usedCouplingName = couplingName

    # If the possibly existing variable was actually a function, this will break...
    couplingVar = self.modelBuilder.out.var(usedCouplingName)

    return couplingVar


@flag_as_method
def MakeWidthExpressions( self ):
    self.chapter( 'Starting model.MakeWidthExpressions()' )

    # ======================================
    # Prepare some variables in expressions in WS

    # Make sure couplings are defined
    kappa_t   = self.GetCouplingOrDefineNew( 'kappa_t', [ 'kappat', 'ct' ] )
    kappa_b   = self.GetCouplingOrDefineNew( 'kappa_b', [ 'kappab', 'cb' ] )
    kappa_c   = self.GetCouplingOrDefineNew( 'kappa_c', [ 'kappac', 'cc' ] )
    kappa_glu = self.GetCouplingOrDefineNew(
        'kappa_glu',
        [ 'kappaglu', 'cglu', 'cg', 'kappag', 'kappa_g' ],
        defaultValue = 0.0
        )

    kappa_V   = self.GetCouplingOrDefineNew( 'kappa_V' )
    kappa_W   = self.GetCouplingOrDefineNew( 'kappa_V' )
    kappa_Z   = self.GetCouplingOrDefineNew( 'kappa_V' )

    kappa_tau = self.GetCouplingOrDefineNew( 'kappa_tau' )
    kappa_mu  = self.GetCouplingOrDefineNew( 'kappa_mu' )

    if self.verbose:
        print '\nPrintout of couplings:'
        kappa_t.Print()
        kappa_b.Print()
        kappa_c.Print()
        kappa_W.Print()
        kappa_Z.Print()
        kappa_tau.Print()
        kappa_mu.Print()
        print ''

    # Other needed expressions
    self.modelBuilder.factory_("expr::kappa_mu_expr(\"@0*@1+(1-@0)*@2\", CMS_use_kmu[0], kappa_mu, kappa_tau)")

    # Invisible BR, needed for the other expressions but fixing to 0. now
    self.modelBuilder.doVar("BRinv[0.,0.,1.]")
    self.modelBuilder.out.var("BRinv").setConstant(True)


    # ======================================
    # Application of SMHiggsBuilder

    if not hasattr( self, 'SMH' ): self.SMH = SMHiggsBuilder(self.modelBuilder)

    # SM BR's, called 'SM_BR_(decayChannel)'
    for decayChannel in self.SM_HIGG_DECAYS:
        self.SMH.makeBR( decayChannel )

    # BR uncertainties
    # doBRU = False
    if self.DoBRUncertainties:
        self.SMH.makePartialWidthUncertainties()
        if self.verbose:
            print '\nPrintout of HiggsDecayWidth_UncertaintyScaling_{decaychannel}:'
            for decayChannel in self.SM_HIGG_DECAYS:
                if decayChannel == 'hss': continue
                print 'HiggsDecayWidth_UncertaintyScaling_{0}:'.format( decayChannel )
                self.modelBuilder.out.function( 'HiggsDecayWidth_UncertaintyScaling_{0}'.format( decayChannel ) ).Print()

            for otherVar in [
                    'HiggsDecayWidthTHU_hvv', 'HiggsDecayWidthTHU_hgg', 'HiggsDecayWidthTHU_hll', 'HiggsDecayWidthTHU_hqq', 'HiggsDecayWidthTHU_hzg', 'HiggsDecayWidthTHU_hll', 'HiggsDecayWidthTHU_hgluglu',
                    'param_mt', 'param_alphaS', 'param_mB', 'param_mC',
                    ]:
                try:
                    self.modelBuilder.out.function( otherVar ).Print()
                except AttributeError:
                    self.modelBuilder.out.var( otherVar ).Print()

    else:
        for decayChannel in self.SM_HIGG_DECAYS: 
            self.modelBuilder.factory_( 'HiggsDecayWidth_UncertaintyScaling_%s[1.0]' % decayChannel )

    # makeScaling functions copied from https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/74x-root6/python/LHCHCGModels.py#L355-L358

    self.modelBuilder.factory_( 'expr::kappa_t_plus_kappa_glu( "@0+@1", {0}, {1} )'.format( kappa_t.GetName(), kappa_glu.GetName() ) )
    self.SMH.makeScaling(
        'hgluglu',
        Cb      = kappa_b.GetName(),
        Ctop    = 'kappa_t_plus_kappa_glu'
        )
    self.SMH.makeScaling(
        'hgg',
        Cb      = kappa_b.GetName(),
        Ctop    = kappa_t.GetName(),
        CW      = kappa_W.GetName(),
        Ctau    = kappa_tau.GetName()
        )
    self.SMH.makeScaling(
        'hzg',
        Cb      = kappa_b.GetName(),
        Ctop    = kappa_t.GetName() ,
        CW      = kappa_W.GetName(),
        Ctau    = kappa_tau.GetName()
        )

    ## partial witdhs, normalized to the SM one
    self.modelBuilder.factory_(
        'expr::c7_Gscal_Z('
        '"@0*@0*@1*@2",'
        ' {0}, SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz)'.format( kappa_Z.GetName() )
        )
    self.modelBuilder.factory_(
        'expr::c7_Gscal_W('
        '"@0*@0*@1*@2",'
        ' {0}, SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww)'.format( kappa_W.GetName() )
        )
    self.modelBuilder.factory_(
        'expr::c7_Gscal_tau('
        '"@0*@0*@1*@4+@2*@2*@3*@5",'
        ' {0}, SM_BR_htt, kappa_mu_expr, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm)'.format( kappa_tau.GetName() )
        )
    self.modelBuilder.factory_(
        'expr::c7_Gscal_bottom('
        '"@0*@0 * (@1*@3+@2)",'
        ' {0}, SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb)'.format( kappa_b.GetName() )
        )
    self.modelBuilder.factory_(
        'expr::c7_Gscal_gluon('
        '"  @0  * @1 * @2",'
        ' Scaling_hgluglu, SM_BR_hgluglu, HiggsDecayWidth_UncertaintyScaling_hgluglu)'
        )
    self.modelBuilder.factory_(
        'expr::c7_Gscal_gamma('
        '"@0*@1*@4 + @2*@3*@5",'
        '  Scaling_hgg, SM_BR_hgg, Scaling_hzg, SM_BR_hzg, HiggsDecayWidth_UncertaintyScaling_hgg, HiggsDecayWidth_UncertaintyScaling_hzg)'
        )
    kappa_c_ForBR = kappa_c if kappa_c.GetName() in self.couplings else kappa_t
    self.modelBuilder.factory_(
        'expr::c7_Gscal_top('
        '"@0*@0 * @1*@2",'
        ' {0}, SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc)'.format( kappa_c_ForBR.GetName() )
        )


    ## fix to have all BRs add up to unity
    self.modelBuilder.factory_(
        'sum::c7_SMBRs({0})'.format(
            ','.join([ 'SM_BR_{0}'.format(decayChannel) for decayChannel in self.SM_HIGG_DECAYS ])
            # "SM_BR_"+X for X in "hzz hww htt hmm hcc hbb hss hgluglu hgg hzg".split())
            ))
    self.modelBuilder.out.function("c7_SMBRs").Print("")      

    ## total witdh, normalized to the SM one
    self.modelBuilder.factory_(
        'expr::c7_Gscal_tot('
        '"(@1+@2+@3+@4+@5+@6+@7)/@8/(1-@0)", BRinv, c7_Gscal_Z, c7_Gscal_W, c7_Gscal_tau, c7_Gscal_top, c7_Gscal_bottom, c7_Gscal_gluon, c7_Gscal_gamma, c7_SMBRs)'
        )

    ## BRs, normalized to the SM ones: they scale as (partial/partial_SM) / (total/total_SM) 
    self.modelBuilder.factory_(
        'expr::c7_BRscal_hww('
        '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hww)'.format( kappa_W.GetName() )
        )
    self.modelBuilder.factory_(
        'expr::c7_BRscal_hzz('
        '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzz)'.format( kappa_Z.GetName() )
        )
    self.modelBuilder.factory_(
        'expr::c7_BRscal_htt('
        '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_htt)'.format( kappa_tau.GetName() )
        )
    self.modelBuilder.factory_(
        'expr::c7_BRscal_hmm('
        '"@0*@0*@2/@1", kappa_mu_expr, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hmm)'
        )
    self.modelBuilder.factory_(
        'expr::c7_BRscal_hbb('
        '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hbb)'.format( kappa_b.GetName() )
        )
    self.modelBuilder.factory_(
        'expr::c7_BRscal_hcc('
        '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hcc)'.format( kappa_c_ForBR.GetName() )
        )
    self.modelBuilder.factory_(
        'expr::c7_BRscal_hgg('
        '"@0*@2/@1", Scaling_hgg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgg)'
        )
    self.modelBuilder.factory_(
        'expr::c7_BRscal_hzg('
        '"@0*@2/@1", Scaling_hzg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzg)'
        )
    self.modelBuilder.factory_(
        'expr::c7_BRscal_hgluglu('
        '"@0*@2/@1", Scaling_hgluglu, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgluglu)'
        )

    if self.verbose:
        print '\nPrintout of c7_Gscal_{decayChannel}:'
        for Gscal in [
                'c7_Gscal_Z',
                'c7_Gscal_W',
                'c7_Gscal_tau',
                'c7_Gscal_bottom',
                'c7_Gscal_gluon',
                'c7_Gscal_gamma',
                'c7_Gscal_top',
                'c7_Gscal_tot',
                ]:
            self.modelBuilder.out.function(Gscal).Print()
        print '\nPrintout of c7_BRscal_{decayChannel}:'
        for BRscal in [
                'c7_BRscal_hww',
                'c7_BRscal_hzz',
                'c7_BRscal_htt',
                'c7_BRscal_hmm',
                'c7_BRscal_hbb',
                'c7_BRscal_hcc',
                'c7_BRscal_hgg',
                'c7_BRscal_hzg',
                'c7_BRscal_hgluglu',
                ]:
            self.modelBuilder.out.function(BRscal).Print()


    # ======================================
    # Create the final scaling parameters

    self.modelBuilder.factory_(
        'expr::scalingBR_hggModifier('
        '"@0", c7_BRscal_hgg )'
        )

    self.modelBuilder.factory_(
        'expr::scalingBR_hzzModifier('
        '"@0", c7_BRscal_hzz )'
        )

    # xH modifier
    self.modelBuilder.factory_(
        'expr::scalingBR_xHModifier('
        '"@0*@0", {0} )'.format( kappa_V.GetName() )
        )

    # Set for convenient access in the ws
    self.modelBuilder.out.defineSet( 'BRvariables', ','.join([
        kappa_t.GetName(), kappa_b.GetName(), kappa_c.GetName(), kappa_W.GetName(), kappa_Z.GetName(), kappa_tau.GetName(), kappa_mu.GetName(),
        'scalingBR_hggModifier', 'scalingBR_hzzModifier', 'scalingBR_xHModifier', 'Scaling_hgg', 'Scaling_hzg', 'Scaling_hgluglu',
        ]))

