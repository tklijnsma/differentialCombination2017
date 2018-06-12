from physicsModels.MethodHandler import flag_as_method
import physicsModels.RooFactoryInterface as RooFactoryInterface

from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder

import ROOT
import sys
import math


@flag_as_method
def list_ws_components(self):
    """Returns list of all components in the workspace self.modelBuilder.out"""
    component_iterator = self.modelBuilder.out.componentIterator()
    components = []
    for i in xrange(self.modelBuilder.out.components().getSize()):
        element = component_iterator.Next()
        components.append(element.GetName())
    return components

@flag_as_method
def implement_freely_floating_BRs(self):
    self.chapter('Defining freely floating BRs')
    self.modelBuilder.doVar('brmodifier_hgg[1.0,-1000.0,1000.0]')
    self.modelBuilder.doVar('brmodifier_hzz[1.0,-1000.0,1000.0]')
    
    if self.constrain_ratio_bb_ZZ:
        self.implement_constrain_ratio_bb_ZZ()
        self.modelBuilder.out.function('brmodifier_hbb').Print()
    else:
        self.modelBuilder.doVar('brmodifier_hbb[1.0,-1000.0,1000.0]')
        self.modelBuilder.out.var('brmodifier_hbb').Print()
    
    self.modelBuilder.out.var('brmodifier_hgg').Print()
    self.modelBuilder.out.var('brmodifier_hzz').Print()
    

@flag_as_method
def implement_constrain_ratio_bb_ZZ(self):
    # Apply a constraint on the ratio bb/ZZ based on HIG-17-031
    mean = 0.84
    up = 0.38
    down = -0.27
    symm = 0.5*(abs(down)+abs(up))
    symm_perc = symm / mean

    self.modelBuilder.doVar('delta_ratio_bb_ZZ[0.0,-10.0,10.0]')
    self.DC.systs.append(
        ( 'delta_ratio_bb_ZZ', False, 'param', ['0', '1'], [] )
        )
    self.modelBuilder.out.var('delta_ratio_bb_ZZ').Print()

    scaler = ROOT.ProcessNormalization('modifier_bb_ZZ', '')
    scaler.addLogNormal(
        math.exp(symm_perc),
        self.modelBuilder.out.var('delta_ratio_bb_ZZ')
        )
    self.modelBuilder.out._import(scaler)
    self.modelBuilder.out.function(scaler.GetName()).Print()

    self.modelBuilder.factory_('prod::brmodifier_hbb(modifier_bb_ZZ, brmodifier_hzz)')


    # scaler = ROOT.ProcessNormalization('ratio_bb_ZZ', '')
    # scaler.addLogNormal(
    #     math.exp(symm_perc),
    #     self.modelBuilder.out.var('delta_ratio_bb_ZZ')
    #     )
    # self.modelBuilder.out._import(scaler)
    # self.modelBuilder.out.function('ratio_bb_ZZ').Print()

    # self.modelBuilder.factory_('prod::brmodifier_hbb(ratio_bb_ZZ, brmodifier_hzz)')


@flag_as_method
def implement_scaling_BRs(self):
    self.chapter('Applying partial width model')
    if not hasattr(self, 'PWB'): self.PWB = PartialWidthBuilder(self)
    self.PWB.get_width_expressions_scaling()

@flag_as_method
def implement_nonscaling_BR_uncertainties(self):
    self.chapter('Applying partial width model')
    if not hasattr(self, 'PWB'): self.PWB = PartialWidthBuilder(self)
    self.PWB.get_width_expressions_nonscaling()


class PartialWidthBuilder(object):
    """docstring for PartialWidthBuilder"""

    SM_HIGG_DECAYS = [ "hww", "hzz", "hgg", "htt", "hbb", 'hzg', 'hmm', 'hcc', 'hgluglu' ]
    SM_HIGG_DECAYS += ['hss']

    def __init__(self, model):
        super(PartialWidthBuilder, self).__init__()
        self.model = model
        if not hasattr(self.model, 'SMH'): self.model.SMH = SMHiggsBuilder(self.model.modelBuilder)
        self.SMH = self.model.SMH # Convenience pointer
        self.modelBuilder = self.model.modelBuilder # Convenience pointer
        self.verbose = self.model.verbose

        self.defined_partial_width_uncertainties = False
        self.defined_SM_BRs = False
        self.defined_BRinv = False
        self.defined_nuis_pars = False

    def get_coupling(self, name, alternatives=None, default_value=1.0 ):
        """
        Gets a RooRealVar coupling from the ws if it exists, or defines a new one if it doesn't.
        Can take various alternatives and will use the first existing alternative that is in the ws
        """
        if alternatives is None: alternatives = []
        if not name in alternatives: alternatives.append(name)

        all_components = self.model.list_ws_components()

        for alternative in alternatives:
            if alternative in all_components:
                print 'Getting coupling for \'{0}\': Using pre-existing variable \'{1}\''.format(name, alternative)
                name = alternative
                break
        else:
            factory_string = '{0}[{1}]'.format(name, default_value)
            print 'Getting coupling for \'{0}\': No alternative, defining new variable {1}'.format(name, factory_string)
            self.model.modelBuilder.doVar(factory_string)

        variable = self.model.modelBuilder.out.var(name)
        variable.setVal(default_value) # Also make sure it is set at SM value
        return variable

    def get_partial_width_uncertainties(self):
        if self.defined_partial_width_uncertainties:
            print 'Partial widths were already defined'
        else:
            # Create expressions for the uncertainties
            self.SMH.makePartialWidthUncertainties()
            self.test_printout_width_uncertainties()
            self.defined_partial_width_uncertainties = True

    def fix_partial_width_uncertainties(self):
        """This disables the partial width uncertainties"""
        if self.defined_partial_width_uncertainties:
            print 'Partial widths were already defined'
        else:
            for decayChannel in self.SM_HIGG_DECAYS: 
                self.modelBuilder.factory_( 'HiggsDecayWidth_UncertaintyScaling_%s[1.0]' % decayChannel )
            self.defined_partial_width_uncertainties = True

    def get_SM_BRs(self):
        if self.defined_SM_BRs:
            print 'SM BRs were already defined'
        else:
            # Get the SM BRs, scaling only with mH (splines)
            # Names will be SM_BR_hgg / SM_BR_hzz / etc.
            for decay_channel in self.SM_HIGG_DECAYS:
                self.SMH.makeBR(decay_channel)
            self.defined_SM_BRs = True

    def get_BRinv(self):
        if self.defined_BRinv:
            print 'BRinv was already defined'
        else:
            # Invisible BR, needed for the other expressions but fixing to 0. now
            # self.modelBuilder.doVar("BRinv[0.,0.,1.]")
            self.modelBuilder.doVar("BRinv[0.]")
            self.modelBuilder.out.var("BRinv").setConstant(True)
            self.defined_BRinv = True


    def implement_nuisance_parameters(self):
        if self.defined_nuis_pars:
            print 'Nuisance parameters were already defined'
        else:
            for syst_name in [
                    'param_alphaS',
                    'param_mB',
                    'param_mC',
                    'param_mt',
                    'HiggsDecayWidthTHU_hqq',
                    'HiggsDecayWidthTHU_hvv',
                    'HiggsDecayWidthTHU_hll',
                    'HiggsDecayWidthTHU_hgg',
                    'HiggsDecayWidthTHU_hzg',
                    'HiggsDecayWidthTHU_hgluglu',
                    ]:

                self.model.DC.systs.append(
                    ( syst_name, False, 'param', ['0', '1'], [] )
                    )
            self.defined_nuis_pars = True

    def get_width_expressions_nonscaling(self):
        self.implement_nuisance_parameters()
        self.get_SM_BRs()
        self.get_partial_width_uncertainties()
        self.get_BRinv()

        # 'Partial widths' (actually scaled BRs, but all will be made relative anyway)
        # Scaling with kappa's removed here!
        self.modelBuilder.factory_(
            'expr::c7_Gscal_Z('
            '"@0*@1",'
            'SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz)'
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_W('
            '"@0*@1",'
            'SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww)'
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_tau('
            '"@0*@2+@1*@3",'
            'SM_BR_htt, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm)'
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_bottom('
            '"@0*@1+@2",'
            'SM_BR_hbb, HiggsDecayWidth_UncertaintyScaling_hbb, SM_BR_hss)'
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_gluon('
            '"@0*@1",'
            'SM_BR_hgluglu, HiggsDecayWidth_UncertaintyScaling_hgluglu)'
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_gamma('
            '"@0*@1 + @2*@3",'
            'SM_BR_hgg, HiggsDecayWidth_UncertaintyScaling_hgg, SM_BR_hzg, HiggsDecayWidth_UncertaintyScaling_hzg)'
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_top('
            '"@0*@1",'
            'SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc)'
            )

        # With this scaling things don't sum up to 1.0, so need to normalize by the sum
        self.modelBuilder.factory_(
            'sum::c7_SMBRs({0})'.format(
                ','.join([ 'SM_BR_{0}'.format(decay_channel) for decay_channel in self.SM_HIGG_DECAYS ])
                ))
        # 'Total width' (only relative, not really), including possible BRinv, 
        # This is 1.0 by defintion if BRinv is 0.0
        self.modelBuilder.factory_(
            'expr::c7_Gscal_tot('
            '"(@1+@2+@3+@4+@5+@6+@7)/@8/(1-@0)",'
            'BRinv, '
            'c7_Gscal_Z, c7_Gscal_W, c7_Gscal_tau, c7_Gscal_top, c7_Gscal_bottom, c7_Gscal_gluon, c7_Gscal_gamma, '
            'c7_SMBRs'
            ')'
            )

        # Non-scaling BRs, normalized to the SM ones: they scale as (partial/partial_SM) / (total/total_SM) 
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hww('
            '"@0/@1", HiggsDecayWidth_UncertaintyScaling_hww, c7_Gscal_tot)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hzz('
            '"@0/@1", HiggsDecayWidth_UncertaintyScaling_hzz, c7_Gscal_tot)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_htt('
            '"@0/@1", HiggsDecayWidth_UncertaintyScaling_htt, c7_Gscal_tot)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hmm('
            '"@0/@1", HiggsDecayWidth_UncertaintyScaling_hmm, c7_Gscal_tot)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hbb('
            '"@0/@1", HiggsDecayWidth_UncertaintyScaling_hbb, c7_Gscal_tot)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hcc('
            '"@0/@1", HiggsDecayWidth_UncertaintyScaling_hcc, c7_Gscal_tot)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hgg('
            '"@0/@1", HiggsDecayWidth_UncertaintyScaling_hgg, c7_Gscal_tot)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hzg('
            '"@0/@1", HiggsDecayWidth_UncertaintyScaling_hzg, c7_Gscal_tot)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hgluglu('
            '"@0/@1", HiggsDecayWidth_UncertaintyScaling_hgluglu, c7_Gscal_tot)'
            )
        self.test_printout_widths()

        # Now get scalers based on SM:
        self.modelBuilder.out.var('MH').setVal(125.)
        self.modelBuilder.factory_(
            'expr::brmodifier_hgg('
            '"@0/{0}", c7_BRscal_hgg)'
            .format( self.modelBuilder.out.function('c7_BRscal_hgg').getVal() )
            )
        self.modelBuilder.factory_(
            'expr::brmodifier_hzz('
            '"@0/{0}", c7_BRscal_hzz)'
            .format( self.modelBuilder.out.function('c7_BRscal_hzz').getVal() )
            )
        self.modelBuilder.factory_(
            'expr::brmodifier_hbb('
            '"@0/{0}", c7_BRscal_hbb)'
            .format( self.modelBuilder.out.function('c7_BRscal_hbb').getVal() )
            )


    def test_printout_width_uncertainties(self):
        # Test printouts
        if self.model.verbose:
            print '\nPrintout of HiggsDecayWidth_UncertaintyScaling_{decaychannel}:'

            for decay_channel in self.SM_HIGG_DECAYS:
                if decay_channel == 'hss': continue
                print 'HiggsDecayWidth_UncertaintyScaling_{0}:'.format(decay_channel)
                self.model.modelBuilder.out.function('HiggsDecayWidth_UncertaintyScaling_{0}'.format(decay_channel)).Print()

            for other_var in [
                    'HiggsDecayWidthTHU_hvv', 'HiggsDecayWidthTHU_hgg', 'HiggsDecayWidthTHU_hll', 'HiggsDecayWidthTHU_hqq', 'HiggsDecayWidthTHU_hzg', 'HiggsDecayWidthTHU_hll', 'HiggsDecayWidthTHU_hgluglu',
                    'param_mt', 'param_alphaS', 'param_mB', 'param_mC',
                    ]:
                try:
                    self.model.modelBuilder.out.function(other_var).Print()
                except AttributeError:
                    self.model.modelBuilder.out.var(other_var).Print()

    def test_printout_widths(self):
        if self.model.verbose:
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

    # Scaling widths
    def import_all_couplings(self):
        # Make sure couplings are defined
        self.kappa_t   = self.get_coupling('kappa_t', [ 'kappat', 'ct' ])
        self.kappa_b   = self.get_coupling('kappa_b', [ 'kappab', 'cb' ])
        self.kappa_c   = self.get_coupling('kappa_c', [ 'kappac', 'cc' ])

        self.cg   = self.get_coupling('cg', default_value=0.0) # Make sure cg is defined, no alternatives
        self.modelBuilder.factory_('expr::kappa_glu("12.*@0", cg)') # 12 times cg, as per the theory ref.
        self.kappa_glu = self.modelBuilder.out.function('kappa_glu')

        self.kappa_V   = self.get_coupling('kappa_V')
        self.kappa_W   = self.get_coupling('kappa_V')
        self.kappa_Z   = self.get_coupling('kappa_V')

        self.kappa_tau = self.get_coupling('kappa_tau')
        self.kappa_mu  = self.get_coupling('kappa_mu')

        # Other needed expressions
        self.modelBuilder.factory_("expr::kappa_mu_expr(\"@0*@1+(1-@0)*@2\", CMS_use_kmu[0], kappa_mu, kappa_tau)")

        if self.verbose:
            print '\nPrintout of couplings:'
            self.kappa_t.Print()
            self.kappa_b.Print()
            self.kappa_c.Print()
            self.kappa_W.Print()
            self.kappa_Z.Print()
            self.kappa_tau.Print()
            self.kappa_mu.Print()
            self.modelBuilder.out.function('kappa_mu_expr').Print()
            print ''


    def get_width_expressions_scaling(self):
        self.import_all_couplings()

        if self.model.do_BR_uncertainties:
            self.implement_nuisance_parameters()
            self.get_partial_width_uncertainties()
        else:
            self.fix_partial_width_uncertainties()

        self.get_SM_BRs()
        self.get_BRinv()

        # Define effective coupling of Higgs to gluons, use it in the scaling
        self.modelBuilder.factory_('sum::kappa_t_plus_kappa_glu({0}, {1})'.format(self.kappa_t.GetName(), self.kappa_glu.GetName()))
        self.SMH.makeScaling(
            'hgluglu',
            Cb      = self.kappa_b.GetName(),
            Ctop    = 'kappa_t_plus_kappa_glu'
            )
        self.SMH.makeScaling(
            'hgg',
            Cb      = self.kappa_b.GetName(),
            Ctop    = self.kappa_t.GetName(),
            CW      = self.kappa_W.GetName(),
            Ctau    = self.kappa_tau.GetName()
            )
        self.SMH.makeScaling(
            'hzg',
            Cb      = self.kappa_b.GetName(),
            Ctop    = self.kappa_t.GetName() ,
            CW      = self.kappa_W.GetName(),
            Ctau    = self.kappa_tau.GetName()
            )

        ## partial witdhs, normalized to the SM one
        self.modelBuilder.factory_(
            'expr::c7_Gscal_Z('
            '"@0*@0*@1*@2",'
            ' {0}, SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz)'.format(self.kappa_Z.GetName())
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_W('
            '"@0*@0*@1*@2",'
            ' {0}, SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww)'.format(self.kappa_W.GetName())
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_tau('
            '"@0*@0*@1*@4+@2*@2*@3*@5",'
            ' {0}, SM_BR_htt, kappa_mu_expr, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm)'.format(self.kappa_tau.GetName())
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_bottom('
            '"@0*@0 * (@1*@3+@2)",'
            ' {0}, SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb)'.format(self.kappa_b.GetName())
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
        kappa_c_ForBR = self.kappa_c if self.kappa_c.GetName() in self.model.couplings else self.kappa_t
        self.modelBuilder.factory_(
            'expr::c7_Gscal_top('
            '"@0*@0 * @1*@2",'
            ' {0}, SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc)'.format(kappa_c_ForBR.GetName())
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
            '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hww)'.format(self.kappa_W.GetName())
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hzz('
            '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzz)'.format(self.kappa_Z.GetName())
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_htt('
            '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_htt)'.format(self.kappa_tau.GetName())
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hmm('
            '"@0*@0*@2/@1", kappa_mu_expr, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hmm)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hbb('
            '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hbb)'.format(self.kappa_b.GetName())
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hcc('
            '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hcc)'.format(kappa_c_ForBR.GetName())
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

        # Now get scalers based on SM:
        self.modelBuilder.out.var('MH').setVal(125.)
        self.modelBuilder.factory_(
            'expr::brmodifier_hgg('
            '"@0/{0}", c7_BRscal_hgg)'
            .format( self.modelBuilder.out.function('c7_BRscal_hgg').getVal() )
            )
        self.modelBuilder.factory_(
            'expr::brmodifier_hzz('
            '"@0/{0}", c7_BRscal_hzz)'
            .format( self.modelBuilder.out.function('c7_BRscal_hzz').getVal() )
            )
        self.modelBuilder.factory_(
            'expr::brmodifier_hbb('
            '"@0/{0}", c7_BRscal_hbb)'
            .format( self.modelBuilder.out.function('c7_BRscal_hbb').getVal() )
            )
        
        self.test_printout_widths()
        print '\nPrintout of brmodifiers'
        self.modelBuilder.out.function('brmodifier_hgg').Print()
        self.modelBuilder.out.function('brmodifier_hzz').Print()
        self.modelBuilder.out.function('brmodifier_hbb').Print()


#____________________________________________________________________
# OLD METHODS

# @flag_as_method
def GetCouplingOrDefineNew( self, couplingName, possibyExistingCouplings=None, defaultValue=1. ):
    # Look for possibly existing couplings that also define the desired coupling

    if possibyExistingCouplings is None:
        possibyExistingCouplings = [ couplingName ]
    if not couplingName in possibyExistingCouplings:
        possibyExistingCouplings.append( couplingName )

    # ======================================
    # Get list of all components in workspace


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


# @flag_as_method
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
    # self.modelBuilder.doVar("BRinv[0.,0.,1.]")
    self.modelBuilder.doVar("BRinv[0.]")
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

    expr_kappa_t_plus_kappa_glu = 'expr::kappa_t_plus_kappa_glu( "@0+12.*@1", {0}, {1} )'.format( kappa_t.GetName(), kappa_glu.GetName() )
    print '\n[WARNING] Adding factor 12. in front of kappa_glu in: {0}'.format(expr_kappa_t_plus_kappa_glu)
    self.modelBuilder.factory_(expr_kappa_t_plus_kappa_glu)
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

    self.modelBuilder.factory_(
        'expr::scalingBR_hbbModifier('
        '"@0", c7_BRscal_hbb )'
        )

    # xH modifier
    self.modelBuilder.factory_(
        'expr::scalingBR_xHModifier('
        '"@0*@0", {0} )'.format( kappa_V.GetName() )
        )

    # Set for convenient access in the ws
    self.modelBuilder.out.defineSet( 'BRvariables', ','.join([
        kappa_t.GetName(), kappa_b.GetName(), kappa_c.GetName(), kappa_W.GetName(), kappa_Z.GetName(), kappa_tau.GetName(), kappa_mu.GetName(),
        'scalingBR_hggModifier', 'scalingBR_hzzModifier', 'scalingBR_hbbModifier', 'scalingBR_xHModifier', 'Scaling_hgg', 'Scaling_hzg', 'Scaling_hgluglu',
        ]))

