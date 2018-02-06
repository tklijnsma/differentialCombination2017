#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys, re
from os.path import *
from glob import glob
from copy import deepcopy

from OptionHandler import flag_as_option

sys.path.append('src')
import Commands
import PhysicsCommands
import TheoryCommands
import LatestPaths
import LatestPathsGetters
import LatestBinning
from Container import Container
import PlotCommands
from differentialTools import *

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

# Test of new bins for hzz and hbb
@flag_as_option
def RenumberHzzProcesses_Jan24(args):
    MergeHGGWDatacards.RenumberProcessesHZZ_Aug21(
        LatestPaths.card_hzz_ggHxH_PTH_newBins_unprocessed,
        )
@flag_as_option
def CombineCards_Jan24_hzz_hbb(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_newBins_hzz_hbb_ggHxH_{0}.txt'.format(datestr),
        'hzz=' + LatestPaths.card_hzz_ggHxH_PTH_newBins,
        'hbb=' + LatestPaths.card_hbb_ggHxH_PTH
        )

#____________________________________________________________________
# Helpers

def get_all_vars_except_POIs_and_envelop(args, postfit, verbose=False):
    with Commands.OpenRootFile(postfit) as postfit_fp:
        w = postfit_fp.Get('w')

    floating_regexps = [
        r'env_pdf_\d+_13TeV_',
        r'n_exp_binhgg_.*_bkg_mass',
        r'shapeBkg_bkg_mass_.*_norm',
        r'CMS_.*_mass',
        ]

    nuisances = Commands.ListSet(w, 'nuisances')
    POIs = Commands.ListSet(w)

    freezing_variables = []
    floating_variables = []
    constant_variables = []

    var_set = w.allVars()
    var_iterable = var_set.createIterator()
    for i_var in xrange(var_set.getSize()):
        variable = var_iterable.Next()
        var_name = variable.GetName()
        if variable.isConstant():
            # No need to freeze
            constant_variables.append(var_name)
            continue
        elif var_name in nuisances or var_name in POIs:
            # Should be handled by 'rgx{.*}'
            continue
        else:
            for pat in floating_regexps:
                if re.match(pat, var_name):
                    floating_variables.append(var_name)
                    break
            else:
                freezing_variables.append(var_name)

    freezing_variables = ['rgx{.*}'] + freezing_variables
    floating_variables.extend(POIs)
    
    if verbose:
        print '\nList of variables to freeze:'
        print ','.join(freezing_variables)
        print '\nList of variables to float:'
        print ','.join(floating_variables)
        print '\nList of variables that are constant:'
        print ','.join(constant_variables)

    if len(','.join(freezing_variables)) > 10000:
        raise ValueError(
            'The string length for the list of frozen variables exceeds 10000, '
            'which will crash when loaded into combine.'
            )
    return freezing_variables

def make_postfit(args, ws, observable_tag=None):
    tag = 'POSTFIT_' + basename(ws).replace('.root','')
    if not(observable_tag is None):
        tag += '_' + observable_tag
    if args.asimov:
        tag += '_asimov'
    cmd = [
        'combineTool.py',
        abspath(ws),
        '-n _{0}'.format(tag),
        '-M MultiDimFit',
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        '--floatOtherPOIs=1',
        # '-P "{0}"'.format( POI ),
        # '--setPhysicsModelParameterRanges "{0}"={1:.3f},{2:.3f} '.format( POI, POIRange[0], POIRange[1] ),
        # '--setPhysicsModelParameters {0}'.format( ','.join([ iterPOI + '=1.0' for iterPOI in allPOIs ]) ),
        '-m 125.00',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--saveWorkspace',
        '--job-mode psi --task-name {0} --sub-opts=\'-q short.q\' '.format(tag),
        ]
    if args.asimov:
        cmd.append( '-t -1' )

    if not 't3' in os.environ['HOSTNAME']:
        raise NotImplementedError('\'t3\' not in $HOSTNAME; Is this executed from the T3?')

    with Commands.EnterDirectory('out/postfitWss_{0}'.format(datestr)):
        Commands.executeCommand(cmd)

class GenericDifferentialObservableScan(object):
    """docstring for GenericDifferentialObservableScan"""

    setPhysicsModelParameters = True
    nPoints = 45
    nPointsPerJob = 3
    POIRange = [ -1.0, 4.0 ]

    def __init__(self, args, obs_name, ws=None):
        self.args = args
        self.obs_name = obs_name
        self.decay_channel = get_decay_channel_tag(args)
        if not(ws is None):
            self.set_ws(ws)

        self.extraOptions = []
        self.physicsModelParameterRanges = []

        self.jobDirectory = 'out/Scan_{0}_{1}_{2}'.format(self.obs_name, datestr, self.decay_channel)
        if args.statonly:
            self.nPointsPerJob = 6
            self.jobDirectory += '_statonly'
            variables_to_freeze = get_all_vars_except_POIs_and_envelop(self.args, self.ws)
            self.extraOptions.extend([
                '--snapshotName MultiDimFit',
                '--skipInitialFit',
                '--freezeNuisances {0}'.format(','.join(variables_to_freeze))
                ])
            self.setPhysicsModelParameters = False

        if args.hzz or args.hbb:
            self.nPointsPerJob = self.nPoints

    def set_ws(self, ws):
        self.ws = ws

    def run(self):
        Commands.BasicCombineTool(
            self.ws,
            POIpattern    = '*',
            POIRange      = self.POIRange,
            nPoints       = self.nPoints,
            nPointsPerJob = self.nPointsPerJob,
            jobDirectory  = self.jobDirectory,
            queue         = 'short.q',
            asimov        = self.args.asimov,
            extraOptions  = self.extraOptions,
            setPhysicsModelParameters = self.setPhysicsModelParameters,
            physicsModelParameterRanges = self.physicsModelParameterRanges
            )

def set_all_false(args):
    args.hgg = False
    args.hzz = False
    args.hbb = False
    args.combWithHbb = False

#____________________________________________________________________
@flag_as_option
def postfit_all(real_args):
    args = deepcopy(real_args)

    for decay_channel in [
            'hgg',
            'hzz',
            # Workspaces are not working yet; Wait for Vs input
            # 'hbb',
            # 'combWithHbb',
            'combination'
            ]:
        set_all_false(args)
        setattr(args, decay_channel, True)
        print '\nMaking postfits for {0}'.format(decay_channel)

        try:
            pth_smH_postfit(args)
        except NotImplementedError:
            pass

        try:
            pth_ggH_postfit(args)
        except NotImplementedError:
            pass

        try:
            ptjet_postfit(args)
        except NotImplementedError:
            pass

        try:
            rapidity_postfit(args)
        except NotImplementedError:
            pass

        try:
            njets_postfit(args)
        except NotImplementedError:
            pass

@flag_as_option
def scan_all(real_args):
    args = deepcopy(real_args)

    for decay_channel in [
            'hgg',
            'hzz',
            # Workspaces are not working yet; Wait for Vs input
            # 'hbb',
            # 'combWithHbb',
            'combination'
            ]:
        set_all_false(args)
        setattr(args, decay_channel, True)
        print '\nSubmitting scans for {0}'.format(decay_channel)

        try:
            pth_smH_scan(args)
        except NotImplementedError:
            pass

        try:
            pth_ggH_scan(args)
        except NotImplementedError:
            pass

        try:
            ptjet_scan(args)
        except NotImplementedError:
            pass

        try:
            rapidity_scan(args)
        except NotImplementedError:
            pass

        try:
            njets_scan(args)
        except NotImplementedError:
            pass


#____________________________________________________________________
# pth_smH

# Preprocessing
@flag_as_option
def RenameHggProcesses_smHcard(args):
    MergeHGGWDatacards.RenameProcesses_Aug21(
        LatestPaths.card_hgg_smH_PTH_unprocessed,
        )
@flag_as_option
def RenumberHzzProcesses_smHcard(args):
    MergeHGGWDatacards.RenumberProcessesHZZ_Aug21(
        LatestPaths.card_hzz_smH_PTH_unprocessed,
        )
@flag_as_option
def CombineCards_smHcard(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_PTH,
        'hzz=' + LatestPaths.card_hzz_smH_PTH
        )

@flag_as_option
def pth_smH_t2ws(args):
    card = LatestPathsGetters.get_card('pth_smH', args)
    if args.hzz:
        Commands.BasicT2WS(
            card,
            manualMaps = [
                '--PO \'map=.*/smH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                ],
            )
    elif args.hbb:
        raise NotImplementedError('Case --hbb in t2ws_smH is not yet implemented')
    elif args.combWithHbb:
        raise NotImplementedError('Case --combWithHbb in t2ws_smH is not yet implemented')
    else:
        Commands.BasicT2WS(
            LatestPaths.card_combined_smH_PTH,
            smartMaps = [
                ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            )

@flag_as_option
def pth_smH_scan(args):
    ws = LatestPathsGetters.get_ws('pth_smH', args) if not args.statonly else LatestPathsGetters.get_postfit(pth_smH, args)
    scan = GenericDifferentialObservableScan(args, 'pth_smH', ws)
    scan.run()

@flag_as_option
def pth_smH_postfit(args):
    make_postfit(args, LatestPathsGetters.get_ws('pth_smH', args), 'pth_smH')


#____________________________________________________________________
# pth_ggH

# Preprocessing
@flag_as_option
def RenameHggProcesses_Aug21(args):
    MergeHGGWDatacards.RenameProcesses_Aug21(
        LatestPaths.card_hgg_ggHxH_PTH_unprocessed,
        )
@flag_as_option
def RenumberHzzProcesses_Aug21(args):
    MergeHGGWDatacards.RenumberProcessesHZZ_Aug21(
        LatestPaths.card_hzz_ggHxH_PTH_unprocessed,
        )
@flag_as_option
def CombineCards_Aug21(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_ggHxH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_ggHxH_PTH,
        'hzz=' + LatestPaths.card_hzz_ggHxH_PTH
        )
@flag_as_option
def CombineCards_Dec15_hbb(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_hgg_hzz_hbb_ggHxH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_ggHxH_PTH,
        'hzz=' + LatestPaths.card_hzz_ggHxH_PTH,
        'hbb=' + LatestPaths.card_hbb_ggHxH_PTH
        )

@flag_as_option
def pth_ggH_t2ws(args):
    datacard = pth_ggH_get_card(args)
    outputWS = basename(datacard).replace( '.txt', '_xHfixed.root' )
    if args.hzz:
        Commands.BasicT2WS(
            datacard,
            manualMaps=[
                '--PO \'map=.*/ggH_PTH_0_15:r_ggH_PTH_0_15[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_15_30:r_ggH_PTH_15_30[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_30_45:r_ggH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_45_85:r_ggH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_85_125:r_ggH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_125_200:r_ggH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_GT350:r_ggH_PTH_GT200[1.0,0.0,3.0]\'',
                ],
            outputWS=outputWS
            )
    elif args.hbb:
        Commands.BasicT2WS(
            datacard,
            manualMaps=[
                # '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_200_350[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_350_600:r_ggH_PTH_350_600[1.0,0.0,10.0]\'',
                '--PO \'map=.*/ggH_PTH_GT600:r_ggH_PTH_GT600[1.0,0.0,10.0]\'',
                ],
            outputWS=outputWS
            )
    elif args.combWithHbb:
        Commands.BasicT2WS(
            datacard,
            manualMaps=[
                '--PO \'map=.*/ggH_PTH_GT350:r_ggH_PTH_350_600[1.0,0.0,10.0]\'',
                '--PO \'map=.*/ggH_PTH_350_600:r_ggH_PTH_350_600[1.0,0.0,10.0]\'',
                '--PO \'map=.*/ggH_PTH_GT600:r_ggH_PTH_GT600[1.0,0.0,10.0]\'',
                ],
            smartMaps = [
                ( r'.*/ggH_PTH_([\d\_GT]+)', r'r_ggH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            outputWS=outputWS
            )
    else:
        Commands.BasicT2WS(
            datacard,
            smartMaps = [
                ( r'.*/ggH_PTH_([\d\_GT]+)', r'r_ggH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            outputWS=outputWS
            )

@flag_as_option
def pth_ggH_scan(args):
    ws = LatestPathsGetters.get_ws('pth_ggH', args) if not args.statonly else LatestPathsGetters.get_postfit(pth_ggH, args)
    scan = GenericDifferentialObservableScan(args, 'pth_ggH', ws)

    if args.hzz:
        scan.POIRange = [ 0.0, 4.0 ]
    elif args.hbb:
        scan.POIRange = [ -10.0, 10.0 ]
        scan.extraOptions.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])
    elif args.combWithHbb:
        scan.POIRange = [ 0.0, 8.5 ]
        scan.nPoints = 150
        scan.nPointsPerJob = 2
        scan.extraOptions.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])
        if not args.statonly:
            scan.physicsModelParameterRanges = [
                [ 'qcdeff', 0.001, 8.0 ],
                [ 'r1p0', 0.0, 8.0 ],
                [ 'r2p0', 0.0, 8.0 ],
                [ 'r3p0', 0.0, 8.0 ],
                [ 'r0p1', 0.0, 8.0 ],
                [ 'r1p1', 0.0, 8.0 ],
                [ 'r2p1', 0.0, 8.0 ],
                [ 'r3p1', 0.0, 8.0 ],
                ]
    scan.run()

@flag_as_option
def pth_ggH_postfit(args):
    make_postfit(args, LatestPathsGetters.get_ws('pth_ggH', args), 'pth_ggH')

#____________________________________________________________________
# njets

@flag_as_option
def rename_hgg_njets(args):
    MergeHGGWDatacards.RenameProcesses_Hgg_nJets(
        LatestPaths.card_hgg_smH_NJ_unprocessed,
        globalReplace = [
            ( 'CMS_hgg_JEC', 'CMS_scale_j' )
            ]
        )

@flag_as_option
def njets_combineCards(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_NJ_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_NJ,
        'hzz=' + LatestPaths.card_hzz_smH_NJ
        )

@flag_as_option
def njets_t2ws(args):
    datacard = LatestPathsGetters.get_card('njets', args)
    if args.hzz:
        Commands.BasicT2WS(
            datacard,
            manualMaps=[
                '--PO \'map=.*/smH_NJ_0:r_smH_NJ_0[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_NJ_1:r_smH_NJ_1[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_NJ_2:r_smH_NJ_2[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_NJ_3:r_smH_NJ_GE3[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_NJ_GE4:r_smH_NJ_GE3[1.0,-1.0,4.0]\'',
                ],
            )
    else:
        Commands.BasicT2WS(
            datacard,
            smartMaps = [
                ( r'.*/smH_NJ_([\d\_GE]+)', r'r_smH_NJ_\1[1.0,-1.0,4.0]' )
                ],
            )

@flag_as_option
def njets_scan(args):
    ws = LatestPathsGetters.get_ws('njets', args) if not args.statonly else LatestPathsGetters.get_postfit(njets, args)
    scan = GenericDifferentialObservableScan(args, 'njets', ws)
    scan.nPoints = 55
    scan.nPointsPerJob = 5
    if args.hzz:
        scan.POIRange = [ 0.0, 4.0 ]
        scan.nPointsPerJob = scan.nPoints
    scan.run()

@flag_as_option
def njets_postfit(args):
    make_postfit(args, LatestPathsGetters.get_ws('njets', args), 'njets')

#____________________________________________________________________
# ptjet

@flag_as_option
def ptjet_rename_hgg(args):
    MergeHGGWDatacards.RenameProcesses_Hgg_differentials(
        LatestPaths.card_hgg_smH_PTJ_unprocessed
        )

@flag_as_option
def ptjet_combineCards(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_PTJ_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_PTJ,
        'hzz=' + LatestPaths.card_hzz_smH_PTJ
        )

@flag_as_option
def ptjet_t2ws(args):
    datacard = LatestPathsGetters.get_card('ptjet', args)
    if args.hzz:
        Commands.BasicT2WS(
            datacard,
            manualMaps=[
                '--PO \'map=.*/smH_PTJ_LT30:r_smH_PTJ_LT30[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTJ_30_55:r_smH_PTJ_30_55[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTJ_55_95:r_smH_PTJ_55_95[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTJ_95_120:r_smH_PTJ_GT95[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTJ_120_200:r_smH_PTJ_GT95[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTJ_GT200:r_smH_PTJ_GT95[1.0,-1.0,4.0]\'',
                ],
            suffix = '_ptjet'
            )
    else:
        Commands.BasicT2WS(
            datacard,
            smartMaps = [
                ( r'.*/smH_PTJ_([pm\d\_GELT]+)', r'r_smH_PTJ_\1[1.0,-1.0,4.0]' )
                ],
            # suffix = '_ptjet'
            )

@flag_as_option
def ptjet_scan(args):
    ws = LatestPathsGetters.get_ws('ptjet', args) if not args.statonly else LatestPathsGetters.get_postfit(ptjet, args)
    scan = GenericDifferentialObservableScan(args, 'ptjet', ws)
    scan.nPoints = 55
    scan.nPointsPerJob = 5
    if args.hzz:
        scan.POIRange = [ 0.0, 4.0 ]
        scan.nPointsPerJob = scan.nPoints
    scan.run()

@flag_as_option
def ptjet_postfit(args):
    make_postfit(args, LatestPathsGetters.get_ws('ptjet', args), 'ptjet')

#____________________________________________________________________
# rapidity

@flag_as_option
def rename_hgg_rapidity(args):
    MergeHGGWDatacards.RenameProcesses_Hgg_differentials(
        LatestPaths.card_hgg_smH_YH_unprocessed,
        )

@flag_as_option
def rapidity_combineCards(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_YH_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_YH,
        'hzz=' + LatestPaths.card_hzz_smH_YH
        )

@flag_as_option
def rapidity_t2ws(args):
    Commands.BasicT2WS(
        LatestPathsGetters.get_card('rapidity', args),
        smartMaps = [
            ( r'.*/smH_YH_([pm\d\_GE]+)', r'r_smH_YH_\1[1.0,-1.0,4.0]' )
            ],
        )

@flag_as_option
def rapidity_scan(args):
    ws = LatestPathsGetters.get_ws('rapidity', args) if not args.statonly else LatestPathsGetters.get_postfit('rapidity', args)
    scan = GenericDifferentialObservableScan(args, 'rapidity', ws)
    scan.nPoints = 55
    scan.nPointsPerJob = 5
    if args.hzz:
        scan.POIRange = [ 0.0, 4.0 ]
        scan.nPointsPerJob = scan.nPoints
    scan.run()

@flag_as_option
def rapidity_postfit(args):
    make_postfit(args, LatestPathsGetters.get_ws('rapidity', args), 'rapidity')
