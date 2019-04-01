#!/usr/bin/env python
"""
Thomas Klijnsma
"""

import re, os


def str_to_float(s):
    return float(s.replace('p','.').replace('m','-'))


class RTransformer(object):

    map_file_base = '/mnt/t3nfs01/data01/shome/vtavolar/FinalFits/CMSSW_7_1_5/src/flashggFinalFit/Plots/FinalResults/'
    map_files = {
        'PtNNLOPS_newBins'                                : 'expectedPrecision_differential_PtNNLOPS_newBins_ub_np40_refined.txt',
        'Njets2p5NNLOPS'                                  : 'expectedPrecision_differential_Njets2p5NNLOPS_ub_np40_refined.txt',
        'AbsRapidityNNLOPS_newBins'                       : 'expectedPrecision_differential_AbsRapidityNNLOPS_newBins_ub_np30_refined.txt',
        'CosThetaStarNNLOPS'                              : 'expectedPrecision_differential_CosThetaStarNNLOPS_ub_np200_refined.txt',
        'Jet2p5Pt0NNLOPS_newBins'                         : 'expectedPrecision_differential_Jet2p5Pt0NNLOPS_newBins_ub_np50_refined.txt',
        'Jet2p5AbsRapidity0NNLOPS_newBins'                : 'expectedPrecision_differential_Jet2p5AbsRapidity0NNLOPS_newBins_ub_np40_refined.txt',
        'AbsDeltaPhiGgJet0Eta2p5NNLOPS_newBins'           : 'expectedPrecision_differential_AbsDeltaPhiGgJet0Eta2p5NNLOPS_newBins_ub_np200_refined.txt',
        'AbsDeltaRapidityGgJet0Eta2p5NNLOPS'              : 'expectedPrecision_differential_AbsDeltaRapidityGgJet0Eta2p5NNLOPS_ub_np30_refined.txt',
        'Jet4p7Pt1NNLOPS_tightJetPuId'                    : 'expectedPrecision_differential_Jet4p7Pt1NNLOPS_tightJetPuId_ub_np50_refined.txt',
        'Jet4p7AbsRapidity1NNLOPS_tightJetPuId'           : 'expectedPrecision_differential_Jet4p7AbsRapidity1NNLOPS_tightJetPuId_ub_np40_refined.txt',
        'AbsDeltaPhiJjEta4p7NNLOPS_tightJetPuId'          : 'expectedPrecision_differential_AbsDeltaPhiJjEta4p7NNLOPS_tightJetPuId_ub_np30_refined.txt',
        'AbsDeltaPhiGgJjEta4p7NNLOPS_tightJetPuId'        : 'expectedPrecision_differential_AbsDeltaPhiGgJjEta4p7NNLOPS_tightJetPuId_ub_np50_refined.txt',
        'ZeppenfeldEta4p7NNLOPS_tightJetPuId'             : 'expectedPrecision_differential_ZeppenfeldEta4p7NNLOPS_tightJetPuId_ub_np40_refined.txt',
        'MjjEta4p7NNLOPS_newBins_tightJetPuId'            : 'expectedPrecision_differential_MjjEta4p7NNLOPS_newBins_tightJetPuId_ub_np50_refined.txt',
        'AbsDeltaEtaJJEta4p7NNLOPS_tightJetPuId'          : 'expectedPrecision_differential_AbsDeltaEtaJJEta4p7NNLOPS_tightJetPuId_ub_np50_refined.txt',
        'Jet4p7Pt1VBFlikeNNLOPS_tightJetPuId'             : 'expectedPrecision_differential_Jet4p7Pt1VBFlikeNNLOPS_tightJetPuId_ub_np80_refined.txt',
        'AbsDeltaPhiJjEta4p7VBFlikeNNLOPS_tightJetPuId'   : 'expectedPrecision_differential_AbsDeltaPhiJjEta4p7VBFlikeNNLOPS_tightJetPuId_ub_np50_refined.txt',
        'AbsDeltaPhiGgJjEta4p7VBFlikeNNLOPS_tightJetPuId' : 'expectedPrecision_differential_AbsDeltaPhiGgJjEta4p7VBFlikeNNLOPS_tightJetPuId_ub_np100_refined.txt',
        'PtNjets2p5NNLOPS_newBins_v2'                     : 'expectedPrecision_differential_PtNjets2p5NNLOPS_newBins_v2_ub_np70_refined.txt',
        'NjetsBflavorTight2p5NNLOPS'                      : 'expectedPrecision_differential_NjetsBflavorTight2p5NNLOPS_ub_np80_refined.txt',
        'NleptonsNNLOPS'                                  : 'expectedPrecision_differential_NleptonsNNLOPS_ub_np70_refined.txt',
        'METNNLOPS_newBins_v2'                            : 'expectedPrecision_differential_METNNLOPS_newBins_v2_ub_np500_refined.txt',
        '1Lepton1BjetNNLOPS'                              : 'expectedPrecision_differential_1Lepton1BjetNNLOPS_ub_np50_refined.txt',
        '1LeptonHighMETNNLOPS'                            : 'expectedPrecision_differential_1LeptonHighMETNNLOPS_ub_np50_refined.txt',
        '1LeptonLowMETNNLOPS'                             : 'expectedPrecision_differential_1LeptonLowMETNNLOPS_ub_np50_refined.txt',
        'InclusiveNNLOPS'                                 : 'expectedPrecision_differential_InclusiveNNLOPS_ub_np20.txt',
        }
    # title_dict = {
    #     'Pt'                           : 'p_{T}^{#gamma#gamma} (GeV)',
    #     'AbsRapidity'                  : '|y_{#gamma#gamma}|',
    #     'CosThetaStar'                 : '|cos(#theta*)|',
    #     'Jet2p5Pt0'                    : 'p_{T}^{j_{1}} (GeV)',
    #     'Jet2p5AbsRapidity0'           : '|y^{j_{1}}|',
    #     'AbsDeltaPhiGgJet0Eta2p5'      : '|#Delta #phi^{#gamma#gamma,j_{1}}|',
    #     'AbsDeltaRapidityGgJet0Eta2p5' : '|#Delta y^{#gamma#gamma,j_{1}}|',
    #     'PtNjets2p5'                   : 'p_{T}^{#gamma#gamma} x n_{jet}',
    #     'PtNjets2p5_0'                 : 'p_{T}^{#gamma#gamma}, n_{jet}=0 (GeV)',
    #     'PtNjets2p5_1'                 : 'p_{T}^{#gamma#gamma}, n_{jet}=1 (GeV)',
    #     'PtNjets2p5_1plus'             : 'p_{T}^{#gamma#gamma}, n_{jet}>1 (GeV)',
    #     'AbsDeltaEtaJJEta4p7'          : '|#Delta #eta^{j_{1},j_{2}}|',
    #     'AbsDeltaPhiGgJjEta4p7'        : '|#Delta #phi^{#gamma#gamma,j_{1}j_{2}}|',
    #     'AbsDeltaPhiJjEta4p7'          : '|#Delta #phi^{j_{1},j_{2}}|',
    #     'Jet4p7Pt1'                    : 'p_{T}^{j_{2}} (GeV)',
    #     'Jet4p7AbsRapidity1'           : '|y^{j_{2}}|',
    #     'ZeppenfeldEta4p7'             : '|#bar{#eta}^{j_{1}j_{2}} - #eta^{#gamma#gamma}|',
    #     'MjjEta4p7'                    : 'm^{j_{1}j_{2}} (GeV)',
    #     'Njets2p5'                     : 'n_{jet}',
    #     'MET'                          : 'p_{T}^{miss} (GeV)',
    #     'NjetsBflavorTight2p5'         : 'n_{jet}^{b}',
    #     'Nleptons'                     : 'n_{lepton}',
    #     '1Lepton1Bjet'                 : 'n_{lepton} = 1, N_{b} = 1',
    #     '1LeptonHighMET'               : 'n_{lepton} = 1, high p_{T}^{miss}',
    #     '1LeptonLowMET'                : 'n_{lepton} = 1, low p_{T}^{miss}',
    #     'Jet4p7Pt1VBFlike'             : 'p_{T}^{j_{2}} (GeV)',
    #     'AbsDeltaPhiGgJjEta4p7VBFlike' : '|#Delta #phi^{#gamma#gamma,j_{1}j_{2}}|',
    #     'AbsDeltaPhiJjEta4p7VBFlike'   : '|#Delta #phi^{j_{1},j_{2}}|',
    #     }
    title_dict = {
        'Pt'                           : 'p_{T}^{#gamma#gamma} (GeV)',
        'Njets2p5'                     : 'N_{jet}',
        'AbsRapidity'                  : '|y^{#gamma#gamma}|',
        'CosThetaStar'                 : '|cos(#theta*)|',
        'Jet2p5Pt0'                    : 'p_{T}^{j_{1}} (GeV)',
        'Jet2p5AbsRapidity0'           : '|y^{j_{1}}|',
        'AbsDeltaPhiGgJet0Eta2p5'      : '|#Delta#phi^{#gamma#gamma,j_{1}}|',
        'AbsDeltaRapidityGgJet0Eta2p5' : '|#Delta y^{#gamma#gamma,j_{1}}|',
        'Jet4p7Pt1'                    : 'p_{T}^{j_{2}} (GeV)',
        'Jet4p7AbsRapidity1'           : '|y^{j_{2}}|',
        'AbsDeltaPhiJjEta4p7'          : '|#Delta#phi^{j_{1},j_{2}}|',
        'AbsDeltaPhiGgJjEta4p7'        : '|#Delta#phi^{#gamma#gamma,j_{1}j_{2}}|',
        'ZeppenfeldEta4p7'             : '|#bar{#eta}^{j_{1}j_{2}} - #eta^{#gamma#gamma}|',
        'MjjEta4p7'                    : 'm^{j_{1}j_{2}} (GeV)',
        'AbsDeltaEtaJJEta4p7'          : '|#Delta#eta^{j_{1},j_{2}}|',
        'Jet4p7Pt1VBFlike'             : 'p_{T}^{j_{2}} (GeV)',
        'AbsDeltaPhiJjEta4p7VBFlike'   : '|#Delta#phi^{j_{1},j_{2}}|',
        'AbsDeltaPhiGgJjEta4p7VBFlike' : '|#Delta#phi^{#gamma#gamma,j_{1}j_{2}}|',
        'PtNjets2p5_0'                 : 'p_{T}^{#gamma#gamma}, N_{jet}=0 (GeV)',
        'PtNjets2p5_1'                 : 'p_{T}^{#gamma#gamma}, N_{jet}=1 (GeV)',
        'PtNjets2p5_1plus'             : 'p_{T}^{#gamma#gamma}, N_{jet}>1 (GeV)',
        'MET'                          : 'p_{T}^{miss} (GeV)',
        'NjetsBflavorTight2p5'         : 'N_{jet}^{b}',
        'Nleptons'                     : 'N_{lepton}',
        '1Lepton1Bjet'                 : '#geq 1 lepton, #geq 1 b',
        '1LeptonHighMET'               : '1-lepton, high-p_{T}^{miss}',
        '1LeptonLowMET'                : '1-lepton, low-p_{T}^{miss}',
        }

    obsnames = map_files.keys()
    obsnames.sort()

    def __init__(self):
        super(RTransformer, self).__init__()
        self.obsname = None
        self.process_str_dict = {}
        self.n_decimals = 1
        self.title = None

    def find_title(self):
        for key in sorted(self.title_dict.keys(), key=len, reverse=True):
            if re.match(key, self.obsname):
                self.title = self.title_dict[key]
                print 'Found title key \'{0}\' for obsname {1}; title: {2}'.format(
                    key, self.obsname, self.title
                    )
                break
        else:
            self.title = self.obsname
            print 'No title key found for {0}; using {0} instead'.format(self.obsname)

        if 'VBFlike' in self.obsname:
            print '  Adding \'(VBF-enriched)\' to the title'
            self.title += ' #scale[0.75]{(VBF-enriched)}'

    def set_obsname_from_corrmat(self, corrmat_file):
        for obsname in self.obsnames:
            if re.search(r'differential_{0}_DataBestFit'.format(obsname), corrmat_file):
            # if obsname in corrmat_file:
                print 'Found \'{0}\' for {1}'.format(obsname, corrmat_file)
                break
        else:
            raise RuntimeError('No map_file for {0}'.format(corrmat_file))
        self.obsname = obsname
        if obsname.startswith('AbsDeltaPhi'):
            self.n_decimals = 2
        self.map_file = os.path.join(self.map_file_base, self.map_files[self.obsname])
        self.get_process_str_dict_from_mapfile()

    def get_process_str_dict_from_mapfile(self):
        with open(self.map_file, 'r') as map_fp:
            text = map_fp.read()
        matches = re.findall(r'\-\-PO \'map=\.\*\/(.*?)\:(r\d+)', text)

        if len(matches) == 0:
            raise RuntimeError(
                'Found no matches for map_file {0}. Contents:\n{1}'
                .format(self.map_file, text)
                )

        for process_str, r in matches:
            self.process_str_dict[r] = process_str

    def get_process_str_from_r(self, r):
        if not r in self.process_str_dict:
            raise RuntimeError(
                '\'{0}\' was not found in the process_str_dict:\n{1}'
                .format(r, self.process_str_dict)
                )
        return self.process_str_dict[r]

    def get_boundaries_from_process_str(self, process_str):
        components = process_str.split('_')
        left = str_to_float(components[2])
        right = str_to_float(components[3])
        return left, right


    def is_overflow_heuristic(self, right):
        return right == 100. or right == 1000. or right == 13000. or right == 10000.

    def is_last_r(self, r):
        return r == self.get_sorted_rs()[-1]


    def format_regular(self, left, right):
        return '{0:.{decimals}f}-{1:.{decimals}f}'.format(
            left, right, decimals=self.n_decimals
            )        

    def format_overflow(self, left):
        return '>{0:.{decimals}f}'.format(left, decimals=self.n_decimals)

    def r_to_binlabel(self, r):
        if len(self.process_str_dict) == 0: self.get_process_str_dict_from_mapfile()
        process_str = self.process_str_dict[r]
        left, right = self.get_boundaries_from_process_str(process_str)

        # Process some exceptions
        if self.obsname == 'AbsRapidityNNLOPS_newBins' and right == 3.0:
            right = 2.5
        if self.obsname.startswith('1Lepton') and left == 0.0 and right == 2.0:
            return 'Selected'
        if self.obsname.startswith('Njets') or self.obsname.startswith('Nleptons'):
            if left == -1000.:
                return 'Underflow'
            elif self.is_overflow_heuristic(right) and self.is_last_r(r):
                return '#geq{0}'.format(int(left)+1)
            else:
                return '{0}'.format(int(right))

        if self.obsname == 'PtNjets2p5NNLOPS_newBins_v2':
            njet_label_dict = {
                0 : 'N_{jet}=0',
                1 : 'N_{jet}=1',
                2 : 'N_{jet}>1',
                }
            njets_left = int(left) / 13000
            njets_right = int(right) / 13000
            left %= 13000.
            right %= 13000.

            njets_str = ' ({0})'.format(njet_label_dict[njets_left])
            if njets_right > njets_left:
                return self.format_overflow(left) + njets_str
            else:
                return self.format_regular(left, right) + njets_str

        # Default procedure
        if left == -1000.:
            return 'Underflow'
        elif self.is_overflow_heuristic(right) and self.is_last_r(r):
            return self.format_overflow(left)
        else:
            return self.format_regular(left, right)

    def rs_to_binlabels(self, rs):
        return [ self.r_to_binlabel(r) for r in rs ]

    def sort_rs(self, rs):
        def get_left(r):
            return self.get_boundaries_from_process_str(self.get_process_str_from_r(r))[0]
        rs.sort(key=get_left)


    def get_sorted_rs(self):
        rs = self.process_str_dict.keys()
        self.sort_rs(rs)
        return rs

    def get_sorted_binlabels(self):
        return self.rs_to_binlabels(self.get_sorted_rs())

    def get_x_title(self):
        if self.title is None:
            self.find_title()
        return self.title


########################################
# Main
########################################

def main():
    """Run some tests on the classes and functions above"""
    r_transformer = RTransformer()

    r_transformer.set_obsname_from_corrmat('../out/corrmats_Mar12/higgsCombine_CORRMAT_differential_PtNNLOPS_newBins_DataBestFit.MultiDimFit.mH125.root')


    print r_transformer.r_to_binlabel('r0')



########################################
# End of Main
########################################
if __name__ == "__main__":
    main()