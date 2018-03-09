import differentials
import differentials.core as core
from differentials.core import AttrDict
import binheuristic
import logging, copy, os, glob, re
import scalecorrelation


class YukawaBinHeuristic(binheuristic.BinHeuristic):
    def __init__(self, *args, **kwargs):
        super(YukawaBinHeuristic, self).__init__(*args, **kwargs)
        self.corr_dict = {
            0.375 : 0.0,
            0.75  : 0.0,
            5.125 : 5.0,
            5.25  : 5.0,
            }

    def get_bin_boundaries(self, bin_centers):
        bin_boundaries = super(YukawaBinHeuristic, self).get_bin_boundaries(bin_centers)
        # For the quark induced histograms there is an irregularity at pt = 5.0
        for val, corr_val in self.corr_dict.iteritems():
            if val in bin_boundaries:
                logging.warning('Correcting bin boundary {0} to {1}'.format(val, corr_val))
                bin_boundaries[bin_boundaries.index(val)] = corr_val
        return bin_boundaries


class KappabKappacInterpreter(object):
    """docstring for KappabKappacInterpreter"""

    SM_qi = AttrDict(
        path = 'suppliedInput/fromPier/13tev-pth_quarkInduced_Aug04/higgs_plus_jet_13tev_1_1_mur050_muf050.pth',
        is_SM = True,
        kappab=1., kappac=1., muR=1., muF=1., Q=1.
        )
    SM_gi = AttrDict(
        path = 'suppliedInput/fromPier/histograms_ggH_May17/H125-LHC13-R04-MSbar-xmur050-xmuf050_1_1-xQ050-NNLO+NNLLmult.1_1.res',
        is_SM = True,
        kappab=1., kappac=1., muR=1., muF=1., Q=1.
        )

    def __init__(self):
        super(KappabKappacInterpreter, self).__init__()
        self.read_attr_dicts()
        self.set_sms()
        
    def get_scales(self, theory_file):
        theory_file = os.path.basename(theory_file)

        mur_match = re.search(r'mur(\d+)', theory_file)
        if mur_match:
            muR = float(mur_match.group(1))/50.
        else:
            muR = 1.0

        muf_match = re.search(r'muf(\d+)', theory_file)
        if muf_match:
            muF = float(muf_match.group(1))/50.
        else:
            muF = 1.0

        xQ_match = re.search(r'xQ(\d+)', theory_file)
        if xQ_match:
            Q = float(xQ_match.group(1))/50.
        else:
            Q = 1.0

        return muR, muF, Q

    def get_couplings(self, theory_file):
        theory_file = os.path.basename(theory_file)

        gI_match = re.search(r'\.([\d\-]+)_([\d\-]+)\.res', theory_file)
        qI_match = re.search(r'_([m\d]+)_([m\d]+)\.pth', theory_file)
        qI_scalevar_match = re.search(r'_([m\d]+)_([m\d]+)_mur\d+_muf\d+\.pth', theory_file)

        if gI_match:
            kappab = float(gI_match.group(1))
            kappac = float(gI_match.group(2))
            return kappab, kappac
        elif qI_match:
            kappab = core.str_to_float(qI_match.group(1))
            kappac = core.str_to_float(qI_match.group(2))
            return kappab, kappac
        elif qI_scalevar_match:
            kappab = core.str_to_float(qI_scalevar_match.group(1))
            kappac = core.str_to_float(qI_scalevar_match.group(2))
            return kappab, kappac
        else:
            raise ValueError('Could not extract couplings from file {0}'.format(theory_file))

    def make_attr_dict(self, theory_file):
        d = AttrDict(path = theory_file)
        muR, muF, Q = self.get_scales(theory_file)
        kappab, kappac = self.get_couplings(theory_file)

        d.muR = muR
        d.muF = muF
        d.Q = Q
        d.kappab = kappab
        d.kappac = kappac
        return d

    def is_scale_variation(self, d):
        if d.muR == 1. and d.muF == 1. and d.Q == 1.:
            return False
        return True

    def read_attr_dicts(self):
        self.coupling_variations_qi = []
        self.scale_variations_qi = []
        for theory_file in glob.glob('suppliedInput/fromPier/13tev-pth_quarkInduced_Aug04/*.pth'):
            if theory_file == self.SM_qi.path: continue
            d = self.make_attr_dict(theory_file)
            if self.is_scale_variation(d):
                self.scale_variations_qi.append(d)
            else:
                self.coupling_variations_qi.append(d)

        self.coupling_variations_gi = []
        self.scale_variations_gi = []
        for theory_file in glob.glob('suppliedInput/fromPier/histograms_ggH_May17/*.res'):
            if theory_file == self.SM_gi.path: continue
            d = self.make_attr_dict(theory_file)
            if self.is_scale_variation(d):
                self.scale_variations_gi.append(d)
            else:
                self.coupling_variations_gi.append(d)


    def set_sms(self):
        SM_gi_theory = YukawaTheoryGluonInduced(self.SM_gi)
        YukawaTheoryGluonInduced.set_sm(SM_gi_theory)
        SM_qi_theory = YukawaTheoryQuarkInduced(self.SM_qi)
        YukawaTheoryQuarkInduced.set_sm(SM_qi_theory)


    def dump_gluon_induced(self):
        for d in self.scale_variations_gi+self.coupling_variations_gi:
            theory = YukawaTheoryGluonInduced(d)
            theory.dump()
        YukawaTheoryGluonInduced.sm.dump()

    def dump_quark_induced(self):
        for d in self.scale_variations_qi+self.coupling_variations_qi:
            theory = YukawaTheoryQuarkInduced(d)
            theory.dump()
        YukawaTheoryQuarkInduced.sm.dump()

    def scale_quark_induced(self):
        theories = []
        for d in self.coupling_variations_qi:
            theories.append(YukawaTheoryQuarkInduced(d))

        parametrization = differentials.parametrization.Parametrization()
        for theory in theories:
            parametrization.add_variation(theory.d.kappab, theory.d.kappac, theory.xs_per_GeV)
        parametrization.parametrize()

        for theory in theories + [YukawaTheoryQuarkInduced.sm]:
            xss_theory = theory.xs_per_GeV
            xss_param  = parametrization.evaluate(theory.d.kappab, theory.d.kappac)

            logging.trace('kappab: {0} / kappac: {1}'.format(theory.d.kappab, theory.d.kappac))
            for xs_theory, xs_param in zip(xss_theory, xss_param):
                logging.trace('{0:+12.8f}  /  {1:+12.8f}'.format(xs_theory, xs_param))





class YukawaTheory(object):
    """docstring for YukawaTheory"""

    bin_heuristic = YukawaBinHeuristic()

    def __init__(self, d):
        super(YukawaTheory, self).__init__()
        self.d = d
        logging.debug('path = {0}'.format(self.d.path))
        logging.debug('d contents:\n{0}'.format(self.d))
        self.tags = ['yukawa']

    def get_outpath(self):
        outpath = 'out/theories_{0}_{1}'.format(
            core.datestr(),
            '_'.join(self.tags),
            )
        return outpath

    def get_outname(self):
        outname = (
            '{tags}'
            '_kappab_{kappab}'
            '_kappac_{kappac}'
            '_muR_{muR}'
            '_muF_{muF}'
            '_Q_{Q}'
            '.txt'
            .format(
                tags='_'.join(self.tags),
                kappab = core.float_to_str(self.d.kappab),
                kappac = core.float_to_str(self.d.kappac),
                muR = core.float_to_str(self.d.muR),
                muF = core.float_to_str(self.d.muF),
                Q = core.float_to_str(self.d.Q),
                )
            )
        return outname

    def dump(self):
        outpath = self.get_outpath()
        logging.debug('Making directory {0}'.format(outpath))
        if not core.is_testmode() and not os.path.isdir(outpath):
            os.makedirs(outpath)

        outname = os.path.join(outpath, self.get_outname())
        logging.info('Dumping contents to {0}'.format(outname))

        contents = [
            'file={0}'.format(self.d.path),
            'kappab={0}'.format(self.d.kappab),
            'kappac={0}'.format(self.d.kappac),
            'muR={0}'.format(self.d.muR),
            'muF={0}'.format(self.d.muF),
            'Q={0}'.format(self.d.Q),
            'binBoundaries={0}'.format(','.join(map(str, self.bin_boundaries))),
            'binCenters={0}'.format(','.join(map(str, self.bin_centers))),
            'crosssection={0}'.format(','.join(map( str, self.xs_per_GeV))),
            'crosssection_integrated={0}'.format(','.join(map( str, self.xs))),
            'ratios={0}'.format(','.join(map(str, self.ratio))),
            ]

        contents = '\n'.join(contents)
        logging.debug('Contents:\n{0}'.format(contents))

        if not core.is_testmode():
            with open(outname, 'w') as out_fp:
                out_fp.write(contents)


    def read(self):
        columns = core.read_data(self.d.path, columns=True, make_float=True)

        self.bin_centers = columns[0]
        self.n_bins = len(self.bin_centers)
        self.bin_boundaries = self.bin_heuristic.get_bin_boundaries(self.bin_centers)

        # Take matched_xs (resummed is columns[2])
        # Default unit is nb, transform to fb
        self.xs_per_GeV = [ 1000.*xs for xs in columns[1] ]

        if getattr(self.d, 'is_SM', False):
            self.ratio = [ 1.0 for i in xrange(self.n_bins) ]
        else:
            self.ratio = [ xs/smxs if smxs!=0.0 else 0.0 for xs, smxs in zip(self.xs_per_GeV, self.sm.xs_per_GeV) ]

        bin_widths = [ r-l for r, l in zip(self.bin_boundaries[1:], self.bin_boundaries[:-1]) ]
        self.xs = [ xs * width for xs, width in zip(self.xs_per_GeV, bin_widths) ]



class YukawaTheoryGluonInduced(YukawaTheory):
    """docstring for YukawaTheoryGluonInduced"""

    sm_is_defined = False
    sm = None

    def __init__(self, d):
        super(YukawaTheoryGluonInduced, self).__init__(d)
        self.tags.append('gluoninduced')
        self.read()

    @staticmethod
    def set_sm(sm):
        YukawaTheoryGluonInduced.sm_is_defined = True
        YukawaTheoryGluonInduced.sm = sm
        


class YukawaTheoryQuarkInduced(YukawaTheory):
    """docstring for YukawaTheoryQuarkInduced"""

    sm_is_defined = False
    sm = None

    mb_old = 4.65
    mb_new = 2.963
    mc_old = 1.275
    mc_new = 0.655

    def __init__(self, d):
        super(YukawaTheoryQuarkInduced, self).__init__(d)
        self.tags.append('quarkinduced')
        self.read()

    @staticmethod
    def set_sm(sm):
        YukawaTheoryQuarkInduced.sm_is_defined = True
        YukawaTheoryQuarkInduced.sm = sm
        
    # def scale(self):






