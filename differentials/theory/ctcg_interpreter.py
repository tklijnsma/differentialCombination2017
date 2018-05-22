import differentials.core as core
from differentials.core import AttrDict
import binheuristic
import logging, copy, os, os.path
import scalecorrelation


def create_all_ctcg():
    create_SM()
    create_scale_variations()
    create_coupling_variations_highpt()

def create_SM():
    if not TopTheory.sm_is_defined:
        sm = TopTheory(SM)
        sm.get_columns()
        TopTheory.set_sm(sm)
        sm.dump()

def create_scale_variations():
    for d in scale_variations:
        theory = TopTheory(d)
        theory.get_columns()
        theory.dump()

def create_coupling_variations_highpt():
    for d in coupling_variations_lowpt:
        theory_lowpt = TopTheory(d)
        theory_lowpt.get_columns()

        # Find corresponding highpt
        for d_highpt in coupling_variations_highpt:
            if theory_lowpt.equals(d_highpt):
                break
        else:
            raise RuntimeError('Could not find a match for low pt coupling variation {0}'.format(theory_lowpt.d))

        theory_highpt = TopTheory(d_highpt)
        theory_highpt.get_columns()

        theory_merged = theory_lowpt.merge(theory_highpt)
        theory_merged.dump()


class TopTheory(object):
    sm_is_defined = False
    sm = None

    @staticmethod
    def set_sm(sm):
        TopTheory.sm_is_defined = True
        TopTheory.sm = sm

    def __init__(self, d):
        logging.debug('Instantiating new TopTheory')
        self.d = self.fill_defaults(d)
        logging.debug('path = {0}'.format(self.d.path))
        logging.debug('d contents:\n{0}'.format(self.d))
        self.tags = ['tophighpt']

    def equals(self, d_to_compare):
        d = self.fill_defaults(d_to_compare)
        for key in self.d.keys():
            if key in ['path', 'is_SM', 'first_column_is_xs']: continue
            try:
                if not(self.d[key] == d[key]):
                    return False
            except KeyError:
                logging.error(
                    'Key {0} was not found in either:\n{1}\nor\n{2}'
                    .format(key, self.d, d)
                    )
                raise
        else:
            return True

    def merge(self, other):
        logging.info(
            'Merging {0} and {1}'
            .format(self.d.path, other.d.path)
            )
        res = copy.deepcopy(self)
        res.other_file = other.d.path
        res.bin_boundaries += other.bin_boundaries[1:]
        res.bin_centers += other.bin_centers
        res.n_bins += other.n_bins
        res.xs_per_GeV += other.xs_per_GeV
        res.xs += other.xs
        res.ratio += other.ratio
        return res

    def fill_defaults(self, d_orig):
        d = copy.copy(d_orig)
        defaults = {
            'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0,
            'cg' : 0.0, 'ct' : 1.0, 'cb' : 1.0,
            'first_column_is_xs' : False, 'is_SM': False
            }
        for key, value in defaults.iteritems():
            if not key in d.keys():
                setattr(d, key, value)
                d[key] = value
        return d

    def get_smxs_per_GeV(self):
        if not self.sm_is_defined:
            raise RuntimeError('Run TopTheory.set_sm(sm) with the SM variation first')

        if self.bin_boundaries[0] == self.sm.bin_boundaries[0]:
            return self.sm.xs_per_GeV
        else:
            # For the high pt variations, the values before the spectrum begins should be dropped
            start_index = self.sm.bin_boundaries.index(self.bin_boundaries[0])
            return self.sm.xs_per_GeV[start_index:]

    def get_columns(self):
        columns = core.read_data(self.d.path, columns=True, make_float=True)
        self.bin_centers = columns[0]
        self.n_bins = len(self.bin_centers)
        self.bin_boundaries = bin_heuristic.get_bin_boundaries(self.bin_centers)

        if self.d.first_column_is_xs:
            logging.warning('Dividing cross section in {0} by factor 2.27'.format(self.d.path))
            self.xs_per_GeV = [ v/2.27 for v in columns[1] ]
            if self.d.is_SM:
                self.ratio = [ 1.0 for xs in self.xs_per_GeV ]
            elif self.sm_is_defined:
                self.ratio = [ xs / sm_xs for xs, sm_xs in zip(self.xs_per_GeV, self.get_smxs_per_GeV()) ]
            else:
                raise RuntimeError('Run TopTheory.set_sm(sm) with the SM variation first')
        else:
            self.ratio = columns[1]
            if self.sm_is_defined:
                self.xs_per_GeV = [ ratio * sm_xs for ratio, sm_xs in zip(self.ratio, self.get_smxs_per_GeV()) ]
            else:
                raise RuntimeError('Run TopTheory.set_sm(sm) with the SM variation first')
        bin_widths = [ r-l for r, l in zip(self.bin_boundaries[1:], self.bin_boundaries[:-1]) ]
        self.xs = [ xs * width for xs, width in zip(self.xs_per_GeV, bin_widths) ]

    def get_outpath(self):
        outpath = 'out/theories_{0}_{1}'.format(
            core.datestr(),
            '_'.join(self.tags),
            )
        return outpath

    def get_outname(self):
        outname = '{0}_ct_{1}_cg_{2}_cb_{3}_muR_{4}_muF_{5}_Q_{6}.txt'.format(
            '_'.join(self.tags),
            core.float_to_str(self.d.ct), core.float_to_str(self.d.cg), core.float_to_str(self.d.cb),
            core.float_to_str(self.d.muR), core.float_to_str(self.d.muF), core.float_to_str(self.d.Q)
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
            'ct={0}'.format(self.d.ct),
            'cg={0}'.format(self.d.cg),
            'cb={0}'.format(self.d.cb),
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


class TopBinHeuristic(binheuristic.BinHeuristic):
    def __init__(self, *args, **kwargs):
        super(TopBinHeuristic, self).__init__(*args, **kwargs)

        self.faulty_bin_centers = [
            401., 410., 420., 430., 440., 450., 460., 470., 480., 490.,
            500., 510., 520., 530., 540., 550., 560., 570., 580., 590.,
            600., 610., 620., 630., 640., 650., 660., 670., 680., 690.,
            700., 710., 720., 730., 740., 750., 760., 770., 780., 790.
            ]
        self.corrected_bin_centers = [
            405., 415., 425., 435., 445., 455., 465., 475., 485., 495.,
            505., 515., 525., 535., 545., 555., 565., 575., 585., 595.,
            605., 615., 625., 635., 645., 655., 665., 675., 685., 695.,
            705., 715., 725., 735., 745., 755., 765., 775., 785., 795.
            ]

    def get_bin_boundaries(self, bin_centers):
        if bin_centers == self.faulty_bin_centers:
            logging.warning('Copying manually corrected bin centers')
            bin_centers = copy.copy(self.corrected_bin_centers)
        bin_boundaries = super(TopBinHeuristic, self).get_bin_boundaries(bin_centers)
        if 51.75 in bin_boundaries:
            logging.warning('Correcting bin boundary 51.75 to 51.0')
            bin_boundaries[bin_boundaries.index(51.75)] = 51
        return bin_boundaries

bin_heuristic = TopBinHeuristic()



SM = AttrDict(
    path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/SM_NNLO',
    first_column_is_xs = True,
    is_SM = True
    )

scale_variations = [
    # This is just the SM
    # AttrDict(
    #     path = 'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/HRes_mR1mF1.top',
    #     first_column_is_xs = True,
    #     ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/HRes_mR1mF1Q2.top',
        Q = 2.0,
        first_column_is_xs = True,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/HRes_mR1mF1Qh.top',
        Q = 0.5,
        first_column_is_xs = True,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/HRes_mR1mF2.top',
        muF = 2.0, 
        first_column_is_xs = True,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/HRes_mR1mFh.top',
        muF = 0.5, 
        first_column_is_xs = True,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/HRes_mR2mF1.top',
        muR = 2.0, 
        first_column_is_xs = True,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/HRes_mR2mF2.top',
        muR = 2.0, muF = 2.0, 
        first_column_is_xs = True,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/HRes_mRhmF1.top',
        muR = 0.5, 
        first_column_is_xs = True,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/HRes_mRhmFh.top',
        muR = 0.5, muF = 0.5, 
        first_column_is_xs = True,
        ),
    ]

coupling_variations_highpt = [
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_ctcg_01_new',
        ct = 0.1, cg = 0.075,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_ctcg_05_new',
        ct = 0.5, cg = 0.042,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_ctcg_15_new',
        ct = 1.5, cg = -0.042,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_ctcg_2_new',
        ct = 2.0, cg = -0.083,
        ),

    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cbct_05_new',
        ct = 0.5, cb = -7.46,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cbct_08_new',
        ct = 0.8, cb = -3.67,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cbct_09_new',
        ct = 0.9, cb = -1.79,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cbct_11_new',
        ct = 1.1, cb = 3.79,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cbct_12_new',
        ct = 1.2, cb = 4.67,
        ),

    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cg003ct12_new',
        ct = 1.2, cb = -2.98 ,  cg = -0.03,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cg003ct13_new',
        ct = 1.3, cb = -0.85 ,  cg = -0.03,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cg003ct14_new',
        ct = 1.4, cb = 3.31 ,   cg = -0.03,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cg004ct12_new',
        ct = 1.2, cb = -4.89 ,  cg = -0.04,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cg004ct13_new',
        ct = 1.3, cb = -3.34 ,  cg = -0.04,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cg004ct15_new',
        ct = 1.5, cb = 1.88 ,   cg = -0.04,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cg005ct14_new',
        ct = 1.4, cb = -3.67 ,  cg = -0.05,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cg005ct15_new',
        ct = 1.5, cb = -1.79 ,  cg = -0.05,
        ),

    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_ctdown_new',
        ct = 0.9,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_ctup_new',
        ct = 1.1,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cgdown_new',
        cg = -0.008,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cgup_new',
        cg = 0.008,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cbdown_new',
        cb = -2.0,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cbup_new',
        cb = 4.0,
        ),
    ]

coupling_variations_lowpt = [
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/ratio_ctcg01_new',
        ct = 0.1, cg = 0.075,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/ratio_ctcg05_new',
        ct = 0.5, cg = 0.042,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/ratio_ctcg15_new',
        ct = 1.5, cg = -0.042,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/ratio_ctcg2_new',
        ct = 2.0, cg = -0.083,
        ),

    AttrDict(
        path = 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/ratio_ctcb05_new',
        ct = 0.5, cb = -7.46,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/ratio_ctcb08_new',
        ct = 0.8, cb = -3.67,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/ratio_ctcb09_new',
        ct = 0.9, cb = -1.79,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/ratio_ctcb11_new',
        ct = 1.1, cb = 3.79,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/ratio_ctcb12_new',
        ct = 1.2, cb = 4.67,
        ),

    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_ctup_new',
        ct = 1.1,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_ctdown_new',
        ct = 0.9,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cgup_new',
        cg = 0.008,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cgdown_new',
        cg = -0.008,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cbup_new',
        cb = 4.0,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cbdown_new',
        cb = -2.0,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cg003ct12_new',
        ct = 1.2, cb = -2.98 ,  cg = -0.03,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cg003ct13sw_new',
        ct = 1.3, cb = -0.85 ,  cg = -0.03,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cg003ct14sw_new',
        ct = 1.4, cb = 3.31 ,   cg = -0.03,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cg004ct12sw_new',
        ct = 1.2, cb = -4.89 ,  cg = -0.04,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cg004ct13sw_new',
        ct = 1.3, cb = -3.34 ,  cg = -0.04,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cg004ct15sw_new',
        ct = 1.5, cb = 1.88 ,   cg = -0.04,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cg005ct14sw_new',
        ct = 1.4, cb = -3.67 ,  cg = -0.05,
        ),
    AttrDict(
        path = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/ratio_cg005ct15sw_new',
        ct = 1.5, cb = -1.79 ,  cg = -0.05,
        ),
    ]



# #____________________________________________________________________
# def DumpContainerToFile_Top( container, prefix, outdir ):

#     # Check if ratios, crosssection and binBoundaries have expected lenghts
#     if not(
#         len(container.crosssection) == len(container.binBoundaries)-1
#         and len(container.ratios) == len(container.binBoundaries)-1
#         ):
#         Commands.ThrowError((
#             'Lists have unexpected lenghts. Found:'
#             '\n  file = {3}'
#             '\n  len(binBoundaries) = {0}'
#             '\n  len(crosssection)  = {1}'
#             '\n  len(ratios)        = {2}'
#             ).format( len(container.binBoundaries), len(container.crosssection), len(container.ratios), container.file )
#             )

#     outname = prefix
#     for coupling in [ 'ct', 'cg', 'cb' ]:
#         if hasattr( container, coupling ):
#             outname += '_{0}_{1}'.format(
#                 coupling, Commands.ConvertFloatToStr( getattr(container, coupling) )
#                 )
    
#     if hasattr( container, 'muR' ) and hasattr( container, 'muF' ):
#         outname += '_muR_{0}_muF_{1}'.format(
#             Commands.ConvertFloatToStr( container.muR ),
#             Commands.ConvertFloatToStr( container.muF )
#             )

#     if hasattr( container, 'Q' ):
#         outname += '_Q_{0}'.format( Commands.ConvertFloatToStr(container.Q) )

#     outname += '.txt'

#     container.derivedTheoryFilePath = abspath( join( outdir, outname ) )

#     if not isdir( outdir ): os.makedirs( outdir )

#     with open( container.derivedTheoryFilePath, 'w' ) as outFp:
#         w = lambda text: outFp.write( text + ' \n' )

#         w( 'file={0}'.format(container.file) )
#         if hasattr( container, 'secondfile' ): w( 'secondfile={0}'.format(container.secondfile) )

#         if hasattr( container, 'ct' ): w( 'ct={0}'.format(container.ct) )
#         if hasattr( container, 'cg' ): w( 'cg={0}'.format(container.cg) )
#         if hasattr( container, 'cb' ): w( 'cb={0}'.format(container.cb) )

#         if hasattr( container, 'muR' ) and hasattr( container, 'muF' ):
#             w( 'muR={0}'.format(container.muR) )
#             w( 'muF={0}'.format(container.muF) )

#         if hasattr( container, 'Q' ):  w( 'Q={0}'.format(container.Q) )

#         w( 'n_binBoundaries={0}'.format(len(container.binBoundaries)) )
#         w( 'n_crosssection={0}'.format(len(container.crosssection)) )
#         w( 'n_ratios={0}'.format(len(container.ratios)) )

#         w( 'binBoundaries={0}'.format( ','.join(map( str, container.binBoundaries )) ) )
#         w( 'crosssection={0}'.format( ','.join(map( str, container.crosssection )) ) )
#         if hasattr( container, 'binCenters' ): w( 'binCenters={0}'.format( ','.join(map( str, container.binCenters )) ) )
#         w( 'ratios={0}'.format( ','.join(map( str, container.ratios )) ) )


# #____________________________________________________________________
# faultyBinCenters = [ 401., 410., 420., 430., 440., 450., 460., 470., 480., 490., 500., 510., 520., 530., 540., 550., 560., 570., 580., 590., 600., 610., 620., 630., 640., 650., 660., 670., 680., 690., 700., 710., 720., 730., 740., 750., 760., 770., 780., 790. ]
# correctedBinCenters = [ 405., 415., 425., 435., 445., 455., 465., 475., 485., 495., 505., 515., 525., 535., 545., 555., 565., 575., 585., 595., 605., 615., 625., 635., 645., 655., 665., 675., 685., 695., 705., 715., 725., 735., 745., 755., 765., 775., 785., 795. ]
# def ReadLinesOfTheoryFile_Top( theoryFile, verbose=False, SM=None, isSMFile=False ):
#     # Read lines
#     with open( theoryFile, 'r' ) as theoryFp:
#         lines = [ line.strip() for line in theoryFp.readlines() ]
#     commentlines = [ line.strip() for line in lines if line.startswith('#') ]
#     lines        = [ line for line in lines if not line.startswith('#') and len(line) > 0 ]

#     if isSMFile:

#         print '\nFound SM file', theoryFile

#         pts  = []
#         xss  = []

#         for line in lines:
#             components = line.split()
#             pt         = float(components[0])
#             xs         = float(components[1]) / 2.27 # 1 over Hgg BR (peculiarity from Agnieszka)

#             if verbose:
#                 if SM is None:
#                     print '    pt = {0:<8.3f} |  xs = {1:<10.6f}'.format( pt, xs )
#                 else:
#                     SMxs = SM.crosssection[ SM.pts.index(pt) ]
#                     print '    pt = {0:<8.3f} |  xs = {1:<10.6f} |  SMxs = {2:<10.6f} |  ratio = {3:<10.6f}'.format(
#                         pt, xs, SMxs, xs/SMxs if not SMxs == 0. else 0. )

#             pts.append( pt )
#             xss.append( xs )

#         return pts, xss

#     else:

#         if SM is None:
#             Commands.ThrowError( 'SM should be specified when reading a non-SM file' )
#         elif theoryFile == SM.file:
#             print 'This file was found to be the SM file; simply return SM values'
#             return SM.pts, SM.crosssection, [ 1.0 for xs in SM.crosssection ]
#         else:
#             print '\nReading', theoryFile


#         pts    = []
#         ratios = []
#         xss    = []

#         # Check for the faulty bin centers
#         correctingForFaultyBinCenters = False
#         binCenters = [ float(line.split()[0]) for line in lines ]
#         if binCenters == faultyBinCenters:
#             Commands.Warning( 'Correcting faulty bin centers to corrected bin centers' )
#             getcorrectedpt = lambda pt: correctedBinCenters[faultyBinCenters.index(pt)]
#             correctingForFaultyBinCenters = True

#         for line in lines:

#             components = line.split()

#             if AgnieszkasFilenameDecoder[ basename(theoryFile).replace('/','') ]['firstColumnIsRatio']:
#                 print '    Assuming the first column is a RATIO'
#                 pt         = float(components[0])
#                 ratio      = float(components[1])
#                 if correctingForFaultyBinCenters: pt = getcorrectedpt(pt)
#                 SMxs = SM.crosssection[ SM.pts.index(pt) ]
#                 xs   = SMxs * ratio
#             else:
#                 print '    Assuming the first column is a CROSS SECTION - ALSO DIVIDING BY 2.27!!'
#                 pt         = float(components[0])
#                 xs         = float(components[1]) / 2.27
#                 if correctingForFaultyBinCenters: pt = getcorrectedpt(pt)
#                 SMxs  = SM.crosssection[ SM.pts.index(pt) ]
#                 ratio = xs / SMxs

#             if verbose:
#                 print '    pt = {0:<8.3f} |  xs = {1:<10.6f} |  SMxs = {2:<10.6f} |  ratio = {3:<10.6f}'.format(
#                     pt, xs, SMxs, ratio )

#             pts.append( pt )
#             ratios.append( ratio )
#             xss.append( xs )

#         return pts, xss, ratios
   


# #____________________________________________________________________
# def CreateDerivedTheoryFiles_Top(
#         theoryDirs = [
#             'suppliedInput/fromAgnieszka/HRes_SMEFT_May16',
#             'suppliedInput/fromAgnieszka/SMEFTscaling_May16',
#             'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17',
#             ],
#         verbose = False,
#         ):

#     outdir = 'derivedTheoryFiles_{0}_Top'.format( datestr )
#     if not isdir( outdir ): os.makedirs( outdir )

#     theoryDirs = [ abspath( theoryDir ) for theoryDir in theoryDirs ]
#     theoryFiles = []
#     for theoryDir in theoryDirs:
#         theoryFiles.extend( glob( join( theoryDir, '*' ) ) )


#     # ======================================
#     # First look for SM file

#     for theoryFile in theoryFiles:
#         if basename(theoryFile).replace('/','') == 'SM_NNLO':
#             smFile = theoryFile
#             break
#     else:
#         Commands.ThrowError( 'Could not find a SM file' )
#         sys.exit()

#     SM = Container()
#     SM.file = smFile

#     binCenters, crosssection = ReadLinesOfTheoryFile_Top( SM.file, isSMFile=True, verbose=verbose )
#     newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic(
#         binCenters, manualSwitchAt50=True, manualSwitchAt5=False )

#     SM.pts           = binCenters
#     SM.binCenters    = deepcopy( newBinCenters )
#     SM.binBoundaries = deepcopy( binBoundaries )
#     SM.crosssection  = deepcopy( crosssection )
#     SM.ratios        = [ 1.0 for xs in SM.crosssection ]

#     for key, value in AgnieszkasFilenameDecoder[basename(SM.file).replace('/','')].iteritems():
#         setattr( SM, key, value )

#     # DumpContainerToFile_Top( SM, prefix='PureSM', outdir=outdir )


#     # ======================================
#     # Process the other theory files

#     for theoryFile in theoryFiles:

#         if verbose:
#             print '\n\n' + '-'*80
#             print 'Processing {0}'.format( theoryFile )
#             print '\n'

#         container = Container()
#         container.file = theoryFile

#         shortFileName = basename(container.file).replace('/','')
#         if not shortFileName in AgnieszkasFilenameDecoder:
#             print 'Skipping \'{0}\''.format(shortFileName)
#             continue

#         binCenters, crosssection, ratios = ReadLinesOfTheoryFile_Top( container.file, verbose, SM=SM )
#         newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic(
#             binCenters, manualSwitchAt50=True, manualSwitchAt5=False )

#         container.binCenters    = newBinCenters
#         container.binBoundaries = binBoundaries
#         container.crosssection  = crosssection
#         container.ratios        = ratios

#         if verbose:
#             print '\nWriting the following to a file:'
#             for i in xrange(len(container.crosssection)):
#                 print '    {0:6.1f} - {1:6.1f}  |  xs = {2:10.6f}  |  ratio = {3:10.6f}'.format(
#                     container.binBoundaries[i], container.binBoundaries[i+1],
#                     container.crosssection[i], container.ratios[i]
#                     )

#         for key, value in AgnieszkasFilenameDecoder[shortFileName].iteritems():
#             setattr( container, key, value )

#         DumpContainerToFile_Top( container, prefix='Top', outdir=outdir )





# #____________________________________________________________________
# def CreateDerivedTheoryFiles_Top_highPt(
#         verbose = True,
#         ):

#     SMfile = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/SM_NNLO'

#     lowRangeFiles = glob( 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/*' )
#     highRangeFiles = (
#         glob( 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cbct_*_new' )
#         + glob( 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_ctcg_*_new' )
#         )


#     # ======================================
#     # Make SM variation

#     SM = Container()
#     SM.file = SMfile

#     binCenters, crosssection = ReadLinesOfTheoryFile_Top( SM.file, isSMFile=True, verbose=verbose )
#     newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic(
#         binCenters, manualSwitchAt50=True, manualSwitchAt5=False )

#     SM.pts           = binCenters
#     SM.binCenters    = deepcopy( newBinCenters )
#     SM.binBoundaries = deepcopy( binBoundaries )
#     SM.crosssection  = deepcopy( crosssection )
#     SM.ratios        = [ 1.0 for xs in SM.crosssection ]

#     for key, value in AgnieszkasFilenameDecoder[basename(SM.file).replace('/','')].iteritems():
#         setattr( SM, key, value )


#     # ======================================
#     # Finding pairs of low and high range files

#     def f(theoryFile):
#         if verbose:
#             print '\n\n' + '-'*80
#             print 'Processing {0}'.format( theoryFile )
#             print '\n'

#         shortFileName = basename(theoryFile).replace('/','')
#         if not shortFileName in AgnieszkasFilenameDecoder:
#             print 'Skipping \'{0}\''.format(shortFileName)
#             return

#         container = Container()
#         container.file = theoryFile

#         binCenters, crosssection, ratios = ReadLinesOfTheoryFile_Top( container.file, verbose, SM=SM )

#         newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic(
#             binCenters, manualSwitchAt50=True, manualSwitchAt5=False )

#         container.binCenters    = newBinCenters
#         container.binBoundaries = binBoundaries
#         container.crosssection  = crosssection
#         container.ratios        = ratios

#         # Load values of coupling into container
#         for key, value in AgnieszkasFilenameDecoder[shortFileName].iteritems():
#             setattr( container, key, value )

#         return container


#     lowRangeContainers = [ f(theoryFile) for theoryFile in lowRangeFiles ]
#     highRangeContainers = [ f(theoryFile) for theoryFile in highRangeFiles ]

#     lowhighPairs = []
#     for lowContainer in lowRangeContainers:
#         for highContainer in highRangeContainers:
#             try:
#                 if (
#                     lowContainer.ct == highContainer.ct
#                     and lowContainer.cg == highContainer.cg
#                     and lowContainer.cb == highContainer.cb
#                     ):
#                     print 'Matched {0} with {1}'.format( lowContainer.file, highContainer.file )
#                     print '   ( ct = {0}, cg = {1}, cb = {2} )'.format( lowContainer.ct, lowContainer.cg, lowContainer.cb )
#                     break
#             except AttributeError:
#                 print 'File {0} or {1} has a coupling problem'.format( lowContainer.file, highContainer.file )
#                 raise
#         else:
#             print 'Could not find matching high range container for {0}'.format( lowContainer.file )
#             continue

#         lowhighPairs.append( [ lowContainer, highContainer ] )


#     outdir = 'derivedTheoryFiles_{0}_TopHighPt'.format( datestr )
#     if not isdir( outdir ): os.makedirs( outdir )

#     print '\n'
#     print 'SM.file:'
#     print SM.file
#     print 'SM.ratios ({0}):'.format( len(SM.ratios) )
#     print SM.ratios
#     print 'SM.crosssection ({0}):'.format( len(SM.crosssection) )
#     print SM.crosssection

#     print '\n' + '='*80
#     for low, high in lowhighPairs:

#         print '\nMerging low {0} with high {1}'.format( low.file, high.file )

#         merged = deepcopy( high )
#         merged.secondfile = low.file

#         assert low.binBoundaries[-1] == high.binBoundaries[0]
#         merged.binBoundaries = low.binBoundaries[:-1] + high.binBoundaries
#         merged.crosssection  = low.crosssection + high.crosssection
#         merged.ratios         = [ xs / SMxs if not SMxs == 0. else 0. for xs, SMxs in zip( merged.crosssection, SM.crosssection ) ]

#         # print len(merged.binBoundaries)
#         # print merged.binBoundaries
#         # print len(merged.crosssection)
#         # print merged.crosssection

#         # print
#         # print len(SM.binBoundaries)
#         # print SM.binBoundaries
#         # print len(SM.crosssection)
#         # print SM.crosssection

#         DumpContainerToFile_Top( merged, prefix='TopHighPt', outdir=outdir )
#     DumpContainerToFile_Top( SM, prefix='TopHighPt', outdir=outdir )

