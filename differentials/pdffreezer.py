import logging
import differentials
import core

import ROOT

class PDFFreezer(object):

    def __init__(self, ws=None):
        super(PDFFreezer, self).__init__()
        self.ws = ws
        self.snapshotname = 'MultiDimFit'
        self._is_read = False
        self.vars_to_freeze = []
        self.vars_to_float = []

    def get_loaded_workspace(self):
        self.w = core.get_ws(self.ws)
        is_loaded = self.w.loadSnapshot(self.snapshotname)
        if not is_loaded:
            raise RuntimeError('Workspace {0} has no snapshot named {1}'.format(self.ws, self.snapshotname))

    def add_categories_to_freeze(self):
        all_cats = self.w.allCats()
        catitr = all_cats.createIterator()
        cat = catitr.Next()
        while cat:
            if cat.GetName().startswith('pdfindex'):
                self.vars_to_freeze.append(cat.GetName())
            cat = catitr.Next()

    def add_pdf_parameters_to_freeze(self):
        all_pdfs = self.w.allPdfs()

        n_pdfs = all_pdfs.getSize()
        logging.info('Looping over {0} pdfs in the workspace'.format(n_pdfs))
        logging.debug('PDF names:')

        pdfitr = all_pdfs.createIterator()
        pdf = pdfitr.Next()
        while pdf:
            if pdf.GetName().startswith("shapeBkg_bkg"):
                # bgks from hzz are RooHistPdfs, not RooMultiPdfs
                if hasattr(pdf, 'getNumPdfs'):
                    # Loop over all shapes in the envelope
                    for ishape in xrange(pdf.getNumPdfs()):
                        shape = pdf.getPdf(ishape)
                        observables = ROOT.RooArgList(shape.getObservables(self.w.allVars()))
                        observables = filter(lambda x: not x.startswith('CMS_hgg_mass'), map(lambda x: observables[x].GetName(), xrange(observables.getSize()) ) )
                        # Freeze all pdf parameters except those from the best fit function
                        if ishape == pdf.getCurrentIndex():
                            self.vars_to_float.extend( observables )
                        else:
                            self.vars_to_freeze.extend( observables )
            pdf = pdfitr.Next()

        # For some reason, the first variable is never frozen; simply append it again at end of list
        self.vars_to_freeze.append( self.vars_to_freeze[0] )

    def read(self):
        self.get_loaded_workspace()
        self.add_categories_to_freeze()
        self.add_pdf_parameters_to_freeze()
        self._is_read = True

    def get_vars_to_freeze(self):
        if not self._is_read:
            self.read()
        return self.vars_to_freeze
