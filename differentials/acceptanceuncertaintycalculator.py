
class AcceptanceUncertaintyCalculator(object):
    """docstring for AcceptanceUncertaintyCalculator"""
    def __init__(self, npz_file):
        super(AcceptanceUncertaintyCalculator, self).__init__()
        self.npz_file = npz_file
        self.npz = numpy.load(self.npz_file)
        self.reinit_lists()

    def reinit_lists(self):
        self.central = []
        self.up = []
        self.down = []
        self.symm = []

    def get_index(self, binstr):
        return int(re.match(r'bin(\d+)', binstr).group(1))

    def get_keys(self):
        r = self.npz.keys()
        r.sort(key=lambda k: self.get_index(k))
        return r

    # def binstr_to_proper_str(self, binstr):
    #     index = self.get_index(binstr)
    #     boundaries = [ 0., 15., 30., 45., 80., 120., 200., 350., 600. ]
    #     if index < 8:
    #         return '[{0:0d},{1:0d})'.format(boundaries[index], boundaries[index+1])
    #     else:
    #         return '[{0:0d},#infty)'.format(boundaries[index])

    def get_unc_for_bin(self, A):
        B = list(A[:].flatten())
        logging.debug(B)
        B.pop(7) # Remove ratio 4 scale variations (higher index first to not mess up the next pop)
        B.pop(5) # Remove ratio 4 scale variations
        central = B[0]
        up   = abs(max(B))/central - 1
        down = 1 - abs(min(B))/central
        symm = 0.5*(abs(down)+abs(up))
        return central, down, up, symm

    def get_cud(self):
        self.reinit_lists()
        for k in self.get_keys():
            logging.debug('Doing key {0}'.format(k))
            central, down, up, symm = self.get_unc_for_bin(self.npz[k])
            self.central.append(central)
            self.up.append(up)
            self.down.append(down)
            self.symm.append(symm)
        
    def get_central_shape(self):
        self.get_cud()
        return self.central

