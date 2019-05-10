import differentials
import os

class Variable(object):
    """docstring for Variable"""
    def __init__(self, title, unit, is_dependent=True):
        super(Variable, self).__init__()
        self.title = title
        self.unit = unit
        self.is_dependent = is_dependent
        self.values = []
        self.errs_up = []
        self.errs_down = []
        self.qualifiers = []

    def add_val(self, val):
        self.values.append(val)

    def add_val_errs(self, val, up, down):
        self.values.append(val)
        self.errs_up.append(up)
        self.errs_down.append(down)

    def add_qualifier(self, name, unit,  val):
        self.qualifiers.append({
            'name' : name, 'units' : unit, 'value': val
            })

    def parse(self):

        out = []

        out.append(
            '- header: {{name: \'{0}\', units: \'{1}\'}}'
            .format(self.title, self.unit)
            )

        if len(self.qualifiers) > 0:
            out.append('  qualifiers:')
            for q in self.qualifiers:

                dictstring = 'name: \'{0}\''.format(q['name'])
                if 'units' in q:
                    dictstring += ', units: \'{0}\''.format(q['units'])

                try:
                    val = float(q['value'])
                    if val == int(val):
                        val = int(val)
                except ValueError:
                    val = "'" + q['value'] + "'"
                dictstring += ', value: {0}'.format(val)

                out.append(
                    '  - {{{dictstring}}}'.format(dictstring=dictstring)
                    )

        out.append('  values:')

        if len(self.errs_down) == len(self.values):
            out.extend(self.parse_values_errors())
        else:
            out.extend(self.parse_values_single())

        return out


    def parse_values_errors(self):
        out = []
        for i in xrange(len(self.values)):
            out.append('  - errors:')
            out.append(
                '    - asymerror: {{minus: {0}, plus: {1}}}'
                .format(self.errs_down[i], self.errs_up[i])
                )
            out.append('      label: 1 s.d.')
            out.append('    value: {0}'.format(self.values[i]))
        return out

    def parse_values_single(self):
        out = []
        for val in self.values:
            out.append('  - {{value: {0}}}'.format(val))
        return out




class HepDataMaker(object):
    """docstring for HepDataMaker"""

    channel_titles = {
        'hgg' : '$H \\to \\gamma\\gamma$',
        'hzz' : '$H \\to ZZ \\to 4 \\ell$',
        'hbb' : '$H \\to b\\overline{b}$',
        'combination' : 'Combination',
        'combWithHbb' : 'Combination'
        }

    obsname_titles = {
        'pth_smH'  : '$p_T^H$',
        'pth_ggH'  : '$p_T^H$ (ggH)',
        'njets'    : '$N_\\text{jets}$',
        'ptjet'    : '$p_T^{j}$',
        'rapidity' : '$\\abs{\\eta_H}$',
        }

    obsname_unit = {
        'pth_smH'  : 'GeV',
        'pth_ggH'  : 'GeV',
        'njets'    : 'NA',
        'ptjet'    : 'GeV',
        'rapidity' : 'NA',
        }


    def __init__(self):
        super(HepDataMaker, self).__init__()

        self.left_var = Variable(
            'Left bin bound', 'GeV', is_dependent = False
            )
        self.right_var = Variable(
            'Right bin bound', 'GeV', is_dependent = False
            )
        self.single_var = Variable(
            'Bin', 'GeV', is_dependent = False
            )

        self.cross_section = Variable(
            'Cross section', 'pb/bin-width'
            )

        self.signal_strength = Variable(
            'Signal strength', 'NA'
            )

        self.qualifiers = [
            {'name': '$\\sqrt{s}$', 'units': 'TeV', 'value': '13'},
            {'name': '$L_{\\mathrm{int}}$', 'units': 'fb$^{-1}$', 'value': '35.9'},
            ]
        self.cross_section.qualifiers = self.qualifiers
        self.signal_strength.qualifiers = self.qualifiers

        self.use_single_var = False


    def from_spectrum(self, spectrum, obsname=None):
        self.spectrum = spectrum
        self.obsname = obsname

        q  = {'name': 'Channel', 'value': self.channel_titles[spectrum.name]}
        self.qualifiers.insert(0, q)

        if obsname:
            q = {
                'name': 'Observable',
                'value': self.obsname_titles[obsname],
                'units': self.obsname_unit[obsname],
                }
            self.qualifiers.insert(0, q)
            self.left_var.unit = self.obsname_unit[obsname]
            self.right_var.unit = self.obsname_unit[obsname]
            self.single_var.unit = self.obsname_unit[obsname]

        bin_boundaries = spectrum.binning()

        for i_scan, scan in enumerate(spectrum.scans):
            left     = bin_boundaries[i_scan]
            right    = bin_boundaries[i_scan+1]

            self.left_var.add_val(left)
            self.right_var.add_val(right)

            center = scan.unc.x_min
            up     = scan.unc.right_error
            down   = scan.unc.left_error

            self.signal_strength.add_val_errs(
                center, up, down
                )

            xs = spectrum.smxs[i_scan]
            bin_width = right-left
            center *= xs
            up *= xs
            down *= xs

            self.cross_section.add_val_errs(
                center, up, down
                )

        return self


    def parse_yaml(self):
        out = []
        out.append('dependent_variables:')
        out.extend(self.cross_section.parse())
        out.extend(self.signal_strength.parse())
        out.append('independent_variables:')
        if self.use_single_var:
            out.extend(self.single_var.parse())
        else:
            out.extend(self.left_var.parse())
            out.extend(self.right_var.parse())
        return '\n'.join(out)


    def parse_yaml_to_file(self, outname):
        dirname = os.path.dirname(outname)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        with open(outname, 'w') as fp:
            fp.write(self.parse_yaml())






