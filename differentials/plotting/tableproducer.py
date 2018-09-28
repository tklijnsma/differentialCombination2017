#!/usr/bin/env python
"""
Thomas Klijnsma
"""

import differentials
import differentials.core as core
import newtables

import logging
import os, sys, re
from math import sqrt



class CellAsymmUncCrossSection(newtables.BaseCell):
    """docstring for CellAsymmUncCrossSection"""

    n_decimals = 2
    scientific_notation = False
    latex_mode = False
    compress_maximally = False

    def __init__(self, center, up, down):
        super(CellAsymmUncCrossSection, self).__init__()
        self.center = center
        self.up = up
        self.down = down

    def num_to_strs(self):
        if self.scientific_notation:
            notation = 'e'
        else:
            notation = 'f'
        center_str = '{0:.{n_decimals}{notation}}'.format(self.center, n_decimals=self.n_decimals, notation=notation)
        up_str = '{0:+.{n_decimals}{notation}}'.format(self.up, n_decimals=self.n_decimals, notation=notation)
        down_str = '{0:+.{n_decimals}{notation}}'.format(self.down, n_decimals=self.n_decimals, notation=notation)

        if self.latex_mode and self.scientific_notation:
            center_str = self.format_10topowerof(self.center)
            up_str     = self.format_10topowerof(self.up, force_sign=True)
            down_str   = self.format_10topowerof(self.down, force_sign=True)
        return center_str, up_str, down_str



    def format_10topowerof(self, number, force_sign=False):
        sci_str = '{0:.4e}'.format(number)
        match = re.search(r'([\d\+\-\.]+)e([\d\+\-\.]+)', sci_str)
        num = float(match.group(1))
        exponent = int(match.group(2))
        return '{0:{sign}.1f} \\cdot 10^{{{1}}}'.format(num, exponent, sign = '+' if force_sign else '')    

    def format_latex(self, number, force_sign=False):
        if not self.compress_maximally:
            return self.format_10topowerof(number, force_sign)
        else:
            if abs(number) < 0.01 or abs(number) >= 1000:
                return self.format_10topowerof(number, force_sign)
            else:
                if abs(number) > 1:
                    n_decimals = 0
                elif abs(number) >= 0.1:
                    n_decimals = 1
                else:
                    n_decimals = 2
                return '{0:{sign}.{n_decimals}f}'.format(
                    number,
                    n_decimals = n_decimals,
                    sign='+' if force_sign else ''
                    )


    def represent_terminal(self):
        center_str, up_str, down_str = self.num_to_strs()
        r = [ '{0}  {1}'.format(center_str, up_str), len(center_str) * ' ' + '  ' + down_str ]
        return r

    def represent_latex(self):
        center_str, up_str, down_str = self.num_to_strs()
        r = '${0} \\, {{}}^{{{1}}}_{{{2}}}$'.format(center_str, up_str, down_str)
        if self.span > 1:
            r = '\\multicolumn{' + str(self.span) + '}{c|}{' + r + '}'
        return [r]

    def represent(self):
        if self.latex_mode:
            return self.represent_latex()
        return self.represent_terminal()



class SpectrumRowProducer(object):
    """docstring for SpectrumRowProducer"""
    def __init__(self, binning, last_bin_is_overflow=False):
        super(SpectrumRowProducer, self).__init__()
        self.binning = binning
        self.do_xs = False
        self.normalize = False
        self.last_bin_is_overflow = last_bin_is_overflow

    def int_binning_where_possible(self):
        self.binning = [ int(b) if int(b)==b else b for b in self.binning ]

    def produce_binning_row(self, title=''):
        self.int_binning_where_possible()
        row = newtables.BaseRow()
        row.append(newtables.CellString(title))
        for left, right in zip(self.binning[:-1], self.binning[1:]):
            if right == 10000. or (self.last_bin_is_overflow and right == self.binning[-1]):
                row.append(newtables.CellString('$>${0}'.format(left)))
            else:
                row.append(newtables.CellString('{0}-{1}'.format(left, right)))
        return row

    def produce_row_given_labels(self, labels):
        row = newtables.BaseRow()
        for label in labels:
            row.append(newtables.CellString(label))
        return row

    def produce(self, spectrum):
        bin_boundaries = spectrum.binning()
        n_bins = len(bin_boundaries)-1

        row = newtables.BaseRow()
        row.append(newtables.CellString(spectrum.latex_title))

        if self.do_xs:
            CellAsymmUncCrossSection.scientific_notation = True
            CellAsymmUncCrossSection.n_decimals = 3

        i_left = self.binning.index(bin_boundaries[0])
        i_right = self.binning.index(bin_boundaries[-1])

        for i in xrange(i_left):
            row.append(newtables.CellString('-'))


        if self.normalize: incl_xs = self.get_incl_xs(spectrum)

        for i_scan, scan in enumerate(spectrum.scans):
            left     = bin_boundaries[i_scan]
            right    = bin_boundaries[i_scan+1]

            center = scan.unc.x_min
            up     = scan.unc.right_error
            down   = scan.unc.left_error

            if self.do_xs:
                xs = spectrum.smxs[i_scan]
                bin_width = right-left
                if self.normalize:
                    center *= xs * bin_width / incl_xs
                    up     *= xs * bin_width / incl_xs
                    down   *= xs * bin_width / incl_xs
                else:
                    center *= xs
                    up *= xs
                    down *= xs

            cell = CellAsymmUncCrossSection(center, up, -abs(down))
            cell.span = self.binning.index(right) - self.binning.index(left)

            row.append(cell)

        for i in xrange(i_right, n_bins):
            row.append(newtables.CellString('-'))

        return row

    def get_incl_xs(self, spectrum):
        incl_xs = 0.0
        bin_widths = [ r-l for l, r in zip(spectrum.binning()[:-1], spectrum.binning()[1:]) ]
        for i_scan, scan in enumerate(spectrum.scans):
            incl_xs += scan.unc.x_min * spectrum.smxs[i_scan] * bin_widths[i_scan]
        return incl_xs





