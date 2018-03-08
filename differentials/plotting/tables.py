#!/usr/bin/env python
"""
Thomas Klijnsma
"""

import differentials
import differentials.core as core

import logging
import os, sys
from math import sqrt


def format_number(number, decimals=2, force_sign=False, verbose_log=False):
    sign = '+' if force_sign else ''
    if isinstance(number, int) or number.is_integer():
        ret =  str(int(number))
    else:
        basic_float_str = '{0:{sign}.{decimals}f}'.format(number, decimals=decimals, sign=sign)
        float_str = "{0:{sign}.{decimals}E}".format(number, decimals=decimals, sign=sign)
        if "E" in float_str:
            base, exponent = float_str.split("E")
            exponent = int(exponent)
            if exponent <= -3:
                if verbose_log:
                    ret = "{{{0} \\times 10^{{{1}}}}}".format(base, exponent)
                else:
                    ret = float_str
            else:
                ret = basic_float_str
        else:
            ret = basic_float_str
    return ret


class Cell(object):
    """docstring for Bin"""
    def __init__(self):
        self.center    = 0.0
        self.err_up    = 0.0
        self.err_down  = 0.0
        self.stat_up   = 0.0
        self.stat_down = 0.0
        self.syst_up   = 0.0
        self.syst_down = 0.0

    def symm_err(self):
        return 0.5*(abs(self.err_up)+abs(self.err_down))

    def absolutize(self):
        absol = lambda number: -abs(number)
        self.err_down = absol(self.err_down)
        self.stat_down = absol(self.stat_down)
        self.syst_down = absol(self.syst_down)

    def repr_terminal(self):
        self.absolutize()
        f = lambda n: '{0:+.2f}'.format(n)
        ret = (
            '{center}  {up}/{down}'
            .format(
                center = f(self.center),
                down   = f(self.err_down),
                up     = f(self.err_up),
                )
            )
        return ret


class SymmetricImprovementCell(Cell):
    def __init__(self):
        super(SymmetricImprovementCell, self).__init__()
        self.symm_improvement = 0.0

    def repr_terminal(self):
        return '{0:.1f}%'.format(100.*self.symm_improvement)





    # def repr_terminal(self):
    #     f = lambda n: '{0:+.2f}'.format(n)
    #     r = '{0}  {1}/{2} (stat.)  {3}/{4} (syst.)'.format(
    #         f(self.center),
    #         f(self.stat_up),
    #         f(-self.stat_down),
    #         f(self.syst_up),
    #         f(-self.syst_down),
    #         )
    #     return r

    # def repr_twiki(self, do_xs=False):
    #     # | One One | One Two | One Three |
    #     # | ^ | Two Two | Two Three |
    #     # | Three One | ^ | Three Three |

    #     # '| center | +stat_up +syst_up |'
    #     # '| ^      | -stat_up -syst_up |'

    #     f = lambda n: format_number(n, force_sign=True)
    #     if do_xs:
    #         f = lambda n: format_number(n*self.SMXS, force_sign=True)

    #     # r1 = '{0:8} | {1:8} {2:8} '.format(f(self.center), f(self.stat_up), f(self.syst_up))
    #     # r2 = '{0:8} | {1:8} {2:8} '.format('^', f(self.stat_down), f(self.syst_down))
    #     r1 = '{0} | {1} {2} '.format(f(self.center), f(self.stat_up), f(self.syst_up))
    #     r2 = '{0} | {1} {2} '.format('^', f(-self.stat_down), f(-self.syst_down))

    #     return r1, r2

    # def repr_twiki_symm(self, do_xs=False):
    #     f = lambda n: format_number(n, force_sign=False)
    #     if do_xs:
    #         f = lambda n: format_number(n*self.SMXS, force_sign=False)
    #     r1 = f(0.5*(abs(self.stat_up)+abs(self.stat_down)))
    #     r2 = f(0.5*(abs(self.syst_up)+abs(self.syst_down)))
    #     return r1, r2


class Row(object):
    def __init__(self, table=None):
        self.name = ''
        self.cells  = []
        self.bin_boundaries = []
        self.table = table  # Pointer to parent

    def append(self, i):
        self.row.append(i)

    def from_spectrum(self, spectrum):
        self.name = spectrum.name
        self.bin_boundaries = spectrum.binning()
        self.n_bins = len(self.bin_boundaries)-1
        for i_scan, scan in enumerate(spectrum.scans):
            cell = Cell()
            cell.left     = self.bin_boundaries[i_scan]
            cell.right    = self.bin_boundaries[i_scan+1]
            cell.center   = scan.unc.x_min
            cell.err_up   = scan.unc.right_error
            cell.err_down = scan.unc.left_error
            if self.table.do_xs:
                xs = spectrum.smxs[i_scan]
                cell.center *= xs
                cell.err_up *= xs
                cell.err_down *= xs
            self.cells.append(cell)


class SpectraTable(object):
    """docstring for SpectraTable"""
    def __init__(self, name, spectra=None, do_xs=False, last_bin_is_overflow=False):
        self.name = name
        self.last_bin_is_overflow = last_bin_is_overflow
        self.do_xs = do_xs
        self.header = []
        self.rows = []
        if not(spectra is None):
            self.spectra = spectra
            self.read_spectra()

    def read_spectra(self):
        for spectrum in self.spectra:
            logging.debug('Reading spectrum {0}'.format(spectrum))
            row = Row(self)
            row.from_spectrum(spectrum)
            self.rows.append(row)
        self.get_bin_boundaries()
        self.make_header()

    def add_symm_improvement_row(self, spectrum_old, spectrum_new):
        row_old = Row(self)
        row_old.from_spectrum(spectrum_old)
        row_new = Row(self)
        row_new.from_spectrum(spectrum_new)
        if not(row_old.n_bins == row_new.n_bins) or row_old.n_bins==0:
            raise ValueError('Inconsistent length of bins or zero')

        row_symm_improv = Row(self)
        row_symm_improv.name = 'Improv. {0} > {1}'.format(spectrum_old.name, spectrum_new.name)
        for i_bin in xrange(row_old.n_bins):
            err_old = row_old.cells[i_bin].symm_err()
            err_new = row_new.cells[i_bin].symm_err()
            symm_cell = SymmetricImprovementCell()
            symm_cell.left = row_old.bin_boundaries[i_bin]
            symm_cell.right = row_old.bin_boundaries[i_bin+1]
            try:
                symm_cell.symm_improvement = 1. - err_new/err_old
            except ZeroDivisionError:
                symm_cell.symm_improvement = -9.
            row_symm_improv.cells.append(symm_cell)
        self.rows.append(row_symm_improv)


    def make_header(self):
        self.header = [ ' ' + self.name + ' ' ]
        for i_bin in xrange(len(self.bin_boundaries)-1):
            left = format_number(self.bin_boundaries[i_bin], decimals=1)
            right = format_number(self.bin_boundaries[i_bin+1], decimals=1)
            self.header.append( ' {0}-{1} '.format(left, right) )

    def get_bin_boundaries(self):
        n_bins_max = 0
        for row in self.rows:
            if row.n_bins > n_bins_max:
                bin_boundaries = row.bin_boundaries
                n_bins_max = row.n_bins
        self.bin_boundaries = bin_boundaries
        logging.info('Set bin_boundaries={0}'.format(self.bin_boundaries))

    def get_degree_of_multispan(self, cell):
        i_left = self.bin_boundaries.index(cell.left)
        i_right = self.bin_boundaries.index(cell.right)
        degree = i_right - i_left
        return degree

    def repr_terminal(self):
        table = [ self.header ]
        for row in self.rows:
            repr_row = [row.name]
            for cell in row.cells:
                repr_row.append(cell.repr_terminal())
                doms = self.get_degree_of_multispan(cell)
                for i_ms in xrange(doms-1):
                    repr_row.append(' ')
            table.append(repr_row)

        return repr_rectangular_table(table)




    # def repr_twiki(self):
    #     rows = []

    #     bin_boundaries = self.get_bin_boundaries()

    #     header = self.make_header(bin_boundaries)
    #     header_row =  '|' + header[0] + '|' + '||'.join(header[1:]) + '||'
    #     rows.append(header_row)

    #     # for row in self.rows:
    #     #     subrow1 = [ row.name ]
    #     #     subrow2 = [ row.name ]
    #     #     for cell in row.row:
    #     #         r1, r2 = cell.repr_twiki()
    #     #         subrow1.append(r1)
    #     #         subrow2.append(r2)
    #     #     rows.append('|' + ' | '.join(subrow1) + '|')
    #     #     rows.append('|' + ' | '.join(subrow2) + '|')

    #     for row in self.rows:
    #         subrow1 = '| ' + row.name + ' | '
    #         subrow2 = '| ' + row.name + ' | '
    #         for cell in row.row:
    #             r1, r2 = cell.repr_twiki(do_xs=self.do_xs)
    #             d = self.get_degree_of_multispan(cell, bin_boundaries)
    #             r1 = r1.replace('|', (2*d-1)*'|')
    #             r2 = r2.replace('|', (2*d-1)*'|')
    #             subrow1 += r1 
    #             subrow2 += r2 
    #             subrow1 += ' | '
    #             subrow2 += ' | '

    #         rows.append(subrow1)
    #         rows.append(subrow2)

    #     return '\n'.join(rows)

    # def repr_twiki_symm(self):
    #     rows = []


    #     bin_boundaries = self.get_bin_boundaries()

    #     header = self.make_header(bin_boundaries)
    #     header_row = '|' + '|'.join(header) + '|'
    #     rows.append(header_row)


    #     for row in self.rows:
    #         subrow1 = '| ' + row.name + ' (stat.) | '
    #         subrow2 = '| ' + row.name + ' (syst.) | '
    #         for cell in row.row:
    #             r1, r2 = cell.repr_twiki_symm(do_xs=self.do_xs)
    #             d = self.get_degree_of_multispan(cell, bin_boundaries)
    #             subrow1 += '   ' + r1 + '   '
    #             subrow2 += '   ' + r2 + '   '
    #             subrow1 += '|'*d
    #             subrow2 += '|'*d
    #         rows.append(subrow1)
    #         rows.append(subrow2)

    #     return '\n'.join(rows)


def assert_table_is_rectangular(table):
    n_cols = len(table[0])
    for i_row, row in enumerate(table):
        if len(row) != n_cols:
            logging.error(
                'Error in row {0}; {1} columns instead of {2}. Printing table per line:\n{3}'
                .format(
                    i_row, len(row), n_cols,
                    '\n'.join([ '{0}'.format(i) for i in table ])
                    )
                )
            raise ValueError('Table has inconsistent columns')


def repr_rectangular_table(table, min_col_width=1, max_col_width=20, sep='  ', newline_sep='\n'):
    assert_table_is_rectangular(table)
    n_rows = len(table)
    n_cols = len(table[0])

    max_col_widths = []
    for i_col in xrange(n_cols):
        max_width = 0
        for i_row in xrange(n_rows):
            entry = table[i_row][i_col]
            if not isinstance(entry, basestring):
                raise TypeError('Entry of a table is not a string: {0}'.format(entry))
            if len(entry) > max_width:
                max_width = len(entry)
        max_col_widths.append( min( max_width, max_col_width ) )

    out = []
    for i_row in xrange(n_rows):
        line = []
        for i_col in xrange(n_cols):
            entry = table[i_row][i_col]
            line.append( format_str_to_width( entry, max_col_widths[i_col] ) )
        out.append( sep.join(line) )
    return newline_sep.join(out)



def format_str_to_width(text, width):
    if len(text) > width:
        text = text[:width-3] + '...'
    ret = '{0:{width}}'.format( text, width=width )
    return ret
