#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys
from math import sqrt

import PlotCommands


########################################
# Main
########################################

# PlotCommands.get_uncs_from_scans



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
        self.stat_up   = 0.0
        self.stat_down = 0.0
        self.syst_up   = 0.0
        self.syst_down = 0.0

    def repr_terminal(self):
        f = lambda n: '{0:+.2f}'.format(n)
        r = '{0}  {1}/{2} (stat.)  {3}/{4} (syst.)'.format(
            f(self.center),
            f(self.stat_up),
            f(-self.stat_down),
            f(self.syst_up),
            f(-self.syst_down),
            )
        return r

    def repr_twiki(self, do_xs=False):
        # | One One | One Two | One Three |
        # | ^ | Two Two | Two Three |
        # | Three One | ^ | Three Three |

        # '| center | +stat_up +syst_up |'
        # '| ^      | -stat_up -syst_up |'

        f = lambda n: format_number(n, force_sign=True)
        if do_xs:
            f = lambda n: format_number(n*self.SMXS, force_sign=True)

        # r1 = '{0:8} | {1:8} {2:8} '.format(f(self.center), f(self.stat_up), f(self.syst_up))
        # r2 = '{0:8} | {1:8} {2:8} '.format('^', f(self.stat_down), f(self.syst_down))
        r1 = '{0} | {1} {2} '.format(f(self.center), f(self.stat_up), f(self.syst_up))
        r2 = '{0} | {1} {2} '.format('^', f(-self.stat_down), f(-self.syst_down))

        return r1, r2

    def repr_twiki_symm(self, do_xs=False):
        f = lambda n: format_number(n, force_sign=False)
        if do_xs:
            f = lambda n: format_number(n*self.SMXS, force_sign=False)
        r1 = f(0.5*(abs(self.stat_up)+abs(self.stat_down)))
        r2 = f(0.5*(abs(self.syst_up)+abs(self.syst_down)))
        return r1, r2



class Row(object):
    def __init__(self):
        self.name = ''
        self.row  = []
        self.bin_boundaries = []

    def append(self, i):
        self.row.append(i)


class DifferentialTable(object):
    """docstring for DifferentialTable"""
    def __init__(self, name=' ', do_xs=False, last_bin_is_overflow=False):
        self.name = name
        self.rows = []
        self.last_bin_is_overflow = last_bin_is_overflow
        self.do_xs = do_xs

    def calculate_stat_syst(self, name, c_statsyst, c_statonly):
        PlotCommands.get_uncs_from_scans(c_statonly)
        PlotCommands.get_uncs_from_scans(c_statsyst)

        row = Row()
        row.n_bins = c_statsyst.nBins
        row.bin_boundaries = c_statsyst.binBoundaries
        row.SMcrosssections = c_statsyst.SMcrosssections
        if self.last_bin_is_overflow:
            row.bin_boundaries[-1] = 10000.

        row.name = name
        for i_bin in xrange(row.n_bins):
            cell = Cell()
            cell.center    = c_statsyst.uncs[i_bin].min
            cell.stat_up   = abs(c_statonly.uncs[i_bin].rightError)
            cell.stat_down = abs(c_statonly.uncs[i_bin].leftError)
            cell.full_up   = abs(c_statsyst.uncs[i_bin].rightError)
            cell.full_down = abs(c_statsyst.uncs[i_bin].leftError)
            cell.syst_up   = sqrt(cell.full_up**2 - cell.stat_up**2) if not cell.stat_up > cell.full_up else 0.0
            cell.syst_down = sqrt(cell.full_down**2 - cell.stat_down**2) if not cell.stat_down > cell.full_down else 0.0
            cell.left = row.bin_boundaries[i_bin]
            cell.right = row.bin_boundaries[i_bin+1]
            cell.SMXS = row.SMcrosssections[i_bin]
            row.append(cell)
        self.rows.append(row)


    def make_header(self, bin_boundaries):
        header = [ ' ' + self.name + ' ' ]
        for i_bin in xrange(len(bin_boundaries)-1):
            left = format_number(bin_boundaries[i_bin], decimals=1)
            right = format_number(bin_boundaries[i_bin+1], decimals=1)
            header.append( ' {0}-{1} '.format(left, right) )
        return header

    def get_bin_boundaries(self):
        n_bins_max = 0
        for row in self.rows:
            if row.n_bins > n_bins_max:
                bin_boundaries = row.bin_boundaries
                n_bins_max = row.n_bins
        return bin_boundaries

    def get_degree_of_multispan(self, cell, bin_boundaries):
        i_left = bin_boundaries.index(cell.left)
        i_right = bin_boundaries.index(cell.right)
        degree = i_right - i_left
        return degree

    def repr_twiki(self):
        rows = []

        bin_boundaries = self.get_bin_boundaries()

        header = self.make_header(bin_boundaries)
        header_row =  '|' + header[0] + '|' + '||'.join(header[1:]) + '||'
        rows.append(header_row)

        # for row in self.rows:
        #     subrow1 = [ row.name ]
        #     subrow2 = [ row.name ]
        #     for cell in row.row:
        #         r1, r2 = cell.repr_twiki()
        #         subrow1.append(r1)
        #         subrow2.append(r2)
        #     rows.append('|' + ' | '.join(subrow1) + '|')
        #     rows.append('|' + ' | '.join(subrow2) + '|')

        for row in self.rows:
            subrow1 = '| ' + row.name + ' | '
            subrow2 = '| ' + row.name + ' | '
            for cell in row.row:
                r1, r2 = cell.repr_twiki(do_xs=self.do_xs)
                d = self.get_degree_of_multispan(cell, bin_boundaries)
                r1 = r1.replace('|', (2*d-1)*'|')
                r2 = r2.replace('|', (2*d-1)*'|')
                subrow1 += r1 
                subrow2 += r2 
                subrow1 += ' | '
                subrow2 += ' | '

            rows.append(subrow1)
            rows.append(subrow2)

        return '\n'.join(rows)

    def repr_twiki_symm(self):
        rows = []


        bin_boundaries = self.get_bin_boundaries()

        header = self.make_header(bin_boundaries)
        header_row = '|' + '|'.join(header) + '|'
        rows.append(header_row)


        for row in self.rows:
            subrow1 = '| ' + row.name + ' (stat.) | '
            subrow2 = '| ' + row.name + ' (syst.) | '
            for cell in row.row:
                r1, r2 = cell.repr_twiki_symm(do_xs=self.do_xs)
                d = self.get_degree_of_multispan(cell, bin_boundaries)
                subrow1 += '   ' + r1 + '   '
                subrow2 += '   ' + r2 + '   '
                subrow1 += '|'*d
                subrow2 += '|'*d
            rows.append(subrow1)
            rows.append(subrow2)

        return '\n'.join(rows)


