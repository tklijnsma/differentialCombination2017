#!/usr/bin/env python
"""
Thomas Klijnsma
"""

import differentials
import differentials.core as core

import logging
import os, sys
from math import sqrt

import collections





class BaseTable(collections.Sequence):
    """docstring for BaseTable"""

    tab_sep = '|'
    line_sep = '\n'

    begin_line_with_tab_sep = False
    end_line_with_tab_sep = False

    def __init__(self):
        super(BaseTable, self).__init__()
        self.rows = []

    def __getitem__(self, i):
        return self.rows[i]

    def __len__(self):
        return len(self.rows)

    def __repr__(self):
        return 'Table(' + ', '.join([ row.__repr__() for row in self.rows ]) + ')'

    def append(self, i):
        self.rows.append(i)

    def extend(self, i):
        self.rows.extend(i)


    def cells(self):
        for row in self.rows:
            for cell in row.cells:
                yield cell
        
    def get_n_cols(self):
        row_lens = [len(row) for row in self.rows]
        n_col_set = set(row_lens)
        if len(n_col_set) > 1:
            raise ValueError(
                'Conflicting column dimensions: {0}'
                .format(row_lens)
                )
        (n_cols,) = n_col_set
        return n_cols

    def get_max_col_widths(self):
        n_rows = len(self.rows)
        n_cols = self.get_n_cols()

        # First get the max colls for only span-1 cells
        max_col_widths = []
        for i_col in xrange(n_cols):
            max_col_width = 0
            for i_row in xrange(n_rows):
                cell = self[i_row][i_col]
                if cell.span != 1: continue
                max_col_width = max(max_col_width, cell.get_max_width())
            max_col_widths.append(max_col_width)

        # Adjust in case multi-span cells are very filled
        for i_row in xrange(n_rows):
            col_iterator = iter(range(n_cols))
            for i_col in col_iterator:
                cell = self[i_row][i_col]
                if cell.span <= 1: continue

                sum_of_cols_spanned = sum(max_col_widths[i_col:i_col+cell.span])
                len_internal_seps = (cell.span-1) * len(self.tab_sep)

                width_of_cell = cell.get_max_width()

                if width_of_cell > (sum_of_cols_spanned + len_internal_seps):
                    # Broaden the rightmost cell to fit the multispan
                    max_col_widths[i_col+cell.span-1] += width_of_cell - (sum_of_cols_spanned + len_internal_seps)

                # Advance multiple columns to not cover the same multispan cell twice
                for i in xrange(cell.span-1):
                    next(col_iterator)

        return max_col_widths


    def latex_mode(self, activate=True):
        self.latex_mode = activate
        differentials.plotting.tableproducer.CellAsymmUncCrossSection.latex_mode = activate
        self.tab_sep = ' & '
        self.line_sep = ' \\\\ \n \\hline \n'


    def produce_formatted_cells(self):
        n_rows = len(self.rows)
        n_cols = self.get_n_cols()
        max_col_widths = self.get_max_col_widths()

        table = []
        for row in self.rows:
            max_n_lines = row.get_max_n_lines()
            lines = [ [] for i in xrange(max_n_lines) ]
            
            col_iterator = iter(range(n_cols))
            for i_col in col_iterator:
                cell = row[i_col]

                width = sum(max_col_widths[i_col:i_col+cell.span])

                if cell.span > 1 and cell.get_max_width() < width:
                    # If multispan but not the longest string, account for tab seps
                    width += (cell.span-1) * len(self.tab_sep)

                output = cell.represent_at_width(width=width, n_rows=max_n_lines)
                for i_line, output_per_line in enumerate(output):
                    lines[i_line].append(output_per_line)

                for i in xrange(cell.span-1):
                    next(col_iterator)

            table.extend(lines)

        return table

    def produce_table_string(self):
        table = self.produce_formatted_cells()
        lines = []
        for components in table:
            line = self.tab_sep.join(components)
            if self.begin_line_with_tab_sep:
                line = self.tab_sep + line
            if self.end_line_with_tab_sep:
                line += self.tab_sep
            lines.append(line)
        return self.line_sep.join(lines)

    def produce_to_file(self, outname):
        dirname = os.path.dirname(outname)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        with open(outname, 'w') as fp:
            fp.write(self.produce_table_string())


class BaseRow(collections.Sequence):
    """docstring for BaseRow"""

    def __init__(self):
        super(BaseRow, self).__init__()
        self.cells = []

    def __getitem__(self, i_col):
        if i_col >= len(self):
            raise IndexError('{0} is outside the range {1}'.format(i_col, len(self)))
        for i_function_call, cell in enumerate(self.cells_per_column()):
            if i_function_call == i_col:
                return cell

    def __len__(self):
        return sum([c.span for c in self.cells])

    def __repr__(self):
        return 'Row' + self.cells.__repr__()

    def cells_per_column(self):
        for cell in self.cells:
            for i in xrange(cell.span):
                yield cell

    def append(self, i):
        self.cells.append(i)

    def extend(self, i):
        self.cells.extend(i)

    def get_max_col_widths(self):
        return [ cell.get_max_width() for cell in self.cells ]

    def get_max_n_lines(self):
        return max([cell.get_n_lines() for cell in self.cells])


class BaseCell(object):
    """docstring for BaseCell"""

    align = '<'

    def __init__(self):
        super(BaseCell, self).__init__()
        self.span = 1
        
    def represent(self):
        return ['']

    def represent_at_width(self, width=None, n_rows=None):
        if width is None:
            width = self.get_max_width()
        r = [ '{0:{align}{width}}'.format(l, width=width, align=self.align) for l in self.represent() ]
        if not(n_rows is None):
            L = len(r)
            if n_rows > L:
                for i in xrange(n_rows-L):
                    r.append(' '*width)
        return r

    def get_max_width(self):
        max_width = 0
        for repr_str in self.represent():
            if len(repr_str) > max_width:
                max_width = len(repr_str)
        return max_width

    def get_min_width(self):
        min_width = 10e5
        for repr_str in self.represent():
            if len(repr_str) < min_width:
                min_width = len(repr_str)
        return min_width

    def get_n_lines(self):
        return len(self.represent())


class CellString(BaseCell):
    """docstring for CellString"""
    def __init__(self, value):
        super(CellString, self).__init__()
        self.value = value
        
    def represent(self):
        if isinstance(self.value, basestring):
            return [ str(self.value) ]
        else:
            return [ str(v) for v in self.value ]


