import logging
import core
from core import AttrDict

import os, glob, re, copy, numpy
# import ROOT



class Rebinner(object):

    def __init__(self, bin_boundaries_old, values_old, bin_boundaries_new):
        super(Rebinner, self).__init__()
        if len(bin_boundaries_old)-1 != len(values_old):
            raise ValueError('Length of given input lists do not match')
        self.bin_boundaries_old = bin_boundaries_old
        self.values_old = values_old
        self.bin_boundaries_new = bin_boundaries_new

    def rebin(self):
        self.n_bins_old = len(self.bin_boundaries_old)-1
        self.n_bins_new = len(self.bin_boundaries_new)-1
        self.integral = Integrator(self.bin_boundaries_old, self.values_old)
        self.values_new = []
        lefts = self.bin_boundaries_new[:-1]
        rights = self.bin_boundaries_new[1:]
        for a, b in zip(lefts, rights):
            self.values_new.append(
                self.integral.integral(a, b) / (b-a)
                )
        return self.values_new


class Integrator(object):

    def __init__(self, bin_boundaries, values):
        super(Integrator, self).__init__()
        if len(bin_boundaries)-1 != len(values):
            raise ValueError('Length of given input lists do not match')
        self.bin_boundaries = bin_boundaries
        self.values = values
        self.n_bins = len(self.bin_boundaries)-1
        self.bin_centers = [0.5*(self.bin_boundaries[i+1] + self.bin_boundaries[i]) for i in xrange(self.n_bins)]

        self.x_min = self.bin_boundaries[0]
        self.x_max = self.bin_boundaries[-1]

        self.set_zero_outside_range()


    def set_zero_outside_range(self):
        self.behavior = 'zero_outside_range'


    def get_underflow(self, a, b):
        if self.behavior == 'zero_outside_range':
            return 0.

    def get_overflow(self, a, b):
        if self.behavior == 'zero_outside_range':
            return 0.

    def get_partial_bin_contribution(self, i, x):
        left = self.bin_boundaries[i]
        right = self.bin_boundaries[i+1]
        value = self.values[i]
        width_left = x - left
        width_right = right - x
        return width_left*value, width_right*value


    def integral(self, a, b):
        ret = 0.0

        logging.debug('Requested integral from a={0} to b={1}'.format(a, b))
        logging.debug('Bin boundaries are {0}'.format(self.bin_boundaries))

        _multiply_by_minus_one = False
        if a == b:
            logging.info('Found a==b; returning zero')
            return 0.
        elif b < a:
            logging.info('Found b < a; will calculate integral from a -> b and multiply by -1 later')
            _ = b
            b = a
            a = _
            _multiply_by_minus_one = True

        logging.debug('Requested behavior is {0}'.format(self.behavior))
        if a < self.x_min:
            logging.debug('Requested left bound {0} < x_min({1})'.format(a, self.x_min))
        if b > self.x_max:
            logging.debug('Requested right bound {0} > x_max({1}); assuming zero above x_max'.format(b, self.x_max))

        ret += self.get_underflow(a, b)

        # New lists to contain actual bin boundaries and values for integration
        _bin_boundaries = []
        _values = []

        for i in xrange(self.n_bins):
            left = self.bin_boundaries[i]
            right = self.bin_boundaries[i+1]
            value = self.values[i]

            _value_added = False
            if left < a and a < right:
                _bin_boundaries.append(a)
                _values.append(value)
                _value_added = True
            elif a <= left and left < b:
                _bin_boundaries.append(left)
                _values.append(value)
                _value_added = True

            if left < b and b <= right:
                _bin_boundaries.append(b)
                if not _value_added:
                    _values.append(value)
                break

        logging.debug('The bin boundaries over which to integrate are: {0}'.format(_bin_boundaries))
        logging.debug('The bin values over which to integrate are: {0}'.format(_values))


        _lefts  = _bin_boundaries[:-1]
        _rights = _bin_boundaries[1:]
        _widths = [ r-l for l, r in zip(_lefts, _rights) ]

        logging.debug('Start evaluating integral')
        logging.debug('underflow contribution:' + ' '*28 + '{0:<+9.3f}'.format(ret))
        for value, width, left, right in zip(_values, _widths, _lefts, _rights):
            contribution = value * width
            logging.debug(
                'val = {0:<+9.3f}, width = {1:<+9.3f}, contribution = {2:<+9.3f} '
                '(l={3}, r={4})'
                .format(value, width, contribution, left, right)
                )
            ret += contribution

        overflow = self.get_overflow(a, b)
        logging.debug('overflow contribution:' + ' '*29 + '{0:<+9.3f}'.format(overflow))
        ret += overflow

        logging.debug('Result: Integral from a={0:<+9.3f} to b={1:<+9.3f} is {2:<+9.3f}'.format(a, b, ret))
        if _multiply_by_minus_one:
            logging.info('Multiplying by -1.0 to correct for the a<>b switch')
            ret *= -1.0
        return ret

