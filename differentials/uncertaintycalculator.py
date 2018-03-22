import core
import logging

import ROOT

from collections import namedtuple
from array import array

def rindex( someList, val ):
    # Regular list.index() finds first instance in list, this function finds the last
    return len(someList) - someList[::-1].index(val) - 1

class UncertaintyCalculator(object):
    """docstring for UncertaintyCalculator"""
    def __init__(self):
        super(UncertaintyCalculator, self).__init__()
        self.cutoff = 0.5

    def create_uncertainties(self, xs, deltaNLLs):
        
        min_deltaNLL   = min(deltaNLLs)
        i_min_deltaNLL = rindex(deltaNLLs, min_deltaNLL)
        x_min          = xs[i_min_deltaNLL]

        logging.debug('Found minimum at index {0}: x={1}, deltaNLL={2}'.format(i_min_deltaNLL, x_min, min_deltaNLL))

        unc_dict = {
            'min_deltaNLL' : min_deltaNLL,
            'i_min' : i_min_deltaNLL,
            'x_min' : x_min,
            'left_bound' : -0,
            'left_error' : -0,
            'right_bound' : 0,
            'right_error' : 0,
            'symm_error' : 0,
            'well_defined_left_bound' : False,
            'well_defined_right_bound' : False,
            'is_hopeless' : True,
            'cutoff_1sigma' : self.cutoff,
            }
        Unc = namedtuple('Unc', unc_dict.keys())

        # Process left uncertainty
        if i_min_deltaNLL < 3:
            logging.debug('Not enough points on the left side for a well defined left bound')
            well_defined_left_bound = False
        else:
            xs_left = xs[:i_min_deltaNLL+1]
            deltaNLLs_left = deltaNLLs[:i_min_deltaNLL+1]
            if min(deltaNLLs_left) > self.cutoff or max(deltaNLLs_left) < self.cutoff:
                logging.debug('Requested dNLL interpolation point is outside the range: min dNLL={0}, max dNLL={1}'.format(min(deltaNLLs_left), max(deltaNLLs_left)))
                well_defined_left_bound = False
            else:
                left_bound = self.interpolate(xs_left, deltaNLLs_left, self.cutoff)
                if left_bound is False:
                    well_defined_left_bound = False
                else:
                    well_defined_left_bound = True

        # Process right uncertainty
        if i_min_deltaNLL > len(xs)-3:
            logging.debug('Not enough points on the right side for a well defined right bound')
            well_defined_right_bound = False
        else:
            xs_right = xs[i_min_deltaNLL:]
            deltaNLLs_right = deltaNLLs[i_min_deltaNLL:]
            if min(deltaNLLs_right) > self.cutoff or max(deltaNLLs_right) < self.cutoff:
                logging.debug('Requested dNLL interpolation point is outside the range: min dNLL={0}, max dNLL={1}'.format(min(deltaNLLs_right), max(deltaNLLs_right)))
                well_defined_right_bound = False
            else:
                right_bound = self.interpolate(xs_right, deltaNLLs_right, self.cutoff)
                if right_bound is False:
                    well_defined_right_bound = False
                else:
                    well_defined_right_bound = True

        is_hopeless = False
        if well_defined_left_bound and well_defined_right_bound:
            pass
        elif well_defined_left_bound and not well_defined_right_bound:
            right_bound = x_min + (x_min - left_bound)
        elif well_defined_right_bound and not well_defined_left_bound:
            left_bound  = x_min - (right_bound - x_min)
        else:
            is_hopeless = True

        if not is_hopeless:
            unc_dict['well_defined_left_bound'] = well_defined_left_bound
            unc_dict['well_defined_right_bound'] = well_defined_right_bound
            left_error = abs(x_min - left_bound)
            right_error = abs(x_min - right_bound)
            unc_dict['left_bound']  = left_bound
            unc_dict['left_error']  = left_error
            unc_dict['right_bound'] = right_bound
            unc_dict['right_error'] = right_error
            unc_dict['symm_error']  = 0.5*(abs(left_error)+abs(right_error))
            unc_dict['is_hopeless'] = False

        unc = Unc(**unc_dict)
        return unc

    def interpolate(self, ys, xs, x_value):
        logging.debug('Interpolating for x_value={0}'.format(x_value))
        logging.trace('  x  /  y:')
        for x, y in zip(xs, ys): logging.trace('    {0:+7.2f}  /  {1:+7.2f}'.format(x, y))

        if min(xs) > x_value or max(xs) < x_value:
            logging.debug('  Requested interpolation {0} is outside the range: {1} to {2}'.format(x_value, min(xs), max(xs)))
            return False
        Tg = ROOT.TGraph(len(xs), array('f', xs), array('f', ys))
        y_value = Tg.Eval(x_value)
        # if y_value < min(ys) or y_value > max(ys):
        #     return False
        if y_value is False:
            logging.debug('  Could not interpolate properly')
        else:
            logging.debug('  Interpolated y_value {0} for x_value {1}'.format(y_value, x_value))
        return y_value

