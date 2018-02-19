import glob, re, copy

import core, logger

from collections import namedtuple


class Scan(object):
    """docstring for Scan"""

    tree_name = 'limit'
    filter_negatives = True
    deltaNLL_threshold = -0.01

    def __init__(self, x_variable, y_variable='deltaNLL', scandir=None, globpat='*'):
        self.x_variable = x_variable
        self.y_variable = y_variable
        self.scandirs = []
        if not(scandir is None): self.scandirs.append(scandir)
        self.globpat = '*'
        self.root_files = []
        self.save_all_variables = False

    def collect_root_files(self):
        root_files = copy.copy(self.root_files)
        for scandir in self.scandirs:
            if not scandir.endswith('/'): scandir += '/'
            root_files.extend(glob.glob(scandir + self.globpat + '.root'))
        return root_files

    def get_list_of_variables_in_tree(self, root_files, accept_pat='*'):
        for root_file in root_files:
            with core.openroot(root_file) as root_fp:
                if not root_fp.GetListOfKeys().Contains('limit'):
                    found_tree = False
                else:
                    variables = []
                    root_array = root_fp.Get(self.tree_name).GetListOfBranches()
                    for i_var in xrange(root_array.GetEntries()):
                        var_name = root_array[i_var].GetTitle()
                        if accept_pat == '*':
                            variables.append( var_name.split('/',1)[0] )
                        else:
                            if re.search( variablePattern, var_name ):
                                variables.append( var_name.split('/',1)[0] )
                    found_tree = True
            if found_tree: break
        else:
            if len(root_files) > 10:
                print_root_files = root_files[:3] + ['...'] + root_files[-4:]
            else:
                print_root_files = root_files
            raise RuntimeError(
                'Not a single root file had a tree called {0}. '
                'List of root files that were checked:\n'
                .format(self.tree_name)
                + '\n'.join(root_files)
                )
        return variables                            


    def read_chain(self, root_files, variables):
        chain = ROOT.TChain(self.tree_name)
        for root_file in root_files:
            chain.Add(root_file)

        Entry = namedtuple('Entry', variables)
        entries = []
        for event in chain:
            entry = {}
            for var_name in variables:
                entry_dict[var_name] = getattr(event, var_name)
                if var_name == self.x_variable:
                    entry_dict['x'] = entry_dict[var_name]
                if var_name == self.y_variable:
                    entry_dict['y'] = entry_dict[var_name]
            entries.append(Entry(**entry_dict))
        return entries


    def read(self):
        root_files = self.collect_root_files()
        if len(root_files) == 0:
            raise RuntimeError(
                'Attemped to retrieve scan for x:{0} y:{1}, '
                'but no .root files were found. Passed list of dirs to look in:\n'
                .format(self.x_variable, self.y_variable)
                + '\n'.join(self.scandirs)
                )

        if self.save_all_variables:
            variables = self.get_list_of_variables_in_tree(self, root_files)
            if not(self.x_variable in variables) or not(self.y_variable in variables):
                raise RuntimeError(
                    'Variables x:{0} and/or y:{1} are not in the found list of variables:\n'
                    .format(self.x_variable, self.y_variable)
                    + '\n'.join(variables)
                    )
        else:
            variables = [ self.x_variable, self.y_variable ]

        self.entries = self.read_chain(root_files, variables)


    def filter_entries(self, inplace=True):
        passed_entries = []
        for entry in self.entries:
            if entry.deltaNLL < self.deltaNLL_threshold:
                if self.filter_negatives:
                    logger.warning(
                        'Dropping entry (deltaNLL<{0}:'.format(self.deltaNLL_threshold)
                        + entry.__repr__()
                        )
                    continue
                else:
                    raise RuntimeError('Not allowed to filter negatives, but found:',entry)
            passed_entries.append(entry)

        if inplace:
            self.entries = passed_entries
        else:
            return passed_entries

    def x(self):
        return [ entry.x for entry in self.entries ]

    def y(self):
        return [ entry.y for entry in self.entries ]

    def deltaNLLs(self):
        return [ entry.deltaNLL for entry in self.entries ]


    def create_uncertainties(self, entries):
        xs = self.x()
        deltaNLLs = self.deltaNLLs()

        min_deltaNLL   = min(deltaNLLs)
        i_min_deltaNLL = rindex(deltaNLLs, min_deltaNLL)
        x_min          = xs[i_min_deltaNLL]

        unc_dict = {
            'min_deltaNLL' : min_deltaNLL,
            'i_min' : i_min_deltaNLL,
            'x_min' : x_min
            }

        # Process left uncertainty
        if i_min < 3:
            well_defined_left_bound = False
        else:
            xs_left = xs[:i_min_deltaNLL+1]
            deltaNLLs_left = deltaNLLs[:i_min_deltaNLL+1]
            if min(deltaNLLs_left) > 0.5 or max(deltaNLLs_left) < 0.5:
                well_defined_left_bound = False
            else:
                left_bound = self.interpolate(xs_left, deltaNLLs_left, 0.5)
                if left_bound is False:
                    well_defined_left_bound = False
                else:
                    well_defined_left_bound = True

        # Process right uncertainty
        if i_min > len(xs)-3:
            well_defined_right_bound = False
        else:
            xs_right = xs[i_min_deltaNLL:]
            deltaNLLs_right = deltaNLLs[i_min_deltaNLL:]
            if min(deltaNLLs_right) > 0.5 or max(deltaNLLs_right) < 0.5:
                well_defined_right_bound = False
            else:
                right_bound = self.interpolate(xs_right, deltaNLLs_right, 0.5)
                if right_bound is False:
                    well_defined_right_bound = False
                else:
                    well_defined_right_bound = True


        if well_defined_left_bound and well_defined_right_bound:
            pass
        elif well_defined_left_bound and not well_defined_right_bound:
            right_bound = x_min + (x_min - left_bound)
        elif well_defined_right_bound and not well_defined_left_bound:
            left_bound  = x_min - (right_bound - x_min)
        else:
            logger.error(
                'Hopeless interpolation case; unable to determine uncertainties for '
                'x = {0}, y = {1}'
                .format(self.x_variable, self.y_variable)
                )
            unc_dict['left_bound'] = -999
            unc_dict['left_error'] = -999
            unc_dict['right_bound'] = 999
            unc_dict['right_error'] = 999



        # Unc = namedtuple('Unc', )


    def interpolate(self, xs, ys, x_value):
        if min(ys) > y_value or max(ys) < y_value:
            return False
        Tg = ROOT.TGraph(len(xs), array('f', xs), array('f', ys))
        y_value = Tg.Eval(x_value)
        # if y_value < min(ys) or y_value > max(ys):
        #     return False
        return y_value








        
def rindex( someList, val ):
    # Regular list.index() finds first instance in list, this function finds the last
    return len(someList) - someList[::-1].index(val) - 1

def FindMinimaAndErrors( POIvals, deltaNLLs, returnContainer=False ):

    minDeltaNLL   = min(deltaNLLs)
    iMin          = rindex( deltaNLLs, minDeltaNLL )
    minimumPOIval = POIvals[iMin]

    # Dict that is returned in case of a total error
    errReturn = {
        'imin'       : iMin,
        'min'        : minimumPOIval,
        'leftError'  : -999,
        'leftBound'  : -999,
        'rightError' : 999,
        'rightBound' : 999,
        'symmError'  : 999,
        'symmBound'  : 999,
        'wellDefinedRightBound' : False,
        'wellDefinedLeftBound'  : False,
        }
    if returnContainer:
        errReturnContainer = TheoryCommands.Container()
        for key, value in errReturn.iteritems():
            setattr( errReturnContainer, key, value )
        errReturn = errReturnContainer

    if iMin > 2:
        # Find left minimum
        POIvalsLeft   = POIvals[:iMin+1]
        deltaNLLsLeft = deltaNLLs[:iMin+1]
        if min(deltaNLLsLeft) > 0.5 or max(deltaNLLsLeft) < 0.5:
            wellDefinedLeftBound = False

        Tg_left = ROOT.TGraph(
            iMin+1,
            array( 'd', deltaNLLsLeft[:iMin+1] ),
            array( 'd', POIvalsLeft[:iMin+1] )
            )
        ROOT.SetOwnership( Tg_left, False )

        leftBound = Tg_left.Eval( 0.5 )
        if leftBound <= POIvals[0]:
            wellDefinedLeftBound = False
        else:
            wellDefinedLeftBound = True
    else:
        wellDefinedLeftBound = False

    if iMin < len(POIvals)-2:
        # Find right minimum
        Tg_right = ROOT.TGraph(
            len(POIvals)-iMin+1,
            array( 'd', deltaNLLs[iMin:] ),
            array( 'd', POIvals[iMin:] )
            )
        ROOT.SetOwnership( Tg_right, False )

        rightBound = Tg_right.Eval( 0.5 )
        if rightBound >= POIvals[-1]:
            wellDefinedRightBound = False
        else:
            wellDefinedRightBound = True
    else:
        wellDefinedRightBound = False


    if wellDefinedLeftBound and wellDefinedRightBound:
        pass
    # Symmetrize if one of the bounds was poorly defined
    elif wellDefinedLeftBound and not wellDefinedRightBound:
        rightBound = minimumPOIval + ( minimumPOIval - leftBound )
    elif wellDefinedRightBound and not wellDefinedLeftBound:
        leftBound  = minimumPOIval - ( rightBound - minimumPOIval )
    else:
        print 'Hopeless interpolation case; unable to determine uncertainties.'
        return errReturn

    leftError = abs(minimumPOIval - leftBound)
    rightError = abs(minimumPOIval - rightBound)


    returnDict = {
        'imin'       : iMin,
        'min'        : minimumPOIval,
        'leftError'  : leftError,
        'leftBound'  : leftBound,
        'rightError' : rightError,
        'rightBound' : rightBound,
        'symmError'  : 0.5*( abs(leftError) + abs(rightError) ),
        'symmBound'  : minimumPOIval + 0.5*( abs(leftError) + abs(rightError) ),
        'wellDefinedRightBound' : wellDefinedRightBound,
        'wellDefinedLeftBound'  : wellDefinedLeftBound,
        }

    if returnContainer:
        container = TheoryCommands.Container()
        for key, value in returnDict.iteritems():
            setattr( container, key, value )
        return container
    else:
        return returnDict


