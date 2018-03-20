import core
import logging, re

class AbstractBin(object):
    """docstring for AbstractBin"""
    def __init__(self, name):
        super(AbstractBin, self).__init__()
        self.name = name

    def __str__(self):
        ret = '{0} {1}: '.format(self.__class__.__name__, self.name)
        if hasattr(self, 'left'):
            ret += 'left={0} '.format(self.left)
        if hasattr(self, 'right'):
            ret += 'right={0} '.format(self.right)
        if self.is_underflow:
            ret += 'is_underflow=True '
        if self.is_overflow:
            ret += 'is_overflow=True '
        if self.is_interval:
            ret += 'is_interval=True '
        return ret

class Process(AbstractBin):
    """docstring for Process"""

    def __init__(self, name):
        super(Process, self).__init__(name)
        self.name = name
        self.prefix = None
        self.is_interval = None
        self.is_overflow = False
        self.is_underflow = False
        self.from_str(name)
        
    def from_str(self, text):
        regular_match = re.search(r'(.*)_([\dpm\.\-]+)_([\dpm\.\-]+)', text)
        overflow_match = re.search(r'(.*)_(GE|GT)([\dpm\.\-]+)', text)
        underflow_match = re.search(r'(.*)_(LE|LT)([\dpm\.\-]+)', text)
        single_match = re.search(r'(.*)_([\dpm\.\-]+)', text)
        if regular_match:
            self.prefix = regular_match.group(1)
            self.left_str = regular_match.group(2)
            self.left = str_to_float(self.left_str)
            self.right_str = regular_match.group(3)
            self.right = str_to_float(self.right_str)
            self.is_interval = True
        elif overflow_match:
            self.prefix = overflow_match.group(1)
            self.overflow_str = overflow_match.group(2)
            self.left_str = overflow_match.group(3)
            self.left = str_to_float(self.left_str)
            self.is_overflow = True
        elif underflow_match:
            self.prefix = underflow_match.group(1)
            self.underflow_str = underflow_match.group(2)
            self.right_str = underflow_match.group(3)
            self.right = str_to_float(self.right_str)
            self.is_underflow = True
        elif single_match:
            self.prefix = single_match.group(1)
            logging.debug('single_match for {0}; matched text is {1}'.format(text, single_match.group(2)))
            self.left_str = single_match.group(2)
            self.left = str_to_float(self.left_str)
            self.is_interval = False
        else:
            raise RuntimeError('Error processing {0}'.format(text))

class POI(AbstractBin):
    """docstring for POI"""
    def __init__(self, *args):
        super(POI, self).__init__(args[0])
        
class YieldParameter(AbstractBin):
    """docstring for YieldParameter"""

    prefix = ''

    def __init__(self, name, left=None, right=None):
        super(YieldParameter, self).__init__(name)
        self.name = self.prefix + name
        self.left = left
        self.right = right
        self.is_interval = None
        self.is_overflow = False
        self.is_underflow = False

        self.scan_left = -1.0
        self.scan_right = 4.0


class ProcessInterpreter(object):
    """docstring for ProcessInterpreter"""
    def __init__(self, process_strs=None, binning=None):
        super(ProcessInterpreter, self).__init__()
        self.processes = []
        self.binning = binning
        self.is_interval = None
        self.last_bin_is_overflow = None
        self.first_bin_is_underflow = None
        self.number_formats = {}
        self.yield_parameters = []

        if not(process_strs is None):
            self.analyze_processes(process_strs)

    def analyze_processes(self, process_strs):
        for process_str in process_strs:
            if 'OutsideAcceptance' in process_str: continue
            self.processes.append(Process(process_str))

        if len(self.processes) == 0:
            raise RuntimeError('len(self.processes) == 0')
        self.prefix = self.processes[0].prefix
        for process in self.processes:
            if not process.prefix == self.prefix:
                raise RuntimeError('Unexpected prefix {0}'.format(process.prefix))

        self.sort_processes()
        logging.debug('Registered the following processes:')
        for process in self.processes:
            logging.debug(process)

        self.non_overflow_processes = [ process for process in self.processes if not(process.is_overflow) and not(process.is_underflow) ]
        self.set_last_bin_is_overflow()
        self.set_first_bin_is_underflow()
        self.set_is_interval()
        self.determine_number_formats()

    def determine_number_formats(self):
        for process in self.processes:
            if hasattr(process, 'left_str'): self.add_number_format(process.left, process.left_str)
            if hasattr(process, 'right_str'): self.add_number_format(process.right, process.right_str)

    def add_number_format(self, val, text):
        if val in self.number_formats:
            if not text == self.number_formats[val]:
                raise RuntimeError(
                    'Inconsistent strings for value {0}: found {1} and {2}'
                    .format(val, text, self.number_formats[val])
                    )
        else:
            self.number_formats[val] = text

    def get_number_format(self, val):
        if not val in self.number_formats:
            return float_to_str(val)
        return self.number_formats[val]

    def sort_processes(self):
        self.processes.sort(key=range_sorter)

    def set_is_interval(self):
        if len(self.non_overflow_processes) == 0:
            self.is_interval = True
            return
        statuses = [p.is_interval for p in self.non_overflow_processes]
        if len(set(statuses)) != 1:
            raise RuntimeError(
                'Inconsistent is_interval:\n{0}'.format('\n'.join([str(p) for p in self.non_overflow_processes]))
                )
        self.is_interval = all(statuses)

    def set_last_bin_is_overflow(self):
        if self.processes[-1].is_overflow: self.last_bin_is_overflow = True

    def set_first_bin_is_underflow(self):
        if self.processes[0].is_underflow: self.first_bin_is_underflow = True

    def make_yield_parameters(self, add_underflow=False, add_overflow=False):
        logging.debug('Making yield yield_parameters')
        if self.binning is None:
            logging.debug('No binning specified; yield parameters will be the process name prefixed with \'r_\'')
            YieldParameter.prefix = 'r_'
            for process in self.processes:
                y = YieldParameter(process.name)
                self.yield_parameters.append(y)
                process.yield_parameter = y
                logging.info('Created {0}'.format(y))
            return

        YieldParameter.prefix = 'r_' + self.prefix + '_'
        logging.debug('Prefix determined from processes is {0}'.format(YieldParameter.prefix))

        self.yield_parameters = []

        if add_underflow:
            underflow_str = self.processes[0].underflow_str
            right = self.binning[0]
            if not(self.is_interval) and underflow_str == 'LT':
                underflow_str = 'LE'
                right -= 1.
            y = YieldParameter(underflow_str + self.get_number_format(right), left=-10e9, right=right)
            y.is_underflow = True
            self.yield_parameters.append(y)
            logging.info('Created {0}'.format(y))

        if self.is_interval:
            for left, right in zip(self.binning[:-1], self.binning[1:]):
                y = YieldParameter(
                    self.get_number_format(left) + '_' + self.get_number_format(right),
                    left = left, right = right
                    )
                self.yield_parameters.append(y)
                logging.info('Created {0}'.format(y))
        else:
            for val in self.binning:
                y = YieldParameter(self.get_number_format(val), left=val)
                self.yield_parameters.append(y)
                logging.info('Created {0}'.format(y))

        if add_overflow:

            if self.is_interval:
                left = self.binning[-1]
                if self.last_bin_is_overflow:
                    overflow_str = self.processes[-1].overflow_str
                else:
                    overflow_str = 'GT'
            else:
                left = self.binning[-1] + 1
                overflow_str = 'GE'

            y = YieldParameter(overflow_str + self.get_number_format(left), left=left, right=10e9)
            y.is_overflow = True
            self.yield_parameters.append(y)
            logging.debug('Created {0}'.format(y))


    def link_processes_to_yield_parameters(self):
        if self.binning is None:
            logging.info('self.binning is None; no need to link')
            return

        if self.first_bin_is_underflow:
            process = self.processes[0]
            for yield_parameter in self.yield_parameters:
                if process.right <= yield_parameter.right:
                    break
            else:
                raise RuntimeError('Error processesing underflow process {0}'.format(process.name))
            if not(yield_parameter is self.yield_parameters[0]):
                logging.warning(
                    'The underflow process {0} is scaled by {1}, not by the underflow yield_parameter {2}!'
                    .format(process.name, yield_parameter.name, self.yield_parameters[0].name)
                    )
            process.yield_parameter = yield_parameter


        for process in self.non_overflow_processes:
            for yield_parameter in self.yield_parameters:
                if self.is_interval:
                    if process.left >= yield_parameter.left and process.right <= yield_parameter.right:
                        break
                else:
                    if process.left == yield_parameter.left:
                        break
            else:
                raise RuntimeError('Could not find a yield_parameter for process {0}'.format(process))
            process.yield_parameter = yield_parameter


        if self.last_bin_is_overflow:
            process = self.processes[-1]
            for yield_parameter in self.yield_parameters[::-1]:
                if process.left >= yield_parameter.left:
                    break
            else:
                raise RuntimeError('Error processesing overflow process {0}'.format(process.name))
            if not(yield_parameter is self.yield_parameters[-1]):
                logging.warning(
                    'The overflow process {0} is scaled by {1}, not by the overflow yield_parameter {2}!'
                    .format(process.name, yield_parameter.name, self.yield_parameters[-1].name)
                    )
            process.yield_parameter = yield_parameter

        for process in self.processes:
            logging.info(
                'Process {0:15} will be scaled by yield parameter {1}'
                .format(process.name, process.yield_parameter.name)
                )

    def make_maps(self):
        maps = []
        for process in self.processes:
            y = process.yield_parameter
            m = (
                '--PO \'map=.*/{process}:{yield_parameter}[1.0,{left},{right}]\''
                .format(process=process.name, yield_parameter=y.name, left=y.scan_left, right=y.scan_right)
                )
            maps.append(m)
        return maps


def float_to_str(number, nDecimals=None):
    number = float(number)
    if not(nDecimals is None):
        string = '{:.{nDecimals}f}'.format(number, nDecimals=nDecimals).replace('-','m').replace('.','p')
        return string
    if number.is_integer():
        number = int(number)
    string = str(number).replace('-','m').replace('.','p')
    return string

def str_to_float(string):
    string = str(string)
    number = float(
        string
        .replace('p','.')
        .replace('m','-')
        .replace('GT','')
        .replace('GE','')
        .replace('LT','')
        .replace('LE','')
        )
    return number

def range_sorter(process):
    if process.is_overflow:
        return 800000
    elif process.is_underflow:
        return -800000
    else:
        return process.left
