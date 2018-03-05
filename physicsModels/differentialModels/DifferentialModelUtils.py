
def str_to_number(number):
    return float(
        number
        .replace('p','.')
        .replace('m','-')
        .replace('GT','')
        .replace('GE','')
        .replace('LT','')
        .replace('LE','')
        )

def left_bound_of_process(process):
    left_str = process.split('_')[2]
    return str_to_number(left_str)

def right_bound_of_process(process):
    components = process.split('_')
    if len(components) < 4:
        return None
    else:
        return str_to_number(components[3])

def boundaries_of_process(process):
    return [ left_bound_of_process(process), right_bound_of_process(process) ]

def is_interval(process):
    return len(process.split('_')) >= 4

def is_overflow(last_process):
    return 'GT' in last_process or 'GE' in last_process


def ensure_format(bin_boundary, processes):
    """
    Returns a string that is formatted as in the processes
    bin_boundary should be a number, processes a list of strings
    """

    # Obtain 1 example of processes that contains this boundary
    for process in processes:
        left, right = boundaries_of_process(process)
        if left == bin_boundary or right == bin_boundary:
            break
    else:
        raise ValueError(
            'The bin boundary {0} could nowhere be found in the passed'
            'list of processes {1}'
            .format(bin_boundary, processes)
            )

    if left == bin_boundary:
        bin_boundary_str = process.split('_')[2]
    else:
        bin_boundary_str = process.split('_')[3]

    return bin_boundary_str



def make_yield_parameters(processes, binning=None, make_underflow=False):
    if binning is None:
        ys = []
        for process in processes:
            ys.append('r_' + process)
        return ys
    else:
        pre = 'r_' + '_'.join(processes[0].split('_')[:2])
        do_intervals = is_interval(processes[0])

        _ensure_format = lambda number: ensure_format(number, processes)
        ys = []

        # Only necessary for ptjet
        if make_underflow:
            ys.append(pre + '_LT' + _ensure_format(binning[0]))

        for i in xrange(len(binning)-1):
            if do_intervals:
                ys.append(pre + '_' + _ensure_format(binning[i]) + '_' + _ensure_format(binning[i+1]))
            else:
                ys.append(pre + '_' + _ensure_format(binning[i]))

        # Process last bin separately
        if is_overflow(processes[-1]):
            if do_intervals:
                ys.append(pre + '_GT' + _ensure_format(binning[-1]))
            else:
                ys.append(pre + '_GE' + _ensure_format(binning[-1]))

        return ys


def find_i_bin(process, binning):
    left, right = boundaries_of_process(process)
    if right is None:
        # Start looking from the right
        reverse_binning = binning[::-1]
        for i in xrange(len(binning)):
            if left >= reverse_binning[i]:
                return len(binning)-i-1
        else:
            raise ValueError(
                'Could not find corresponding bin for {0} in {1}'
                .format(process, binning)
                )
    else:
        if left >= binning[-1]:
            # This will have to be the overflow; may raise IndexError later if no overflow was made
            return len(binning)-1
        for i in xrange(len(binning)-1):
            if left >= binning[i] and right <= binning[i+1]:
                return i
        else:
            # raise ValueError(
            #     'Could not find corresponding bin for {0} in {1}'
            #     .format(process, binning)
            #     )
            print 'Could not find corresponding bin for {0} in {1}; Assuming it is meant to be unscaled'.format(process, binning)
            return -1



