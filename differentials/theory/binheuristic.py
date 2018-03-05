#____________________________________________________________________
# Theory files contain a somewhat non-consistent binning, turn it into a well-defined binning
def BinningHeuristic(
    binCenters,
    manualSwitchAt50 = True,
    manualSwitchAt5  = False,
    ):

    nBins = len(binCenters)
    binBoundaries = []

    # Assume first bin center is in the middle; use second bin to determine first boundary
    binBoundaries.append( binCenters[0] - 0.5*(binCenters[1]-binCenters[0]) )

    # Bin boundaries are defined by being the middle of bin centers
    for iBin in xrange(0,nBins-1):

        # Manually overwrite for the cross-over at pt = 50 (which is irregular)
        if manualSwitchAt50 and binCenters[iBin] == 48.5:
            binBoundaries.append( 51 )

        # For the quark induced histograms there is an irregularity at pt = 5.0
        elif manualSwitchAt5 and ( binCenters[iBin] == 3.0 or binCenters[iBin] == 2.75 ):
            binBoundaries.append( 5.0 )

        else:
            binBoundaries.append(
            binCenters[iBin] + 0.5*(binCenters[iBin+1]-binCenters[iBin])
            )

    binBoundaries.sort()

    # Last point gets the same width as previous bin
    binBoundaries.append( binCenters[-1] + 0.5*(binCenters[-1]-binCenters[-2]) )

    # Bin centers may have changed because of this definition
    newBinCenters = [ 0.5*(binBoundaries[i+1]+binBoundaries[i]) for i in xrange(nBins-1) ]
    binWidths = [ binBoundaries[i+1]-binBoundaries[i] for i in xrange(nBins-1) ]

    return newBinCenters, binBoundaries, binWidths



class BinHeuristic(object):
    """docstring for BinHeuristic"""
    def __init__(self):
        super(BinHeuristic, self).__init__()

    def get_bin_boundaries(self, bin_centers):
        n_bins = len(bin_centers)

        leftmost = bin_centers[0] - 0.5*(bin_centers[1]-bin_centers[0])
        rightmost = bin_centers[-1] + 0.5*(bin_centers[-1]-bin_centers[-2])

        bin_boundaries = []
        for i_bin in xrange(n_bins-1):
            boundary = bin_centers[i_bin] + 0.5*(bin_centers[i_bin+1]-bin_centers[i_bin])
            bin_boundaries.append(boundary)

        return [leftmost] + bin_boundaries + [rightmost]