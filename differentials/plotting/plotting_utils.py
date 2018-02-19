import differentials.core
import ROOT

ROOTCOUNTER = 1000
def get_unique_rootname():
    global ROOTCOUNTER
    name = 'root{0}_{1}'.format( ROOTCOUNTER, differentials.core.__uniqueid__().next() )
    ROOTCOUNTER += 1
    return name

def get_plot_base(
        x_min = 0, x_max = 1,
        y_min = 0, y_max = 1,
        x_title = 'x', y_title = 'y',
        set_title_sizes = True,
        ):
    base = ROOT.TH1F()
    ROOT.SetOwnership( base, False )
    base.SetName( get_unique_rootname() )
    base.GetXaxis().SetLimits( x_min, x_max )
    base.SetMinimum( y_min )
    base.SetMaximum( y_max )
    base.SetMarkerColor(0)
    base.GetXaxis().SetTitle( x_title )
    base.GetYaxis().SetTitle( y_title )
    if set_title_sizes:
        base.GetXaxis().SetTitleSize( 0.06 )
        base.GetYaxis().SetTitleSize( 0.06 )
    return base

def set_color_palette(option=None):
    # n_stops = 3
    # stops  = [ 0.0, 0.5, 1.0 ]
    # reds   = [ 0.0, 1.0, 1.0 ]
    # blues  = [ 1.0, 1.0, 0.0 ]
    # greens = [ 0.0, 1.0, 0.0 ]

    # n_stops = 2
    # stops  = [ 0.0, 1.0 ]
    # reds   = [ 55./255.,  1.0 ]
    # greens = [ 138./255., 1.0 ]
    # blues  = [ 221./255., 1.0 ]

    if option == 'twocolor':
        n_stops = 3
        stops  = [ 0.0, 0.3, 1.0 ]
        reds   = [ 55./255.,  1.0, 1.0 ]
        greens = [ 138./255., 1.0, 26./255. ]
        blues  = [ 221./255., 1.0, 26./255. ]
    else:
        n_stops = 3
        stops  = [ 0.0, 0.3, 1.0 ]
        reds   = [ 55./255.,  166./255., 1.0 ]
        greens = [ 138./255., 203./255., 1.0 ]
        blues  = [ 221./255., 238./255., 1.0 ]

    ROOT.TColor.CreateGradientColorTable(
        n_stops,
        array('d', stops ),
        array('d', reds ),
        array('d', greens ),
        array('d', blues ),
        255 )
