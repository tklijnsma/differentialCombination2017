import os
import logging
import differentials
import ROOT

class CanvasReader(object):
    """docstring for CanvasReader"""

    outdir = 'fermilabcode/input'

    def __init__(self, root_file):
        super(CanvasReader, self).__init__()
        self.root_file = root_file
        self.canvas_name = 'c'
        self.primitives = []

    def read_c(self):
        with differentials.core.openroot(self.root_file) as root_fp:
            self.c = root_fp.Get(self.canvas_name)
            self.primitives.extend(self.read_primitives(self.c))

    def read_primitives(self, c):
        primitives_TList = c.GetListOfPrimitives()
        n_primitives = primitives_TList.GetSize()
        primitives = []
        for i in xrange(n_primitives):
            primitives.append(primitives_TList.At(i))

        logging.debug('List of found primitives in {0}:'.format(c.GetName()))
        for primitive in primitives:
            name = primitive.GetName().strip()
            if name == '': name = '(unnamed)'
            title = primitive.GetTitle().strip()
            if title == '': title = '(unnamed)'
            logging.debug('{0:60}; name: {1:15}; title: {2:15}'.format(primitive, name, title))
        return primitives

    def read_tgraph(self, tg):
        if not isinstance(tg, ROOT.TGraph):
            logging.error('Passed object {0} is not an (subclass) instance of ROOT.TGraph'.format(tg))
            return

        logging.trace('Found points in TGraph {0}:'.format(tg))
        if isinstance(tg, ROOT.TGraphAsymmErrors):
            points = differentials.plotting.plotting_utils.get_x_y_from_TGraphAsymmErrors(tg, per_point=True)
        else:
            points = differentials.plotting.plotting_utils.get_x_y_from_TGraph(tg, per_point=True)
        self.log_points(points)
        return points

    def log_points(self, points, x_name='x', y_name='y'):
        for i_point, point in enumerate(points):
            if hasattr(point, 'x_err_down'):
                logging.trace(
                    'Point {0:3}: '
                    '{1} = {2:+8.3f}  {3:+8.3f}/{4:+8.3f}; '
                    '{5} = {6:+8.3f}  {7:+8.3f}/{8:+8.3f}'
                    .format(i_point,
                        x_name, point.x, point.x_err_down, point.x_err_up,
                        y_name, point.y, point.y_err_down, point.y_err_up
                        )
                    )
            else:
                logging.trace(
                    'Point {0:3}: '
                    '{1} = {2:+8.3f}; '
                    '{3} = {4:+8.3f}'
                    .format(i_point, x_name, point.x, y_name, point.y)
                    )

    def get_graphs(self, primitives):
        graphs = []
        for primitive in primitives:
            if not isinstance(primitive, ROOT.TGraphAsymmErrors):
                primitive.is_graph = False
                continue
            primitive.is_graph = True
            primitive.points = self.read_tgraph(primitive)
            graphs.append(primitive)
        return graphs

    def dump(self, name):
        out_file = os.path.join(self.outdir, name + '_{0}.py'.format(differentials.core.datestr()))
        contents = [
            '# ' + differentials.core.gittag(),
            'from differentials.core import AttrDict',
            'data = AttrDict(',
            '    binning = {0},'.format(self.binning),
            '    mu      = {0},'.format([ p.y for p in self.mu_data.points ]),
            '    mu_up   = {0},'.format([ p.y_err_up for p in self.mu_data.points ]),
            '    mu_down = {0},'.format([ p.y_err_down for p in self.mu_data.points ]),
            '    xs      = {0},'.format([ p.y for p in self.xs_data.points ]),
            '    xs_up   = {0},'.format([ p.y_err_up for p in self.xs_data.points ]),
            '    xs_down = {0},'.format([ p.y_err_down for p in self.xs_data.points ]),
            '    )'
            ]
        contents = '\n'.join(contents)
        logging.info(
            'Writing following contents to {0}:\n{1}'
            .format(out_file, contents)
            )
        with open(out_file, 'w') as out_fp:
            out_fp.write(contents)


class CanvasReaderHgg(CanvasReader):
    """docstring for CanvasReaderHgg"""
    def __init__(self, root_file):
        super(CanvasReaderHgg, self).__init__(root_file)

    def read(self):
        self.read_c()
        self.top_primitives = self.read_primitives(self.primitives[0])
        self.bottom_primitives = self.read_primitives(self.primitives[1])

        # Only 1 graph in canvas, makes things easier
        self.xs_data = self.get_graphs(self.top_primitives)[0]
        self.mu_data = self.get_graphs(self.bottom_primitives)[0]

        logging.info('Found data; xs:')
        self.log_points(self.xs_data.points)
        logging.info('Found data; mu:')
        self.log_points(self.mu_data.points)

        self.binning = [ 0. ]
        for i_point, point in enumerate(self.xs_data.points):
            self.binning.append(
                self.binning[-1] + 2.*(point.x - self.binning[-1])
                )
        self.binning[-1] = 10000.
        logging.info('Extracted binning: {0}'.format(self.binning))


class CanvasReaderHzz(CanvasReader):
    """docstring for CanvasReaderHzz"""
    def __init__(self, root_file):
        super(CanvasReaderHzz, self).__init__(root_file)

    def read(self):
        self.read_c()
        self.subpad = [ p for p in self.primitives if isinstance(p, ROOT.TPad) ][0]
        self.ratio_primitives = self.read_primitives(self.subpad)

        xs_graphs = self.get_graphs(self.primitives)
        self.xs_data = self.get_data(xs_graphs)

        mu_graphs = self.get_graphs(self.ratio_primitives)
        self.mu_data = self.get_data(mu_graphs)

        logging.info('Found data; xs:')
        self.log_points(self.xs_data.points)
        logging.info('Found data; mu:')
        self.log_points(self.mu_data.points)

        self.binning = [ 0. ]
        for i_point, point in enumerate(self.xs_data.points):
            self.binning.append(
                self.binning[-1] + 2.*(point.x - self.binning[-1])
                )
        self.binning[-1] = 10000.
        logging.info('Extracted binning: {0}'.format(self.binning))


    def get_data(self, graphs):
        for graph in graphs:
            marker_style = graph.GetMarkerStyle()
            line_color   = graph.GetLineColor()
            if marker_style >= 8 and marker_style <= 20 and line_color == 1:
                break
        else:
            raise RuntimeError('Could not find data; no graph fits the requirements')
        return graph




