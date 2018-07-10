import differentials
import copy

class OneDimScanFilter(object):
    """docstring for OneDimScanFilter"""
    def __init__(self, xs, ys):
        super(OneDimScanFilter, self).__init__()
        self.xs = xs
        self.ys = ys
        self.filter_duplicates_and_sort()
        self.set_minimum()
        self.deleted_points = []

        self.multiply_by_two_to_graph = True


    def set_minimum(self):
        self.y_min = min(self.ys)
        self.i_min = self.ys.index(self.y_min)
        self.x_min = self.xs[self.i_min]

    def left_branch(self):
        return self.xs[:self.i_min+1], self.ys[:self.i_min+1]

    def right_branch(self):
        return self.xs[self.i_min:], self.ys[self.i_min:]


    def filter_duplicates_and_sort(self):
        d = {}
        for x, y in zip(self.xs, self.ys):
            d[x] = y
        self.xs = d.keys()
        self.xs.sort()
        self.ys = [ d[x] for x in self.xs ]


    def filter_clear_nonsense(self):
        xs_left, ys_left = self.left_branch()
        xs_left, ys_left, xs_deleted, ys_deleted = self.filter_branch(xs_left, ys_left, is_left=True)
        for x, y in zip(xs_deleted, ys_deleted): self.deleted_points.append((x, y))

        xs_right, ys_right = self.right_branch()
        xs_right, ys_right, xs_deleted, ys_deleted = self.filter_branch(xs_right, ys_right)
        for x, y in zip(xs_deleted, ys_deleted): self.deleted_points.append((x, y))

        self.xs = xs_left + xs_right
        self.ys = ys_left + ys_right
        self.set_minimum()

    def raise_by_minimum(self):
        self.set_minimum()
        self.ys = [ y - self.y_min for y in self.ys ]
        self.set_minimum()


    def filter_branch(self, xs, ys, is_left=False):
        if is_left:
            xs = xs[::-1]
            ys = ys[::-1]
        else:
            xs = xs[:]
            ys = ys[:]
        n = len(xs)
        xs_deleted = []
        ys_deleted = []
        # Filter n times (cannot be more often), but break when nothing to delete
        for i in xrange(n):
            to_delete = self.filter_branch_once(xs, ys)
            if len(to_delete) == 0: break
            xs, ys, xs_del, ys_del = self.delete_points(to_delete, xs, ys)
            xs_deleted += xs_del
            ys_deleted += ys_del
        if is_left:
            xs = xs[::-1]
            ys = ys[::-1]
            xs_deleted = xs_deleted[::-1]
            ys_deleted = ys_deleted[::-1]
        return xs, ys, xs_deleted, ys_deleted

    def filter_branch_once(self, xs, ys):
        """Built for the right branch; invert ys to use for left branch"""
        n = len(xs)
        to_delete = []
        for i in xrange(n-1):
            if ys[i+1] < ys[i]:
                # y[i] is non-sensical if the next point is lower
                # For the last point, nothing can be said
                to_delete.append(i)
        return to_delete

    def delete_points(self, indices, xs, ys):
        xs_kept = []
        ys_kept = []
        xs_deleted = []
        ys_deleted = []
        for i in xrange(len(xs)):
            if i in indices:
                xs_deleted.append(xs[i])
                ys_deleted.append(ys[i])
            else:
                xs_kept.append(xs[i])
                ys_kept.append(ys[i])
        return xs_kept, ys_kept, xs_deleted, ys_deleted


    def to_graph(self, style=None):
        if self.multiply_by_two_to_graph:
            ys = [ 2. * y for y in self.ys ]
        else:
            ys = self.ys
        graph = differentials.plotting.pywrappers.Graph(
            'auto', 'title', self.xs, ys
            )
        if not style is None: style.apply(graph)
        return graph

    def show(self):
        print '    {0:<9s}    {1:>9s}   {2:>9s}'.format('i', 'x', 'y')
        for i, (x, y) in enumerate(zip(self.xs, self.ys)):
            print '    {0:<9}    {1:9.4f}   {2:9.4f}'.format(i, x, y)
        if len(self.deleted_points) > 0:
            print '\n    Deleted points: {0}'.format(self.deleted_points)







