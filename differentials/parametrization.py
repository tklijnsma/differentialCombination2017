import core
import ROOT
import logging

class WSParametrization(object):
    """docstring for WSParametrization"""
    def __init__(self, ws_file):
        logging.debug('Initializing parametrization with file {0}'.format(ws_file))
        self.ws_file = ws_file
        with core.openroot(ws_file) as ws_fp:
            self.w = ws_fp.Get('w')
        self.old_style = False
        self.yield_parameters = []
        self.smxs = []

    def get_yield_parameters(self):
        yp_list = self.get_yield_parameter_arglist()
        for i in xrange(yp_list.getSize()):
            yp = yp_list[i]
            yp_name = yp.GetName()
            logging.debug('Found yield parameter: '.format(yp))

            if not(self.old_style) and yp_name.startswith('r_'):
                new_yp_name = yp_name.replace('r_', 'parametrization_')
                logging.info('Taking {0} instead of {1}'.format(new_yp_name, yp_name))
                new_yp = self.w.function(new_yp_name)
                if new_yp == None:
                    logging.error('Variable {0} does not exist in ws; taking {1} instead'.format(new_yp_name, yp_name))
                else:
                    yp = new_yp
            self.yield_parameters.append(yp)

    def get_yield_parameter_arglist(self):
        if self.set_exists('all_ggH_yieldParameters'):
            logging.debug('Found set called all_ggH_yieldParameters')
            argset = self.w.set('all_ggH_yieldParameters')
        elif self.set_exists('yieldParameters'):
            logging.debug('Found set called yieldParameters')
            self.old_style = True
            argset = self.w.set('yieldParameters')
        else:
            raise RuntimeError(
                'Sets \'{0}\' and \'{1}\' do not exist in {2}'
                .format('all_ggH_yieldParameters', 'yieldParameters', self.ws_file)
                )
            # raise RuntimeError(
            #     'Set \'{0}\' does not exist in {1}'
            #     .format(set_name, self.ws_file)
            #     )
        ROOT.SetOwnership(argset, False)
        arglist = ROOT.RooArgList(argset)
        ROOT.SetOwnership(arglist, False)
        return arglist

    def set_exists(self, set_name):
        logging.debug('Checking if set {0} exists'.format(set_name))
        s = self.w.set(set_name)
        logging.debug('  Raw repr of w.set: {0}'.format(s))
        if s == None:
            return False
        else:
            return True

    def set(self, name, value):
        logging.debug('Setting {0} to {1}'.format(name, value))
        roovar = self.w.var(name)
        if roovar == None:
            raise RuntimeError(
                'Variable \'{0}\' does not exist in {1}'
                .format(name, self.ws_file)
                )
        roovar.setVal(value)

    def set_smxs(self, smxs):
        self.smxs = smxs

    def get_mus_exp(self, **kwargs):
        if len(self.yield_parameters) == 0:
            self.get_yield_parameters()

        for name, value in kwargs.iteritems():
            self.set(name, value)
        mus = [ yp.getVal() for yp in self.yield_parameters ]
        return mus

    def get_xs_exp(self, **kwargs):
        if len(self.smxs)==0:
            raise RuntimeError('Need to set smxs to a list of xs first')
        mus = self.get_mus_exp(**kwargs)
        if not len(mus) == len(self.smxs):
            raise ValueError(
                'Found {0} yield parameters, and {1} smxs'
                .format(len(mus), len(self.smxs))
                )
        return [ mu * xs for mu, xs in zip(mus, self.smxs) ]


