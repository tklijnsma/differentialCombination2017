import glob
import differentials.core as core

def read_theory_file(theory_file):
    with open(theory_file, 'r') as theory_fp:
        text = theory_fp.read()

    d = {}
    for line in text.split('\n'):
        line = line.strip()
        if len(line)==0 or line.startswith('#'): continue
        key, value = line.split('=',1)
        if key in ['file']:
            pass
        elif ',' in value:
            value = [ core.str_to_float(v) for v in value.split(',') ]
        else:
            value = core.str_to_float(value)
        d[key] = value

    d['theory_file'] = theory_file
    return core.AttrDict(**d)


class FileFinder(object):
    def __init__(self, **kwargs):
        super(FileFinder, self).__init__()
        option_keys = ['directory', 'expect_only_one']
        for option_key in option_keys:
            if option_key in kwargs:
                setattr(self, option_key, kwargs.pop(option_key))
        self.par_dict = kwargs

    def get(self):
        theory_files = glob.glob(self.directory + '/*.txt')
        theories = [ read_theory_file(f) for f in theory_files ]

        accepted_theories = []
        for theory in theories:
            for key, value in self.par_dict.iteritems():
                if getattr(theory, key, None) != value:
                    break
            else:
                accepted_theories.append(theory)

        return accepted_theories
