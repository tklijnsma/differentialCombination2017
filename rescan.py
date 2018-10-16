#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, argparse, differentials, logging, re, shutil, tempfile
import os.path as osp
from time import strftime


########################################
# Main
########################################


class JobSplitter(object):
    """docstring for JobSplitter"""
    def __init__(self, accountant):
        super(JobSplitter, self).__init__()
        self.accountant = accountant
        self.sh_files = []
        self._example_sh_text = False

    def make_new_scandir(self):
        self.old_scandir = self.accountant.scandir
        if self.old_scandir.endswith('/'): self.old_scandir = self.old_scandir[:-1]
        self.scandir = differentials.core.make_unique_directory(
            self.old_scandir + '_rescan_{0}'.format(differentials.core.datestr())
            )
        if not(differentials.core.is_testmode()):
            os.makedirs(self.scandir)

    def copy_postfit_and_fastscan(self):
        src = osp.join(self.old_scandir, 'postfit_and_fastscan')
        dst = osp.join(self.scandir, 'postfit_and_fastscan')
        if osp.isdir(src):
            if differentials.core.is_testmode():
                logging.info('Would now copy {0} to {1}'.format(src, dst))
            else:
                logging.info('Copying {0} to {1}'.format(src, dst))
                shutil.copytree(src, dst)


    def make_new_sh_files(self, job):
        with open(job.sh_file, 'r') as fp:
            base_sh_text = fp.read()

        points = job.points
        for i_chunk, chunk in enumerate(self.split_points(points)):
            sh_text = self.get_new_shfile_text(base_sh_text, i_chunk, chunk)
            outfile = osp.join(
                self.scandir,
                osp.basename(job.sh_file).replace('.sh', '_chunk{0}.sh'.format(i_chunk))
                )
            self.dump_sh_text_to_file(outfile, sh_text)
            if self._example_sh_text is False: self._example_sh_text = sh_text

    def split_points(self, points):
        i_half = int(0.5*len(points))
        return [ points[:i_half], points[i_half:] ]

    def get_new_shfile_text(self, base_sh_text, i_chunk, points):
        sh_text = base_sh_text
        sh_text = sh_text.replace(self.old_scandir, self.scandir)
        sh_text = re.sub(r'--doPoints ([\d,]+)', '', sh_text)
        sh_text = re.sub(r'--firstPoint (\d+)', '', sh_text)
        sh_text = re.sub(r'--lastPoint (\d+)', '', sh_text)
        sh_text = sh_text.strip()
        sh_text += ' --doPoints {0}\n'.format(','.join([str(p) for p in points]))
        sh_text = re.sub(r'\s-n\s+(\w+)', r' -n \1_chunk{0}'.format(i_chunk), sh_text)
        return sh_text

    def dump_sh_text_to_file(self, outfile, sh_text):
        if differentials.core.is_testmode():
            logging.info(
                'Would now dump new sh_text to {0}:\n{1}'
                .format(outfile, sh_text)
                )
        else:
            logging.warning('Writing {0}'.format(outfile))
            with open(outfile, 'w') as f:
                f.write(sh_text)
        self.sh_files.append(outfile)

    def submit(self):
        logging.warning('Submitting jobs')
        full_submission_output = ''
        with differentials.core.enterdirectory(self.scandir):
            for sh_file in self.sh_files:
                submit_cmd = 'qsub -q {0} {1}'.format(
                    self.accountant.queue,
                    osp.basename(sh_file)
                    )
                if differentials.core.is_testmode():
                    logging.info(submit_cmd)
                    output = 'TESTMODEOUTPUT'
                else:
                    output = differentials.core.execute(submit_cmd, py_capture_output=True)
                full_submission_output += output
        self.register_jobids_in_jobmanager(full_submission_output)

    def register_jobids_in_jobmanager(self, submission_output):
        if differentials.core.is_testmode():
            logging.info('Not writing any jobmanager files')
            return

        jobids = re.findall(r'Your job (\d+)', submission_output)
        if len(jobids) == 0:
            logging.error('No jobids were found in the passed submission output; nothing to register for the jobmanager')

        header = [
            osp.basename(self.scandir),
            'scandir: {0}'.format(self.scandir),
            'registration time: {0}'.format(strftime('%y-%m-%d %H:%M:%S')),
            'example cmd:\n\n{0}\n'.format(self._example_sh_text)
            ]
        contents = '\n'.join(header) + '\n' + '\n'.join(jobids) + '\n'

        _, jobman_file = tempfile.mkstemp(
            prefix = 'tklijnsm_queuegroup_',
            suffix = '.jobman',
            dir    = '/tmp'
            )
        logging.warning('Dumping following jobmanager contents to {0}:\n'.format(jobman_file))
        logging.warning(contents + '\n')
        with open(jobman_file, 'w') as jobman_fp:
            jobman_fp.write(contents)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( 'scandirs', metavar='N', type=str, nargs='+', help='list of strings' )
    parser.add_argument( '--test', action='store_true', help='boolean')
    parser.add_argument( '--dry', action='store_true', help='boolean')
    # parser.add_argument( '--allq', action='store_true', help='boolean')
    # parser.add_argument( '--shortq', action='store_true', help='boolean')
    # parser.add_argument( '--longq', action='store_true', help='boolean')
    args = parser.parse_args()

    logging.getLogger().setLevel(logging.WARNING)
    
    scandir = args.scandirs[0]
    accountant = differentials.scan_accounting.ScanAccountant(scandir)
    jobs = accountant.get_failed_jobs()

    if args.test:
        differentials.core.testmode()
        logging.getLogger().setLevel(logging.INFO)

    splitter = JobSplitter(accountant)
    splitter.make_new_scandir()
    splitter.copy_postfit_and_fastscan()

    for job in jobs:
        splitter.make_new_sh_files(job)

    splitter.submit()










########################################
# End of Main
########################################
if __name__ == "__main__":
    main()