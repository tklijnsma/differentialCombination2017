
import differentials
import os, re, logging
import os.path as osp
from glob import glob

import math
from math import sqrt
from subprocess import CalledProcessError




class Job(object):
    """docstring for Job"""
    def __init__(self, sh_file):
        super(Job, self).__init__()
        self.sh_file = sh_file
        self._n_points_set = False
        self.parameter_ranges = {}
        self._via_dopoints = False
    
    def short_filename(self, filename):
        return fixed_length_string(osp.basename(filename), 30)

    def get_o_file(self):
        output_files = glob(self.sh_file + '.o*')
        if len(output_files) == 0:
            logging.error(
                'Job {0}: No .o files'
                .format(fixed_length_basename(self.sh_file))
                )
            return False
        if len(output_files) > 1:
            output_files.sort(key = lambda o: int(o.rsplit('.o', 1)[1]))
            self.o_file = output_files[-1]
            logging.warning(
                'Job {0}: Found {1}, picking: {2}'
                .format(
                    fixed_length_basename(self.sh_file),
                    len(output_files),
                    fixed_length_basename(self.o_file)
                    )
                )
        else:
            self.o_file = output_files[0]


    def get_combine_settings(self):
        self.get_n_points()

    def get_n_points(self):
        if self._n_points_set: return
        match = re.search(r'--doPoints ([\d,]+)', self.combine_command)
        if match:
            self.points = [ int(i) for i in match.group(1).split(',') ]
            self.n_points = len(self.points)
            self._n_points_set = True
            self._via_dopoints = True
            return
        else:
            match_first = re.search(r'--firstPoint (\d+)', self.combine_command)
            match_last = re.search(r'--lastPoint (\d+)', self.combine_command)
            if match_first and match_last:
                self.first_point = int(match_first.group(1))
                self.last_point = int(match_last.group(1))
                self.n_points = self.last_point - self.first_point + 1
                self.points = range(self.first_point, self.last_point+1)
                self._n_points_set = True
                self._via_dopoints = False
                return
        logging.error('Could not determine points for {0}'.format(self.sh_file))

    def get_n_points_total(self):
        match = re.search(r'--points (\d+)', self.combine_command)
        if match:
            return int(match.group(1))
        else:
            logging.error('Could not determine points from {0}'.format(self.sh_file))
            return

    def get_parameter_ranges(self):
        match = re.search(r'--setPhysicsModelParameterRanges ([^ ]+)', self.combine_command)
        parameter_ranges = {}
        if match:
            for rangestr in match.group(1).split(':'):
                name, ranges = rangestr.split('=',1)
                left, right = [ float(i) for i in ranges.split(',') ]
                parameter_ranges[name] = (left, right)
        return parameter_ranges

    def cpu_per_point(self):
        self.get_n_points()
        return float(self.cpu_time.number_of_hours()) / self.n_points

    def get_jobid(self):
        if hasattr(self, 'jobid'): return self.jobid
        if not hasattr(self, 'o_file'): self.get_o_file()
        self.jobid = int(self.o_file.rsplit('.o', 1)[1])
        return self.jobid

    def get_qacct(self):
        if not hasattr(self, 'jobid'): self.get_jobid()
        cmd = 'qacct -j {0}'.format(self.jobid)
        try:
            self.qacct = differentials.core.execute(cmd, py_capture_output=True)
        except CalledProcessError:
            logging.error('Could not qacct job {0}; too old, or not yet finished?'.format(self.jobid))
            self.qacct = ''


    def is_failed(self):
        if hasattr(self, '_is_failed'):
            return self._is_failed
        else:
            return self.is_failed_using_qacct()

    def is_failed_using_qacct(self):
        if hasattr(self, '_is_failed'): return self._is_failed
        if not hasattr(self, 'qacct'): self.get_qacct()
        match = re.search(r'failed\s+(\d+)', self.qacct)
        if not match:
            logging.error('Problem obtaining is_failed')
            return
        self.status = int(match.group(1))
        if self.status != 0:
            self._is_failed = True
        else:
            self._is_failed = False
        return self._is_failed


class CPUTime(object):
    """Simple time object"""
    def __init__(self):
        super(CPUTime, self).__init__()
        self.hours = 0.
        self.minutes = 0.
        self.seconds = 0.
        self._initialized = False

    def from_combine_output(self, output):
        match = re.search(r'cpu=([\d\:]+),', output)
        if not match:
            logging.error('Problem processing: {0}'.format(output))
            return self
        components = [ float(i) for i in match.group(1).split(':') ]
        self.hours = components[0]
        self.minutes = components[1]
        self.seconds = components[2]
        self._initialized = True
        return self

    def from_number_of_seconds(self, n_seconds):
        self.hours = int(n_seconds) / 3600
        n_seconds -= 3600. * self.hours
        self.minutes = int(n_seconds) / 60
        n_seconds -= 60. * self.minutes
        self.seconds = n_seconds

    def from_number_of_hours(self, n_hours):
        self.hours = int(n_hours)
        n_hours -= self.hours
        self.minutes = int(n_hours*60.)
        n_hours -= self.minutes/60.
        self.seconds = n_hours * 3600.

    def number_of_seconds(self):
        return 3600. * self.hours + 60. * self.minutes + self.seconds

    def number_of_hours(self):
        return self.hours + self.minutes/60. + self.seconds/3600.




class ScanAccountant(object):
    """docstring for ScanAccountant"""

    scandir_printlength = 30
    n_decimals = 2
    get_is_failed_via_done_str = True
    only_good_jobs = True

    def __init__(self, scandir):
        super(ScanAccountant, self).__init__()
        self.scandir = osp.abspath(scandir)
        self.jobs = []
        self._ready = False
    
    def exists(self):
        return osp.isdir(scandir)

    def process(self):
        status = self.get_jobs()
        if status is False: return False
        self._ready = True
        self.produce_statistics()

    def get_jobs(self):
        status = self.build_ofile_cpuline_dict()
        if status is False: return False
        status = self.build_shfile_combineline_dict()
        if status is False: return False

        for sh_file in glob(osp.join(self.scandir, '*.sh')):
            job = Job(sh_file)
            
            key_sh = osp.basename(job.sh_file)
            if not(key_sh in self.shfile_combineline_dict):
                logging.warning(
                    '{0} does not have an eval combine line'
                    .format(fixed_length_basename(sh_file))
                    )
                continue
            else:
                job.combine_command = self.shfile_combineline_dict[key_sh]

            status = job.get_o_file()
            if status is False:
                continue

            key = osp.basename(job.o_file)

            if not(key in self.ofile_cputime_dict):
                logging.warning(
                    '{0} does not have an extracted cpu time'
                    .format(fixed_length_basename(job.o_file))
                    )
                continue

            job.cpu_time = self.ofile_cputime_dict[key]
            if self.get_is_failed_via_done_str: job._is_failed = not(self.ofile_done_dict[key])

            job.get_combine_settings()
            self.jobs.append(job)


    def get_failed_jobs(self):
        if not self._ready: self.process()
        if hasattr(self, 'failed_jobs'): return self.failed_jobs

        self.failed_jobs = []
        # logging.warning('Doing qacct for all jobs to determine exit status')
        for i, job in enumerate(self.jobs):
            # differentials.core.print_progress_bar(i, len(self.jobs)-1)
            if job.is_failed():
                self.failed_jobs.append(job)
        return self.failed_jobs


    def build_ofile_cpuline_dict(self):
        self.ofile_cputime_dict = {}
        self.ofile_done_dict = {}

        globpat = osp.join(self.scandir, '*.o*')
        cpuline_cmd = 'tail -n 15 {0}'.format(globpat)

        try:
            output = differentials.core.execute(cpuline_cmd, py_capture_output=True)
        except CalledProcessError:
            logging.error('Could not call tail; probably no .o files in the directory.')
            return False

        _got_queue = False
        for block in output.split('==> '):
            block = block.strip()
            if len(block) == 0: continue
            o_file = block.split('<==')[0].strip()

            key = osp.basename(o_file)

            cputime = CPUTime().from_combine_output(block)
            if cputime._initialized:
                self.ofile_cputime_dict[key] = cputime

            match_done = re.search(r'Done in', block)
            if match_done:
                self.ofile_done_dict[key] = True
            else:
                self.ofile_done_dict[key] = False

            if not _got_queue:
                match = re.search(r'queue@host = (\w+?\.q)', block)
                if not match:
                    continue
                self.queue = match.group(1)
                _got_queue = True

        if not _got_queue:
            self.queue = 'unknwn'
            logging.error('Could not determine the queue')
        return True

    def build_shfile_combineline_dict(self):
        globpat = osp.join(self.scandir, '*.sh')
        cpuline_cmd = 'tail -n 2 {0}'.format(globpat)
        output = differentials.core.execute(cpuline_cmd, py_capture_output=True)
        self.shfile_combineline_dict = {}
        for block in output.split('==> '):
            sh_file = block.split('<==')[0].strip()
            for line in block.split('\n'):
                if line.startswith('eval combine'):
                    self.shfile_combineline_dict[osp.basename(sh_file)] = line

    def produce_statistics(self):
        if self.only_good_jobs:
            jobs = [ j for j in self.jobs if not j in self.get_failed_jobs() ]
            if len(jobs) == 0:
                logging.error('All jobs failed; using the failed jobs instead')
                jobs = self.jobs
        else:
            jobs = self.jobs

        self.hours = [ job.cpu_time.number_of_hours() for job in jobs ]
        self.mean = sum(self.hours)/len(self.hours)
        self.std  = sqrt(sum([ (h - self.mean)**2 for h in self.hours]) / (len(self.hours)-1))
        self.perc_10 = percentile(self.hours, 0.1)
        self.perc_90 = percentile(self.hours, 0.9)
        self.min = min(self.hours)
        self.max = max(self.hours)

        self.n_points_per_job = list(set([ job.n_points for job in jobs ]))
        self.n_points_per_job.sort(reverse=True)
        if len(self.n_points_per_job) == 1: self.n_points_per_job = self.n_points_per_job[0]

        # Do for only one job, since they should all have the same settings anyway
        job = jobs[0]
        self.n_points = job.get_n_points_total()
        self.parameter_ranges = job.get_parameter_ranges()

        self.real_n_points = sum([ job.n_points for job in jobs ])
        self.total_cpu_time = float(sum(self.hours))


        cpu_per_point_per_job = [ job.cpu_per_point() for job in jobs ]
        cpu_per_point_per_job.sort()
        self.mean_cpu_per_point = self.total_cpu_time / self.real_n_points
        self.std_cpu_per_point = std(cpu_per_point_per_job, self.mean_cpu_per_point)
        self.perc90_cpu_per_point = percentile(cpu_per_point_per_job, 0.9)
        self.perc98_cpu_per_point = percentile(cpu_per_point_per_job, 0.98)


    def overview_queue_advice(self, queuename, max_cpu, perc90, perc98):
        d = self.get_optimal_n_points_per_job_for_queue(perc90, perc98, max_cpu)
        r = (
            '{0:<7}: {1:<4} points/job (90% pass, {2:4} jobs), '
            '{3:<4} points/job, (98% pass, {4} jobs)'
            .format(
                queuename,
                d.optimistic_n_points_per_job, d.optimistic_n_jobs,
                d.pessimistic_n_points_per_job, d.pessimistic_n_jobs,
                )
            )
        return r

    def get_optimal_n_points_per_job_for_queue(self, perc90, perc98, max_cpu):
        d = differentials.core.AttrDict()
        d.optimistic_n_points_per_job = int(float(max_cpu) / perc90)
        d.optimistic_n_jobs           = int(math.ceil(float(self.real_n_points) / d.optimistic_n_points_per_job)) if d.optimistic_n_points_per_job != 0. else 0.
        d.pessimistic_n_points_per_job = int(float(max_cpu) / perc98)
        d.pessimistic_n_jobs           = int(math.ceil(float(self.real_n_points) / d.pessimistic_n_points_per_job)) if d.pessimistic_n_points_per_job != 0. else 0.
        return d


    def only_allq_recommendation(self):
        self.only_one_queue_recommendation(10.0)

    def only_shortq_recommendation(self):
        self.only_one_queue_recommendation(1.5)

    def only_longq_recommendation(self):
        self.only_one_queue_recommendation(96.0)

    def only_one_queue_recommendation(self, max_cpu):
        if not(self._ready):
            status = self.process()
            if status is False: return False
        queue = self.get_optimal_n_points_per_job_for_queue(
            self.perc90_cpu_per_point, self.perc98_cpu_per_point, max_cpu
            )
        print '{0}: {1:<4} pnts/job (98%)'.format(
            fixed_length_basename(self.scandir, 68),
            queue.pessimistic_n_points_per_job
            )


    def overview(self):
        if not(self._ready):
            status = self.process()
            if status is False: return False

        l = [ fixed_length_basename(self.scandir, 70) ]

        l.append(
            'cpu/pnt = {0:.{ndec}f} h +/- {1:.{ndec}f}'
            '  (={2:3.{ndec}f} min +/- {3:.{ndec}f})'
            '  (total pnts = {4:<4})'
            .format(
                self.mean_cpu_per_point, self.std_cpu_per_point,
                self.mean_cpu_per_point*60., self.std_cpu_per_point*60.,
                self.real_n_points,
                ndec = self.n_decimals                
                )
            )

        l.append(
            'cpu/pnt < {0:.{ndec}f}h ({1:.{ndec}f}min) (90%), cpu/pnt < {2:.{ndec}f}h ({3:.{ndec}f}min) (98%)'
            .format(
                self.perc90_cpu_per_point, self.perc90_cpu_per_point*60.,
                self.perc98_cpu_per_point, self.perc98_cpu_per_point*60.,
                ndec = self.n_decimals                
                )
            )

        l.append(
            'cpu/job = {0:.{ndec}f} h +/- {1:.{ndec}f}'
            '  (={2:3.{ndec}f} min +/- {3:.{ndec}f})'
            '  (total jobs = {4}, min cpu = {5:.{ndec}f}, max cpu = {6:.{ndec}f})'
            .format(
                self.mean, self.std,
                self.mean*60., self.std*60.,
                len(self.jobs),
                self.min, self.max,
                ndec = self.n_decimals
                )
            )

        l.append(
            'points/job: {0}  ({1}/{2} points scanned, {3} jobs on {4}, total cpu = {5:.{ndec}f})'
            .format(
                self.n_points_per_job,
                self.real_n_points,
                self.n_points,
                len(self.jobs),
                self.queue,
                self.total_cpu_time,
                ndec = self.n_decimals
                )
            )

        l.append('Parameter ranges:')
        par_str_length = max(map(len, self.parameter_ranges.keys()))
        for par in sorted(self.parameter_ranges.keys()):
            l.append(
                '  {0:{parlen}} = {1:.{ndec}f} to {2:.{ndec}f}'
                .format(
                    par, self.parameter_ranges[par][0], self.parameter_ranges[par][1],
                    parlen=par_str_length, ndec=self.n_decimals
                    )
                )

        l.append(self.overview_queue_advice('shortq', 1.5, self.perc90_cpu_per_point, self.perc98_cpu_per_point))
        l.append(self.overview_queue_advice('allq', 10.0,  self.perc90_cpu_per_point, self.perc98_cpu_per_point))
        l.append(self.overview_queue_advice('longq', 96.0, self.perc90_cpu_per_point, self.perc98_cpu_per_point))

        failed_jobs = self.get_failed_jobs()
        frac_failed = float(len(failed_jobs)) / len(self.jobs)
        l.append(
            '{0:.1f}% failed ({1}/{2})'
            .format(
                100. * frac_failed,
                len(failed_jobs), len(self.jobs)
                )
            )

        print '\n    '.join(l)



def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1


def fixed_length_string(s, max_length=30):
    if len(s) > max_length:
        n_head = int(0.5*max_length) - 1
        n_tail = int(0.5*max_length) - 1
        if max_length % 2 == 0: n_head -= 1
        s = s[:n_head] + '...' + s[-n_tail:]
    else:
        s = '{0:<{L}s}'.format(s, L=max_length)
    return s

def fixed_length_basename(s, max_length=30):
    if osp.isdir(s):
        s = osp.split(osp.abspath(s))[1]
    elif osp.isfile(s):
        s = osp.basename(s)
    else:
        raise ValueError('{0} is not a file or a dir'.format(s))
    return fixed_length_string(s, max_length)


def tail(f, n):
    stdin,stdout = os.popen2("tail -n "+str(n)+" "+f)
    stdin.close()
    lines = stdout.readlines()
    stdout.close()
    return lines

def std(X, mean):
    N = len(X)
    return sqrt(
        sum([ (x-mean)**2 for x in X ]) / (N-1)
        )
