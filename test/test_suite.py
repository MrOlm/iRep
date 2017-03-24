#!/usr/bin/env python

'''
Testing suite for iRep
'''

from subprocess import call
import os
import shutil
import glob

import iRep.iRep_filter

def load_data(datum):
    '''
    return the system path to the datum requested
    '''

    if datum == 'test_outdir':
        loc = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), \
        'tmp/test_outdir')
        return loc

    elif datum =='test_genome':
        loc = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), \
        '../sample_data/l_gasseri.fna')
        return loc

    elif datum =='test_stb':
        loc = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), \
        '../sample_data/l_gasseri.stb')
        return loc

    elif datum == 'test_sams':
        loc = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), \
        '../sample_data/*.sam')
        return glob.glob(loc)

    elif datum == 'solution_irep':
        loc = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), \
        '../sample_output/test.iRep.tsv')
        return loc

    else:
        raise AttributeError("datum {0} is not found".format(datum))

def execute_cmd(cmd, dry=False, shell=False, stdout=None, stderr=None):
    '''
    just a wrapper to call commands
    '''

    devnull = open(os.devnull, 'w')
    if stdout == None:
        stdout = devnull

    if stderr == None:
        stderr = devnull

    if not shell:
        print(' '.join(cmd))
        #if not dry: call(cmd, stderr=stderr, stdout=stdout)
        if not dry: call(cmd)

    else:
        cmd = ' '.join(cmd)
        print(cmd)
        #if not dry: call(cmd, shell=True, stderr=stderr, stdout=stdout)
        if not dry: call(cmd, shell=True)

    return

class TestiRep():
    '''
    class to do regression tests of iRep
    '''

    def run_all(self):
        self.setUp()

        self.pipe_test2()
        self.tearDown()

        self.pipe_test()
        self.tearDown()

        self.stb_test()
        self.tearDown()

        self.regression_test()

    def setUp(self):
        self.outdir = load_data('test_outdir')
        if os.path.exists(self.outdir):
            shutil.rmtree(self.outdir)
        os.makedirs(self.outdir)

        self.genome = load_data('test_genome')
        self.stb = load_data('test_stb')
        self.sams = load_data('test_sams')
        self.irep_solution = load_data('solution_irep')

    def tearDown(self):
        if os.path.exists(self.outdir):
            shutil.rmtree(self.outdir)
        os.makedirs(self.outdir)

    def pipe_test(self):
        '''
        test the functioning of pipeing in .sam files
        '''

        # generate command
        test_out = os.path.join(self.outdir, 'test.iRep')
        cmd = ['cat', self.sams[0], '|', 'iRep', '-f', self.genome, '-o', test_out, \
            '-s', '-']

        # run command
        execute_cmd(cmd, shell=True)

        # verify output
        out_pdf = test_out + '.pdf'
        assert os.stat(out_pdf).st_size > 0

        out_tsv = test_out + '.tsv'
        assert os.stat(out_tsv).st_size > 0

        sol_tsv = self.irep_solution

        sol_ireps = self.extract_ireps(sol_tsv)
        test_ireps = self.extract_ireps(out_tsv)

        assert sorted(sol_ireps)[0] == sorted(test_ireps)[0]

    def pipe_test2(self):
        '''
        test the functioning of pipeing in .sam files, with -b instead of -f
        '''

        # generate command
        test_out = os.path.join(self.outdir, 'test.iRep')
        cmd = ['cat', self.sams[0], '|', 'iRep', '--no-gc-correction', '-b', self.stb, '-o', test_out, '-s', '-']

        # run command
        execute_cmd(cmd, shell=True)

        # verify output
        out_pdf = test_out + '.pdf'
        assert os.stat(out_pdf).st_size > 0

        out_tsv = test_out + '.tsv'
        assert os.stat(out_tsv).st_size > 0

        sol_tsv = self.irep_solution

        sol_ireps = self.extract_ireps(sol_tsv)
        test_ireps = self.extract_ireps(out_tsv)

        assert abs(sorted(sol_ireps)[0] - sorted(test_ireps)[0]) <= .05

    def regression_test(self):
        '''
        just call iRep as an external script, compare output to sample output
        '''

        # generate command
        test_out = os.path.join(self.outdir, 'test.iRep')
        cmd = ['iRep', '-f', self.genome, '-o', test_out, '-s'] + self.sams

        # run command
        execute_cmd(cmd)

        # verify output
        out_pdf = test_out + '.pdf'
        assert os.stat(out_pdf).st_size > 0

        out_tsv = test_out + '.tsv'
        assert os.stat(out_tsv).st_size > 0

        sol_tsv = self.irep_solution

        sol_ireps = self.extract_ireps(sol_tsv)
        test_ireps = self.extract_ireps(out_tsv)

        assert sorted(sol_ireps) == sorted(test_ireps)

    def stb_test(self):
        '''
        test the functioning of .stb input instead of .fasta
        '''

        # generate command
        test_out = os.path.join(self.outdir, 'test.iRep')
        cmd = ['iRep', '--no-gc-correction', '-b', self.stb, '-o', test_out, '-s'] + self.sams

        # run command
        execute_cmd(cmd)

        # verify output
        out_pdf = test_out + '.pdf'
        assert os.stat(out_pdf).st_size > 0

        out_tsv = test_out + '.tsv'
        assert os.stat(out_tsv).st_size > 0

        sol_tsv = self.irep_solution

        sol_ireps = self.extract_ireps(sol_tsv)
        test_ireps = self.extract_ireps(out_tsv)

        # Allow some liniency because it's not GC corrected, and the solution is
        for solution_value, test_value in zip(sorted(sol_ireps), sorted(test_ireps)):
            assert abs(solution_value - test_value) <= .05, "{0} vs {1}".format(\
                solution_value, test_value)

    def extract_ireps(self, tsv):
        '''
        from the path to the .tsv file, return a list of iRep values
        '''
        values = []
        irep = iRep.iRep_filter.parse_tables([tsv])
        for genome in irep.keys():
            for sam in irep[genome].keys():
                val = irep[genome][sam]['iRep']
                values.append(float("{0:.2f}".format(val)))
        return values

def test_short():
    '''
    tests that shouldn't take very long
    '''
    pass

def test_long():
    '''
    full test suite (could take a while)
    '''

    t = TestiRep()
    t.run_all()
    return

if __name__ == "__main__":
    test_short()
    test_long()
    print("All tests passed!")
