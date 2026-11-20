#!/usr/bin/env python3
"""General ED2 combinatorial counting and overflow regression (no BLAS needed)."""
import argparse
import math
from pathlib import Path
import subprocess
import tempfile

DRIVER = """program cpc_count_probe
  use state_lists
  implicit none
  integer :: n
  integer, allocatable :: caps(:)
  read(*,*) n
  allocate(caps(n))
  caps=1
  call get_cum_ltd_rep_combination(n,3,3,caps)
  print *, 'COUNT',n,cum_ltd_rep_combination(n,3)
end program cpc_count_probe
"""

def run(cmd, work, stdin=None):
    return subprocess.run(cmd, cwd=str(work), input=stdin, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, universal_newlines=True, timeout=45)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source', type=Path, default=Path('.'))
    p.add_argument('--compiler', default='gfortran')
    a=p.parse_args()
    source=a.source.resolve()/'source'
    with tempfile.TemporaryDirectory(prefix='ed2-count-') as name:
        w=Path(name)
        (w/'probe.f90').write_text(DRIVER)
        cp=run([a.compiler, '-O0', '-fopenmp', '-c',
                str(source/'input_param.f90'),str(source/'state_lists.f90')],w)
        assert cp.returncode==0,cp.stderr
        cp=run([a.compiler, '-fopenmp', 'probe.f90', 'state_lists.o',
                'input_param.o', '-o', 'probe'],w)
        assert cp.returncode==0,cp.stderr
        for n in (64,216,512,1000):
            r=run([str(w/'probe')],w,str(n)+'\n')
            expected=n*(n-1)*(n-2)//6
            assert r.returncode==0,(n,r.stderr)
            assert r.stdout.split()==['COUNT',str(n),str(expected)],(n,r.stdout)
            print('PASS',n,expected)
        r=run([str(w/'probe')],w,'2500\n')
        assert r.returncode!=0
        assert 'Hilbert count exceeds supported integer range' in r.stderr
        print('PASS overflow rejection')
    print('GENERIC_COUNT_GUARD_PASS')

if __name__=='__main__':
    main()
