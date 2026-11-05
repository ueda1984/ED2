#!/usr/bin/env python3
"""Compare full current-run ED2 outputs with unchanged historical references.

Run after make check-long in a fresh build directory. NumPy is test-only.
Energy tolerances: atol=rtol=1e-9; observables: atol=rtol=1e-6.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='new JSON report path')
    args = parser.parse_args()
    build = args.build.resolve(strict=True)
    report = {'build':str(build), 'comparisons':[], 'examples':[]}
    expected_names = {'chain','j1j2_chain','mixed_spin_chain','square','triangular',
                      'honeycomb','kagome','cubic','cubic_sp_HB'}
    refroot = build/'examples/ref_data'
    if {p.name for p in refroot.iterdir() if p.is_dir()} != expected_names:
        parser.error('historical reference example membership differs from the nine-case baseline')
    for name in sorted(expected_names):
        actual = build/'examples'/name
        reference = refroot/name
        item = {'name':name,'input_sha256':{}}
        for p in sorted(actual.glob('*.dat')):
            if p.name.startswith('list_') or p.name in ['input.dat','input_param.dat']:
                item['input_sha256'][p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        text = (actual/'output.dat').read_text()
        item['dimensions'] = {k:int(v) for k,v in re.findall(r'\b(THS(?:\(k\))?)\s*=\s*(\d+)',text)}
        report['examples'].append(item)
        for filename, indices, tolerance in [('eigenvalues.dat',1,1e-9),
                                              ('local_mag.dat',1,1e-6),
                                              ('two_body_cf_xyz.dat',2,1e-6)]:
            result = {'example':name,'file':filename,'atol':tolerance,'rtol':tolerance}
            try:
                a = np.loadtxt(actual/filename,ndmin=2)
                b = np.loadtxt(reference/filename,ndmin=2)
                if a.shape != b.shape or a.shape[1] <= indices:
                    raise ValueError('output shape differs or has no numerical columns')
                if not np.isfinite(a).all() or not np.isfinite(b).all():
                    raise ValueError('nonfinite output/reference values')
                np.testing.assert_array_equal(a[:,:indices],b[:,:indices])
                np.testing.assert_allclose(a[:,indices:],b[:,indices:],atol=tolerance,rtol=tolerance)
                result.update(status='passed', max_absolute_error=float(np.max(np.abs(a[:,indices:]-b[:,indices:]))))
                if filename=='eigenvalues.dat':
                    item['eigenvalues'] = a[:,indices:].tolist()
            except (OSError,ValueError,AssertionError) as error:
                result.update(status='failed',error=str(error))
            report['comparisons'].append(result)
    report['passed'] = sum(c['status']=='passed' for c in report['comparisons'])
    report['failed'] = sum(c['status']=='failed' for c in report['comparisons'])
    # Refuse to overwrite an earlier run's report.
    with args.output.open('x') as stream:
        stream.write(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(f"Historical comparisons: {report['passed']} passed; {report['failed']} failed")
    return int(report['failed']!=0)


if __name__=='__main__':
    sys.exit(main())
