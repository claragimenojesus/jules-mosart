#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import xarray as xr

OUTDIR = '/home/clara/project/rahu/mosart-output/rahu'
OUTLET_X=-72.6271
OUTLET_Y=-13.0045

def main():
    fs = [f for f in os.listdir(os.path.join(OUTDIR)) if re.match(r'.*_[0-9]{4}_[0-9]{2}\.nc', f)]
    fs.sort()
    fs = [os.path.join(OUTDIR, f) for f in fs]
    x = xr.open_mfdataset(fs)
    xi = x.sel(lat=OUTLET_Y, lon=OUTLET_X, method='nearest')
    xi = xi.to_dataframe()
    xi.to_csv(os.path.join(OUTDIR, 'output_pt.csv'))

if __name__ == '__main__':
    main()
