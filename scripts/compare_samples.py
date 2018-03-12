#!/usr/bin/env python

import os
import sys
from csv import DictReader


samples=dict()
with open(sys.argv[1], 'r') as fin:
    reader = DictReader(fin)
    fnames = []
    for r in reader:
        fname, patient, stype = r['Filename'], r['Patient'], r['Type']
        if patient in samples:
            if stype in samples[patient]:
                samples[patient][stype].append(fname)
            else:
                samples[patient][stype] = [fname]
        else:
            samples[patient] = dict()
            samples[patient][stype] = [fname]

print("PATH,TUMOR,NORMAL")
for patient in samples.keys():
    stypes = [s for s in samples[patient].keys()]
    for i in range(0, len(stypes)):
        for j in range(1, len(stypes)):
            if i > j: break
            st1 = stypes[i]
            st2 = stypes[j]
            if st1 is st2: continue
            fn1 = samples[patient][st1]
            fn2 = samples[patient][st2]
            for f1 in fn1:
                for f2 in fn2:
                    #print("%s\t%s: %s\t%s: %s" % (patient, st1, f1, st2, f2))
                    path = os.path.dirname(f1)
                    fname1 = os.path.splitext(os.path.basename(f1))[0]
                    fname2 = os.path.splitext(os.path.basename(f2))[0]
                    print("%s,%s,%s" % (path, fname2, fname1))
