# -*- coding: utf-8 -*-

import sys
sys.path.append('../CoCoPy/modules/')
import analysis.autofit as afit

fname = sys.argv[1]
N = int(sys.argv[2])
scratch = sys.argv[3]

task = afit.load_task(fname, N)[0][0]
task.report = task.prefix
task.prefix = scratch + '/' + task.prefix

task.run_autofit_GWDG()
task.write_report()
task.del_files()
