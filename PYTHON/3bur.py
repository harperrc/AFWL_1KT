#!/usr/bin/python3

import math
import numpy as np
import sys
from afwl_1kt import *

#  N.B. afwl_1kt uses nukelib (see SearchForOverpressure repo)

#  distance to burst (Miles)

dist = [10,20,30]

#  height of burst (meters)

hobm = 6000.0

#  yield (Kt)

yld  = 300.0

#  time step (sec)

dt   = 0.01

for d in dist:
   a = afwl1kt()

   grm    = d * 5280.0 / 3.2808

   oname = '%d.out' % (d)
   fo = open(oname,'w')
   t = 0.0
   while (t <= 200.0):
      a.shock(yld,hobm,t,grm)
      mph = a.vr * 2.23694
      fo.write('%15.5e %15.5e %15.5e %15.5e\n' % (t,a.opr*a.cp2psi,a.odr,mph))
      t = t + dt

   fo.close()
