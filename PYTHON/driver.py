#!/usr/bin/python3

import math
import numpy as np
import sys
from afwl_1kt import *

a = afwl1kt()

yldkt  = 1.0
hobm   = 0.0
grm    = 100.0

grft  = grm * 3.2808
hobft = hobm * 3.2808

tarv = timar(grft,hobft,yldkt)

fo = open('out','w')
t = 0.0
dt = 0.001
while (t <= 10.0):
   a.shock(yldkt,hobm,t,grm)
   fo.write('%15.5e %15.5e %15.5e %15.5e\n' % (t,a.opr*a.cp2psi,a.odr,a.vr))
   t = t + dt

fo.close()
