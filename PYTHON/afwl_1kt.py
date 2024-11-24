#!/usr/bin/python3

# AFWL 1kt standard after c. needham

import math
import numpy as np
import sys

class matm62:
   def __init__(self):

#  pad fortran arrays with 0 to keep original logic (i'm lazy)

      self.nz    = 21
      self.rhoz  = 1.225e-3
      self.tabat = [0.0,8.31440000e+07, 6.36748800e+08, 9.80665000e+02, 2.89644000e+01]

      self.tabz  = [0.0,
                    0.            , 1.10190000e+06, 2.00630000e+06, 3.21620000e+06,
                    4.73500000e+06, 5.24290000e+06, 6.15910000e+06, 7.99940000e+06,
                    9.00000000e+06, 1.00000000e+07, 1.10000000e+07, 1.20000000e+07,
                    1.50000000e+07, 1.60000000e+07, 1.70000000e+07, 1.90000000e+07,
                    2.20000000e+07, 3.00000000e+07, 4.00000000e+07, 5.00000000e+07,
                    6.00000000e+07, 7.00000000e+07]

      self.tabl  = [0.0,
                    -6.49291769e-05, 9.28018576e-08, 9.86255062e-06, 2.77080373e-05,
                    -1.72248474e-07,-1.96000240e-05,-3.91696897e-05, 1.60821507e-07,
                     2.98166740e-05, 5.02020080e-05, 9.97762300e-05, 2.00108809e-04,
                     1.49589024e-04, 1.00407490e-04, 6.97598500e-05, 6.68801467e-05,
                     3.49035000e-05, 3.31099360e-05, 2.58868500e-05, 1.71252960e-05,
                     1.09162420e-05]

      self.tabt = [0.0,
                   2.88150000e+02, 2.16604540e+02, 2.16688470e+02, 2.28621170e+02,
                   2.70704137e+02, 2.70616652e+02, 2.52659110e+02, 1.80575130e+02,
                   1.80736048e+02, 2.10552722e+02, 2.60754730e+02, 3.60530960e+02,
                   9.60857386e+02, 1.11044641e+03, 1.21085390e+03, 1.35037360e+03,
                   1.55101404e+03, 1.83024204e+03, 2.16134140e+03, 2.42020990e+03,
                   2.59146286e+03, 2.70062528e+03]

      self.tabp = [0.0,
                   1.01325000e+06, 2.26320000e+05, 5.47486994e+04, 8.68013979e+03,
                   1.10904998e+03, 5.90004987e+02, 1.82098959e+02, 1.03769924e+01,
                   1.64379881e+00, 3.00749781e-01, 7.35439451e-02, 2.52169805e-02,
                   5.06169601e-03, 3.69429709e-03, 2.79259780e-03, 1.68519867e-03,
                   8.67381898e-04, 1.94317430e-04, 4.15743157e-05, 1.13023464e-05,
                   3.55894454e-06, 1.22936355e-06]

      self.init   = True

      self.re    = 0.0
      self.cons  = 0.0
      self.trho  = 0.0

   def returnConditions(self,ty):

      if (self.init):
         self.init = False

         self.re   = self.tabat[2]
         self.cons = self.re * self.re * self.tabat[3] * self.tabat[4] / self.tabat[1]
         self.trho = self.tabt[1] * self.rhoz * 1.0e3

      tty = ty * 100.0
      j   = 1

      if (tty > 0.0):
         j = self.nz
         for k in range(2,self.nz+1):
            dif = tty - self.tabz[k]
            if (dif <= 0.0):
               if (abs(dif) <= 1.0e-3):
                  j = k
                  break
               elif (dif < 0.0):
                  j = k - 1
                  break

      dum2  = (self.tabz[j] - tty) / ((self.re + tty) * (self.re + self.tabz[j]))
      dum3  = (self.re + self.tabz[j]) / (self.re + tty)
      var1  = (self.re + self.tabz[j]) * self.tabl[j] - self.tabt[j]
      rvar1 = 1.0 / var1
      var2  = ((tty - self.tabz[j]) * self.tabl[j] + self.tabt[j]) / self.tabt[j]
      fs    = -self.cons * (((self.tabl[j] * math.log(dum3 * var2)) * rvar1 + \
               dum2) * rvar1)

      wsp   = self.tabp[j] * math.exp(fs)
      wst   = self.tabt[j] + self.tabl[j] * (tty - self.tabz[j])
      wsr   = self.trho * wsp / (wst * self.tabp[1])
      wsp   = 0.10 * wsp
      cs    = math.sqrt(1.4 * wsp / wsr)

      return (wsp,cs,wsr,wst)

class afwl1kt:
   def __init__(self):

      self.atm       = matm62()

      self.cp2psi = 0.000145038

#  air constants

      self.rhoz      = 1.225e-3
      self.opz       = 1.01325e6
      self.gamm1     = 0.404574

#  wfpkod constants

      self.told      = 0.0
      self.psca      = 0.10
      self.vsca      = 0.010
      self.rhosca    = 1000.00

#  scalkt constants

      self.p1        = self.opz / 10.0
      self.c1        = 3.4029399e2
      self.r1        = self.rhoz * 1000.0
      self.t1        = 288.15

      self.tscale    = 0.0
      self.vscale    = 0.0
      self.pscale    = 0.0
      self.dscale    = 0.0
      self.cscale    = 0.0

#  common /wfrt/

      self.prad      = 0.0
      self.oppk      = 0.0
      self.odpk      = 0.0
      self.vpk       = 0.0
      self.opr       = 0.0
      self.odr       = 0.0
      self.vr        = 0.0
      self.rzp       = 0.0
      self.rzd       = 0.0
      self.rzv       = 0.0
      self.opmn      = 0.0 
      self.odmn      = 0.0
      self.vmn       = 0.0

#  peak

      self.prado = 0.0
      self.oppko = 0.0
      self.odpko = 0.0
      self.vpko  = 0.0

# wfprmt 

      self.ttold_1 = 0.0

# wfdrmt

      self.ttold_2 = 0.0

# wfvrmt

      self.ttold_3 = 0.0

   def air(self,eee,rrr):

      e     = eee * 1.0e-10
      rholn = math.log(773.39520495 * rrr)

#  original code multiplied 2 constants together
#  for e1, did it once

      e1    =  0.008717948717600

#  original code had check if (abs(e1).lt.5.0)goto 20
#  e1 is ALWAYS < 5 since its a fixed value
#  so we dropped the 2 if checks out

      ws  = (8.5 + 0.15504314 * rholn - e) * \
            math.exp(-0.05 * rholn + 0.02531780798)
      ws  = max(-60.0,min(60.0,ws))
      ws  = 1.0 / (math.exp(-ws) + 1.0)
      fo  = math.exp(-0.22421524664 * e) * ws
      fon = math.exp(-0.15082856259 * e) * (1.0 - ws)

      beta = 0.0
      if (e > 1.0):
         beta = (6.9487e-3 * ws + 1.38974e-2) * math.log(e)

      e2  = (e - 40.0) * 0.33333333333333333
      fn  = 0.0

      if (abs(e2) < 5.0):
         ws = (e - math.exp(0.0157 * rholn + 3.806662489)) * \
               math.exp(-0.085 * rholn - 1.38629436)
         ws = 1.0 / (math.exp(-ws) + 1.0)
         fn  = math.exp(-0.039215686275 * e) * ws
      elif (e2 > 0.0):
         fn  = math.exp(-0.039215686275 * e)
         ws  = 1.0
      else:
         fn = 0.0
         ws = 0.0

      beta = beta + ws * (0.045 - beta)
      ws   = (e - 160.0) / max(30.0 + 3.474356 * rholn,6.0)
      fe   = 0.0

      if (ws > -5.0):
         fe = 1.0 / (math.exp(-ws) + 1.0)

      gm   = (0.161 + 0.225 * fo + 0.280 * fon +  \
             0.137 * fn + 0.050 * fe) *  \
             math.exp(beta * rholn)

      return (gm)

   def cubrt(self,x):
      val = math.exp(math.log(abs(x)) / 3.0)
      if (x < 0.0):
        val = -1.0 * val
      return (val)

   def peak(self,t,ra):
   
      r = ra * 100.0

      if (t != self.told):
         self.rzp  = self.wfzr(t)
         self.rzd  = self.wfdzr(t)
         self.rzv  = self.wfvzr(t)
         self.prad = self.wfpr(t)

         self.prado = self.prad * self.vsca
         self.oppk  = self.wfpkop(self.prad)

         self.oppko = self.oppk * self.psca
         self.odpk  = self.wfpkod(self.prad)
         self.odpko = self.odpk * self.rhosca

         self.vpk   = self.wfpkv(self.prad)
         self.vpko  = self.vpk * self.vsca

         self.told = t

      self.opr  = 0.0
      self.odr  = 0.0
      self.vr   = 0.0
      self.opro = 0.0
      self.odro = 0.0
      self.vro  = 0.0

      if (r <= self.prad):
         self.wfprmt(t,r)
         self.opro = self.opr * self.psca

         self.wfdrmt(t,r)

         odw = self.well(t,r)
         self.odro = odw * self.rhosca

         self.wfvrmt(t,r)
         self.vro = self.vr * self.vsca

   def scalkt(self,hfpt,wb):

# hfpt   height of point (m)
# wb     yield (kt)

      (p3,c3,r3,t3) = self.atm.returnConditions(hfpt)

      rp3 = self.cubrt(self.p1 / p3)
      wsw = self.cubrt(wb)

      self.vscale = c3 / self.c1
      self.dscale = r3 / self.r1
      self.tscale = c3 / (wsw * rp3 * self.c1)
      self.cscale = rp3 * wsw
      self.pscale = p3 / self.p1

   def shock(self,yld,height,tim,ra):

      alt = height

      self.scalkt(alt,yld)

      t = tim * self.tscale
      r = ra / self.cscale

      self.peak(t,r)

      self.opr = self.opro * self.pscale
      self.odr = self.odro * self.dscale
      self.vr  = self.vro * self.vscale

   def speed(self,t):

      r1 = 24210.0 * t**0.371 * (1.0 + (1.23 * t + 0.123) * \
           (1.0 - exp(-26.25 * t**0.79)))

      r2 = (1.0 - 0.03291 * t**(-1.086)) * \
           (33897.0 * t + 8490.0) + 8.36e3 + 2.5e3 * \
           math.log(t) + 800.0 * t**(-0.21)

      v1 = 10086.685 * t**(-0.629) + 40826.049 * t**0.371 + \
           math.exp(-26.25 * t**0.79) * \
           (617527.5 * t**1.161 - 40826.048 * t**0.371 - \
            1104.7749 * t**(-0.629) +  \
            61752.75 * t**0.161)

      v2 = 33897.0 + 95.937326 * t**(-1.086) +  \
           303.43481 * t**(-2.086) + 2.5e3 / t -  \
           168.0 * t**(-1.21)

      v3 = (v2 * (t - 0.21) + v1 * (0.28 - t)) / 0.07

      if (t < 0.21):
         spd = v1
      elif (t > 0.28):
         spd = v2
      else:
         spd = v3

      return (spd)

   def well(self,t,r):

      depth = self.odr

      rmsw  = min(0.7 * self.rzd,1.55e4)

      if (t <= 1.2 and r <= rmsw):
         depth = max(-1.21e-3,-1.5e-3 * math.exp(-8.0e-13 * r * r * r))
         if (t > 1.0):
            depth = (1.2 - t) * depth * 5.0
         trm = 0.90 * rmsw

         if (r >= trm):
            depth = 10.0 * (depth * (rmsw - r) + self.odr * (r - trm)) / rmsw

      return (depth)

   def wfdrmt(self,t,r):

      rpk = self.prad

      if (t >= 0.20):
         rpk = rpk * 1.0e-5
         r   = r * 1.0e-5

         if (t != self.ttold_2):
            rz         = self.rzd
            rmn        = rz - 9.7163e3 * t**0.12115
            rz         = rz * 1.0e-5
            rmn        = rmn * 1.0e-5
            self.odmn  = -0.50 * self.odpk + 2.2e-5*t**(-1.6026)
            rneg       = rz - rmn
            rpls       = rpk - rz
            rnp        = rpk - rmn
            a1n        = self.odpk / rpls
            b1n        = self.odpk - a1n * rpk
            odmn1n     = a1n * rmn + b1n
            fmlt       = abs(self.odmn / self.odpk)
            odmhy      = self.odmn + fmlt * (odmn1n - self.odmn)
            alpha      = (rnp / (odmhy - self.odpk) + rpls / self.odpk) / rneg
            beta       = -rpls * (alpha + 1.0 / self.odpk)
            fngz       = rnp / rpls
            dnom       = alpha * rnp + beta
            bcrmn      = 1.0 - self.odmn / (rnp / dnom + self.odpk)
            crmn1b     = math.log(bcrmn)
            cgz1       = (beta * (1.0 / bcrmn - 1.0) / dnom) / \
                         ((fngz * crmn1b * rnp**(fngz - 1.0)) *  \
                         (rnp + self.odpk * dnom))
            cgz        = math.exp(cgz1)
            bgz1       = crmn1b / cgz**(rnp**fngz)
            bgz        = math.exp(bgz1)
         
         rbr  = rpk - r
         gr   = 1.0 - bgz**(cgz**(rbr**fngz))

         if (r > rz):
            gr = (rpk - r) / rpls * gr + (r - rz) / rpls

         hr       = rbr / (alpha * rbr + beta) + self.odpk
         self.odr = gr * hr
         rpk      = rpk * 1.0e5
         r        = r * 1.0e5

      if (t <= 1.0):
         if (t != self.ttold_2):
            a = -1.2e-3
            c = math.log(max(1.0e-20,-a / (self.odpk - a))) / (self.rzd - rpk)
            b = (self.odpk - a) * math.exp(-c * rpk)

         wflt = a + b * math.exp(c * r)

         if (t <= 0.20):
            self.odr = wflt
         else:
            self.odr = (wflt * (1.0 - t) + self.odr * (t - 0.20)) * 1.25

      self.ttold_2 = t
   
   def wfdzr(self,t):

      rwfdzr = 0.0

      if (t < 0.0):
         print('wfdzr t < 0 %15.5f' % (t))
         sys.exit()

      if (t > 0.0 and t < 0.265):
         rwfdzr = 2.568e4 * t**0.395
      else:
         rwfdzr = (1.0 - 0.03499 *t**(-1.068)) * (33897.0 * t + 8490.0) + 500.0

      return (rwfdzr)

   def wfpkod(self,dummy):

      op    = self.oppk
      rtio  = op / self.opz
      p     = op + self.opz
      gmone = self.gamm1
      gamma = 1.0 + self.gamm1
      gamra = gamma

#  iterate to solution

      loop = 1
      while (loop):

         rho1 = self.rhoz * ((2.0 * gamra + (gamra + 1.0) * rtio) / \
                             (2.0 * gamra + (gamra - 1.0) * rtio))

         ee = p / (gmone * rho1)

         gmone = self.air(ee,rho1)

         gamrao = gamra

         gamra  = 2.0 * gmone / (2.5 * gmone + 1.0) + 1.0

         if (abs(gamra - gamrao) < 1.0e-5):
            return (rho1 - self.rhoz)

         loop = loop + 1

         if (loop > 20):
            break

      return (rho1 - self.rhoz)

   def wfpkop(self,r):

      rr   = 1.0 / r

      rtio = 2.24517e-5 * r
      cf   = math.sqrt(math.log(rtio + 3.0 * math.exp(-(math.sqrt(rtio) / 3.0))))

      rwfpkop = ((3.04e18 * rr + 1.13e14) * rr + 7.9e9 / cf) * rr

      return (rwfpkop)
   
   def wfpkv(self,dummy):

      rwfpkv = math.sqrt(self.oppk * self.odpk / (self.rhoz * (self.rhoz + self.odpk)))

      return (rwfpkv)

   def wfpr(self,t):

      r1 = 24210.0 * t**0.371 * (1.0 + (1.23 * t + 0.123) * (1.0 - math.exp(-26.25 * t**0.79)))

      r2 = (1.0 - 0.03291 * t**(-1.086)) * \
          (33897.0 * t + 8490.0) + \
          8.36e3 + 2.5e3 * \
          math.log(t) + 800.0 * t**(-0.21)

      if (t < 0.21):
         rwfpr = r1
      elif (t > 0.28):
         rwfpr = r2
      else:
         rwfpr = (r2 * (t - 0.21) + r1 * (0.28 - t)) / 0.07

      return (rwfpr)

   def wfprmt(self,t,r):

      rpk = self.prad

      if (t >= 0.10):
         rpk = rpk * 1.0e-5
         r   = r * 1.0e-5

         if (t != self.ttold_1):
            rz        = self.rzp
            rmn       = rz - 9.7163e3 * t**0.12115
            rz        = rz * 1.0e-5
            rmn       = rmn * 1.0e-5
            self.opmn = 2.2e4 / (t * math.sqrt(t)) - 0.50 * self.oppk
            rneg      = rz - rmn
            rpls      = rpk - rz
            rnp       = rpk - rmn
            a1n       = self.oppk / rpls
            b1n       = self.oppk - a1n * rpk
            opmn1n    = a1n * rmn + b1n
            fmlt      = abs(self.opmn / self.oppk)
            opmhy     = self.opmn + fmlt * (opmn1n - self.opmn)
            alpha     = (rnp / (opmhy - self.oppk) + rpls / self.oppk) / rneg
            beta      = -rpls * (alpha + 1.0 / self.oppk)
            fngz      = rnp / rpls
            dnom      = alpha * rnp + beta
            bcrmn     = 1.0 - self.opmn / (rnp / dnom + self.oppk)
            crmn1b    = math.log(bcrmn)
            cgz1      = (beta * (1.00 / bcrmn - 1.00) / dnom) / \
                         ((fngz * crmn1b * rnp**(fngz - 1.00)) *  \
                         (rnp + self.oppk * dnom))
            cgz       = math.exp(cgz1)
            bgz1      = crmn1b / cgz**(rnp**fngz)
            bgz       = math.exp(bgz1)

         rbr   = rpk - r
         gr    = 1.00 - bgz**(cgz**(rbr**fngz))

         if (r > rz):
            gr = (rpk - r) / rpls * gr + (r - rz) / rpls

         hr       = rbr / (alpha * rbr + beta) + self.oppk
         self.opr = gr * hr
         rpk      = rpk * 1.0e5
         r        = r * 1.0e5

      if (t < 0.95):
         if (t != self.ttold_1):
            a = 5.446e4 * t**(-1.22) - 7.135e5 * (1.0 - t**2)
            if (t < 0.20):
               c = 7.41e-5 * t**(-0.885)
            else:
               c = math.log(max(1.0e-20,-a / (self.oppk - a))) / (self.rzp - rpk)

            b = (self.oppk - a) * math.exp(-c * rpk)

         arg = c * r
         arg = max(-60.0,min(60.0,arg))

         wflt = a + b * math.exp(arg)

         if (t < 0.10):
            self.opr = wflt
         else:
            self.opr = (wflt * (0.95 - t) + self.opr * (t - 0.10)) / 0.85

      self.ttold_1 = t

   def wfvrmt(self,t,r):

      rpk = self.prad

      if (t >= 0.150):
         rpk = rpk * 1.0e-5
         r   = r * 1.0e-5

         if (t != self.ttold_3):
            rz     = self.rzv
            rmn    = rz - 9.31e3 * t**0.12115
            rz     = rz * 1.0e-5
            rmn    = rmn * 1.0e-5
            vmn    = 650.0 / (t * math.sqrt(t)) - 0.5 * self.vpk
            rneg   = rz - rmn
            rpls   = rpk - rz
            rnp    = rpk - rmn
            a1n    = self.vpk / rpls
            b1n    = self.vpk - a1n * rpk
            ovmn1n = a1n * rmn + b1n
            fmlt   = abs(vmn / self.vpk)
            ovmhy  = vmn + fmlt * (ovmn1n - vmn)
            alpha  = (rnp / (ovmhy - self.vpk) + rpls / self.vpk) / rneg
            beta   = -rpls * (alpha + 1.00 / self.vpk)
            fngz   = rnp / rpls
            dnom   = alpha * rnp + beta
            bcrmn  = 1.00 - vmn / (rnp / dnom + self.vpk)
            crmn1b = math.log(bcrmn)
            cgz1   = (beta * (1.00 / bcrmn - 1.00) / dnom) / \
                     ((fngz * crmn1b * rnp**(fngz - 1.00)) * \
                     (rnp + self.vpk * dnom))
            cgz    = math.exp(cgz1)
            bgz1   = crmn1b / cgz**(rnp**fngz)
            bgz    = math.exp(bgz1)

         rbr   = rpk - r
         gr    = 1.00 - bgz**(cgz**(rbr**fngz))

         if (r > rz):
            gr = (rpk - r) / rpls * gr + (r - rz) / rpls

         hr      = rbr / (alpha * rbr + beta) + self.vpk
         self.vr = gr * hr
         rpk     = rpk * 1.0e5
         r       = r * 1.0e5

      if (t <= 0.20):
         rtio = r / self.prad
         wflt = self.vpk * rtio * math.sqrt(rtio)

         if (t < 0.150):
            self.vr = wflt
         else:
            self.vr = (wflt * (0.200 - t) + self.vr * (t - 0.150)) * 20.0

      self.ttold_3 = t

   def wfvzr(self,t):

      rwfvzr = 0.0

      if (t < 0.0):
         print('wfvzr t < 0 %15.5f' % (t))
         sys.exit()

      if (t > 0.0 and t < 0.263):
         rwfvzr = 2.5e4 * t**0.80
      else:
         rwfvzr = (1.0 - 0.08459 * t**(-1.340)) * (33897.0 * t + 8490.0)

      return (rwfvzr)

   def wfzr(self,t):

      rwfzr = 0.0

      if (t < 0.0):
         print('wfzr t < 0 %15.5f' % (t))
         sys.exit()

      rwfzr = (1.0 - 0.03291 * t**(-1.086)) * (33897.0 * t + 8490.0)

      return (rwfzr)
