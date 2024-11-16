      subroutine wfprmt(t,r)

      include 'real8.h'

      save

      include "wfrt.h"

      data ttold /0.0d0/

      rpk = prad

      if (t.ge.0.1) then
         rpk = rpk * 1.0d-5
         r   = r * 1.0d-5

         if (t.ne.ttold) then
            rz     = rzp
            rmn    = rz - 9.7163d3 * t**0.12115d0
            rz     = rz * 1.0d-5
            rmn    = rmn * 1.0d-5
            opmn   = 2.2d4 / (t * sqrt(t)) - 0.50d0 * oppk
            rneg   = rz - rmn
            rp1s   = rpk - rz
            rnp    = rpk - rmn
            a1n    = oppk / rp1s
            b1n    = oppk - a1n * rpk
            opmn1n = a1n * rmn + b1n
            fmlt   = abs(opmn / oppk)
            opmhy  = opmn + fmlt * (opmn1n - opmn)
            alpha  = (rnp / (opmhy - oppk) + rp1s / oppk) / rneg
            beta   = -rp1s * (alpha + 1.0 / oppk)
            fngz   = rnp / rp1s
            dnom   = alpha * rnp + beta
            bcrmn  = 1.0 - opmn / (rnp / dnom + oppk)
            crmn1b = log(bcrmn)
            cgz1   = (beta * (1.0d0 / bcrmn - 1.0d0) / dnom) /
     *               ((fngz * crmn1b * rnp**(fngz - 1.0d0)) * 
     *               (rnp + oppk * dnom))
            cgz    = exp(cgz1)
            bgz1   = crmn1b / cgz**(rnp**fngz)
            bgz    = exp(bgz1)
         endif

         rbr   = rpk - r
         gr    = 1.0d0 - bgz**(cgz**(rbr**fngz))
         if (r.gt.rz)gr = (rpk - r) / rp1s * gr + (r - rz) / rp1s
         hr  = rbr / (alpha * rbr + beta) + oppk
         opr = gr * hr
         rpk = rpk * 1.0d5
         r   = r * 1.0d5
      endif

      if (t.lt.0.95d0) then
         if (t.ne.ttold) then
            a = 5.446d4 * t**(-1.22d0) - 7.135d5 * (1.0d0 - t**2)
            if (t.lt.0.2d0)c = 7.41d-5 * t**(-0.885d0)
            if (t.ge.0.2d0) then
               c = log(max(1.0d-20,-a / (oppk - a))) / (rzp - rpk)
            endif
            b = (oppk - a) * exp(-c * rpk)
         endif
         arg = c * r
         arg = max(-60.0d0,min(60.0d0,arg))
         wflt = a + b * exp(arg)
         if (t.lt.0.1d0) then
            opr = wflt
         else
            opr = (wflt * (0.95d0 - t) + opr * (t - 0.1d0)) / 0.85d0
         endif
      endif

      ttold = t

      return
      end
