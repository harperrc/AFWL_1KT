      subroutine wfvrmt(t,r)

      include 'real8.h'

      save

      include "wfrt.h"

      data ttold/0.0d0/

      rpk = prad

      if (t.ge.0.15d0) then
         rpk = rpk * 1.0d-5
         r   = r * 1.0d-5

         if (t.ne.ttold) then
            rz     = rzv
            rmn    = rz - 9.31d3 * t**0.12115d0
            rz     = rz * 1.0d-5
            rmn    = rmn * 1.0d-5
            vmn    = 650.0d0 / (t * sqrt(t)) - 0.5d0 * vpk
            rneg   = rz - rmn
            rp1s   = rpk - rz
            rnp    = rpk - rmn
            a1n    = vpk / rp1s
            b1n    = vpk - a1n * rpk
            ovmn1n = a1n * rmn + b1n
            fmlt   = abs(vmn / vpk)
            ovmhy  = vmn + fmlt * (ovmn1n - vmn)
            alpha  = (rnp / (ovmhy - vpk) + rp1s / vpk) / rneg
            beta   = -rp1s * (alpha + 1.0d0 / vpk)
            fngz   = rnp / rp1s
            dnom   = alpha * rnp + beta
            bcrmn  = 1.0d0 - vmn / (rnp / dnom + vpk)
            crmn1b = log(bcrmn)
            cgz1   = (beta * (1.0d0 / bcrmn - 1.0d0) / dnom) /
     *               ((fngz * crmn1b * rnp**(fngz - 1.0d0)) *
     *               (rnp + vpk * dnom))
            cgz    = exp(cgz1)
            bgz1   = crmn1b / cgz**(rnp**fngz)
            bgz    = exp(bgz1)
         endif

         rbr   = rpk - r
         gr    = 1.0d0 - bgz**(cgz**(rbr**fngz))
         if (r.gt.rz)gr = (rpk - r) / rp1s * 
     *               gr + (r - rz) / rp1s
         hr  = rbr / (alpha * rbr + beta) + vpk
         vr  = gr * hr
         rpk = rpk * 1.0d5
         r   = r* 1.0d5
      endif

      if (t.le.0.2d0) then
         rtio = r / prad
         wflt = vpk * rtio * sqrt(rtio)
         if (t.lt.0.15d0) then
            vr = wflt
         else
            vr = (wflt * (0.20d0 - t) + vr * (t - 0.15d0)) * 20.0d0
         endif
      endif

      ttold = t

      return
      end
