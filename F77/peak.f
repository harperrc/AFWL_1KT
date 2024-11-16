      subroutine peak(t,ra,prado,oppko,odpko,opro,odro,vpko,vro)

      include 'real8.h'

      save

      include "wfrt.h"

      data told   /0.0d0/
     *     psca   /0.1d0/,
     *     vsca   /0.01d0/,
     *     rhosca /1000.0d0/

      r = ra * 100.0d0

      if (t.ne.told) then
         rzp   = wfzr(t)
         rzd   = wfdzr(t)
         rzv   = wfvzr(t)
         prad  = wfpr(t)

         prado = prad * vsca
         oppk  = wfpkop(prad)
         oppko = oppk * psca
         odpk  = wfpkod(prad)
         odpko = odpk * rhosca
         vpk   = wfpkv(prad)
         vpko  = vpk * vsca
         told  = t
      endif

      write(16,10)t,ra,rzp,rzd,rzv,prad,oppk,odpk,vpk
 10   format(9(1x,1pe15.5))

      opr   = 0.0d0
      odr   = 0.0d0
      vr    = 0.0d0
      opro  = 0.0d0
      odro  = 0.0d0
      vro   = 0.0d0

      if (r.le.prad) then
         call wfprmt(t,r)
         opro = opr * psca

         call wfdrmt(t,r)

         call well(t,r,odw)
         odro = odw * rhosca

         call wfvrmt(t,r)
         vro = vr * vsca
      endif

      return
      end
