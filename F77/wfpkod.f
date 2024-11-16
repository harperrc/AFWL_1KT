      function wfpkod(dummy)

      include 'real8.h'

      save

      include "wfrt.h"

      data rhoz  /1.225d-3/,
     *     opz   /1.01325d6/,
     *     gamm1 /0.404574d0/

      op    = oppk
      rtio  = op / opz
      p     = op + opz
      gmone = gamm1
      gamma = 1.0d0 + gamm1
      gamra = gamma

c  iterate to solution

      do 2 n=1,20
      rho1 = rhoz * ((2.0d0 * gamra + (gamra + 1.0d0) * rtio) /
     *       (2.0d0 * gamra + (gamra - 1.0d0) * rtio))

      ee = p / (gmone * rho1)

      call air(ee,rho1,gmone,dum1,dum2)

      gamrao = gamra

      gamra  = 2.0d0 * gmone / (2.5d0 * gmone + 1.0d0) + 1.0d0

      if (abs(gamra - gamrao).lt.1.0d-5) goto 3

 2    continue

 3    wfpkod = rho1 - rhoz

      return
      end
