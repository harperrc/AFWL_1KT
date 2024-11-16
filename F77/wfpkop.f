      function wfpkop(r)

      include 'real8.h'

      save

      include "wfrt.h"

      data ac    /3.04d18/,
     *     aq    /1.13d14/,
     *     astar /7.9d9/,
     *     rstar /4.454d4/

      rr   = 1.0d0 / r

      rtio = 2.24517d-5 * r

      cf   = sqrt(log(rtio + 3.0d0 * exp(-(sqrt(rtio) / 3.0d0))))

      wfpkop = ((ac * rr + aq) * rr + astar / cf) * rr

      return
      end
