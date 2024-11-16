      function wfpkv(dummy)

      include 'real8.h'

      save

      include "wfrt.h"

      data rhoz /1.225d-3/

      wfpkv = sqrt(oppk * odpk / (rhoz * (rhoz + odpk)))

      return
      end
