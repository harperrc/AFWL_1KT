      subroutine shock(yield,height,tim,ra,opr,odr,vr)

      include 'real8.h'

      alt = height

      call scalkt(alt,yield,vs1,ds1,ts1,cs1,ps1)

      t = tim * ts1
      r = ra / cs1

      call peak(t,r,prado,oppko,odpko,opro,odro,vpko,vro)

      opr = opro * ps1
      odr = odro * ds1
      vr  = vro * vs1

      return
      end
