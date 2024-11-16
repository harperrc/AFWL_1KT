      subroutine scalkt(hfpt,wb,vscale,dscale,tscale,cscale,pscale)

      include 'real8.h'

c  input:
c     hfpt   height of field point (m)
c     wb     yield of burst (Kt)

c  output:
c     vscale scale 1Kt sea level velocity to actual
c     dscale scale 1Kt sea level densities to actual
c     tscale scale actuval time to 1Kt at sea level 
c     cscale scale 1Kt sea level dimensions to actual
c     pscale scale 1Kt sea level pressures to actual

      save

c  0 (m)

      data p1 /1.01325d5/,c1/3.4029399d2/,r1/1.225d0/,t1/288.15d0/

      cubrt(x) = sign(exp(log(abs(x))/3.0d0),x)

      call matm62(hfpt,p3,c3,r3,t3)

      rp3    = cubrt(p1/p3)
      
      wsw    = cubrt(wb)

      tscale = c3 / (wsw * rp3 * c1)
      vscale = c3 / c1
      pscale = p3 / p1
      dscale = r3 / r1
      cscale = rp3 * wsw

      return
      end
