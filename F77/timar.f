      double precision function timar(grm,hobm,yield)  
      include 'real8.h'
c   
c  this routine computes the time of arrival of a nuclear   
c  blast given yield, ground range and height of burst of leaker
c   
      data third/0.3333333333333/   
c   
      gr  = grm * 3.2808 / 1000.0d0
      hob = hobm * 3.2808 / 1000.0d0

      r = sqrt(gr*gr+hob*hob)   
      r2 = r*r  
      r3 = r2*r 
      top = 0.54291-21.1858*r+361.81*r2+2383.0*r3   
      bot = 1.00 + 2.04797*r+ 2.68717*r2
      u = top/bot   
      xm=170.0*hob/(1.0+337.0*hob**0.25)+0.914*hob**2.5 
      if(gr.ge.xm) then 
         top=1.0858-33.629*r+455.84*r2+2383*r3  
         bot=1.5874+2.58*r+2.6872*r2
         tim=u*xm/gr+top/bot*(1.0-xm/gr)
      else  
         tim=u  
      endif 
      timar=tim*yield**third  
ccc      timar=tim*0.001*yield**third  
      return
      end   
