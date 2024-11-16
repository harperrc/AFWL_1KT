      subroutine terror(t,name)
      include 'real8.h'
      character*6 name
      write(6,2)t,name
      stop
 2    format(1x,'Time Error t = ',1pe13.5,' seconds in routine ',a6)
      end
