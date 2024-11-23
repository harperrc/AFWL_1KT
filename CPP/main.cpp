#include "AFWL.h"

double timar(double grft,
             double hobft,
             double yld)
{

//  return time of shock arrival in seconds

   double y13  = pow(yld,0.333333333333333);
   double x    = MAX(grft / y13,1.0e-9);
   double y    = MAX(hobft / y13,1.0e-9);
   double capr = sqrt(x * x + y * y);
   double r    = capr / 1000.0;
   double z    = y / x;

   if (z > 100.0)z = 100.0;

   double r2   = r * r;
   double r3   = r * r2;
   double r4   = r * r3;
   double r6   = r4 * r2;
   double r8   = pow(r,8.0);

//  time of shock arrival

   double xm   = 170.0 * y / (1.0 + 60.0 * pow(y,0.25)) + 2.89 * pow(y / 100.0,2.5);
   double top  = (0.543 - 21.8 * r + 386.0 * r2 + 2383.0 * r3) * r8;
   double bot  = 2.99e-14 - 1.91e-10 * r2 + 1.032e-6 * r4 -
                 4.43e-6 * r6 + (1.632 + 2.629 * r + 2.69 * r2) * r8;

   double tau  = top / bot;

   if (x < xm)return (0.001 * tau * y13);

//  if in mach region

   top         = (1.086 - 34.605 * r + 436.30 * r2 + 2383 * r3) * r8;
   bot         = 3.0137e-13 - 1.2128e-9 * r2 + 4.128e-6 * r4 -
                 1.116e-5 * r6 + (1.632 + 2.629 * r + 2.69 * r2) * r8;

   double w    = top / bot;

   tau = tau * xm / x + w * ( 1.0 - xm / x);

   return (0.001 * tau * y13);
}
int main(int argc,char *argv[])
{
   cout.setf(ios::fixed);

   AFWL afwl;

   double xodr;
   double xopr;
   double xvr;

   double yld    = 1.0;
   double grft   = 18.896;
   double hobft  = 0.0;

   double grm    = grft / 3.2808;
   double hobm   = hobft / 3.2808;

   double tau = timar(grft,hobft,yld);

   double tm = 5.0 * tau;
   double dt = tm / 500.0;

   ofstream outp("out");
   outp.setf(ios::fixed);

   double t = 0.0;

   while (t <= tm){
      afwl.shock(yld,hobm,t,grm,xopr,xodr,xvr);

      outp <<setw(20)<<setprecision(10)
           <<t<<" " 
           <<setw(20)<<setprecision(10)
           <<xopr * 0.000145038<<" "
           <<setw(20)<<setprecision(10)
           <<xodr<<" "
           <<setw(20)<<setprecision(10)
           <<xvr<<" "
           <<1<<" "
           <<endl;

      t += dt;
   }

   tm = 1.0;
   dt = (tm - t) / 5000.0;

   while (t <= tm){
      afwl.shock(yld,hobm,t,grm,xopr,xodr,xvr);

      outp <<setw(20)<<setprecision(10)
           <<t<<" " 
           <<setw(20)<<setprecision(10)
           <<xopr * 0.000145038<<" "
           <<setw(20)<<setprecision(10)
           <<xodr<<" "
           <<setw(20)<<setprecision(10)
           <<xvr<<" "
           <<2<<" "
           <<endl;

      t += dt;
   }

   outp.close();

   return (0);
}
