#include "AFWL.h"

AFWL::AFWL()
{
   rhoz    = 1.225e-3;
   opz     = 1.01325e6;
   gamm1   = 0.404574;

// wfpkod

   told    =  0.0;
   psca    = 0.10;
   vsca    = 0.010;
   rhosca  = 1000.0;

//  scalkt

   p1      = opz / 10.0;
   c1      = 3.4029399e2;
   r1      = rhoz * 1000.0;
   t1      = 288.15;

   tscale  = 0.0;
   vscale  = 0.0;
   pscale  = 0.0;
   dscale  = 0.0;
   cscale  = 0.0;

//  wfrt

   prad    = 0.0;
   oppk    = 0.0;
   odpk    = 0.0;
   vpk     = 0.0;
   opr     = 0.0;
   odr     = 0.0;
   vr      = 0.0;
   rzp     = 0.0;
   rzd     = 0.0;
   rzv     = 0.0;
   opmn    = 0.0;
   odmn    = 0.0;
   vmn     = 0.0;

//  wfprmt

   ttold_1 =  0.0;

//  wfvrmt

   ttold_3 =  0.0;

//  wfdrmt

   ttold_2 =  0.0;
}
AFWL::~AFWL()
{
}
double AFWL::air(double eee,
                 double rrr)
{

//  modified to only compute gmone, p and temp never used in 1KT standard

   double e     = eee * 1.0e-10;
   double rholn = log(773.39520495 * rrr);

// original code had check if (abs(e1).lt.5.0)goto 20
// e1 is ALWAYS < 5 since its a fixed value
// so we dropped the 2 if checks out and the calculation of e1

   double fn  = 0.0;
   double ws  = (8.5 + 0.15504314 * rholn - e) * exp(-0.05 * rholn + 0.02531780798);
          ws  = MAX(-60.0,MIN(60.0,ws));
          ws  = 1.0 / (exp(-ws) + 1.0);
   double fo  = exp(-0.22421524664 * e) * ws;
   double fon = exp(-0.15082856259 * e) * (1.0 - ws);

   double beta = 0.0;
   if (e > 1.0){
      beta = (6.9487e-3 * ws + 1.38974e-2) * log(e);
   }

   double e2  = (e - 40.0) * 0.33333333333333333;

   if (abs(e2) < 5.0){
      ws = (e - exp(0.0157 * rholn + 3.806662489)) * exp(-0.085 * rholn - 1.38629436);
      ws = 1.0 / (exp(-ws) + 1.0);
      fn  = exp(-0.039215686275 * e) * ws;
   }else if (e2 > 0.0){
      fn  = exp(-0.039215686275 * e);
      ws  = 1.0;
   }else{
      fn = 0.0;
      ws = 0.0;
   }

   beta = beta + ws * (0.045 - beta);
   ws   = (e - 160.0) / MAX(30.0 + 3.474356 * rholn,6.0);
   double fe   = 0.0;

   if (ws > -5.0){
      fe = 1.0 / (exp(-ws) + 1.0);
   }

   double gm   = (0.161 + 0.225 * fo + 0.280 * fon +  
                  0.137 * fn + 0.050 * fe) *  
                  exp(beta * rholn);

   return (gm);
}
double AFWL::cubrt(double x)
{
   return (SIGN(exp(log(abs(x)) / 3.0),x));
}
void AFWL::peak(double t,
                double ra)
{
   double r = ra * 100.0;

   if (t != told){
      rzp   = wfzr(t);
      rzd   = wfdzr(t);
      rzv   = wfvzr(t);
      prad  = wfpr(t);
      prado = prad * vsca;
      oppk  = wfpkop(prad);
      oppko = oppk * psca;
      odpk  = wfpkod(prad);
      odpko = odpk * rhosca;
      vpk   = wfpkv();
      vpko  = vpk * vsca;
      told  = t;
   }

   opr  = 0.0;
   odr  = 0.0;
   vr   = 0.0;
   opro = 0.0;
   odro = 0.0;
   vro  = 0.0;

   if (r <= prad){
      wfprmt(t,r);
      opro = opr * psca;

      wfdrmt(t,r);

      double odw = well(t,r);
      odro = odw * rhosca;

      wfvrmt(t,r);
      vro = vr * vsca;
   }
}
void AFWL::scalkt(double hfpt,
                  double wb)
{
   double p3 = 101325.000000000000;
   double c3 =    340.293990543471;  
   double r3 =      1.225000000000;

   double rp3 = cubrt(p1 / p3);
   double wsw = cubrt(wb);

   vscale     = c3 / c1;
   dscale     = r3 / r1;
   tscale     = c3 / (wsw * rp3 * c1);
   cscale     = rp3 / wsw;
   pscale     = p3 / p1;
}
void AFWL::shock(double yld,
                 double height,
                 double tim,
                 double ra,
                 double &xopr,
                 double &xodr,
                 double &xvr)
{
   double alt = height;

   scalkt(alt,yld);

   double t = tim * tscale;
   double r = ra / cscale;

   peak(t,r);

   xopr = opro * pscale;
   xodr = odro * dscale;
   xvr  = vro * vscale;
}
double AFWL::well(double t,
                  double r)
{
   double depth = odr;

   double rmsw = MIN(0.7 * rzd,1.55e4);

   if (t <= 1.2 && r <= rmsw){
      depth = MAX(-1.21e-3,-1.5e-3 * exp(-8.0e-13 * r * r * r));

      if (t > 1.0){
         depth = (1.2 - t) * depth * 5.0;
      }

      double trm = 0.90 * rmsw;

      if (r >= trm){
         depth = 10.0 * (depth * (rmsw - r) + odr * (r - trm)) / rmsw;
      }
   }

   return (depth);
}
void AFWL::wfdrmt(double t,
                  double r)
{
   double a;
   double a1n;
   double alpha;
   double b;
   double b1n;
   double bcrmn;
   double beta;
   double bgz;
   double bgz1;
   double c;
   double cgz;
   double cgz1;
   double crmn1b;
   double dnom;
   double fmlt;
   double fngz;
   double gr;
   double hr;
   double odmhy;
   double odmn1n;
   double rbr;
   double rmn;
   double rneg;
   double rnp;
   double rpk;
   double rpls;
   double rz;
   double wflt;

   rpk = prad;

   if (t >= 0.20){
      rpk = rpk * 1.0e-5;
      r   = r * 1.0e-5;

      if (t != ttold_2){
         rz         = rzd;
         rmn        = rz - 9.7163e3 * pow(t,0.12115);
         rz         = rz * 1.0e-5;
         rmn        = rmn * 1.0e-5;
         odmn  = -0.50 * odpk + 2.2e-5 * pow(t,-1.6026);
         rneg       = rz - rmn;
         rpls       = rpk - rz;
         rnp        = rpk - rmn;
         a1n        = odpk / rpls;
         b1n        = odpk - a1n * rpk;
         odmn1n     = a1n * rmn + b1n;
         fmlt       = abs(odmn / odpk);
         odmhy      = odmn + fmlt * (odmn1n - odmn);
         alpha      = (rnp / (odmhy - odpk) + rpls / odpk) / rneg;
         beta       = -rpls * (alpha + 1.0 / odpk);
         fngz       = rnp / rpls;
         dnom       = alpha * rnp + beta;
         bcrmn      = 1.0 - odmn / (rnp / dnom + odpk);
         crmn1b     = log(bcrmn);
         cgz1       = (beta * (1.0 / bcrmn - 1.0) / dnom) / 
                      ((fngz * crmn1b * pow(rnp,fngz - 1.0)) *  
                      (rnp + odpk * dnom));
         cgz        = exp(cgz1);
         bgz1       = crmn1b / pow(cgz,pow(rnp,fngz));
         bgz        = exp(bgz1);
      }
      
      rbr  = rpk - r;
      gr   = 1.0 - pow(bgz,pow(cgz,pow(rbr,fngz)));

      if (r > rz) {
         gr = (rpk - r) / rpls * gr + (r - rz) / rpls;
      }

      hr       = rbr / (alpha * rbr + beta) + odpk;
      odr = gr * hr;
      rpk      = rpk * 1.0e5;
      r        = r * 1.0e5;
   }

   if (t <= 1.0){
      if (t != ttold_2){
         a = -1.2e-3;
         c = log(MAX(1.0e-20,-a / (odpk - a))) / (rzd - rpk);
         b = (odpk - a) * exp(-c * rpk);
      }

      wflt = a + b * exp(c * r);

      if (t <= 0.20) {
         odr = wflt;
      }else{
         odr = (wflt * (1.0 - t) + odr * (t - 0.20)) * 1.25;
      }
   }

   ttold_2 = t;
}
double AFWL::wfdzr(double t)
{
   double val = 0.0;

   if (t < 0.0) {
      cout <<"wfdzr t < 0 "<<t<<endl;
      exit(-1);
   }

   if (t > 0.0 && t < 0.265){
      val = 2.568e4 * pow(t,0.395);
   }else{
      val = (1.0 - 0.03499 * pow(t,-1.068)) * (33897.0 * t + 8490.0) + 500.0;
   }

   return (val);
}
double AFWL::wfpkod(double r)
{
   int loop;

   double ee;
   double gamma;
   double gamra;
   double gamrao;
   double gmone;
   double op;
   double p;
   double rtio;
   double rho1;

   op    = oppk;
   rtio  = op / opz;
   p     = op + opz;
   gmone = gamm1;
   gamma = 1.0 + gamm1;
   gamra = gamma;

// iterate to solution

   loop = 1;
   while (loop){

      rho1 = rhoz * ((2.0 * gamra + (gamra + 1.0) * rtio) / 
                          (2.0 * gamra + (gamra - 1.0) * rtio));

      ee = p / (gmone * rho1);

      gmone = air(ee,rho1);

      gamrao = gamra;

      gamra  = 2.0 * gmone / (2.5 * gmone + 1.0) + 1.0;

      if (abs(gamra - gamrao) < 1.0e-5){
         return (rho1 - rhoz);
      }

      loop = loop + 1;

      if (loop > 20) {
         break;
      }
   }

   return (rho1 - rhoz);
}
double AFWL::wfpkop(double r)
{
   double val  = 0.0;

   double rr   = 1.0 / r;

   double rtio = 2.24517e-5 * r;
   double cf   = sqrt(log(rtio + 3.0 * exp(-(sqrt(rtio) / 3.0))));

   val = ((3.04e18 * rr + 1.13e14) * rr + 7.9e9 / cf) * rr;

   return (val);
}
double AFWL::wfpkv(void)
{
   double val = sqrt(oppk * odpk / (rhoz * (rhoz + odpk)));

   return (val);
}
void AFWL::wfprmt(double t,
                  double r)
{
   double a;
   double a1n;
   double alpha;
   double arg;
   double b;
   double b1n;
   double bcrmn;
   double beta;
   double bgz;
   double bgz1;
   double c;
   double cgz;
   double cgz1;
   double crmn1b;
   double dnom;
   double fmlt;
   double fngz;
   double gr;
   double hr;
   double opmn;
   double opmhy;
   double opmn1n;
   double rbr;
   double rmn;
   double rneg;
   double rnp;
   double rpk;
   double rpls;
   double rz;
   double wflt;

   rpk = prad;

   if (t >= 0.10){
      rpk = rpk * 1.0e-5;
      r   = r * 1.0e-5;

      if (t != ttold_1){
         rz        = rzp;
         rmn       = rz - 9.7163e3 * pow(t,0.12115);
         rz        = rz * 1.0e-5;
         rmn       = rmn * 1.0e-5;
         opmn      = 2.2e4 / (t * sqrt(t)) - 0.50 * oppk;
         rneg      = rz - rmn;
         rpls      = rpk - rz;
         rnp       = rpk - rmn;
         a1n       = oppk / rpls;
         b1n       = oppk - a1n * rpk;
         opmn1n    = a1n * rmn + b1n;
         fmlt      = abs(opmn / oppk);
         opmhy     = opmn + fmlt * (opmn1n - opmn);
         alpha     = (rnp / (opmhy - oppk) + rpls / oppk) / rneg;
         beta      = -rpls * (alpha + 1.0 / oppk);
         fngz      = rnp / rpls;
         dnom      = alpha * rnp + beta;
         bcrmn     = 1.0 - opmn / (rnp / dnom + oppk);
         crmn1b    = log(bcrmn);
         cgz1      = (beta * (1.00 / bcrmn - 1.00) / dnom) / 
                      ((fngz * crmn1b * pow(rnp,fngz - 1.00)) * 
                      (rnp + oppk * dnom));
         cgz       = exp(cgz1);
         bgz1      = crmn1b / pow(cgz,pow(rnp,fngz));
         bgz       = exp(bgz1);
      }

      rbr   = rpk - r;
      gr    = 1.00 - pow(bgz,pow(cgz,pow(rbr,fngz)));

      if (r > rz){
         gr = (rpk - r) / rpls * gr + (r - rz) / rpls;
      }

      hr       = rbr / (alpha * rbr + beta) + oppk;
      opr      = gr * hr;
      rpk      = rpk * 1.0e5;
      r        = r * 1.0e5;
   }

   if (t < 0.95) {
      if (t != ttold_1){
         a = 5.446e4 * pow(t,-1.22) - 7.135e5 * (1.0 - pow(t,2));

         if (t < 0.20){
            c = 7.41e-5 * pow(t,-0.885);
         }else{
            c = log(max(1.0e-20,-a / (oppk - a))) / (rzp - rpk);
         }

         b = (oppk - a) * exp(-c * rpk);
      }

      arg = c * r;
      arg = MAX(-60.0,MIN(60.0,arg));

      wflt = a + b * exp(arg);

      if (t < 0.10){
         opr = wflt;
      }else{
         opr = (wflt * (0.95 - t) + opr * (t - 0.10)) / 0.85;
      }
   }

   ttold_1 = t;
}
double AFWL::wfpr(double t)
{
   double rwfpr;
   double xr1;
   double xr2;

   xr1 = 24210.0 * pow(t,0.371) * (1.0 + (1.23 * t + 0.123) * (1.0 - exp(-26.25 * pow(t,0.79))));

   xr2 = (1.0 - 0.03291 * pow(t,-1.086)) * (33897.0 * t + 8490.0) + 
         8.36e3 + 2.5e3 * 
         log(t) + 800.0 * pow(t,-0.21);

   if (t < 0.21){
      rwfpr = xr1;
   }else if (t > 0.28){
      rwfpr = xr2;
   }else{
      rwfpr = (xr2 * (t - 0.21) + xr1 * (0.28 - t)) / 0.07;
   }

   return (rwfpr);
}
void AFWL::wfvrmt(double t,
                  double r)
{
   double a1n;
   double alpha;
   double b1n;
   double bcrmn;
   double beta;
   double bgz;
   double bgz1;
   double cgz;
   double cgz1;
   double crmn1b;
   double dnom;
   double fmlt;
   double fngz;
   double gr;
   double hr;
   double ovmhy;
   double ovmn1n;
   double rbr;
   double rmn;
   double rneg;
   double rnp;
   double rpk;
   double rpls;
   double rtio;
   double rz;
   double wflt;

   rpk = prad;

   if (t >= 0.150){
      rpk = rpk * 1.0e-5;
      r   = r * 1.0e-5;

      if (t != ttold_3){
         rz     = rzv;
         rmn    = rz - 9.31e3 * pow(t,0.12115);
         rz     = rz * 1.0e-5;
         rmn    = rmn * 1.0e-5;
         vmn    = 650.0 / (t * sqrt(t)) - 0.5 * vpk;
         rneg   = rz - rmn;
         rpls   = rpk - rz;
         rnp    = rpk - rmn;
         a1n    = vpk / rpls;
         b1n    = vpk - a1n * rpk;
         ovmn1n = a1n * rmn + b1n;
         fmlt   = abs(vmn / vpk);
         ovmhy  = vmn + fmlt * (ovmn1n - vmn);
         alpha  = (rnp / (ovmhy - vpk) + rpls / vpk) / rneg;
         beta   = -rpls * (alpha + 1.00 / vpk);
         fngz   = rnp / rpls;
         dnom   = alpha * rnp + beta;
         bcrmn  = 1.00 - vmn / (rnp / dnom + vpk);
         crmn1b = log(bcrmn);
         cgz1   = (beta * (1.00 / bcrmn - 1.00) / dnom) / 
                  ((fngz * crmn1b * pow(rnp,fngz - 1.00)) * 
                  (rnp + vpk * dnom));
         cgz    = exp(cgz1);
         bgz1   = crmn1b / pow(cgz,pow(rnp,fngz));
         bgz    = exp(bgz1);
      }

      rbr   = rpk - r;
      gr    = 1.00 - pow(bgz,pow(cgz,pow(rbr,fngz)));

      if (r > rz){
         gr = (rpk - r) / rpls * gr + (r - rz) / rpls;
      }

      hr      = rbr / (alpha * rbr + beta) + vpk;
      vr      = gr * hr;
      rpk     = rpk * 1.0e5;
      r       = r * 1.0e5;
   }

   if (t <= 0.20){
      rtio = r / prad;
      wflt = vpk * rtio * sqrt(rtio);

      if (t < 0.150){
         vr = wflt;
      }else{
         vr = (wflt * (0.200 - t) + vr * (t - 0.150)) * 20.0;
      }
   }

   ttold_3 = t;
}
double AFWL::wfvzr(double t)
{
   if (t < 0.0){
      cout <<"wfvzr t < 0 "<<t<<endl;
      exit(-1);
   }

   double val = 0.0;

   if (t > 0.0 && t < 0.263){
      val = 2.5e4 * pow(t,0.80);
   }else{
      val = (1.0 - 0.09459 * pow(t,-1.340)) * (33987.0 * t + 8490.0);
   }

   return (val);
}
double AFWL::wfzr(double t)
{
   if (t < 0.0){
      cout <<"wfzr t < 0 "<<t<<endl;
      exit(-1);
   }

   double val = (1.0 - 0.03291 * pow(t,-1.086)) * (33897.0 * t + 8490.0);

   return (val);
}
