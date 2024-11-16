#ifndef AFWL_H
#define AFWL_H

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>

using namespace std;

template <class T>
inline T MAX(const T &a,
             const T &b)
{
   return ((a > b) ? a : b);
}

template <class T>
inline T MIN(const T &a,
             const T &b)
{
   return ((a > b) ? b : a);
}

template <class T>
inline T SIGN(const T &a,const T &b)
{
   return (b >= 0 ? (a>=0 ? a : -a) : (a>=0 ? -a : a));
}

class AFWL
{
   public:

      AFWL();

      ~AFWL();

      double air(double eee,
                 double rrr);

      double cubrt(double x);

      void peak(double t,
                double ra);

      void scalkt(double hfpt,
                  double wb);

      void shock(double yld,
                 double height,
                 double tim,
                 double ra,
                 double &xopr,
                 double &xodr,
                 double &xvr);

      double well(double t,
                  double r);

      void wfdrmt(double t,
                  double r);

      double wfdzr(double t);

      double wfpkod(double r);

      double wfpkop(double r);

      double wfpkv(void);

      void wfprmt(double t,
                  double r);

      void wfvrmt(double t,
                  double r);

      double wfpr(double t);

      double wfvzr(double t);

      double wfzr(double t);

   private:

      double c1;
      double cscale;
      double dscale;
      double gamm1;
      double odmn;
      double odpk;
      double odpko;
      double opmn;
      double oppk;
      double oppko;
      double odr;
      double odro;
      double opr;
      double opro;
      double opz;
      double p1;
      double prad;
      double prado;
      double psca;
      double pscale;
      double r1;
      double rhosca;
      double rhoz;
      double rzd;
      double rzv;
      double rzp;
      double t1;
      double told;
      double ttold_1;
      double ttold_2;
      double ttold_3;
      double tscale;
      double vmn;
      double vpk;
      double vpko;
      double vr;
      double vro;
      double vsca;
      double vscale;
};

#endif
