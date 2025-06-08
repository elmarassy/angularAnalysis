//
// Created by Mero Elmarassy on 6/6/25.
//

#ifndef CERN_TIMEDEPENDENTFIT_H
#define CERN_TIMEDEPENDENTFIT_H

#include "TMath.h"

#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooRealProxy.h"



class timeDependentB: public RooAbsPdf {
public:
    timeDependentB(const char *name, const char *title, RooAbsReal& costhetal, RooAbsReal& costhetak, RooAbsReal& phi,
               RooAbsReal& t, RooAbsReal& x, RooAbsReal& y, RooAbsReal& J1s, RooAbsReal& J1c, RooAbsReal& J2s, RooAbsReal& J2c,
               RooAbsReal& J3, RooAbsReal& J4, RooAbsReal& J5, RooAbsReal& J6s, RooAbsReal& J7, RooAbsReal& J8,
               RooAbsReal& J9, RooAbsReal &J1sBar, RooAbsReal &J1cBar,
                               RooAbsReal &J2sBar, RooAbsReal &J2cBar, RooAbsReal &J3Bar, RooAbsReal &J4Bar, RooAbsReal &J5Bar,
                               RooAbsReal &J6sBar, RooAbsReal &J7Bar, RooAbsReal &J8Bar, RooAbsReal &J9Bar, RooAbsReal& h1s, RooAbsReal& h1c, RooAbsReal& h2s, RooAbsReal& h2c, RooAbsReal& h3,
               RooAbsReal& h4, RooAbsReal& h5, RooAbsReal& h6, RooAbsReal& h7, RooAbsReal& h8, RooAbsReal& h9,
               RooAbsReal& s1s, RooAbsReal& s1c, RooAbsReal& s2s, RooAbsReal& s2c, RooAbsReal& s3,
               RooAbsReal& s4, RooAbsReal& s5, RooAbsReal& s6, RooAbsReal& s7, RooAbsReal& s8, RooAbsReal& s9);

    timeDependentB(timeDependentB const &other, const char *name=nullptr);

    TObject* clone(const char *newname) const override {
        return new timeDependentB(*this, newname);
    }

protected:
    RooRealProxy cosThetaL_;
    RooRealProxy cosThetaK_;
    RooRealProxy phi_;
    RooRealProxy t_;
    RooRealProxy x_;
    RooRealProxy y_;
    RooRealProxy J1s_;
    RooRealProxy J1c_;
    RooRealProxy J2s_;
    RooRealProxy J2c_;
    RooRealProxy J3_;
    RooRealProxy J4_;
    RooRealProxy J5_;
    RooRealProxy J6s_;
    RooRealProxy J7_;
    RooRealProxy J8_;
    RooRealProxy J9_;
    RooRealProxy J1sBar_;
    RooRealProxy J1cBar_;
    RooRealProxy J2sBar_;
    RooRealProxy J2cBar_;
    RooRealProxy J3Bar_;
    RooRealProxy J4Bar_;
    RooRealProxy J5Bar_;
    RooRealProxy J6sBar_;
    RooRealProxy J7Bar_;
    RooRealProxy J8Bar_;
    RooRealProxy J9Bar_;
    RooRealProxy h1s_;
    RooRealProxy h1c_;
    RooRealProxy h2s_;
    RooRealProxy h2c_;
    RooRealProxy h3_;
    RooRealProxy h4_;
    RooRealProxy h5_;
    RooRealProxy h6_;
    RooRealProxy h7_;
    RooRealProxy h8_;
    RooRealProxy h9_;
    RooRealProxy s1s_;
    RooRealProxy s1c_;
    RooRealProxy s2s_;
    RooRealProxy s2c_;
    RooRealProxy s3_;
    RooRealProxy s4_;
    RooRealProxy s5_;
    RooRealProxy s6_;
    RooRealProxy s7_;
    RooRealProxy s8_;
    RooRealProxy s9_;

    inline double evaluate_prob(double cosThetaL, double cosThetaK, double phi, double t, double x, double y,
                double J1s, double J1c, double J2s, double J2c, double J3, double J4, double J5, double J6s, double J7, double J8, double J9,
                double J1sBar, double J1cBar, double J2sBar, double J2cBar, double J3Bar, double J4Bar, double J5Bar, double J6sBar, double J7Bar, double J8Bar, double J9Bar,
                double h1s, double h1c, double h2s, double h2c, double h3, double h4, double h5, double h6, double h7, double h8, double h9,
                double s1s, double s1c, double s2s, double s2c, double s3, double s4, double s5, double s6, double s7, double s8, double s9) const;
public:
    inline double evaluate() const override;

    inline void doEval(RooFit::EvalContext &ctx) const override;

     int getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char */*rangeName*/) const override {
         if (matchArgs(allVars, analVars, cosThetaL_, cosThetaK_, phi_, t_)) return 1 ;
         return 0 ;
     }

     double analyticalIntegral(int code, const char *rangeName) const override {
         if (code==1) {
             double x = x_;
             double y = y_;
             double J1s = J1s_;
             double J1c = J1c_;
             double J2s = J2s_;
             double J2c = J2c_;
             double J1sBar = J1sBar_;
             double J1cBar = J1cBar_;
             double J2sBar = J2sBar_;
             double J2cBar = J2cBar_;
             double h1s = h1s_;
             double h1c = h1c_;
             double h2s = h2s_;
             double h2c = h2c_;
             double s1s = s1s_;
             double s1c = s1c_;
             double s2s = s2s_;
             double s2c = s2c_;

             return 1/(8*(1 + x*x)*(-1 + y*y)) * (
                     -6*J1c - 12*J1s + 2*J2c + 4*J2s + 3*s1c*x + 6*s1s*x - s2c*x - 2*s2s*x +
                     (-3*J1c - 3*J1cBar - 6*J1s - 6*J1sBar + J2c + J2cBar + 2*(J2s + J2sBar))*x*x +
                     (3*h1c + 6*h1s - h2c - 2*h2s)*(1 + x*x)*y +
                     (3*J1c - 3*J1cBar + 6*J1s - 6*J1sBar - J2c + J2cBar - 2*J2s + 2*J2sBar - 3*s1c*x - 6*s1s*x + s2c*x + 2*s2s*x)*y*y);




//             double factor = 1/((1+x*x)*(1-y*y));
//             double term1 = 12*J1s - 2*J2c - 4*J2s - (3*s1c + 6*s1s - s2c - 2*s2s)*x;
//             double term2 = (3*J1cBar + 6*J1s + 6*J1sBar - J2c - J2cBar - 2*(J2s + J2sBar)) * x*x;
//             double term3 = (3*J1cBar - 6*J1s + 6*J1sBar + J2c - J2cBar + 2*J2s - 2*J2sBar +(3*s1c + 6*s1s - s2c - 2*s2s)*x)*y*y;
//             double term4 = 3 * J1c * (2 + x*x - y*y);
//             double term5 = (3*h1c + 6*h1s - h2c - 2*h2s) * y / (y*y-1);
//             return (factor*(term1 + term2 + term3 + term4) + term5)/8;

         }
         return 0 ;
     }
};

class timeDependentBBar: public RooAbsPdf {
public:
    timeDependentBBar(const char *name, const char *title, RooAbsReal& costhetal, RooAbsReal& costhetak, RooAbsReal& phi, RooAbsReal& t,
                RooAbsReal& x, RooAbsReal& y, RooAbsReal& J1s, RooAbsReal& J1c, RooAbsReal& J2s, RooAbsReal& J2c,
                     RooAbsReal& J3, RooAbsReal& J4, RooAbsReal& J5, RooAbsReal& J6s, RooAbsReal& J7, RooAbsReal& J8,
                     RooAbsReal& J9,
                  RooAbsReal& J1sBar, RooAbsReal& J1cBar, RooAbsReal& J2sBar, RooAbsReal& J2cBar, RooAbsReal& J3Bar, RooAbsReal& J4Bar,
                  RooAbsReal& J5Bar, RooAbsReal& J6sBar, RooAbsReal& J7Bar, RooAbsReal& J8Bar, RooAbsReal& J9Bar,
                  RooAbsReal& h1s, RooAbsReal& h1c, RooAbsReal& h2s, RooAbsReal& h2c, RooAbsReal& h3,
                  RooAbsReal& h4, RooAbsReal& h5, RooAbsReal& h6, RooAbsReal& h7, RooAbsReal& h8, RooAbsReal& h9,
                  RooAbsReal& s1s, RooAbsReal& s1c, RooAbsReal& s2s, RooAbsReal& s2c, RooAbsReal& s3,
                  RooAbsReal& s4, RooAbsReal& s5, RooAbsReal& s6, RooAbsReal& s7, RooAbsReal& s8, RooAbsReal& s9);

    timeDependentBBar(timeDependentBBar const &other, const char *name=nullptr);

    TObject* clone(const char *newname) const override {
        return new timeDependentBBar(*this, newname);
    }

protected:
    RooRealProxy cosThetaL_;
    RooRealProxy cosThetaK_;
    RooRealProxy phi_;
    RooRealProxy t_;
    RooRealProxy x_;
    RooRealProxy y_;
    RooRealProxy J1s_;
    RooRealProxy J1c_;
    RooRealProxy J2s_;
    RooRealProxy J2c_;
    RooRealProxy J3_;
    RooRealProxy J4_;
    RooRealProxy J5_;
    RooRealProxy J6s_;
    RooRealProxy J7_;
    RooRealProxy J8_;
    RooRealProxy J9_;
    RooRealProxy J1sBar_;
    RooRealProxy J1cBar_;
    RooRealProxy J2sBar_;
    RooRealProxy J2cBar_;
    RooRealProxy J3Bar_;
    RooRealProxy J4Bar_;
    RooRealProxy J5Bar_;
    RooRealProxy J6sBar_;
    RooRealProxy J7Bar_;
    RooRealProxy J8Bar_;
    RooRealProxy J9Bar_;
    RooRealProxy h1s_;
    RooRealProxy h1c_;
    RooRealProxy h2s_;
    RooRealProxy h2c_;
    RooRealProxy h3_;
    RooRealProxy h4_;
    RooRealProxy h5_;
    RooRealProxy h6_;
    RooRealProxy h7_;
    RooRealProxy h8_;
    RooRealProxy h9_;
    RooRealProxy s1s_;
    RooRealProxy s1c_;
    RooRealProxy s2s_;
    RooRealProxy s2c_;
    RooRealProxy s3_;
    RooRealProxy s4_;
    RooRealProxy s5_;
    RooRealProxy s6_;
    RooRealProxy s7_;
    RooRealProxy s8_;
    RooRealProxy s9_;

    inline double evaluate_prob(double cosThetaL, double cosThetaK, double phi, double t, double x, double y,
                                double J1s, double J1c, double J2s, double J2c, double J3, double J4, double J5, double J6s, double J7, double J8, double J9,
                                double J1sBar, double J1cBar, double J2sBar, double J2cBar, double J3Bar, double J4Bar, double J5Bar, double J6sBar, double J7Bar, double J8Bar, double J9Bar,
                                double h1s, double h1c, double h2s, double h2c, double h3, double h4, double h5, double h6, double h7, double h8, double h9,
                                double s1s, double s1c, double s2s, double s2c, double s3, double s4, double s5, double s6, double s7, double s8, double s9) const;
public:
    inline double evaluate() const override;

    inline void doEval(RooFit::EvalContext &ctx) const override;

    double analyticalIntegral(int code, const char *rangeName) const override {
        if (code==1) {
            double x = x_;
            double y = y_;
            double J1s = J1s_;
            double J1c = J1c_;
            double J2s = J2s_;
            double J2c = J2c_;
            double J1sBar = J1sBar_;
            double J1cBar = J1cBar_;
            double J2sBar = J2sBar_;
            double J2cBar = J2cBar_;
            double h1s = h1s_;
            double h1c = h1c_;
            double h2s = h2s_;
            double h2c = h2c_;
            double s1s = s1s_;
            double s1c = s1c_;
            double s2s = s2s_;
            double s2c = s2c_;

            return 1/(8*(1 + x*x)*(-1 + y*y)) * (
                    -6*J1c - 12*J1s + 2*J2c + 4*J2s + 3*s1c*x + 6*s1s*x - s2c*x - 2*s2s*x +
                    (-3*J1c - 3*J1cBar - 6*J1s - 6*J1sBar + J2c + J2cBar + 2*(J2s + J2sBar))*x*x +
                    (3*h1c + 6*h1s - h2c - 2*h2s)*(1 + x*x)*y +
                    (3*J1c - 3*J1cBar + 6*J1s - 6*J1sBar - J2c + J2cBar - 2*J2s + 2*J2sBar - 3*s1c*x - 6*s1s*x + s2c*x + 2*s2s*x)*y*y);

        }
        return 0 ;
    }
};

class timeDependent: public RooAbsPdf {
public:
    timeDependent(const char *name, const char *title, RooAbsReal& costhetal, RooAbsReal& costhetak, RooAbsReal& phi,
                  RooAbsReal& t, RooAbsReal& x, RooAbsReal& y,
                  RooAbsReal& J1s, RooAbsReal& J1c, RooAbsReal& J2s, RooAbsReal& J2c, RooAbsReal& J3, RooAbsReal& J4,
                  RooAbsReal& J5, RooAbsReal& J6s, RooAbsReal& J7, RooAbsReal& J8, RooAbsReal& J9,
                  RooAbsReal& J1sBar, RooAbsReal& J1cBar, RooAbsReal& J2sBar, RooAbsReal& J2cBar, RooAbsReal& J3Bar, RooAbsReal& J4Bar,
                  RooAbsReal& J5Bar, RooAbsReal& J6sBar, RooAbsReal& J7Bar, RooAbsReal& J8Bar, RooAbsReal& J9Bar,
                  RooAbsReal& h1s, RooAbsReal& h1c, RooAbsReal& h2s, RooAbsReal& h2c, RooAbsReal& h3,
                  RooAbsReal& h4, RooAbsReal& h5, RooAbsReal& h6, RooAbsReal& h7, RooAbsReal& h8, RooAbsReal& h9,
                  RooAbsReal& s1s, RooAbsReal& s1c, RooAbsReal& s2s, RooAbsReal& s2c, RooAbsReal& s3,
                  RooAbsReal& s4, RooAbsReal& s5, RooAbsReal& s6, RooAbsReal& s7, RooAbsReal& s8, RooAbsReal& s9);

    timeDependent(timeDependent const &other, const char *name=nullptr);

    TObject* clone(const char *newname) const override {
        return new timeDependent(*this, newname);
    }

protected:
    RooRealProxy cosThetaL_;
    RooRealProxy cosThetaK_;
    RooRealProxy phi_;
    RooRealProxy t_;
    RooRealProxy x_;
    RooRealProxy y_;
    RooRealProxy J1s_;
    RooRealProxy J1c_;
    RooRealProxy J2s_;
    RooRealProxy J2c_;
    RooRealProxy J3_;
    RooRealProxy J4_;
    RooRealProxy J5_;
    RooRealProxy J6s_;
    RooRealProxy J7_;
    RooRealProxy J8_;
    RooRealProxy J9_;
    RooRealProxy J1sBar_;
    RooRealProxy J1cBar_;
    RooRealProxy J2sBar_;
    RooRealProxy J2cBar_;
    RooRealProxy J3Bar_;
    RooRealProxy J4Bar_;
    RooRealProxy J5Bar_;
    RooRealProxy J6sBar_;
    RooRealProxy J7Bar_;
    RooRealProxy J8Bar_;
    RooRealProxy J9Bar_;
    RooRealProxy h1s_;
    RooRealProxy h1c_;
    RooRealProxy h2s_;
    RooRealProxy h2c_;
    RooRealProxy h3_;
    RooRealProxy h4_;
    RooRealProxy h5_;
    RooRealProxy h6_;
    RooRealProxy h7_;
    RooRealProxy h8_;
    RooRealProxy h9_;
    RooRealProxy s1s_;
    RooRealProxy s1c_;
    RooRealProxy s2s_;
    RooRealProxy s2c_;
    RooRealProxy s3_;
    RooRealProxy s4_;
    RooRealProxy s5_;
    RooRealProxy s6_;
    RooRealProxy s7_;
    RooRealProxy s8_;
    RooRealProxy s9_;

    inline double evaluate_prob(double cosThetaL, double cosThetaK, double phi, double t, double x, double y,
                                double J1s, double J1c, double J2s, double J2c, double J3, double J4, double J5, double J6s, double J7, double J8, double J9,
                                double J1sBar, double J1cBar, double J2sBar, double J2cBar, double J3Bar, double J4Bar, double J5Bar, double J6sBar, double J7Bar, double J8Bar, double J9Bar,
                                double h1s, double h1c, double h2s, double h2c, double h3, double h4, double h5, double h6, double h7, double h8, double h9,
                                double s1s, double s1c, double s2s, double s2c, double s3, double s4, double s5, double s6, double s7, double s8, double s9) const;
public:
    inline double evaluate() const override;

    inline void doEval(RooFit::EvalContext &ctx) const override;

    // inline int getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char */*rangeName*/) const override {
    //     if (matchArgs(allVars, analVars, cosThetaL_, cosThetaK_, phi_, t_)) return 1 ;
    //     return 0 ;
    // }
    //
    // inline double analyticalIntegral(int code, const char *rangeName) const override {
    //     if (code==1) {
    //         return 1.0;
    //     }
    //     return 0 ;
    // }
};

void fitTimeDependent(int nruns, int ngen);
#endif //CERN_TIMEDEPENDENTFIT_H
