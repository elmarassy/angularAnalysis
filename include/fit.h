//
// Created by Mero Elmarassy on 6/5/25.
//

#ifndef CERN_FIT_H
#define CERN_FIT_H

#include "TMath.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooMCStudy.h"
#include "RooPlot.h"

class timeIndependent: public RooAbsPdf {
public:
    timeIndependent(const char *name, const char *title, RooAbsReal& costhetal, RooAbsReal& costhetak, RooAbsReal& phi,
                    RooAbsReal& J1s, RooAbsReal& J1c, RooAbsReal& J2s, RooAbsReal& J2c, RooAbsReal& J3, RooAbsReal& J4,
                    RooAbsReal& J5, RooAbsReal& J6, RooAbsReal& J7, RooAbsReal& J8, RooAbsReal& J9);

    timeIndependent(timeIndependent const &other, const char *name=nullptr);

    TObject* clone(const char *newname) const override {
        return new timeIndependent(*this, newname);
    }

protected:
    RooRealProxy cosThetaL_;
    RooRealProxy cosThetaK_;
    RooRealProxy phi_;
    RooRealProxy J1s_;
    RooRealProxy J1c_;
    RooRealProxy J2s_;
    RooRealProxy J2c_;
    RooRealProxy J3_;
    RooRealProxy J4_;
    RooRealProxy J5_;
    RooRealProxy J6_;
    RooRealProxy J7_;
    RooRealProxy J8_;
    RooRealProxy J9_;

    inline double evaluate_prob(double cosThetaL, double cosThetaK, double phi, double J1s, double J1c, double J2s, double J2c, double J3, double J4, double J5, double J6, double J7, double J8, double J9) const;
public:
    inline double evaluate() const override;

    inline void doEval(RooFit::EvalContext &ctx) const override;

    inline int getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char */*rangeName*/) const override {
        if (matchArgs(allVars, analVars, cosThetaL_, cosThetaK_, phi_)) return 1 ;
        return 0 ;
    }

    inline double analyticalIntegral(int code, const char *rangeName) const override {
        if (code==1) {
            return 1.0;
        }
        return 0 ;
    }
};

class timeDependentUntagged: public RooAbsPdf {
public:
    timeDependentUntagged(const char *name, const char *title, RooAbsReal& costhetal, RooAbsReal& costhetak, RooAbsReal& phi,
                          RooAbsReal& t, RooAbsReal& gamma, RooAbsReal& y,
                          RooAbsReal& SJ1s, RooAbsReal& SJ1c, RooAbsReal& SJ2s, RooAbsReal& SJ2c, RooAbsReal& SJ3,
                          RooAbsReal& SJ4, RooAbsReal& DJ5, RooAbsReal& DJ6s, RooAbsReal& SJ7, RooAbsReal& DJ8, RooAbsReal& DJ9,
                          RooAbsReal& h1s, RooAbsReal& h1c, RooAbsReal& h2s, RooAbsReal& h2c, RooAbsReal& h3,
                          RooAbsReal& h4, RooAbsReal& h5, RooAbsReal& h6, RooAbsReal& h7, RooAbsReal& h8, RooAbsReal& h9);

    timeDependentUntagged(timeDependentUntagged const &other, const char *name=nullptr);

    TObject* clone(const char *newname) const override {
        return new timeDependentUntagged(*this, newname);
    }

protected:
    RooRealProxy cosThetaL_;
    RooRealProxy cosThetaK_;
    RooRealProxy phi_;
    RooRealProxy t_;
    RooRealProxy y_;
    RooRealProxy gamma_;
    RooRealProxy SJ1s_;
    RooRealProxy SJ1c_;
    RooRealProxy SJ2s_;
    RooRealProxy SJ2c_;
    RooRealProxy SJ3_;
    RooRealProxy SJ4_;
    RooRealProxy DJ5_;
    RooRealProxy DJ6s_;
    RooRealProxy SJ7_;
    RooRealProxy DJ8_;
    RooRealProxy DJ9_;
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

    inline double evaluate_prob(double cosThetaL, double cosThetaK, double phi, double t, double gamma, double y,
                                double SJ1s, double SJ1c, double SJ2s, double SJ2c, double SJ3, double SJ4, double DJ5, double DJ6s, double SJ7, double DJ8, double DJ9,
                                double h1s, double h1c, double h2s, double h2c, double h3, double h4, double h5, double h6, double h7, double h8, double h9) const;
public:
    inline double evaluate() const override;

    inline void doEval(RooFit::EvalContext &ctx) const override;

    inline int getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char */*rangeName*/) const override {
        if (matchArgs(allVars, analVars, cosThetaL_, cosThetaK_, phi_, t_)) return 1 ;
        return 0 ;
    }

    inline double analyticalIntegral(int code, const char *rangeName) const {
        if (code==1) {
            return 1.0;
        }
        return 0 ;
    }
};

class timeDependentTagged: public RooAbsPdf {
public:
    timeDependentTagged(const char *name, const char *title, RooAbsReal& costhetal, RooAbsReal& costhetak, RooAbsReal& phi,
                          RooAbsReal& t, RooAbsReal& gamma, RooAbsReal& x,
                          RooAbsReal& DJ1s, RooAbsReal& DJ1c, RooAbsReal& DJ2s, RooAbsReal& DJ2c, RooAbsReal& DJ3,
                          RooAbsReal& DJ4, RooAbsReal& SJ5, RooAbsReal& SJ6s, RooAbsReal& DJ7, RooAbsReal& SJ8, RooAbsReal& SJ9,
                          RooAbsReal& s1s, RooAbsReal& s1c, RooAbsReal& s2s, RooAbsReal& s2c, RooAbsReal& s3,
                          RooAbsReal& s4, RooAbsReal& s5, RooAbsReal& s6, RooAbsReal& s7, RooAbsReal& s8, RooAbsReal& s9);

    timeDependentTagged(timeDependentTagged const &other, const char *name=nullptr);

    TObject* clone(const char *newname) const override {
        return new timeDependentTagged(*this, newname);
    }

protected:
    RooRealProxy cosThetaL_;
    RooRealProxy cosThetaK_;
    RooRealProxy phi_;
    RooRealProxy t_;
    RooRealProxy x_;
    RooRealProxy gamma_;
    RooRealProxy DJ1s_;
    RooRealProxy DJ1c_;
    RooRealProxy DJ2s_;
    RooRealProxy DJ2c_;
    RooRealProxy DJ3_;
    RooRealProxy DJ4_;
    RooRealProxy SJ5_;
    RooRealProxy SJ6s_;
    RooRealProxy DJ7_;
    RooRealProxy SJ8_;
    RooRealProxy SJ9_;
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

    inline double evaluate_prob(double cosThetaL, double cosThetaK, double phi, double t, double gamma, double x,
                                double DJ1s, double DJ1c, double DJ2s, double DJ2c, double DJ3, double DJ4, double SJ5, double SJ6s, double DJ7, double SJ8, double SJ9,
                                double s1s, double s1c, double s2s, double s2c, double s3, double s4, double s5, double s6, double s7, double s8, double s9) const;
public:
    inline double evaluate() const override;

    inline void doEval(RooFit::EvalContext &ctx) const override;

    inline int getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char */*rangeName*/) const override {
        if (matchArgs(allVars, analVars, cosThetaL_, cosThetaK_, phi_, t_)) return 1 ;
        return 0 ;
    }

    inline double analyticalIntegral(int code, const char *rangeName) const {
        if (code==1) {
            return 1.0;
        }
        return 0 ;
    }
};

#endif //CERN_FIT_H