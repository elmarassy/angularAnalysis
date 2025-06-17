//
// Created by Mero Elmarassy on 6/10/25.
//

#ifndef CERN_NEWTIMEDEPENDENTFIT_H
#define CERN_NEWTIMEDEPENDENTFIT_H


#include "TMath.h"

#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooRealProxy.h"


class timeDependent: public RooAbsPdf {
public:
    timeDependent(const char *name, const char *title,
                  RooAbsReal &cosThetaL, RooAbsReal &cosThetaK, RooAbsReal &phi, RooAbsReal &t,
                  int sign, RooAbsReal &x, RooAbsReal &y,
                  RooAbsReal &K1s, RooAbsReal &K1c, RooAbsReal &K2s, RooAbsReal &K2c, RooAbsReal &K3,
                  RooAbsReal &K4, RooAbsReal &K5, RooAbsReal &K6s, RooAbsReal &K7, RooAbsReal &K8, RooAbsReal &K9,
                  RooAbsReal &W1s, RooAbsReal &W1c, RooAbsReal &W2s, RooAbsReal &W2c, RooAbsReal &W3,
                  RooAbsReal &W4, RooAbsReal &W5, RooAbsReal &W6s, RooAbsReal &W7, RooAbsReal &W8, RooAbsReal &W9,
                  RooAbsReal &H1s, RooAbsReal &H1c, RooAbsReal &H2s, RooAbsReal &H2c, RooAbsReal &H3,
                  RooAbsReal &H4, RooAbsReal &H5, RooAbsReal &H6s, RooAbsReal &H7, RooAbsReal &H8, RooAbsReal &H9,
                  RooAbsReal &Z1s, RooAbsReal &Z1c, RooAbsReal &Z2s, RooAbsReal &Z2c, RooAbsReal &Z3,
                  RooAbsReal &Z4, RooAbsReal &Z5, RooAbsReal &Z6s, RooAbsReal &Z7, RooAbsReal &Z8, RooAbsReal &Z9);

    timeDependent(timeDependent const &other, const char *name = nullptr);

    TObject *clone(const char *newname) const override {
        return new timeDependent(*this, newname);
    }

protected:
    int sign_;
    RooRealProxy cosThetaL_;
    RooRealProxy cosThetaK_;
    RooRealProxy phi_;
    RooRealProxy t_;
    RooRealProxy x_;
    RooRealProxy y_;
    RooRealProxy K1s_;
    RooRealProxy K1c_;
    RooRealProxy K2s_;
    RooRealProxy K2c_;
    RooRealProxy K3_;
    RooRealProxy K4_;
    RooRealProxy K5_;
    RooRealProxy K6s_;
    RooRealProxy K7_;
    RooRealProxy K8_;
    RooRealProxy K9_;
    RooRealProxy W1s_;
    RooRealProxy W1c_;
    RooRealProxy W2s_;
    RooRealProxy W2c_;
    RooRealProxy W3_;
    RooRealProxy W4_;
    RooRealProxy W5_;
    RooRealProxy W6s_;
    RooRealProxy W7_;
    RooRealProxy W8_;
    RooRealProxy W9_;
    RooRealProxy H1s_;
    RooRealProxy H1c_;
    RooRealProxy H2s_;
    RooRealProxy H2c_;
    RooRealProxy H3_;
    RooRealProxy H4_;
    RooRealProxy H5_;
    RooRealProxy H6s_;
    RooRealProxy H7_;
    RooRealProxy H8_;
    RooRealProxy H9_;
    RooRealProxy Z1s_;
    RooRealProxy Z1c_;
    RooRealProxy Z2s_;
    RooRealProxy Z2c_;
    RooRealProxy Z3_;
    RooRealProxy Z4_;
    RooRealProxy Z5_;
    RooRealProxy Z6s_;
    RooRealProxy Z7_;
    RooRealProxy Z8_;
    RooRealProxy Z9_;

    inline double
    evaluate_prob(double cosThetaL, double cosThetaK, double phi, double t, int sign, double x, double y,
                  double K1s, double K1c, double K2s, double K2c, double K3, double K4, double K5, double K6s,
                  double K7, double K8, double K9,
                  double W1s, double W1c, double W2s, double W2c, double W3, double W4, double W5, double W6s,
                  double W7, double W8, double W9,
                  double H1s, double H1c, double H2s, double H2c, double H3, double H4, double H5, double H6s,
                  double H7, double H8, double H9,
                  double Z1s, double Z1c, double Z2s, double Z2c, double Z3, double Z4, double Z5, double Z6s,
                  double Z7, double Z8, double Z9) const;

public:
    inline double evaluate() const override {
        return evaluate_prob(cosThetaL_, cosThetaK_, phi_, t_, sign_, x_, y_,
                             K1s_, K1c_, K2s_, K2c_, K3_, K4_, K5_, K6s_, K7_, K8_, K9_,
                             W1s_, W1c_, W2s_, W2c_, W3_, W4_, W5_, W6s_, W7_, W8_, W9_,
                             H1s_, H1c_, H2s_, H2c_, H3_, H4_, H5_, H6s_, H7_, H8_, H9_,
                             Z1s_, Z1c_, Z2s_, Z2c_, Z3_, Z4_, Z5_, Z6s_, Z7_, Z8_, Z9_
        );
    }

    inline void doEval(RooFit::EvalContext &ctx) const override {
        std::span<const double> cosThetaLSpan = ctx.at(cosThetaL_);
        std::span<const double> cosThetaKSpan = ctx.at(cosThetaK_);
        std::span<const double> tSpan = ctx.at(t_);
        std::span<const double> phiSpan = ctx.at(phi_);
        std::span<const double> xSpan = ctx.at(x_);
        std::span<const double> ySpan = ctx.at(y_);

        std::span<const double> K1sSpan = ctx.at(K1s_);
        std::span<const double> K1cSpan = ctx.at(K1c_);
        std::span<const double> K2sSpan = ctx.at(K2s_);
        std::span<const double> K2cSpan = ctx.at(K2c_);
        std::span<const double> K3Span = ctx.at(K3_);
        std::span<const double> K4Span = ctx.at(K4_);
        std::span<const double> K5Span = ctx.at(K5_);
        std::span<const double> K6sSpan = ctx.at(K6s_);
        std::span<const double> K7Span = ctx.at(K7_);
        std::span<const double> K8Span = ctx.at(K8_);
        std::span<const double> K9Span = ctx.at(K9_);
        std::span<const double> W1sSpan = ctx.at(W1s_);
        std::span<const double> W1cSpan = ctx.at(W1c_);
        std::span<const double> W2sSpan = ctx.at(W2s_);
        std::span<const double> W2cSpan = ctx.at(W2c_);
        std::span<const double> W3Span = ctx.at(W3_);
        std::span<const double> W4Span = ctx.at(W4_);
        std::span<const double> W5Span = ctx.at(W5_);
        std::span<const double> W6sSpan = ctx.at(W6s_);
        std::span<const double> W7Span = ctx.at(W7_);
        std::span<const double> W8Span = ctx.at(W8_);
        std::span<const double> W9Span = ctx.at(W9_);
        std::span<const double> H1sSpan = ctx.at(H1s_);
        std::span<const double> H1cSpan = ctx.at(H1c_);
        std::span<const double> H2sSpan = ctx.at(H2s_);
        std::span<const double> H2cSpan = ctx.at(H2c_);
        std::span<const double> H3Span = ctx.at(H3_);
        std::span<const double> H4Span = ctx.at(H4_);
        std::span<const double> H5Span = ctx.at(H5_);
        std::span<const double> H6sSpan = ctx.at(H6s_);
        std::span<const double> H7Span = ctx.at(H7_);
        std::span<const double> H8Span = ctx.at(H8_);
        std::span<const double> H9Span = ctx.at(H9_);
        std::span<const double> Z1sSpan = ctx.at(Z1s_);
        std::span<const double> Z1cSpan = ctx.at(Z1c_);
        std::span<const double> Z2sSpan = ctx.at(Z2s_);
        std::span<const double> Z2cSpan = ctx.at(Z2c_);
        std::span<const double> Z3Span = ctx.at(Z3_);
        std::span<const double> Z4Span = ctx.at(Z4_);
        std::span<const double> Z5Span = ctx.at(Z5_);
        std::span<const double> Z6sSpan = ctx.at(Z6s_);
        std::span<const double> Z7Span = ctx.at(Z7_);
        std::span<const double> Z8Span = ctx.at(Z8_);
        std::span<const double> Z9Span = ctx.at(Z9_);

        std::size_t n = ctx.output().size();
        for (std::size_t i = 0; i < n; ++i) {
            ctx.output()[i] = evaluate_prob(cosThetaLSpan.size() > 1 ? cosThetaLSpan[i] : cosThetaLSpan[0],
                                            cosThetaKSpan.size() > 1 ? cosThetaKSpan[i] : cosThetaKSpan[0],
                                            phiSpan.size() > 1 ? phiSpan[i] : phiSpan[0],
                                            tSpan.size() > 1 ? tSpan[i] : tSpan[0],
                                            sign_,
                                            xSpan.size() > 1 ? xSpan[i] : xSpan[0],
                                            ySpan.size() > 1 ? ySpan[i] : ySpan[0],

                                            K1sSpan.size() > 1 ? K1sSpan[i] : K1sSpan[0],
                                            K1cSpan.size() > 1 ? K1cSpan[i] : K1cSpan[0],
                                            K2sSpan.size() > 1 ? K2sSpan[i] : K2sSpan[0],
                                            K2cSpan.size() > 1 ? K2cSpan[i] : K2cSpan[0],
                                            K3Span.size() > 1 ? K3Span[i] : K3Span[0],
                                            K4Span.size() > 1 ? K4Span[i] : K4Span[0],
                                            K5Span.size() > 1 ? K5Span[i] : K5Span[0],
                                            K6sSpan.size() > 1 ? K6sSpan[i] : K6sSpan[0],
                                            K7Span.size() > 1 ? K7Span[i] : K7Span[0],
                                            K8Span.size() > 1 ? K8Span[i] : K8Span[0],
                                            K9Span.size() > 1 ? K9Span[i] : K9Span[0],

                                            W1sSpan.size() > 1 ? W1sSpan[i] : W1sSpan[0],
                                            W1cSpan.size() > 1 ? W1cSpan[i] : W1cSpan[0],
                                            W2sSpan.size() > 1 ? W2sSpan[i] : W2sSpan[0],
                                            W2cSpan.size() > 1 ? W2cSpan[i] : W2cSpan[0],
                                            W3Span.size() > 1 ? W3Span[i] : W3Span[0],
                                            W4Span.size() > 1 ? W4Span[i] : W4Span[0],
                                            W5Span.size() > 1 ? W5Span[i] : W5Span[0],
                                            W6sSpan.size() > 1 ? W6sSpan[i] : W6sSpan[0],
                                            W7Span.size() > 1 ? W7Span[i] : W7Span[0],
                                            W8Span.size() > 1 ? W8Span[i] : W8Span[0],
                                            W9Span.size() > 1 ? W9Span[i] : W9Span[0],

                                            H1sSpan.size() > 1 ? H1sSpan[i] : H1sSpan[0],
                                            H1cSpan.size() > 1 ? H1cSpan[i] : H1cSpan[0],
                                            H2sSpan.size() > 1 ? H2sSpan[i] : H2sSpan[0],
                                            H2cSpan.size() > 1 ? H2cSpan[i] : H2cSpan[0],
                                            H3Span.size() > 1 ? H3Span[i] : H3Span[0],
                                            H4Span.size() > 1 ? H4Span[i] : H4Span[0],
                                            H5Span.size() > 1 ? H5Span[i] : H5Span[0],
                                            H6sSpan.size() > 1 ? H6sSpan[i] : H6sSpan[0],
                                            H7Span.size() > 1 ? H7Span[i] : H7Span[0],
                                            H8Span.size() > 1 ? H8Span[i] : H8Span[0],
                                            H9Span.size() > 1 ? H9Span[i] : H9Span[0],

                                            Z1sSpan.size() > 1 ? Z1sSpan[i] : Z1sSpan[0],
                                            Z1cSpan.size() > 1 ? Z1cSpan[i] : Z1cSpan[0],
                                            Z2sSpan.size() > 1 ? Z2sSpan[i] : Z2sSpan[0],
                                            Z2cSpan.size() > 1 ? Z2cSpan[i] : Z2cSpan[0],
                                            Z3Span.size() > 1 ? Z3Span[i] : Z3Span[0],
                                            Z4Span.size() > 1 ? Z4Span[i] : Z4Span[0],
                                            Z5Span.size() > 1 ? Z5Span[i] : Z5Span[0],
                                            Z6sSpan.size() > 1 ? Z6sSpan[i] : Z6sSpan[0],
                                            Z7Span.size() > 1 ? Z7Span[i] : Z7Span[0],
                                            Z8Span.size() > 1 ? Z8Span[i] : Z8Span[0],
                                            Z9Span.size() > 1 ? Z9Span[i] : Z9Span[0]);
        }
    }

    int getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars, const char */*rangeName*/) const override {
        if (matchArgs(allVars, analVars, cosThetaL_, cosThetaK_, phi_, t_)) return 1;
        if (matchArgs(allVars, analVars, cosThetaK_, phi_, t_)) return 6;
        if (matchArgs(allVars, analVars, cosThetaL_, phi_, t_)) return 7;
        if (matchArgs(allVars, analVars, cosThetaL_, cosThetaK_, t_)) return 8;
        if (matchArgs(allVars, analVars, cosThetaL_, cosThetaK_, phi_)) return 9;
        return 0;
    }

    double analyticalIntegral(int code, const char *rangeName) const override {
        double cosThetaL = cosThetaL_;
        double cosThetaK = cosThetaK_;
        double phi = phi_;
        double t = t_;

        double x = x_;
        double y = y_;
        double K1s = K1s_;
        double K1c = K1c_;
        double K2s = K2s_;
        double K2c = K2c_;
        double K3 = K3_;
        double K6s = K6s_ * sign_;
        double K9 = K9_ * sign_;
        double W1s = W1s_ * sign_;
        double W1c = W1c_ * sign_;
        double W2s = W2s_ * sign_;
        double W2c = W2c_ * sign_;
        double W3 = W3_ * sign_;
        double W6s = W6s_;
        double W9 = W9_;
        double H1s = H1s_;
        double H1c = H1c_;
        double H2s = H2s_;
        double H2c = H2c_;
        double H3 = H3_;
        double H6s = H6s_;
        double H9 = H9_;
        double Z1s = Z1s_ * sign_;
        double Z1c = Z1c_ * sign_;
        double Z2s = Z2s_ * sign_;
        double Z2c = Z2c_ * sign_;
        double Z3 = Z3_ * sign_;
        double Z6s = Z6s_ * sign_;
        double Z9 = Z9_ * sign_;

        if (code > 0 && code <= 5) {

            return (1./8)*(
                    (-3*K1c - 6*K1s + K2c + 2*K2s + (3*H1c + 6*H1s - H2c - 2*H2s)*y)/(-1 + y*y) +
                    sign_ * (3*W1c + 6*W1s - W2c - 2*W2s - (3*Z1c + 6*Z1s - Z2c - 2*Z2s)*x)/(1 + x*x));
        }

        if (code == 6) {
            double cosThetaL2 = cosThetaL * cosThetaL;
            return -(3 / (16 * (1 + x * x) * (-1 + y * y))) * (
                    K1c + 2 * K1s + 2 * K6s * cosThetaL + W1c + 2 * W1s + 2 * cosThetaL * W6s - Z1c * x -
                    2 * Z1s * x - 2 * Z6s * cosThetaL * x + K1c * x * x + 2 * K1s * x * x +
                    2 * cosThetaL * W6s * x * x +
                    K2c * (-1 + 2 * cosThetaL2) * (1 + x * x) + 2 * K2s * (-1 + 2 * cosThetaL2) * (1 + x * x) -
                    H1c * y - 2 * H1s * y - 2 * H6s * cosThetaL * y - H1c * x * x * y - 2 * H1s * x * x * y -
                    2 * H6s * cosThetaL * x * x * y + H2c * (1 - 2 * cosThetaL2) * (1 + x * x) * y -
                    2 * H2s * (-1 + 2 * cosThetaL2) * (1 + x * x) * y - 2 * K6s * cosThetaL * y * y - W1c * y * y -
                    2 * W1s * y * y + Z1c * x * y * y + 2 * Z1s * x * y * y +
                    2 * Z6s * cosThetaL * x * y * y - (-1 + 2 * cosThetaL2) * W2c * (-1 + y * y) -
                    2 * (-1 + 2 * cosThetaL2) * W2s * (-1 + y * y) + Z2c * (-1 + 2 * cosThetaL2) * x * (-1 + y * y) +
                    2 * Z2s * (-1 + 2 * cosThetaL2) * x * (-1 + y * y));
        }
        if (code == 7) {
            double cosThetaK2 = cosThetaK * cosThetaK;

            return (1 / (16 * (1 + x * x) * (-1 + y * y))) * 3 * (
                    -3 * W1s + W2s + 3 * Z1s * x - Z2s * x + 3 * K1s * (-1 + cosThetaK2) * (1 + x * x) -
                    K2s * (-1 + cosThetaK2) * (1 + x * x) +
                    cosThetaK2 * (K2c - 3 * W1c + 3 * W1s + W2c - W2s + (3 * Z1c - 3 * Z1s - Z2c + Z2s) * x +
                                  K2c * x * x - 3 * K1c * (1 + x * x)) + 3 * H1s * y - H2s * y +
                    ((3 * H1c - 3 * H1s - H2c + H2s) * cosThetaK2 -
                     (-3 * H1s + H2s + (-3 * H1c + 3 * H1s + H2c - H2s) * cosThetaK2) * x * x) * y +
                    (3 * W1s - W2s - 3 * Z1s * x + Z2s * x + cosThetaK2 * (3 * W1c - 3 * W1s - W2c + W2s +
                                                                           (-3 * Z1c + 3 * Z1s + Z2c - Z2s) * x)) * y * y);
        }
        if (code == 8) {
            return 1 / (16 * TMath::Pi() * (1 + x * x) * (-1 + y * y)) *
                   (-3 * K1c - 6 * K1s + K2c + 2 * K2s - 3 * W1c - 6 * W1s + W2c + 2 * W2s +
                    3 * Z1c * x + 6 * Z1s * x - Z2c * x - 2 * Z2s * x - 3 * K1c * x * x - 6 * K1s * x * x +
                    K2c * x * x + 2 * K2s * x * x + 3 * H1c * y + 6 * H1s * y -
                    H2c * y - 2 * H2s * y + 3 * H1c * x * x * y + 6 * H1s * x * x * y - H2c * x * x * y -
                    2 * H2s * x * x * y + 3 * W1c * y * y + 6 * W1s * y * y - W2c * y * y -
                    2 * W2s * y * y - 3 * Z1c * x * y * y - 6 * Z1s * x * y * y + Z2c * x * y * y +
                    2 * Z2s * x * y * y -
                    4 * (K3 + W3 - Z3 * x + K3 * x * x - H3 * y - H3 * x * x * y - W3 * y * y * Z3 * x * y * y) *
                    cos(2 * phi) -
                    4 * (K9 + W9 - Z9 * x + W9 * x * x - H9 * y - H9 * x * x * y - K9 * y * y + Z9 * x * y * y) *
                    sin(2 * phi));
        }
        if (code == 9) {
            return exp(-t) * (
                    (3 * W1c + 6 * W1s - W2c - 2 * W2s) * cos(t * x) +
                    (3 * K1c + 6 * K1s - K2c - 2 * K2s) * cosh(t * y) +
                    (-3 * Z1c - 6 * Z1s + Z2c + 2 * Z2s) * sin(t * x) +
                    (-3 * H1c - 6 * H1s + H2c + 2 * H2s) * sinh(t * y)) / 8;
        }
        return 0;
    }
};

void runFit(int nEvents, int nToys, bool useApproximation);


#endif //CERN_NEWTIMEDEPENDENTFIT_H
