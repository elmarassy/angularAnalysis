//
// Created by Mero Elmarassy on 6/6/25.
//

#include <RooCategory.h>
#include <RooSimultaneous.h>
#include <RooMCStudy.h>
#include "../include/timeDependentFit.h"

#include <RooFitResult.h>

using namespace RooFit;

timeDependentB::timeDependentB(const char *name, const char *title, RooAbsReal &costhetal, RooAbsReal &costhetak,
                               RooAbsReal& t,
                               RooAbsReal &phi, RooAbsReal &x, RooAbsReal &y, RooAbsReal &J1s, RooAbsReal &J1c,
                               RooAbsReal &J2s, RooAbsReal &J2c, RooAbsReal &J3, RooAbsReal &J4, RooAbsReal &J5,
                               RooAbsReal &J6s, RooAbsReal &J7, RooAbsReal &J8, RooAbsReal &J9, RooAbsReal &J1sBar, RooAbsReal &J1cBar,
                               RooAbsReal &J2sBar, RooAbsReal &J2cBar, RooAbsReal &J3Bar, RooAbsReal &J4Bar, RooAbsReal &J5Bar,
                               RooAbsReal &J6sBar, RooAbsReal &J7Bar, RooAbsReal &J8Bar, RooAbsReal &J9Bar, RooAbsReal &h1s,
                               RooAbsReal &h1c, RooAbsReal &h2s, RooAbsReal &h2c, RooAbsReal &h3, RooAbsReal &h4,
                               RooAbsReal &h5, RooAbsReal &h6, RooAbsReal &h7, RooAbsReal &h8, RooAbsReal &h9,
                               RooAbsReal &s1s, RooAbsReal &s1c, RooAbsReal &s2s, RooAbsReal &s2c, RooAbsReal &s3,
                               RooAbsReal &s4, RooAbsReal &s5, RooAbsReal &s6, RooAbsReal &s7, RooAbsReal &s8,
                               RooAbsReal &s9):
        RooAbsPdf(name, title),
        cosThetaL_("costhetal", "costhetal", this, costhetal),
        cosThetaK_("costhetak", "costhetak", this, costhetak),
        phi_("phi", "phi", this, phi),
        t_("t", "t", this, t),
        x_("x", "x", this, x),
        y_("y", "y", this, y),
        J1s_("J1s", "J1s", this, J1s),
        J1c_("J1c", "J1c", this, J1c),
        J2s_("J2s", "J2s", this, J2s),
        J2c_("J2c", "J2c", this, J2c),
        J3_("J3", "J3", this, J3),
        J4_("J4", "J4", this, J4),
        J5_("J5", "J5", this, J5),
        J6s_("J6s", "J6s", this, J6s),
        J7_("J7", "J7", this, J7),
        J8_("J8", "J8", this, J8),
        J9_("J9", "J9", this, J9),
        J1sBar_("J1s", "J1s", this, J1sBar),
        J1cBar_("J1c", "J1c", this, J1cBar),
        J2sBar_("J2s", "J2s", this, J2sBar),
        J2cBar_("J2c", "J2c", this, J2cBar),
        J3Bar_("J3", "J3", this, J3Bar),
        J4Bar_("J4", "J4", this, J4Bar),
        J5Bar_("J5", "J5", this, J5Bar),
        J6sBar_("J6s", "J6s", this, J6sBar),
        J7Bar_("J7", "J7", this, J7Bar),
        J8Bar_("J8", "J8", this, J8Bar),
        J9Bar_("J9", "J9", this, J9Bar),
        h1s_("h1s", "h1s", this, h1s),
        h1c_("h1c", "h1c", this, h1c),
        h2s_("h2s", "h2s", this, h2s),
        h2c_("h2c", "h2c", this, h2c),
        h3_("h3", "h3", this, h3),
        h4_("h4", "h4", this, h4),
        h5_("h5", "h5", this, h5),
        h6_("h6", "h6", this, h6),
        h7_("h7", "h7", this, h7),
        h8_("h8", "h8", this, h8),
        h9_("h9", "h9", this, h9),
        s1s_("s1s", "s1s", this, s1s),
        s1c_("s1c", "s1c", this, s1c),
        s2s_("s2s", "s2s", this, s2s),
        s2c_("s2c", "s2c", this, s2c),
        s3_("s3", "s3", this, s3),
        s4_("s4", "s4", this, s4),
        s5_("s5", "s5", this, s5),
        s6_("s6", "s6", this, s6),
        s7_("s7", "s7", this, s7),
        s8_("s8", "s8", this, s8),
        s9_("s9", "s9", this, s9) {}


timeDependentB::timeDependentB(timeDependentB const &other, const char *name):
        RooAbsPdf(other, name),
        cosThetaL_("costhetal", this, other.cosThetaL_),
        cosThetaK_("costhetak", this, other.cosThetaK_),
        phi_("phi", this, other.phi_),
        t_("t", this, other.t_),
        x_("x", this, other.x_),
        y_("y", this, other.y_),
        J1s_("J1s", this, other.J1s_),
        J1c_("J1c", this, other.J1c_),
        J2s_("J2s", this, other.J2s_),
        J2c_("J2c", this, other.J2c_),
        J3_("J3", this, other.J3_),
        J4_("J4", this, other.J4_),
        J5_("J5", this, other.J5_),
        J6s_("J6s", this, other.J6s_),
        J7_("J7", this, other.J7_),
        J8_("J8", this, other.J8_),
        J9_("J9", this, other.J9_),
        J1sBar_("J1sBar", this, other.J1sBar_),
        J1cBar_("J1cBar", this, other.J1cBar_),
        J2sBar_("J2sBar", this, other.J2sBar_),
        J2cBar_("J2cBar", this, other.J2cBar_),
        J3Bar_("J3Bar", this, other.J3Bar_),
        J4Bar_("J4Bar", this, other.J4Bar_),
        J5Bar_("J5Bar", this, other.J5Bar_),
        J6sBar_("J6sBar", this, other.J6sBar_),
        J7Bar_("J7Bar", this, other.J7Bar_),
        J8Bar_("J8Bar", this, other.J8Bar_),
        J9Bar_("J9Bar", this, other.J9Bar_),
        h1s_("h1s", this, other.h1s_),
        h1c_("h1c", this, other.h1c_),
        h2s_("h2s", this, other.h2s_),
        h2c_("h2c", this, other.h2c_),
        h3_("h3", this, other.h3_),
        h4_("h4", this, other.h4_),
        h5_("h5", this, other.h5_),
        h6_("h6", this, other.h6_),
        h7_("h7", this, other.h7_),
        h8_("h8", this, other.h8_),
        h9_("h9", this, other.h9_),
        s1s_("s1s", this, other.s1s_),
        s1c_("s1c", this, other.s1c_),
        s2s_("s2s", this, other.s2s_),
        s2c_("s2c", this, other.s2c_),
        s3_("s3", this, other.s3_),
        s4_("s4", this, other.s4_),
        s5_("s5", this, other.s5_),
        s6_("s6", this, other.s6_),
        s7_("s7", this, other.s7_),
        s8_("s8", this, other.s8_),
        s9_("s9", this, other.s9_) {}

inline double timeDependentB::evaluate_prob(double cosThetaL, double cosThetaK, double phi, double t, double x, double y,
                            double J1s, double J1c, double J2s, double J2c, double J3, double J4, double J5, double J6s, double J7, double J8, double J9,
                            double J1sBar, double J1cBar, double J2sBar, double J2cBar, double J3Bar, double J4Bar, double J5Bar, double J6sBar, double J7Bar, double J8Bar, double J9Bar,
                            double h1s, double h1c, double h2s, double h2c, double h3, double h4, double h5, double h6, double h7, double h8, double h9,
                            double s1s, double s1c, double s2s, double s2c, double s3, double s4, double s5, double s6, double s7, double s8, double s9) const {

    const double c = 9.0/32.0/TMath::Pi();
    const double cosThetaL2 = cosThetaL * cosThetaL;
    const double cosThetaK2 = cosThetaK * cosThetaK;
    const double cos2ThetaL = 2.0 * cosThetaL2 - 1.0;
    const double sinThetaK2 = 1.0 - cosThetaK2;
    const double sinThetaL2 = 1.0 - cosThetaL2;
    const double sinThetaL = sqrt(sinThetaL2);
    const double sinThetaK = sqrt(sinThetaK2);

    const double sin2ThetaL = 2.0 * sinThetaL * cosThetaL;
    const double sin2ThetaK = 2.0 * sinThetaK * cosThetaK;

    const double cosh_yt = cosh(y * t);
    const double sinh_yt = sinh(y * t);
    const double cos_xt = cos(x * t);
    const double sin_xt = sin(x * t);
    const double decay = exp( -t);

    auto helper = [=](double Ji, double JiBar, double si, double hi) {
        return (Ji + JiBar) * cosh_yt + (Ji-JiBar) * cos_xt - (hi * sinh_yt + si * sin_xt);
    };

    return 0.5 * c * decay * (
            helper(J1s, J1sBar, s1s, h1s) * sinThetaK2 +
            helper(J1c, J1cBar, s1c, h1c) * cosThetaK2 +
            helper(J2s, J2sBar, s2s, h2s) * sinThetaK2 * cos2ThetaL +
            helper(J2c, J2cBar, s2c, h2c) * cosThetaK2 * cos2ThetaL +
            helper(J3, J3Bar, s3, h3) * sinThetaK2 * sinThetaL2 * cos(2 * phi) +
            helper(J4, J4Bar, s4, h4) * sin2ThetaK * sin2ThetaL * cos(phi) +
            helper(J5, -J5Bar, s5, h5) * sin2ThetaK * sinThetaL * cos(phi) +
            helper(J6s, -J6sBar, s6, h6) * sinThetaK2 * cosThetaL +
            helper(J7, J7Bar, s7, h7) * sin2ThetaK * sinThetaL * sin(phi) +
            helper(J8, -J8Bar, s8, h8) * sin2ThetaK * sin2ThetaL * sin(phi) +
            helper(J9, -J9Bar, s9, h9) * sinThetaK2 * sinThetaL2 * sin(2 * phi)
    );
}

double timeDependentB::evaluate() const {
    return evaluate_prob( cosThetaL_,  cosThetaK_,  phi_,  t_,  x_,  y_,
             J1s_,  J1c_,  J2s_,  J2c_,  J3_,  J4_,  J5_,  J6s_,  J7_,  J8_,  J9_,
             J1sBar_,  J1cBar_,  J2sBar_,  J2cBar_,  J3Bar_,  J4Bar_,  J5Bar_,  J6sBar_,  J7Bar_,  J8Bar_,  J9Bar_,
             h1s_,  h1c_,  h2s_,  h2c_,  h3_,  h4_,  h5_,  h6_,  h7_,  h8_,  h9_,
             s1s_,  s1c_,  s2s_,  s2c_,  s3_,  s4_,  s5_,  s6_,  s7_,  s8_,  s9_);
}


void timeDependentB::doEval(RooFit::EvalContext &ctx) const {
    std::span<const double> cosThetaLSpan = ctx.at(cosThetaL_);
    std::span<const double> cosThetaKSpan = ctx.at(cosThetaK_);
    std::span<const double> tSpan = ctx.at(t_);
    std::span<const double> phiSpan = ctx.at(phi_);
    std::span<const double> xSpan = ctx.at(x_);
    std::span<const double> ySpan = ctx.at(y_);
    std::span<const double> J1sSpan = ctx.at(J1s_);
    std::span<const double> J1cSpan = ctx.at(J1c_);
    std::span<const double> J2sSpan = ctx.at(J2s_);
    std::span<const double> J2cSpan = ctx.at(J2c_);
    std::span<const double> J3Span = ctx.at(J3_);
    std::span<const double> J4Span = ctx.at(J4_);
    std::span<const double> J5Span = ctx.at(J5_);
    std::span<const double> J6sSpan = ctx.at(J6s_);
    std::span<const double> J7Span = ctx.at(J7_);
    std::span<const double> J8Span = ctx.at(J8_);
    std::span<const double> J9Span = ctx.at(J9_);
    std::span<const double> J1sBarSpan = ctx.at(J1sBar_);
    std::span<const double> J1cBarSpan = ctx.at(J1cBar_);
    std::span<const double> J2sBarSpan = ctx.at(J2sBar_);
    std::span<const double> J2cBarSpan = ctx.at(J2cBar_);
    std::span<const double> J3BarSpan = ctx.at(J3Bar_);
    std::span<const double> J4BarSpan = ctx.at(J4Bar_);
    std::span<const double> J5BarSpan = ctx.at(J5Bar_);
    std::span<const double> J6sBarSpan = ctx.at(J6sBar_);
    std::span<const double> J7BarSpan = ctx.at(J7Bar_);
    std::span<const double> J8BarSpan = ctx.at(J8Bar_);
    std::span<const double> J9BarSpan = ctx.at(J9Bar_);
    std::span<const double> h1sSpan = ctx.at(h1s_);
    std::span<const double> h1cSpan = ctx.at(h1c_);
    std::span<const double> h2sSpan = ctx.at(h2s_);
    std::span<const double> h2cSpan = ctx.at(h2c_);
    std::span<const double> h3Span = ctx.at(h3_);
    std::span<const double> h4Span = ctx.at(h4_);
    std::span<const double> h5Span = ctx.at(h5_);
    std::span<const double> h6Span = ctx.at(h6_);
    std::span<const double> h7Span = ctx.at(h7_);
    std::span<const double> h8Span = ctx.at(h8_);
    std::span<const double> h9Span = ctx.at(h9_);
    std::span<const double> s1sSpan = ctx.at(s1s_);
    std::span<const double> s1cSpan = ctx.at(s1c_);
    std::span<const double> s2sSpan = ctx.at(s2s_);
    std::span<const double> s2cSpan = ctx.at(s2c_);
    std::span<const double> s3Span = ctx.at(s3_);
    std::span<const double> s4Span = ctx.at(s4_);
    std::span<const double> s5Span = ctx.at(s5_);
    std::span<const double> s6Span = ctx.at(s6_);
    std::span<const double> s7Span = ctx.at(s7_);
    std::span<const double> s8Span = ctx.at(s8_);
    std::span<const double> s9Span = ctx.at(s9_);
    std::size_t n = ctx.output().size();
    for (std::size_t i = 0; i < n; ++i) {
        ctx.output()[i] = evaluate_prob(cosThetaLSpan.size() > 1 ? cosThetaLSpan[i] : cosThetaLSpan[0],
                                        cosThetaKSpan.size() > 1 ? cosThetaKSpan[i] : cosThetaKSpan[0],
                                        phiSpan.size() > 1 ? phiSpan[i] : phiSpan[0],
                                        tSpan.size() > 1 ? tSpan[i] : tSpan[0],
                                        xSpan.size() > 1 ? xSpan[i] : xSpan[0],
                                        ySpan.size() > 1 ? ySpan[i] : ySpan[0],
                                        J1sSpan.size() > 1 ? J1sSpan[i] : J1sSpan[0],
                                        J1cSpan.size() > 1 ? J1cSpan[i] : J1cSpan[0],
                                        J2sSpan.size() > 1 ? J2sSpan[i] : J2sSpan[0],
                                        J2cSpan.size() > 1 ? J2cSpan[i] : J2cSpan[0],
                                        J3Span.size() > 1 ? J3Span[i] : J3Span[0],
                                        J4Span.size() > 1 ? J4Span[i] : J4Span[0],
                                        J5Span.size() > 1 ? J5Span[i] : J5Span[0],
                                        J6sSpan.size() > 1 ? J6sSpan[i] : J6sSpan[0],
                                        J7Span.size() > 1 ? J7Span[i] : J7Span[0],
                                        J8Span.size() > 1 ? J8Span[i] : J8Span[0],
                                        J9Span.size() > 1 ? J9Span[i] : J9Span[0],
                                        J1sBarSpan.size() > 1 ? J1sBarSpan[i] : J1sBarSpan[0],
                                        J1cBarSpan.size() > 1 ? J1cBarSpan[i] : J1cBarSpan[0],
                                        J2sBarSpan.size() > 1 ? J2sBarSpan[i] : J2sBarSpan[0],
                                        J2cBarSpan.size() > 1 ? J2cBarSpan[i] : J2cBarSpan[0],
                                        J3BarSpan.size() > 1 ? J3BarSpan[i] : J3BarSpan[0],
                                        J4BarSpan.size() > 1 ? J4BarSpan[i] : J4BarSpan[0],
                                        J5BarSpan.size() > 1 ? J5BarSpan[i] : J5BarSpan[0],
                                        J6sBarSpan.size() > 1 ? J6sBarSpan[i] : J6sBarSpan[0],
                                        J7BarSpan.size() > 1 ? J7BarSpan[i] : J7BarSpan[0],
                                        J8BarSpan.size() > 1 ? J8BarSpan[i] : J8BarSpan[0],
                                        J9BarSpan.size() > 1 ? J9BarSpan[i] : J9BarSpan[0],
                                        h1sSpan.size() > 1 ? h1sSpan[i] : h1sSpan[0],
                                        h1cSpan.size() > 1 ? h1cSpan[i] : h1cSpan[0],
                                        h2sSpan.size() > 1 ? h2sSpan[i] : h2sSpan[0],
                                        h2cSpan.size() > 1 ? h2cSpan[i] : h2cSpan[0],
                                        h3Span.size() > 1 ? h3Span[i] : h3Span[0],
                                        h4Span.size() > 1 ? h4Span[i] : h4Span[0],
                                        h5Span.size() > 1 ? h5Span[i] : h5Span[0],
                                        h6Span.size() > 1 ? h6Span[i] : h6Span[0],
                                        h7Span.size() > 1 ? h7Span[i] : h7Span[0],
                                        h8Span.size() > 1 ? h8Span[i] : h8Span[0],
                                        h9Span.size() > 1 ? h9Span[i] : h9Span[0],
                                        s1sSpan.size() > 1 ? s1sSpan[i] : s1sSpan[0],
                                        s1cSpan.size() > 1 ? s1cSpan[i] : s1cSpan[0],
                                        s2sSpan.size() > 1 ? s2sSpan[i] : s2sSpan[0],
                                        s2cSpan.size() > 1 ? s2cSpan[i] : s2cSpan[0],
                                        s3Span.size() > 1 ? s3Span[i] : s3Span[0],
                                        s4Span.size() > 1 ? s4Span[i] : s4Span[0],
                                        s5Span.size() > 1 ? s5Span[i] : s5Span[0],
                                        s6Span.size() > 1 ? s6Span[i] : s6Span[0],
                                        s7Span.size() > 1 ? s7Span[i] : s7Span[0],
                                        s8Span.size() > 1 ? s8Span[i] : s8Span[0],
                                        s9Span.size() > 1 ? s9Span[i] : s9Span[0]);
    }
}

timeDependentBBar::timeDependentBBar(const char *name, const char *title, RooAbsReal &costhetal, RooAbsReal &costhetak,
                               RooAbsReal& t,
                               RooAbsReal &phi, RooAbsReal &x, RooAbsReal &y, RooAbsReal &J1s, RooAbsReal &J1c,
                               RooAbsReal &J2s, RooAbsReal &J2c, RooAbsReal &J3, RooAbsReal &J4, RooAbsReal &J5,
                               RooAbsReal &J6s, RooAbsReal &J7, RooAbsReal &J8, RooAbsReal &J9, RooAbsReal &J1sBar, RooAbsReal &J1cBar,
                               RooAbsReal &J2sBar, RooAbsReal &J2cBar, RooAbsReal &J3Bar, RooAbsReal &J4Bar, RooAbsReal &J5Bar,
                               RooAbsReal &J6sBar, RooAbsReal &J7Bar, RooAbsReal &J8Bar, RooAbsReal &J9Bar, RooAbsReal &h1s,
                               RooAbsReal &h1c, RooAbsReal &h2s, RooAbsReal &h2c, RooAbsReal &h3, RooAbsReal &h4,
                               RooAbsReal &h5, RooAbsReal &h6, RooAbsReal &h7, RooAbsReal &h8, RooAbsReal &h9,
                               RooAbsReal &s1s, RooAbsReal &s1c, RooAbsReal &s2s, RooAbsReal &s2c, RooAbsReal &s3,
                               RooAbsReal &s4, RooAbsReal &s5, RooAbsReal &s6, RooAbsReal &s7, RooAbsReal &s8,
                               RooAbsReal &s9):
        RooAbsPdf(name, title),
        cosThetaL_("costhetal", "costhetal", this, costhetal),
        cosThetaK_("costhetak", "costhetak", this, costhetak),
        phi_("phi", "phi", this, phi),
        t_("t", "t", this, t),
        x_("x", "x", this, x),
        y_("y", "y", this, y),
        J1s_("J1s", "J1s", this, J1s),
        J1c_("J1c", "J1c", this, J1c),
        J2s_("J2s", "J2s", this, J2s),
        J2c_("J2c", "J2c", this, J2c),
        J3_("J3", "J3", this, J3),
        J4_("J4", "J4", this, J4),
        J5_("J5", "J5", this, J5),
        J6s_("J6s", "J6s", this, J6s),
        J7_("J7", "J7", this, J7),
        J8_("J8", "J8", this, J8),
        J9_("J9", "J9", this, J9),
        J1sBar_("J1s", "J1s", this, J1sBar),
        J1cBar_("J1c", "J1c", this, J1cBar),
        J2sBar_("J2s", "J2s", this, J2sBar),
        J2cBar_("J2c", "J2c", this, J2cBar),
        J3Bar_("J3", "J3", this, J3Bar),
        J4Bar_("J4", "J4", this, J4Bar),
        J5Bar_("J5", "J5", this, J5Bar),
        J6sBar_("J6s", "J6s", this, J6sBar),
        J7Bar_("J7", "J7", this, J7Bar),
        J8Bar_("J8", "J8", this, J8Bar),
        J9Bar_("J9", "J9", this, J9Bar),
        h1s_("h1s", "h1s", this, h1s),
        h1c_("h1c", "h1c", this, h1c),
        h2s_("h2s", "h2s", this, h2s),
        h2c_("h2c", "h2c", this, h2c),
        h3_("h3", "h3", this, h3),
        h4_("h4", "h4", this, h4),
        h5_("h5", "h5", this, h5),
        h6_("h6", "h6", this, h6),
        h7_("h7", "h7", this, h7),
        h8_("h8", "h8", this, h8),
        h9_("h9", "h9", this, h9),
        s1s_("s1s", "s1s", this, s1s),
        s1c_("s1c", "s1c", this, s1c),
        s2s_("s2s", "s2s", this, s2s),
        s2c_("s2c", "s2c", this, s2c),
        s3_("s3", "s3", this, s3),
        s4_("s4", "s4", this, s4),
        s5_("s5", "s5", this, s5),
        s6_("s6", "s6", this, s6),
        s7_("s7", "s7", this, s7),
        s8_("s8", "s8", this, s8),
        s9_("s9", "s9", this, s9) {}

timeDependentBBar::timeDependentBBar(timeDependentBBar const &other, const char *name):
        RooAbsPdf(other, name),
        cosThetaL_("costhetal", this, other.cosThetaL_),
        cosThetaK_("costhetak", this, other.cosThetaK_),
        phi_("phi", this, other.phi_),
        t_("t", this, other.t_),
        x_("x", this, other.x_),
        y_("y", this, other.y_),
        J1s_("J1s", this, other.J1s_),
        J1c_("J1c", this, other.J1c_),
        J2s_("J2s", this, other.J2s_),
        J2c_("J2c", this, other.J2c_),
        J3_("J3", this, other.J3_),
        J4_("J4", this, other.J4_),
        J5_("J5", this, other.J5_),
        J6s_("J6s", this, other.J6s_),
        J7_("J7", this, other.J7_),
        J8_("J8", this, other.J8_),
        J9_("J9", this, other.J9_),
        J1sBar_("J1sBar", this, other.J1sBar_),
        J1cBar_("J1cBar", this, other.J1cBar_),
        J2sBar_("J2sBar", this, other.J2sBar_),
        J2cBar_("J2cBar", this, other.J2cBar_),
        J3Bar_("J3Bar", this, other.J3Bar_),
        J4Bar_("J4Bar", this, other.J4Bar_),
        J5Bar_("J5Bar", this, other.J5Bar_),
        J6sBar_("J6sBar", this, other.J6sBar_),
        J7Bar_("J7Bar", this, other.J7Bar_),
        J8Bar_("J8Bar", this, other.J8Bar_),
        J9Bar_("J9Bar", this, other.J9Bar_),
        h1s_("h1s", this, other.h1s_),
        h1c_("h1c", this, other.h1c_),
        h2s_("h2s", this, other.h2s_),
        h2c_("h2c", this, other.h2c_),
        h3_("h3", this, other.h3_),
        h4_("h4", this, other.h4_),
        h5_("h5", this, other.h5_),
        h6_("h6", this, other.h6_),
        h7_("h7", this, other.h7_),
        h8_("h8", this, other.h8_),
        h9_("h9", this, other.h9_),
        s1s_("s1s", this, other.s1s_),
        s1c_("s1c", this, other.s1c_),
        s2s_("s2s", this, other.s2s_),
        s2c_("s2c", this, other.s2c_),
        s3_("s3", this, other.s3_),
        s4_("s4", this, other.s4_),
        s5_("s5", this, other.s5_),
        s6_("s6", this, other.s6_),
        s7_("s7", this, other.s7_),
        s8_("s8", this, other.s8_),
        s9_("s9", this, other.s9_) {}

inline double timeDependentBBar::evaluate_prob(double cosThetaL, double cosThetaK, double phi, double t, double x, double y,
                            double J1s, double J1c, double J2s, double J2c, double J3, double J4, double J5, double J6s, double J7, double J8, double J9,
                            double J1sBar, double J1cBar, double J2sBar, double J2cBar, double J3Bar, double J4Bar, double J5Bar, double J6sBar, double J7Bar, double J8Bar, double J9Bar,
                            double h1s, double h1c, double h2s, double h2c, double h3, double h4, double h5, double h6, double h7, double h8, double h9,
                            double s1s, double s1c, double s2s, double s2c, double s3, double s4, double s5, double s6, double s7, double s8, double s9) const {

    const double c = 9.0/32.0/TMath::Pi();
    const double cosThetaL2 = cosThetaL * cosThetaL;
    const double cosThetaK2 = cosThetaK * cosThetaK;
    const double cos2ThetaL = 2.0 * cosThetaL2 - 1.0;
    const double sinThetaK2 = 1.0 - cosThetaK2;
    const double sinThetaL2 = 1.0 - cosThetaL2;
    const double sinThetaL = sqrt(sinThetaL2);
    const double sinThetaK = sqrt(sinThetaK2);

    const double sin2ThetaL = 2.0 * sinThetaL * cosThetaL;
    const double sin2ThetaK = 2.0 * sinThetaK * cosThetaK;

    const double cosh_yt = cosh(y * t);
    const double sinh_yt = sinh(y * t);
    const double cos_xt = cos(x * t);
    const double sin_xt = sin(x * t);
    const double decay = exp( -t);

    auto helper = [=](double Ji, double JiBar, double si, double hi) {
        return (Ji + JiBar) * cosh_yt - (Ji-JiBar) * cos_xt - (hi * sinh_yt - si * sin_xt);
    };

    return 0.5 * c * decay * (
            helper(J1s, J1sBar, s1s, h1s) * sinThetaK2 +
            helper(J1c, J1cBar, s1c, h1c) * cosThetaK2 +
            helper(J2s, J2sBar, s2s, h2s) * sinThetaK2 * cos2ThetaL +
            helper(J2c, J2cBar, s2c, h2c) * cosThetaK2 * cos2ThetaL +
            helper(J3, J3Bar, s3, h3) * sinThetaK2 * sinThetaL2 * cos(2 * phi) +
            helper(J4, J4Bar, s4, h4) * sin2ThetaK * sin2ThetaL * cos(phi) +
            helper(J5, -J5Bar, s5, h5) * sin2ThetaK * sinThetaL * cos(phi) +
            helper(J6s, -J6sBar, s6, h6) * sinThetaK2 * cosThetaL +
            helper(J7, J7Bar, s7, h7) * sin2ThetaK * sinThetaL * sin(phi) +
            helper(J8, -J8Bar, s8, h8) * sin2ThetaK * sin2ThetaL * sin(phi) +
            helper(J9, -J9Bar, s9, h9) * sinThetaK2 * sinThetaL2 * sin(2 * phi)
    );
}

double timeDependentBBar::evaluate() const {
    return evaluate_prob( cosThetaL_,  cosThetaK_,  phi_,  t_,  x_,  y_,
             J1s_,  J1c_,  J2s_,  J2c_,  J3_,  J4_,  J5_,  J6s_,  J7_,  J8_,  J9_,
             J1sBar_,  J1cBar_,  J2sBar_,  J2cBar_,  J3Bar_,  J4Bar_,  J5Bar_,  J6sBar_,  J7Bar_,  J8Bar_,  J9Bar_,
             h1s_,  h1c_,  h2s_,  h2c_,  h3_,  h4_,  h5_,  h6_,  h7_,  h8_,  h9_,
             s1s_,  s1c_,  s2s_,  s2c_,  s3_,  s4_,  s5_,  s6_,  s7_,  s8_,  s9_);
}

void timeDependentBBar::doEval(RooFit::EvalContext &ctx) const {
    std::span<const double> cosThetaLSpan = ctx.at(cosThetaL_);
    std::span<const double> cosThetaKSpan = ctx.at(cosThetaK_);
    std::span<const double> tSpan = ctx.at(t_);
    std::span<const double> phiSpan = ctx.at(phi_);
    std::span<const double> xSpan = ctx.at(x_);
    std::span<const double> ySpan = ctx.at(y_);
    std::span<const double> J1sSpan = ctx.at(J1s_);
    std::span<const double> J1cSpan = ctx.at(J1c_);
    std::span<const double> J2sSpan = ctx.at(J2s_);
    std::span<const double> J2cSpan = ctx.at(J2c_);
    std::span<const double> J3Span = ctx.at(J3_);
    std::span<const double> J4Span = ctx.at(J4_);
    std::span<const double> J5Span = ctx.at(J5_);
    std::span<const double> J6sSpan = ctx.at(J6s_);
    std::span<const double> J7Span = ctx.at(J7_);
    std::span<const double> J8Span = ctx.at(J8_);
    std::span<const double> J9Span = ctx.at(J9_);
    std::span<const double> J1sBarSpan = ctx.at(J1sBar_);
    std::span<const double> J1cBarSpan = ctx.at(J1cBar_);
    std::span<const double> J2sBarSpan = ctx.at(J2sBar_);
    std::span<const double> J2cBarSpan = ctx.at(J2cBar_);
    std::span<const double> J3BarSpan = ctx.at(J3Bar_);
    std::span<const double> J4BarSpan = ctx.at(J4Bar_);
    std::span<const double> J5BarSpan = ctx.at(J5Bar_);
    std::span<const double> J6sBarSpan = ctx.at(J6sBar_);
    std::span<const double> J7BarSpan = ctx.at(J7Bar_);
    std::span<const double> J8BarSpan = ctx.at(J8Bar_);
    std::span<const double> J9BarSpan = ctx.at(J9Bar_);
    std::span<const double> h1sSpan = ctx.at(h1s_);
    std::span<const double> h1cSpan = ctx.at(h1c_);
    std::span<const double> h2sSpan = ctx.at(h2s_);
    std::span<const double> h2cSpan = ctx.at(h2c_);
    std::span<const double> h3Span = ctx.at(h3_);
    std::span<const double> h4Span = ctx.at(h4_);
    std::span<const double> h5Span = ctx.at(h5_);
    std::span<const double> h6Span = ctx.at(h6_);
    std::span<const double> h7Span = ctx.at(h7_);
    std::span<const double> h8Span = ctx.at(h8_);
    std::span<const double> h9Span = ctx.at(h9_);
    std::span<const double> s1sSpan = ctx.at(s1s_);
    std::span<const double> s1cSpan = ctx.at(s1c_);
    std::span<const double> s2sSpan = ctx.at(s2s_);
    std::span<const double> s2cSpan = ctx.at(s2c_);
    std::span<const double> s3Span = ctx.at(s3_);
    std::span<const double> s4Span = ctx.at(s4_);
    std::span<const double> s5Span = ctx.at(s5_);
    std::span<const double> s6Span = ctx.at(s6_);
    std::span<const double> s7Span = ctx.at(s7_);
    std::span<const double> s8Span = ctx.at(s8_);
    std::span<const double> s9Span = ctx.at(s9_);
    std::size_t n = ctx.output().size();
    for (std::size_t i = 0; i < n; ++i) {
        ctx.output()[i] = evaluate_prob(cosThetaLSpan.size() > 1 ? cosThetaLSpan[i] : cosThetaLSpan[0],
                                        cosThetaKSpan.size() > 1 ? cosThetaKSpan[i] : cosThetaKSpan[0],
                                        phiSpan.size() > 1 ? phiSpan[i] : phiSpan[0],
                                        tSpan.size() > 1 ? tSpan[i] : tSpan[0],
                                        xSpan.size() > 1 ? xSpan[i] : xSpan[0],
                                        ySpan.size() > 1 ? ySpan[i] : ySpan[0],
                                        J1sSpan.size() > 1 ? J1sSpan[i] : J1sSpan[0],
                                        J1cSpan.size() > 1 ? J1cSpan[i] : J1cSpan[0],
                                        J2sSpan.size() > 1 ? J2sSpan[i] : J2sSpan[0],
                                        J2cSpan.size() > 1 ? J2cSpan[i] : J2cSpan[0],
                                        J3Span.size() > 1 ? J3Span[i] : J3Span[0],
                                        J4Span.size() > 1 ? J4Span[i] : J4Span[0],
                                        J5Span.size() > 1 ? J5Span[i] : J5Span[0],
                                        J6sSpan.size() > 1 ? J6sSpan[i] : J6sSpan[0],
                                        J7Span.size() > 1 ? J7Span[i] : J7Span[0],
                                        J8Span.size() > 1 ? J8Span[i] : J8Span[0],
                                        J9Span.size() > 1 ? J9Span[i] : J9Span[0],
                                        J1sBarSpan.size() > 1 ? J1sBarSpan[i] : J1sBarSpan[0],
                                        J1cBarSpan.size() > 1 ? J1cBarSpan[i] : J1cBarSpan[0],
                                        J2sBarSpan.size() > 1 ? J2sBarSpan[i] : J2sBarSpan[0],
                                        J2cBarSpan.size() > 1 ? J2cBarSpan[i] : J2cBarSpan[0],
                                        J3BarSpan.size() > 1 ? J3BarSpan[i] : J3BarSpan[0],
                                        J4BarSpan.size() > 1 ? J4BarSpan[i] : J4BarSpan[0],
                                        J5BarSpan.size() > 1 ? J5BarSpan[i] : J5BarSpan[0],
                                        J6sBarSpan.size() > 1 ? J6sBarSpan[i] : J6sBarSpan[0],
                                        J7BarSpan.size() > 1 ? J7BarSpan[i] : J7BarSpan[0],
                                        J8BarSpan.size() > 1 ? J8BarSpan[i] : J8BarSpan[0],
                                        J9BarSpan.size() > 1 ? J9BarSpan[i] : J9BarSpan[0],
                                        h1sSpan.size() > 1 ? h1sSpan[i] : h1sSpan[0],
                                        h1cSpan.size() > 1 ? h1cSpan[i] : h1cSpan[0],
                                        h2sSpan.size() > 1 ? h2sSpan[i] : h2sSpan[0],
                                        h2cSpan.size() > 1 ? h2cSpan[i] : h2cSpan[0],
                                        h3Span.size() > 1 ? h3Span[i] : h3Span[0],
                                        h4Span.size() > 1 ? h4Span[i] : h4Span[0],
                                        h5Span.size() > 1 ? h5Span[i] : h5Span[0],
                                        h6Span.size() > 1 ? h6Span[i] : h6Span[0],
                                        h7Span.size() > 1 ? h7Span[i] : h7Span[0],
                                        h8Span.size() > 1 ? h8Span[i] : h8Span[0],
                                        h9Span.size() > 1 ? h9Span[i] : h9Span[0],
                                        s1sSpan.size() > 1 ? s1sSpan[i] : s1sSpan[0],
                                        s1cSpan.size() > 1 ? s1cSpan[i] : s1cSpan[0],
                                        s2sSpan.size() > 1 ? s2sSpan[i] : s2sSpan[0],
                                        s2cSpan.size() > 1 ? s2cSpan[i] : s2cSpan[0],
                                        s3Span.size() > 1 ? s3Span[i] : s3Span[0],
                                        s4Span.size() > 1 ? s4Span[i] : s4Span[0],
                                        s5Span.size() > 1 ? s5Span[i] : s5Span[0],
                                        s6Span.size() > 1 ? s6Span[i] : s6Span[0],
                                        s7Span.size() > 1 ? s7Span[i] : s7Span[0],
                                        s8Span.size() > 1 ? s8Span[i] : s8Span[0],
                                        s9Span.size() > 1 ? s9Span[i] : s9Span[0]);
    }
}

void fitTimeDependent(int nruns, int ngen) {

    RooRealVar x("x", "x", 26.93);
    RooRealVar y("y", "y", 0.124);

    RooRealVar J1s("J1s", "J1s", 0.7, 0, 1);
    RooRealVar J1c("J1c", "J1c", 0.4, -1, 1);
    RooRealVar J2s("J2s", "J2s", 0, -1, 1);
    RooRealVar J2c("J2c", "J2c", 0, -1, 1);
    RooRealVar J3("J3", "J3", 0, -1, 1);
    RooRealVar J4("J4", "J4", 0, -1, 1);
    RooRealVar J5("J5", "J5", 0, -1, 1);
    RooRealVar J6s("J6s", "J6s", 0, -1, 1);
    RooRealVar J7("J7", "J7", 0, -1, 1);
    RooRealVar J8("J8", "J8", 0, -1, 1);
    RooRealVar J9("J9", "J9", 0, -1, 1);

    RooRealVar J1sBar("J1sBar", "J1sBar", 0.7, 0, 1);
    RooRealVar J1cBar("J1cBar", "J1cBar", 0.4, -1, 1);
    RooRealVar J2sBar("J2sBar", "J2sBar", 0, -1, 1);
    RooRealVar J2cBar("J2cBar", "J2cBar", 0, -1, 1);
    RooRealVar J3Bar("J3Bar", "J3Bar", 0, -1, 1);
    RooRealVar J4Bar("J4Bar", "J4Bar", 0, -1, 1);
    RooRealVar J5Bar("J5Bar", "J5Bar", 0, -1, 1);
    RooRealVar J6sBar("J6sBar", "J6sBar", 0, -1, 1);
    RooRealVar J7Bar("J7Bar", "J7Bar", 0, -1, 1);
    RooRealVar J8Bar("J8Bar", "J8Bar", 0, -1, 1);
    RooRealVar J9Bar("J9Bar", "J9Bar", 0, -1, 1);

    RooRealVar h1s("h1s", "h1s", 0, -1, 1);
    RooRealVar h1c("h1c", "h1c", 0, -1, 1);
    RooRealVar h2s("h2s", "h2s", 0, -1, 1);
    RooRealVar h2c("h2c", "h2c", 0, -1, 1);
    RooRealVar h3("h3", "h3", 0, -1, 1);
    RooRealVar h4("h4", "h4", 0, -1, 1);
    RooRealVar h5("h5", "h5", 0, -1, 1);
    RooRealVar h6("h6s", "h6s", 0, -1, 1);
    RooRealVar h7("h7", "h7", 0, -1, 1);
    RooRealVar h8("h8", "h8", 0, -1, 1);
    RooRealVar h9("h9", "h9", 0, -1, 1);

    RooRealVar s1s("s1s", "s1s", 0, -1, 1);
    RooRealVar s1c("s1c", "s1c", 0, -1, 1);
    RooRealVar s2s("s2s", "s2s", 0, -1, 1);
    RooRealVar s2c("s2c", "s2c", 0, -1, 1);
    RooRealVar s3("s3", "s3", 0, -1, 1);
    RooRealVar s4("s4", "s4", 0, -1, 1);
    RooRealVar s5("s5", "s5", 0, -1, 1);
    RooRealVar s6("s6s", "s6s", 0, -1, 1);
    RooRealVar s7("s7", "s7", 0, -1, 1);
    RooRealVar s8("s8", "s8", 0, -1, 1);
    RooRealVar s9("s9", "s9", 0, -1, 1);

    // RooFormulaVar J1s("J1s", "("
    //                          "(16*(1+x**2)*(1-y**2)) - "
    //                          "(12*J1c-4*J2c-8*J2s)*(1-y**2)*(2+x**2-y**2) - "
    //                          "2*x*(2*s2s+s2c-6*s1s-3*s1c)*(1-y**2) - "
    //                          "y*(2*h2c+4*h2s-6*h1c-12*h1s)*(1+x**2)*(1-y**2)"
    //                          ") / (24*(1-y**2)*(2+x**2))",
    //                          RooArgList(x, y, J1c, J2s, J2c, h1s, h1c, h2s, h2c, s1s, s1c, s2s, s2c));
    //
    // RooFormulaVar J1sBar("J1sBar", "("
    //                          "(16*(1+x**2)*(1-y**2)) - "
    //                          "(12*J1cBar-4*J2cBar-8*J2sBar)*(1-y**2)*(2+x**2-y**2) - "
    //                          "2*x*(2*s2s+s2c-6*s1s-3*s1c)*(1-y**2) + "
    //                          "y*(2*h2c+4*h2s-6*h1c-12*h1s)*(1+x**2)*(1-y**2)"
    //                          ") / (24*(1-y**2)*(2+x**2))",
    //                   RooArgList(x, y, J1cBar, J2sBar, J2cBar, h1s, h1c, h2s, h2c, s1s, s1c, s2s, s2c));

    RooRealVar cosThetaK("cosThetaK", "cosThetaK", -1.0, 1.0);
    RooRealVar cosThetaL("cosThetaL", "cosThetaL", -1, 1);
    RooRealVar phi("phi", "phi", -TMath::Pi(), TMath::Pi());
    RooRealVar t("gamma*t", "gamma*t", 0.0, 10.0);

    timeDependentB b("b", "b", cosThetaL, cosThetaK, phi, t, x, y, J1s, J1c, J2s,J2c,J3, J4,
    J5, J6s, J7, J8,J9, J1sBar, J1cBar, J2sBar,J2cBar,J3Bar, J4Bar,
    J5Bar, J6sBar, J7Bar, J8Bar,J9Bar, h1s, h1c, h2s, h2c, h3,h4, h5, h6, h7, h8, h9,
                     s1s, s1c, s2s, s2c, s3,s4, s5, s6, s7, s8, s9);

    timeDependentB bBar("bBar", "bBar", cosThetaL, cosThetaK, phi, t, x, y, J1s, J1c, J2s,J2c,J3, J4,
    J5, J6s, J7, J8,J9, J1sBar, J1cBar, J2sBar,J2cBar,J3Bar, J4Bar,
    J5Bar, J6sBar, J7Bar, J8Bar,J9Bar, h1s, h1c, h2s, h2c, h3,h4, h5, h6, h7, h8, h9,
                     s1s, s1c, s2s, s2c, s3,s4, s5, s6, s7, s8, s9);


    RooCategory type("type", "type");
    type.defineType("B0");
    type.defineType("B0Bar");
    std::cout << "creating data" << std::endl;

    std::unique_ptr<RooDataSet> bData{b.generate(RooArgSet(cosThetaL, cosThetaK, phi, t), ngen)};
    std::unique_ptr<RooDataSet> bBarData{bBar.generate(RooArgSet(cosThetaL, cosThetaK, phi, t), 0)};
    RooFitResult* result = b.fitTo(*bData, Save(), FitOptions(Strategy(2), Hesse(true), Save(true), PrintEvalErrors(0), EvalBackend("legacy"), NumCPU(16)));
    //
    // RooDataSet combData("combData", "combined data", RooArgSet(cosThetaL, cosThetaK, phi, t), Index(type),
    //     Import({{"B0", bData.get()}, {"B0Bar", bBarData.get()}}));
    //
    // RooSimultaneous simPdf("simPdf", "Simultaneous PDF", {{"B0", &b}, {"B0Bar", &bBar}}, type);
    //
    //
    // std::cout << "created data" << std::endl;
    //
    // RooFitResult* result = simPdf.fitTo(combData, Save(), FitOptions(Strategy(2), Hesse(true), Save(true), PrintEvalErrors(0), EvalBackend("legacy"), NumCPU(16)));

    // RooDataSet test("test", "test", RooArgSet(type));
    // type.setLabel("B0");
    // test.add(RooArgSet(type), ngen);
    // type.setLabel("B0Bar");
    // test.add(RooArgSet(type), ngen);


//     RooMCStudy *mc = new RooMCStudy(simPdf, RooArgSet(cosThetaL, cosThetaK, phi, t, type),
//                   RooFit::Binned(false),
// //                  RooFit::Silence(),
// //                  RooFit::Verbose(false),
//                   RooFit::FitModel(simPdf),
//                   RooFit::ProtoData(test),
//                   FitOptions(Strategy(2), Hesse(true), Save(true), PrintEvalErrors(0), EvalBackend("legacy"), NumCPU(16))
//     );

    // mc->generateAndFit(nruns, ngen, true);

}

