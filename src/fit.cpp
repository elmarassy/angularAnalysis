//
// Created by Mero Elmarassy on 6/5/25.
//

#include "../include/fit.h"

timeIndependent::timeIndependent(const char *name, const char *title, RooAbsReal& costhetal, RooAbsReal& costhetak, RooAbsReal& phi,
                                 RooAbsReal& J1s, RooAbsReal& J1c, RooAbsReal& J2s, RooAbsReal& J2c, RooAbsReal& J3, RooAbsReal& J4,
                                 RooAbsReal& J5, RooAbsReal& J6, RooAbsReal& J7, RooAbsReal& J8, RooAbsReal& J9):
    RooAbsPdf(name, title),
    cosThetaL_("costhetal", "costhetal", this, costhetal),
    cosThetaK_("costhetak", "costhetak", this, costhetak),
    phi_("phi", "phi", this, phi),
    J1s_("J1s", "J1s", this, J1s),
    J1c_("J1c", "J1c", this, J1c),
    J2s_("J2s", "J2s", this, J2s),
    J2c_("J2c", "J2c", this, J2c),
    J3_("J3", "J3", this, J3),
    J4_("J4", "J4", this, J4),
    J5_("J5", "J5", this, J5),
    J6_("J6", "J6", this, J6),
    J7_("J7", "J7", this, J7),
    J8_("J8", "J8", this, J8),
    J9_("J9", "J9", this, J9) {}

timeIndependent::timeIndependent(timeIndependent const &other, const char *name):
    RooAbsPdf(other, name),
    cosThetaL_("costhetal", this, other.cosThetaL_),
    cosThetaK_("costhetak", this, other.cosThetaK_),
    phi_("phi", this, other.phi_),
    J1s_("J1s", this, other.J1s_),
    J1c_("J1c", this, other.J1c_),
    J2s_("J2s", this, other.J2s_),
    J2c_("J2c", this, other.J2c_),
    J3_("J3", this, other.J3_),
    J4_("J4", this, other.J4_),
    J5_("J5", this, other.J5_),
    J6_("J6", this, other.J6_),
    J7_("J7", this, other.J7_),
    J8_("J8", this, other.J8_),
    J9_("J9", this, other.J9_) {}

double timeIndependent::evaluate_prob(double cosThetaL, double cosThetaK, double phi, double J1s, double J1c, double J2s,
                                      double J2c, double J3, double J4, double J5, double J6, double J7, double J8,
                                      double J9) const {
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
    return c * (
            J1s * sinThetaK2 +
            J1c * cosThetaK2 +
            J2s * sinThetaK2 * cos2ThetaL +
            J2c * cosThetaK2 * cos2ThetaL +
            J3 * sinThetaK2 * cosThetaK2 * cos(2 * phi) +
            J4 * sin2ThetaK * sin2ThetaL * cos(phi) +
            J5 * sin2ThetaK * sinThetaL * cos(phi) +
            J6 * sinThetaK2 * cosThetaL +
            J7 * sin2ThetaK * sinThetaL * sin(phi) +
            J8 * sin2ThetaK * sin2ThetaL * sin(phi) +
            J9 * sinThetaK2 * sinThetaL2 * sin(2 * phi)
    );
}

double timeIndependent::evaluate() const {
    return evaluate_prob(cosThetaL_, cosThetaK_, phi_, J1s_, J1c_, J2s_, J2c_, J3_, J4_, J5_, J6_, J7_, J8_, J9_);
}

void timeIndependent::doEval(RooFit::EvalContext &ctx) const {
    std::span<const double> costhetalSpan = ctx.at(cosThetaL_);
    std::span<const double> costhetakSpan = ctx.at(cosThetaK_);
    std::span<const double> phiSpan = ctx.at(phi_);

    std::span<const double> J1sSpan = ctx.at(J1s_);
    std::span<const double> J1cSpan = ctx.at(J1c_);
    std::span<const double> J2sSpan = ctx.at(J2s_);
    std::span<const double> J2cSpan = ctx.at(J2c_);
    std::span<const double> J3Span = ctx.at(J3_);
    std::span<const double> J4Span = ctx.at(J4_);
    std::span<const double> J5Span = ctx.at(J5_);
    std::span<const double> J6Span = ctx.at(J6_);
    std::span<const double> J7Span = ctx.at(J7_);
    std::span<const double> J8Span = ctx.at(J8_);
    std::span<const double> J9Span = ctx.at(J9_);

    std::size_t n = ctx.output().size();
    for (std::size_t i = 0; i < n; ++i) {
    ctx.output()[i] = evaluate_prob(costhetalSpan.size() > 1 ? costhetalSpan[i] : costhetalSpan[0],
                                    costhetakSpan.size() > 1 ? costhetakSpan[i] : costhetakSpan[0],
                                    phiSpan.size() > 1 ? phiSpan[i] : phiSpan[0],
                                    J1sSpan.size() > 1 ? J1sSpan[i] : J1sSpan[0],
                                    J1cSpan.size() > 1 ? J1cSpan[i] : J1cSpan[0],
                                    J2sSpan.size() > 1 ? J2sSpan[i] : J2sSpan[0],
                                    J2cSpan.size() > 1 ? J2cSpan[i] : J2cSpan[0],
                                    J3Span.size() > 1 ? J3Span[i] : J3Span[0],
                                    J4Span.size() > 1 ? J4Span[i] : J4Span[0],
                                    J5Span.size() > 1 ? J5Span[i] : J5Span[0],
                                    J6Span.size() > 1 ? J6Span[i] : J6Span[0],
                                    J7Span.size() > 1 ? J7Span[i] : J7Span[0],
                                    J8Span.size() > 1 ? J8Span[i] : J8Span[0],
                                    J9Span.size() > 1 ? J9Span[i] : J9Span[0]);
    }
}


timeDependentUntagged::timeDependentUntagged(const char *name, const char *title, RooAbsReal& costhetal, RooAbsReal& costhetak, RooAbsReal& phi,
                                             RooAbsReal& t, RooAbsReal& gamma, RooAbsReal& y,
                                             RooAbsReal& SJ1s, RooAbsReal& SJ1c, RooAbsReal& SJ2s, RooAbsReal& SJ2c, RooAbsReal& SJ3,
                                             RooAbsReal& SJ4, RooAbsReal& DJ5, RooAbsReal& DJ6s, RooAbsReal& SJ7, RooAbsReal& DJ8, RooAbsReal& DJ9,
                                             RooAbsReal& h1s, RooAbsReal& h1c, RooAbsReal& h2s, RooAbsReal& h2c, RooAbsReal& h3,
                                             RooAbsReal& h4, RooAbsReal& h5, RooAbsReal& h6, RooAbsReal& h7, RooAbsReal& h8, RooAbsReal& h9):
        RooAbsPdf(name, title),
        cosThetaL_("costhetal", "costhetal", this, costhetal),
        cosThetaK_("costhetak", "costhetak", this, costhetak),
        phi_("phi", "phi", this, phi),
        t_("t", "t", this, t),
        gamma_("gamma", "gamma", this, gamma),
        y_("y", "y", this, y),
        SJ1s_("SJ1s", "SJ1s", this, SJ1s),
        SJ1c_("SJ1c", "SJ1c", this, SJ1c),
        SJ2s_("SJ2s", "SJ2s", this, SJ2s),
        SJ2c_("SJ2c", "SJ2c", this, SJ2c),
        SJ3_("SJ3", "SJ3", this, SJ3),
        SJ4_("SJ4", "SJ4", this, SJ4),
        DJ5_("SJ5", "DJ5", this, DJ5),
        DJ6s_("DJ6s", "DJ6s", this, DJ6s),
        SJ7_("SJ7", "SJ7", this, SJ7),
        DJ8_("DJ8", "DJ8", this, DJ8),
        DJ9_("DJ9", "DJ9", this, DJ9),
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
        h9_("h9", "h9", this, h9) {}

timeDependentUntagged::timeDependentUntagged(timeDependentUntagged const &other, const char *name):
        RooAbsPdf(other, name),
        cosThetaL_("costhetal", this, other.cosThetaL_),
        cosThetaK_("costhetak", this, other.cosThetaK_),
        phi_("phi", this, other.phi_),
        t_("t", this, other.t_),
        gamma_("gamma", this, other.gamma_),
        y_("y", this, other.y_),
        SJ1s_("SJ1s", this, other.SJ1s_),
        SJ1c_("SJ1c", this, other.SJ1c_),
        SJ2s_("SJ2s", this, other.SJ2s_),
        SJ2c_("SJ2c", this, other.SJ2c_),
        SJ3_("SJ3", this, other.SJ3_),
        SJ4_("SJ4", this, other.SJ4_),
        DJ5_("DJ5", this, other.DJ5_),
        DJ6s_("DJ6s", this, other.DJ6s_),
        SJ7_("SJ7", this, other.SJ7_),
        DJ8_("DJ8", this, other.DJ8_),
        DJ9_("DJ9", this, other.DJ9_),
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
        h9_("h9", this, other.h9_){}

double timeDependentUntagged::evaluate_prob(double cosThetaL, double cosThetaK, double phi, double t, double gamma, double y,
                                            double SJ1s, double SJ1c, double SJ2s, double SJ2c, double SJ3, double SJ4, double DJ5, double DJ6s, double SJ7, double DJ8, double DJ9,
                                            double h1s, double h1c, double h2s, double h2c, double h3, double h4, double h5, double h6, double h7, double h8, double h9) const {
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

    const double cosh_yGt = cosh(y * gamma * t);
    const double sinh_yGt = sinh(y * gamma * t);
    const double decay = gamma * exp( - gamma * t);

    return c * decay * (
            (SJ1s * cosh_yGt - h1s * sinh_yGt) * sinThetaK2 +
            (SJ1c * cosh_yGt - h1c * sinh_yGt) * cosThetaK2 +
            (SJ2s * cosh_yGt - h2s * sinh_yGt) * sinThetaK2 * cos2ThetaL +
            (SJ2c * cosh_yGt - h2c * sinh_yGt) * cosThetaK2 * cos2ThetaL +
            (SJ3 * cosh_yGt - h3 * sinh_yGt) * sinThetaK2 * sinThetaL2 * cos(2 * phi) +
            (SJ4 * cosh_yGt - h4 * sinh_yGt) * sin2ThetaK * sin2ThetaL * cos(phi) +
            (DJ5 * cosh_yGt - h5 * sinh_yGt) * sin2ThetaK * sinThetaL * cos(phi) +
            (DJ6s * cosh_yGt - h6 * sinh_yGt) * sinThetaK2 * cosThetaL +
            (SJ7 * cosh_yGt - h7 * sinh_yGt) * sin2ThetaK * sinThetaL * sin(phi) +
            (DJ8 * cosh_yGt - h8 * sinh_yGt) * sin2ThetaK * sin2ThetaL * sin(phi) +
            (DJ9 * cosh_yGt - h9 * sinh_yGt) * sinThetaK2 * sinThetaL2 * sin(2 * phi)
    );
}

double timeDependentUntagged::evaluate() const {
    return evaluate_prob(cosThetaL_, cosThetaK_, phi_, t_,gamma_, y_,
                        SJ1s_, SJ1c_, SJ2s_, SJ2c_, SJ3_, SJ4_, DJ5_, DJ6s_, SJ7_, DJ8_, DJ9_,
                        h1s_, h1c_, h2s_, h2c_, h3_, h4_, h5_, h6_, h7_, h8_, h9_);
}

void timeDependentUntagged::doEval(RooFit::EvalContext &ctx) const {
    std::span<const double> costhetalSpan = ctx.at(cosThetaL_);
    std::span<const double> costhetakSpan = ctx.at(cosThetaK_);
    std::span<const double> tSpan = ctx.at(t_);
    std::span<const double> phiSpan = ctx.at(phi_);
    std::span<const double> gammaSpan = ctx.at(gamma_);
    std::span<const double> ySpan = ctx.at(y_);
    std::span<const double> SJ1sSpan = ctx.at(SJ1s_);
    std::span<const double> SJ1cSpan = ctx.at(SJ1c_);
    std::span<const double> SJ2sSpan = ctx.at(SJ2s_);
    std::span<const double> SJ2cSpan = ctx.at(SJ2c_);
    std::span<const double> SJ3Span = ctx.at(SJ3_);
    std::span<const double> SJ4Span = ctx.at(SJ4_);
    std::span<const double> DJ5Span = ctx.at(DJ5_);
    std::span<const double> DJ6sSpan = ctx.at(DJ6s_);
    std::span<const double> SJ7Span = ctx.at(SJ7_);
    std::span<const double> DJ8Span = ctx.at(DJ8_);
    std::span<const double> DJ9Span = ctx.at(DJ9_);
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
    std::size_t n = ctx.output().size();
    for (std::size_t i = 0; i < n; ++i) {
        ctx.output()[i] = evaluate_prob(costhetalSpan.size() > 1 ? costhetalSpan[i] : costhetalSpan[0],
                                        costhetakSpan.size() > 1 ? costhetakSpan[i] : costhetakSpan[0],
                                        phiSpan.size() > 1 ? phiSpan[i] : phiSpan[0],
                                        tSpan.size() > 1 ? tSpan[i] : tSpan[0],
                                        gammaSpan.size() > 1 ? gammaSpan[i] : gammaSpan[0],
                                        ySpan.size() > 1 ? ySpan[i] : ySpan[0],
                                        SJ1sSpan.size() > 1 ? SJ1sSpan[i] : SJ1sSpan[0],
                                        SJ1cSpan.size() > 1 ? SJ1cSpan[i] : SJ1cSpan[0],
                                        SJ2sSpan.size() > 1 ? SJ2sSpan[i] : SJ2sSpan[0],
                                        SJ2cSpan.size() > 1 ? SJ2cSpan[i] : SJ2cSpan[0],
                                        SJ3Span.size() > 1 ? SJ3Span[i] : SJ3Span[0],
                                        SJ4Span.size() > 1 ? SJ4Span[i] : SJ4Span[0],
                                        DJ5Span.size() > 1 ? DJ5Span[i] : DJ5Span[0],
                                        DJ6sSpan.size() > 1 ? DJ6sSpan[i] : DJ6sSpan[0],
                                        SJ7Span.size() > 1 ? SJ7Span[i] : SJ7Span[0],
                                        DJ8Span.size() > 1 ? DJ8Span[i] : DJ8Span[0],
                                        DJ9Span.size() > 1 ? DJ9Span[i] : DJ9Span[0],
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
                                        h9Span.size() > 1 ? h9Span[i] : h9Span[0]);
    }
}

timeDependentTagged::timeDependentTagged(const char *name, const char *title, RooAbsReal& costhetal, RooAbsReal& costhetak, RooAbsReal& phi,
                                         RooAbsReal& t, RooAbsReal& gamma, RooAbsReal& x,
                                         RooAbsReal& DJ1s, RooAbsReal& DJ1c, RooAbsReal& DJ2s, RooAbsReal& DJ2c, RooAbsReal& DJ3,
                                         RooAbsReal& DJ4, RooAbsReal& SJ5, RooAbsReal& SJ6s, RooAbsReal& DJ7, RooAbsReal& SJ8, RooAbsReal& SJ9,
                                         RooAbsReal& s1s, RooAbsReal& s1c, RooAbsReal& s2s, RooAbsReal& s2c, RooAbsReal& s3,
                                         RooAbsReal& s4, RooAbsReal& s5, RooAbsReal& s6, RooAbsReal& s7, RooAbsReal& s8, RooAbsReal& s9):
        RooAbsPdf(name, title),
        cosThetaL_("costhetal", "costhetal", this, costhetal),
        cosThetaK_("costhetak", "costhetak", this, costhetak),
        phi_("phi", "phi", this, phi),
        t_("t", "t", this, t),
        gamma_("gamma", "gamma", this, gamma),
        x_("x", "x", this, x),
        DJ1s_("DJ1s", "DJ1s", this, DJ1s),
        DJ1c_("DJ1c", "DJ1c", this, DJ1c),
        DJ2s_("DJ2s", "DJ2s", this, DJ2s),
        DJ2c_("DJ2c", "DJ2c", this, DJ2c),
        DJ3_("DJ3", "DJ3", this, DJ3),
        DJ4_("DJ4", "DJ4", this, DJ4),
        SJ5_("SJ5", "SJ5", this, SJ5),
        SJ6s_("SJ6s", "SJ6s", this, SJ6s),
        DJ7_("DJ7", "DJ7", this, DJ7),
        SJ8_("SJ8", "SJ8", this, SJ8),
        SJ9_("SJ9", "SJ9", this, SJ9),
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

timeDependentTagged::timeDependentTagged(timeDependentTagged const &other, const char *name):
        RooAbsPdf(other, name),
        cosThetaL_("costhetal", this, other.cosThetaL_),
        cosThetaK_("costhetak", this, other.cosThetaK_),
        phi_("phi", this, other.phi_),
        t_("t", this, other.t_),
        gamma_("gamma", this, other.gamma_),
        x_("x", this, other.x_),
        DJ1s_("DJ1s", this, other.DJ1s_),
        DJ1c_("DJ1c", this, other.DJ1c_),
        DJ2s_("DJ2s", this, other.DJ2s_),
        DJ2c_("DJ2c", this, other.DJ2c_),
        DJ3_("DJ3", this, other.DJ3_),
        DJ4_("DJ4", this, other.DJ4_),
        SJ5_("SJ5", this, other.SJ5_),
        SJ6s_("SJ6s", this, other.SJ6s_),
        DJ7_("DJ7", this, other.DJ7_),
        SJ8_("DJ8", this, other.SJ8_),
        SJ9_("SJ9", this, other.SJ9_),
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
        s9_("s9", this, other.s9_){}

double timeDependentTagged::evaluate_prob(double cosThetaL, double cosThetaK, double phi, double t, double gamma, double x,
                                            double DJ1s, double DJ1c, double DJ2s, double DJ2c, double DJ3, double DJ4, double SJ5, double SJ6s, double DJ7, double SJ8, double SJ9,
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

    const double cos_xGt = cos(x * gamma * t);
    const double sin_xGt = sin(x * gamma * t);
    const double decay = gamma * exp( - gamma * t);

    return c * decay * (
                               (DJ1s * cos_xGt - s1s * sin_xGt) * sinThetaK2 +
                               (DJ1c * cos_xGt - s1c * sin_xGt) * cosThetaK2 +
                               (DJ2s * cos_xGt - s2s * sin_xGt) * sinThetaK2 * cos2ThetaL +
                               (DJ2c * cos_xGt - s2c * sin_xGt) * cosThetaK2 * cos2ThetaL +
                               (DJ3 * cos_xGt - s3 * sin_xGt) * sinThetaK2 * sinThetaL2 * cos(2 * phi) +
                               (DJ4 * cos_xGt - s4 * sin_xGt) * sin2ThetaK * sin2ThetaL * cos(phi) +
                               (SJ5 * cos_xGt - s5 * sin_xGt) * sin2ThetaK * sinThetaL * cos(phi) +
                               (SJ6s * cos_xGt - s6 * sin_xGt) * sinThetaK2 * cosThetaL +
                               (DJ7 * cos_xGt - s7 * sin_xGt) * sin2ThetaK * sinThetaL * sin(phi) +
                               (SJ8 * cos_xGt - s8 * sin_xGt) * sin2ThetaK * sin2ThetaL * sin(phi) +
                               (SJ9 * cos_xGt - s9 * sin_xGt) * sinThetaK2 * sinThetaL2 * sin(2 * phi)
    );
}

double timeDependentTagged::evaluate() const {
    return evaluate_prob(cosThetaL_, cosThetaK_, phi_, t_,gamma_, x_,
                         DJ1s_, DJ1c_, DJ2s_, DJ2c_, DJ3_, DJ4_, SJ5_, SJ6s_, DJ7_, SJ8_, SJ9_,
                         s1s_, s1c_, s2s_, s2c_, s3_, s4_, s5_, s6_, s7_, s8_, s9_);
}

void timeDependentTagged::doEval(RooFit::EvalContext &ctx) const {
    std::span<const double> cosThetaLSpan = ctx.at(cosThetaL_);
    std::span<const double> cosThetaKSpan = ctx.at(cosThetaK_);
    std::span<const double> tSpan = ctx.at(t_);
    std::span<const double> phiSpan = ctx.at(phi_);
    std::span<const double> gammaSpan = ctx.at(gamma_);
    std::span<const double> xSpan = ctx.at(x_);
    std::span<const double> DJ1sSpan = ctx.at(DJ1s_);
    std::span<const double> DJ1cSpan = ctx.at(DJ1c_);
    std::span<const double> DJ2sSpan = ctx.at(DJ2s_);
    std::span<const double> DJ2cSpan = ctx.at(DJ2c_);
    std::span<const double> DJ3Span = ctx.at(DJ3_);
    std::span<const double> DJ4Span = ctx.at(DJ4_);
    std::span<const double> SJ5Span = ctx.at(SJ5_);
    std::span<const double> SJ6sSpan = ctx.at(SJ6s_);
    std::span<const double> DJ7Span = ctx.at(DJ7_);
    std::span<const double> SJ8Span = ctx.at(SJ8_);
    std::span<const double> SJ9Span = ctx.at(SJ9_);
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
                                        gammaSpan.size() > 1 ? gammaSpan[i] : gammaSpan[0],
                                        xSpan.size() > 1 ? xSpan[i] : xSpan[0],
                                        DJ1sSpan.size() > 1 ? DJ1sSpan[i] : DJ1sSpan[0],
                                        DJ1cSpan.size() > 1 ? DJ1cSpan[i] : DJ1cSpan[0],
                                        DJ2sSpan.size() > 1 ? DJ2sSpan[i] : DJ2sSpan[0],
                                        DJ2cSpan.size() > 1 ? DJ2cSpan[i] : DJ2cSpan[0],
                                        DJ3Span.size() > 1 ? DJ3Span[i] : DJ3Span[0],
                                        DJ4Span.size() > 1 ? DJ4Span[i] : DJ4Span[0],
                                        SJ5Span.size() > 1 ? SJ5Span[i] : SJ5Span[0],
                                        SJ6sSpan.size() > 1 ? SJ6sSpan[i] : SJ6sSpan[0],
                                        DJ7Span.size() > 1 ? DJ7Span[i] : DJ7Span[0],
                                        SJ8Span.size() > 1 ? SJ8Span[i] : SJ8Span[0],
                                        SJ9Span.size() > 1 ? SJ9Span[i] : SJ9Span[0],
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
