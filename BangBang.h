#pragma once 

#include <cmath>
#include <limits>
//#include <algorithm>

/// <summary>
/// Computes 1-switch bang-bang jerk trajectories between arbitrary
/// initial and final (v, a) boundary conditions.
///
/// The trajectory consists of two phases:
///
///   Phase 1:  j =  j0   for 0 ≤ t ≤ t1
///   Phase 2:  j = -j0   for t1 ≤ t ≤ T
///
/// Two solving modes are supported:
///
///   solveForMinJerk(v0, a0, vf, af, T)
///       Finds the minimal jerk magnitude j0 and switch time t1
///       that satisfy the boundary conditions exactly over a fixed
///       duration T.
///
///   solveForMinTime(v0, a0, vf, af, Jmax)
///       Finds the minimum feasible duration T using maximal jerk
///       magnitude |j0| = Jmax, solving for both T and t1.
///
/// After solving, evaluate(t, deriv) returns position, velocity,
/// acceleration, or jerk at time t (clamped to [0, T]).
///
/// All math is analytic (no integration), and the solver enforces the
/// boundary conditions exactly within numerical precision.
///
/// </summary>
class BangBang
{
public:
    BangBang() :
        v0_(0), a0_(0), vf_(0), af_(0),
        T_(0), t1_(0), j0_(0), j1_(0),
        valid_(false)
    {
    }

    // -------------------------------------------------------------
    //  solveForMinJerk
    //
    //  Given initial state (v0, a0),
    //        final state (vf, af),
    //        and duration T,
    //  find the symmetric bang-bang jerk profile that satisfies
    //        a(t1) = a1,
    //        a(T) = af,
    //        v(T) = vf,
    //  using minimal jerk magnitude j0.
    // -------------------------------------------------------------
    void solveForMinJerk(double v0, double a0,
        double vf, double af,
        double T)
    {
        v0_ = v0;  a0_ = a0;
        vf_ = vf;  af_ = af;
        T_ = T;
        valid_ = false;
        j0_ = j1_ = 0.0;
        t1_ = 0.0;

        const double eps = 1e-6;

        if (T_ <= 0)
            return;

        // Quadratic coefficients for t1:
        // A t1^2 + B t1 + C = 0
        double A = (a0_ - af_);
        double B = 2.0 * T_ * af_ + 2.0 * (v0_ - vf_);
        double C = -0.5 * T_ * T_ * (a0_ + af_) - T_ * (v0_ - vf_);

        double t1_candidate = std::numeric_limits<double>::quiet_NaN();

        if (std::abs(A) > eps)
        {
            double D = B * B - 4.0 * A * C;
            if (D < 0.0)
            {
                valid_ = false;
                return;
            }

            double sqrtD = std::sqrt(D);
            double t1a = (-B + sqrtD) / (2.0 * A);
            double t1b = (-B - sqrtD) / (2.0 * A);

            auto accept = [&](double t)->bool
                {
                    return (t >= -eps && t <= T_ + eps);
                };

            if (accept(t1a)) t1_candidate = t1a;
            else if (accept(t1b)) t1_candidate = t1b;
            else
            {
                valid_ = false;
                return;
            }
        }
        else
        {
            // Linear case: B t1 + C = 0
            if (std::abs(B) < eps)
            {
                if (std::abs(C) < eps)
                    t1_candidate = 0.0;   // degenerate family, pick something
                else
                {
                    valid_ = false;
                    return;
                }
            }
            else
            {
                double t1 = -C / B;
                if (t1 < -eps || t1 > T_ + eps)
                {
                    valid_ = false;
                    return;
                }
                t1_candidate = t1;
            }
        }

        t1_ = std::clamp(t1_candidate, 0.0, T_);

        // Compute j0.
        // In general: af = a0 + j0*(2*t1 - T)
        // j0 = (af - a0) / (2*t1 - T)
        double denom = 2.0 * t1_ - T_;
        if (std::abs(denom) >= eps)
        {
            j0_ = (af_ - a0_) / denom;
        }
        else
        {
            // a0 == af and 2*t1 ~= T:
            // acceleration boundary is automatically satisfied; solve j0
            // from the velocity boundary:
            //
            // vf = v0 + a0*T - 0.5*j0*T*T + 2*j0*T*t1 - j0*t1*t1
            // => 0 = v0 + a0*T - vf + j0 * (-0.5*T*T + 2*T*t1 - t1*t1)
            //
            double coeff = -0.5 * T_ * T_ + 2.0 * T_ * t1_ - t1_ * t1_;
            // don't worry for coeff ~= 0, because t1 ~= T/2, and in that case coeff ~= T^2/4
            j0_ = (vf_ - v0_ - a0_ * T_) / coeff;
        }

        j1_ = -j0_;
        valid_ = true;
    }

    void solveForMinTime(double v0, double a0,
        double vf, double af,
        double Jmax)
    {
        v0_ = v0;  a0_ = a0;
        vf_ = vf;  af_ = af;

        T_ = 0.0;
        t1_ = 0.0;
        j0_ = 0.0;
        j1_ = 0.0;
        valid_ = false;

        const double eps = 1e-12;

        if (Jmax <= 0.0)
            return;

        // special case where no motion is needed
        if (v0 == 0 && a0 == 0 && vf == 0 && af == 0) {
            valid_ = true;
            return;
        }

        double bestT = std::numeric_limits<double>::infinity();
        double best_t1 = 0.0;
        double best_j0 = 0.0;

        // Try both jerk signs: +Jmax and -Jmax
        for (int s = -1; s <= 1; s += 2)
        {
            double j0 = s * Jmax;

            // Quadratic in t1:
            //   A t1^2 + B t1 + C = 0
            double A = j0;
            double B = 2.0 * a0_;
            double C = (a0_ * a0_ - af_ * af_) / (2.0 * j0) + (v0_ - vf_);

            double D = B * B - 4.0 * A * C;
            if (D < -eps)
                continue;
            if (D < 0.0)
                D = 0.0;

            double sqrtD = std::sqrt(D);

            auto try_t1 = [&](double t1_val)
                {
                    if (t1_val < -eps)
                        return;

                    double Tval = (2.0 * j0 * t1_val + a0_ - af_) / j0;
                    if (Tval <= eps)
                        return;

                    if (t1_val - Tval > eps)
                        return; // t1 must be <= T

                    // Optional sanity check (can be disabled in release):
                    {
                        double a1 = a0_ + j0 * t1_val;
                        double tau = Tval - t1_val;
                        double af_calc = a1 - j0 * tau;

                        double v1 = v0_ + a0_ * t1_val + 0.5 * j0 * t1_val * t1_val;
                        double vf_calc = v1 + a1 * tau - 0.5 * j0 * tau * tau;

                        if (std::abs(af_calc - af_) > 1e-6 ||
                            std::abs(vf_calc - vf_) > 1e-6)
                            return;
                    }

                    if (Tval < bestT)
                    {
                        bestT = Tval;
                        best_t1 = t1_val;
                        best_j0 = j0;
                    }
                };

            // Two possible roots
            if (std::abs(A) > eps)
            {
                double t1a = (-B + sqrtD) / (2.0 * A);
                double t1b = (-B - sqrtD) / (2.0 * A);
                try_t1(t1a);
                try_t1(t1b);
            }
            else
            {
                // Degenerate (A ~ 0): linear B t1 + C = 0
                if (std::abs(B) > eps)
                {
                    double t1_lin = -C / B;
                    try_t1(t1_lin);
                }
                // else: no t1 dependence (unlikely here); ignore
            }
        }

        if (!std::isfinite(bestT))
        {
            valid_ = false;
            return;
        }

        // Store final solution
        valid_ = true;
        T_ = bestT;
        t1_ = best_t1;
        j0_ = best_j0;
        j1_ = -best_j0;
    }

    // -------------------------------------------------------------
    double switchTime() const
    {
        return valid_ ? t1_ : std::numeric_limits<double>::quiet_NaN();
    }
    double endTime() const { return T_; }
    // -------------------------------------------------------------
    //   evaluate:
    //     deriv = 0 → position x(t)
    //     deriv = 1 → velocity v(t)
    //     deriv = 2 → acceleration a(t)
    //     deriv = 3 → jerk j(t)
    // -------------------------------------------------------------
    double evaluate(double t, int deriv = 0) const
    {
        if (!valid_)
            return std::numeric_limits<double>::quiet_NaN();

        t = std::clamp(t, 0.0, T_);
        int d = std::clamp(deriv, 0, 3);

        double tt1 = t1_;
        double j0 = j0_;
        double j1 = j1_;

        // Phase 1 values:
        double a1 = a0_ + j0 * tt1;
        double v1 = v0_ + a0_ * tt1 + 0.5 * j0 * tt1 * tt1;
        double x1 = v0_ * tt1
            + 0.5 * a0_ * tt1 * tt1
            + (1.0 / 6.0) * j0 * tt1 * tt1 * tt1;

        // ---------------------------------------------------------
        // Phase 1 (t <= t1)
        // ---------------------------------------------------------
        if (t <= tt1)
        {
            switch (d)
            {
            case 3: return j0;
            case 2: return a0_ + j0 * t;
            case 1: return v0_ + a0_ * t + 0.5 * j0 * t * t;
            default:
                return v0_ * t
                    + 0.5 * a0_ * t * t
                    + (1.0 / 6.0) * j0 * t * t * t;
            }
        }

        // ---------------------------------------------------------
        // Phase 2 (t > t1)
        // ---------------------------------------------------------
        double dt = t - tt1;

        switch (d)
        {
        case 3: return j1;
        case 2: return a1 + j1 * dt;
        case 1: return v1 + a1 * dt + 0.5 * j1 * dt * dt;
        default:
            return x1
                + v1 * dt
                + 0.5 * a1 * dt * dt
                + (1.0 / 6.0) * j1 * dt * dt * dt;
        }
    }

protected:
    double v0_, a0_, vf_, af_;
    double T_;
    double t1_, j0_, j1_;
    bool valid_;
    friend class BangBangLimitedAccel;
};

class BangBangLimitedAccel {
public:
    void solveForMinTime(double v0, double a0,
        double vf, double af,
        double aMax, double Jmax)
    {
        if (std::abs(a0) > aMax) throw std::invalid_argument("|a0| > aMax");
        if (std::abs(af) > aMax) throw std::invalid_argument("|af| > aMax");

        v0_ = v0; a0_ = a0;
        vf_ = vf; af_ = af;
        aMax_ = aMax; jMax_ = Jmax;

        // --------------------------------------------------------
        // 1. Solve unconstrained jerk bang-bang
        // --------------------------------------------------------
        BangBang bb;
        bb.solveForMinTime(v0, a0, vf, af, Jmax);
        if (!bb.valid_) { valid_ = false; return; }

        double t1_bb = bb.switchTime();
        double T_bb = bb.endTime();

        auto aBB = [&](double t) { return bb.evaluate(t, 2); };

        // peak accel of unconstrained profile
        double a_peak = std::abs(aBB(0));
        a_peak = std::max(a_peak, std::abs(aBB(t1_bb)));
        a_peak = std::max(a_peak, std::abs(aBB(T_bb)));

        // --------------------------------------------------------
        // 2. Trivial case: no clipping needed
        // --------------------------------------------------------
        if (std::abs(a_peak) <= aMax + 1e-12)
        {
            isTrivial_ = true;
            valid_ = true;

            // Extract the values needed to evaluate later
            t1_ = t1_bb;
            T_ = T_bb;

            // jerks of the unconstrained solution
            double j0 = bb.j0_;
            double j1 = bb.j1_;

            j0_ = j0;
            j1_ = j1;

            t2_.reset();   // no 2nd switch
            return;
        }

        // --------------------------------------------------------
        // 3. S-lobe clipping needed
        // --------------------------------------------------------
        isTrivial_ = false;
        bool peak_accel_positive = aBB(t1_bb) > 0;
        a_peak = (peak_accel_positive ? aMax_ : -aMax_);
        double j_head = bb.evaluate(0, 3);
        double j_tail = -j_head;

        // SEG1: jerk to a_max
        t1_ = (a_peak - a0_) / j_head;

        // SEG3: jerk-down to a_f
        double t3 = (af_ - a_peak) / j_tail;

        // Required total Δv
        double dv_req = vf_ - v0_;

        // Δv from SEG1
        double dv1 = a0_ * t1_ + 0.5 * j_head * t1_ * t1_;

        // Δv from SEG3
        double dv3 = a_peak * t3 + 0.5 * j_tail * t3 * t3;

        // Flat accel segment: dv2 = aMax * t_flat
        double dv2 = dv_req - (dv1 + dv3);
        double t_flat = dv2 / a_peak;
        if (t_flat < 0) t_flat = 0.0;

        t2_ = t1_ + t_flat;
        T_ = *t2_ + t3;

        // Store resulting jerk pattern
        j0_ = j_head;     // segment 1
        j1_ = 0.0;      // flat acceleration
        j2_ = j_tail;    // segment 3
        t_flat_ = t_flat;
        t3_ = t3;

        valid_ = true;
    }

    std::vector<double> switchingTimes() const {
        if (t2_) return std::vector<double>{ t1_, *t2_ };
        return std::vector<double>{ t1_ };
    }

    double endTime() const { return T_; }

    // ----------------------------------------------------------
    // EVALUATE
    // ----------------------------------------------------------
    double evaluate(double t, int deriv = 0) const
    {
        if (!valid_) return std::numeric_limits<double>::quiet_NaN();
        t = std::clamp(t, 0.0, T_);
        int d = std::clamp(deriv, 0, 3);

        // ----------------------------------------
        // CASE A: trivial – this is plain BangBang
        // ----------------------------------------
        if (isTrivial_)
        {
            double tt1 = t1_;
            double j0 = j0_;
            double j1 = j1_;

            // Precompute end of phase 1
            double a1 = a0_ + j0 * tt1;
            double v1 = v0_ + a0_ * tt1 + 0.5 * j0 * tt1 * tt1;
            double x1 = v0_ * tt1
                + 0.5 * a0_ * tt1 * tt1
                + (1.0 / 6.0) * j0 * tt1 * tt1 * tt1;

            if (t <= tt1)
            {
                if (d == 3) return j0;
                if (d == 2) return a0_ + j0 * t;
                if (d == 1) return v0_ + a0_ * t + 0.5 * j0 * t * t;
                return v0_ * t + 0.5 * a0_ * t * t + (1.0 / 6.0) * j0 * t * t * t;
            }

            // second phase
            double dt = t - tt1;
            if (d == 3) return j1;
            if (d == 2) return a1 + j1 * dt;
            if (d == 1) return v1 + a1 * dt + 0.5 * j1 * dt * dt;
            return x1 + v1 * dt + 0.5 * a1 * dt * dt + (1.0 / 6.0) * j1 * dt * dt * dt;
        }

        // ----------------------------------------
        // CASE B: clipped S-lobe (3 segments)
        // ----------------------------------------
        double t1 = t1_;
        double t2 = *t2_;

        // SEG1: jerk up
        if (t <= t1)
        {
            double tt = t;
            if (d == 3) return j0_;
            if (d == 2) return a0_ + j0_ * tt;
            if (d == 1) return v0_ + a0_ * tt + 0.5 * j0_ * tt * tt;
            return v0_ * tt + 0.5 * a0_ * tt * tt + (1.0 / 6.0) * j0_ * tt * tt * tt;
        }

        // Precompute end of SEG1
        double a1 = a0_ + j0_ * t1;
        double v1 = v0_ + a0_ * t1 + 0.5 * j0_ * t1 * t1;
        double x1 = v0_ * t1 + 0.5 * a0_ * t1 * t1 + (1.0 / 6.0) * j0_ * t1 * t1 * t1;

        // SEG2: flat accel
        if (t <= t2)
        {
            double dt = t - t1;
            if (d == 3) return j1_; // which is 0
            if (d == 2) return a1;
            if (d == 1) return v1 + a1 * dt;
            return x1 + v1 * dt + 0.5 * a1 * dt * dt;
        }

        // Precompute end of SEG2
        double dt12 = t2 - t1;
        double a2 = a1;
        double v2 = v1 + a1 * dt12;
        double x2 = x1 + v1 * dt12 + 0.5 * a1 * dt12 * dt12;

        // SEG3: jerk down
        double t3 = T_ - t2;
        double dt3 = t - t2;

        if (d == 3) return j2_;
        if (d == 2) return a2 + j2_ * dt3;
        if (d == 1) return v2 + a2 * dt3 +0.5 * j2_ * dt3 * dt3;

        return x2 + v2 * dt3 + 0.5 * a2 * dt3 * dt3 + (1.0 / 6.0) * j2_ * dt3 * dt3 * dt3;
    }

private:
    // inputs
    double v0_, a0_, vf_, af_;
    double aMax_, jMax_;

    // outputs
    double t1_;
    std::optional<double> t2_;
    double T_;

    // jerk values
    double j0_, j1_, j2_;

    double t_flat_ = 0;
    double t3_ = 0;

    bool valid_ = false;
    bool isTrivial_ = false;
};
