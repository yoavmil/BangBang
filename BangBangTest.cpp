#include "LineClusterPch.h"
#pragma optimize("", off)
#include "BangBang.h"

using std::vector;
TEST(BangBang, min_time)
{
    struct Test {
        double v0, a0, vf, af, J;
        double expected_switch_time;
        double expected_T;
    };

    std::vector<Test> tests;

    // ---------------------------------------------------------
    // Test 1: Zero motion (no movement required)
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 0;
        test.vf = 0;
        test.af = 0;
        test.J = 10;
        test.expected_switch_time = 0;
        test.expected_T = 0;
    }

    // ---------------------------------------------------------
    // Test 2: Acceleration reversal: +1 to -1 with J=2
    // Expected: t1=0.5, T=1.0
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = +1;
        test.vf = 0;
        test.af = -1;
        test.J = 2;
        test.expected_switch_time = 1;
        test.expected_T = 1.0;
    }

    // ---------------------------------------------------------
    // Test 3: Decelerate from v0=10 to v=0 using J=20
    // Expected:
    //   t1 = sqrt(10/20) = sqrt(0.5)
    //   T  = 2*t1
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 10;
        test.a0 = 0;
        test.vf = 0;
        test.af = 0;
        test.J = 20;
        test.expected_switch_time = std::sqrt(10.0 / 20.0);
        test.expected_T = 2.0 * std::sqrt(10.0 / 20.0);
    }

    // ---------------------------------------------------------
    // Test 4: a:0 → +2 → 0 with J=4 (symmetric)
    // We know vf = 1.0 and T = 1.0, t1 = 0.5
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 0;
        test.vf = 1;   // computed analytically
        test.af = 0;
        test.J = 4;
        test.expected_switch_time = 0.5;
        test.expected_T = 1.0;
    }

    // ---------------------------------------------------------
    // Test 5: Decelerate from v0=4 to 0 with J=8
    // Expected:
    //   t1 = sqrt(4/8)
    //   T  = 2*t1
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 4;
        test.a0 = 0;
        test.vf = 0;
        test.af = 0;
        test.J = 8;
        test.expected_switch_time = std::sqrt(4.0 / 8.0);
        test.expected_T = 2.0 * std::sqrt(4.0 / 8.0);
    }

    {
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 1;
        test.vf = 0;
        test.af = 1;
        test.J = 2;
        test.expected_switch_time = 1.0;
        test.expected_T = 2.0;
    }

    { // measured data at ACS
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 0;
        test.vf = 500;
        test.af = 0;
        test.J = 100'000;
        test.expected_switch_time = 0.0707;
        test.expected_T = 0.0707 * 2;
    }

    for (long i = 0; i < tests.size(); i++) {
        auto& t = tests[i];

        BangBang bb;
        bb.solveForMinTime(t.v0, t.a0, t.vf, t.af, t.J);

        EXPECT_NEAR(bb.switchTime(), t.expected_switch_time, 1e-1) << i;
        EXPECT_NEAR(bb.endTime(), t.expected_T, 1e-1) << i;
    }
}

TEST(BangBang, min_time_limited_accel)
{
    struct Test {
        double v0, a0, vf, af, J, A;
        vector<double> expected_switch_times;
        double expected_T;
    };

    std::vector<Test> tests;

    // first test cases are copied from `min_time` test, just to see there is no regression
    // ---------------------------------------------------------
    // Test 1: Zero motion (no movement required)
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 0;
        test.vf = 0;
        test.af = 0;
        test.A = 1;
        test.J = 10;
        test.expected_switch_times = { 0 };
        test.expected_T = 0;
    }

    // ---------------------------------------------------------
    // Test 2: Acceleration reversal: +1 to -1 with J=2
    // Expected: t1=0.5, T=1.0
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = +1;
        test.vf = 0;
        test.af = -1;
        test.A = 1;
        test.J = 2;
        test.expected_switch_times = { 1 };
        test.expected_T = 1.0;
    }

    // ---------------------------------------------------------
    // Test 3: Decelerate from v0=10 to v=0 using J=20
    // Expected:
    //   t1 = sqrt(10/20) = sqrt(0.5)
    //   T  = 2*t1
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 10;
        test.a0 = 0;
        test.vf = 0;
        test.af = 0;
        test.A = 20; // enough acceleration
        test.J = 20;
        test.expected_switch_times = { std::sqrt(10.0 / 20.0) };
        test.expected_T = 2.0 * std::sqrt(10.0 / 20.0);
    }

    // ---------------------------------------------------------
    // Test 4: a:0 → +2 → 0 with J=4 (symmetric)
    // We know vf = 1.0 and T = 1.0, t1 = 0.5
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 0;
        test.vf = 1;   // computed analytically
        test.af = 0;
        test.A = 2; // just enough
        test.J = 4;
        test.expected_switch_times = { 0.5 };
        test.expected_T = 1.0;
    }

    // ---------------------------------------------------------
    // Test 5: Decelerate from v0=4 to 0 with J=8
    // Expected:
    //   t1 = sqrt(4/8)
    //   T  = 2*t1
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 4;
        test.a0 = 0;
        test.vf = 0;
        test.af = 0;
        test.A = 20;
        test.J = 8;
        test.expected_switch_times = { std::sqrt(4.0 / 8.0) };
        test.expected_T = 2.0 * std::sqrt(4.0 / 8.0);
    }

    { // measured data at ACS, enough acceleration
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 0;
        test.vf = 500;
        test.af = 0;
        test.A = 10'000;
        test.J = 100'000;
        test.expected_switch_times = { 0.0707 };
        test.expected_T = 0.0707 * 2;
    }

    { // measured data at ACS, not enough acceleration
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 0;
        test.vf = 500;
        test.af = 0;
        test.A = 5'000;
        test.J = 100'000;
        test.expected_switch_times = { 0.05, 0.1 }; // 0.05 = 5000/100000
        test.expected_T = 0.1 + 0.05;
    }

    { // decel measured data at ACS, the opposite of previos test
        Test& test = tests.emplace_back();
        test.v0 = 500;
        test.a0 = 0;
        test.vf = 0;
        test.af = 0;
        test.A = 5'000;
        test.J = 100'000;
        test.expected_switch_times = { 0.05, 0.1 }; // 0.05 = 5000/100000
        test.expected_T = 0.1 + 0.05;
    }

    { // v0>0, a0>0, just one half of the acceleration ramp
        Test& test = tests.emplace_back();
        test.v0 = 5;
        test.a0 = 14.135; // measured at ACS
        test.vf = 10;
        test.af = 0;
        test.A = 20;
        test.J = 20;
        test.expected_switch_times = { 0 }; // no switchs, just acceleration all along
        test.expected_T = 0.707;
    }

    { // v0>0, a0>0, 1/4 of the acceleration ramp
        Test& test = tests.emplace_back();
        test.v0 = 2.5;
        test.a0 = 10; // measured at ACS
        test.vf = 10;
        test.af = 0;
        test.A = 20;
        test.J = 20;
        test.expected_switch_times = { 0.207 };
        test.expected_T = 0.914;
    }

    { // reverse motion
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 0; // measured at ACS
        test.vf = -10;
        test.af = 0;
        test.A = 10;
        test.J = 20;
        test.expected_switch_times = { 0.5, 1.0 };
        test.expected_T = 1.5;
    }

    { // reverse motion, v0<0, a0<0, vf<0, af<0
        Test& test = tests.emplace_back();
        test.v0 = -0.637;
        test.a0 = -5.06;
        test.vf = -9.387;
        test.af = -4.94;
        test.A = 10;
        test.J = 20;
        test.expected_switch_times = { 0.247, 0.747 };
        test.expected_T = 1.0;
    }

    { // reverse through zero, limited accel
        Test& test = tests.emplace_back();
        test.v0 = 5;
        test.a0 = 0;
        test.vf = -5;
        test.af = 0;
        test.A = 2;      // aMax
        test.J = 10;     // Jmax
        test.expected_switch_times = { 0.2, 5.0 };
        test.expected_T = 5.2;
    }

    { // asymmetric braking from negative speed, nonzero initial accel
        Test& test = tests.emplace_back();
        test.v0 = -4;    // start moving backward
        test.a0 = +1;    // already accelerating forward
        test.vf = 0;     // come to a stop
        test.af = 0;     // with zero acceleration
        test.A = 2;     // max accel
        test.J = 10;    // max jerk

        // Expected: t1 = 0.1, t2 = 1.925, T = 2.125
        test.expected_switch_times = { 0.1, 1.925 };
        test.expected_T = 2.125;
    }

    { // asymmetric braking from negative speed, starting at aMax
        Test& test = tests.emplace_back();
        test.v0 = -4;
        test.a0 = +1;    // already at accel limit
        test.vf = 0;
        test.af = 0;
        test.A = 1;     // lower accel limit
        test.J = 10;    // jerk limit

        // Expected: t1 = 0, t2 = 3.95, T = 4.05
        test.expected_switch_times = { 0.0, 3.95 };
        test.expected_T = 4.05;
    }

    for (long i = 0; i < tests.size(); i++) {
        auto& t = tests[i];

        BangBangLimitedAccel bb;
        bb.solveForMinTime(t.v0, t.a0, t.vf, t.af, t.A, t.J);
        auto switchs = bb.switchingTimes();
        ASSERT_EQ(switchs.size(), t.expected_switch_times.size()) << i;
        for (int st_idx = 0; st_idx < switchs.size(); st_idx++) {
            EXPECT_NEAR(switchs[st_idx], t.expected_switch_times[st_idx], 1e-3) << i << ":" << st_idx;
            if (switchs.size() == 2) {
                auto acc = bb.evaluate(switchs[st_idx], 2);
                ASSERT_NEAR(std::abs(acc), t.A, 1e-5) << i << ":" << st_idx;
            }
        }
        EXPECT_NEAR(bb.endTime(), t.expected_T, 1e-3) << i;
    }
}

TEST(BangBang, min_jerk)
{
    struct Test {
        double v0, a0, vf, af, T;
        double expected_t1;
        double expected_j0;
    };

    std::vector<Test> tests;

    // ---------------------------------------------------------
    // Test 1: trivial motion — no change, no jerk needed
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 0;
        test.vf = 0;
        test.af = 0;
        test.T = 1.0;
        test.expected_t1 = 0.0;
        test.expected_j0 = 0.0;
    }

    // ---------------------------------------------------------
    // Test 2: symmetric stop with fixed time
    //
    // Case: start with v0 = 4, stop at vf = 0 in T = 2
    // Using symmetric bang-bang jerk:
    //   t1 = T/2 = 1
    //   a1 = a0 + j0*t1 = j0
    // 
    //   v(T) = 0 ==> j0 = -v0/(0.5*T*T - T*t1 + 0.5*t1*t1)
    // But for t1 = T/2, the denominator simplifies, giving:
    //   j0 = - 2*v0 / T^2 = -2*4/4 = -2
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 4;
        test.a0 = 0;
        test.vf = 0;
        test.af = 0;
        test.T = 2.0;
        test.expected_t1 = 1.0;
        test.expected_j0 = -4.0;
    }

    // ---------------------------------------------------------
    // Test 4: nontrivial symmetric case
    //
    // v0 = 6, a0 = 0, vf = 0, af = 0, T = 3
    //
    // For symmetric bang-bang jerk (t1=T/2=1.5):
    //
    // Velocity constraint gives:
    //   j0 = -4*v0 / T^2 = -4*6/9 = -24/9 = -2.6666666667
    //
    // ---------------------------------------------------------
    {
        Test& test = tests.emplace_back();
        test.v0 = 6;
        test.a0 = 0;
        test.vf = 0;
        test.af = 0;
        test.T = 3.0;
        test.expected_t1 = 1.5;
        test.expected_j0 = -2.6666666666666667;
    }

    // this test data failed at real data, so I added it. 
    // it was numerically instable because of the very small a0 value
    {
        Test& test = tests.emplace_back();
        test.v0 = -188.75957428404541;
        test.a0 = -5.1781495757907692e-12;
        test.vf = 0;
        test.af = 0;
        test.T = 0.12649110640673519;
        test.expected_t1 = test.T / 2; // more or less symmetric
        test.expected_j0 = -4 * test.v0 / (test.T * test.T);
    }

    // if a0 or af is very small, expect t1 ~= T/2
    {
        Test& test = tests.emplace_back();
        test.v0 = 0;
        test.a0 = 0;
        test.vf = 4;
        test.af = 1e-6;
        test.T = 2;
        test.expected_t1 = test.T / 2; // more or less symmetric
        test.expected_j0 = 4 * test.vf / (test.T * test.T);
    }
   
    for (long i = 0; i < tests.size(); i++) {
        auto& t = tests[i];

        BangBang bb;
        bb.solveForMinJerk(t.v0, t.a0, t.vf, t.af, t.T);

        EXPECT_NEAR(bb.switchTime(), t.expected_t1, 1e-5) << i;
        EXPECT_NEAR(bb.evaluate(0, 3), t.expected_j0, 1e-5) << i;

        // Optional: verify boundary conditions
        EXPECT_NEAR(bb.evaluate(t.T, 1), t.vf, 1e-5) << i; // velocity(T)
        EXPECT_NEAR(bb.evaluate(t.T, 2), t.af, 1e-5) << i; // accel(T)
    }
}
