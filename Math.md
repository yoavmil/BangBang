## Mathematical Model – `solveForMinTime`

The goal of `solveForMinTime` is to find the **shortest possible trajectory duration**
that satisfies given boundary conditions under a **maximum jerk constraint**.

### Boundary Conditions

Given:

- initial velocity and acceleration:  
  $v(0) = v_0,\; a(0) = a_0$

- final velocity and acceleration:  
  $v(T) = v_f,\; a(T) = a_f$

- maximum jerk magnitude:  
  $|j(t)| \le J_{\max}$

The total trajectory duration $T$ is **unknown** and must be minimized.

---

### Jerk Profile Assumption

The jerk is assumed to be piecewise constant with a single switching time:

- $j(t) = +j_0 \quad \text{for } 0 \le t \le t_1$
- $j(t) = -j_0 \quad \text{for } t_1 < t \le T$

where $|j_0| = J_{\max}$.

Both $t_1$ and $T$ are unknowns.

---

### Kinematic Equations

Integrating jerk analytically gives the following expressions.

#### Acceleration

For the first phase:  
$a(t) = a_0 + j_0 t$

Acceleration at the switching time:  
$a_1 = a_0 + j_0 t_1$

For the second phase:  
$a(t) = a_1 - j_0 (t - t_1)$

---

#### Velocity

First phase:  
$v(t) = v_0 + a_0 t + \frac{1}{2} j_0 t^2$

Velocity at the switching time:  
$v_1 = v_0 + a_0 t_1 + \frac{1}{2} j_0 t_1^2$

Second phase:  
$v(t) = v_1 + a_1 (t - t_1) - \frac{1}{2} j_0 (t - t_1)^2$

---

### Boundary Condition Equations

Applying the final acceleration constraint:

$a(T) = a_f = a_1 - j_0 (T - t_1)$

Substituting $a_1$ gives:

$a_f = a_0 + j_0 (2 t_1 - T)$

This provides a linear relation between $t_1$ and $T$.

---

Applying the final velocity constraint:

$v(T) = v_f$

Substituting the expressions above yields:

$v_f = v_0 + a_0 T
      + \frac{1}{2} j_0 t_1^2
      + a_1 (T - t_1)
      - \frac{1}{2} j_0 (T - t_1)^2$

After simplification, this results in a **quadratic equation in $t_1$**.

---

### Solving Strategy

1. Fix the jerk magnitude:  
   $j_0 = \pm J_{\max}$

2. For each jerk sign:
   - solve the quadratic equation for valid $t_1$
   - compute the corresponding $T$
   - reject non-physical solutions

3. Select the solution with the **minimum positive duration $T$**.

---

### Notes

- The solution is fully analytic.
- No numerical integration is used.
- Degenerate and symmetric cases are handled explicitly in the implementation.

## Mathematical Model – `solveForMinJerk` (Fixed Time)

This solver uses the same bang-bang jerk model and kinematic equations
described in `solveForMinTime`.

The key difference is that the **total duration is fixed**, and the solver
finds the **minimum jerk magnitude** required to satisfy the boundary conditions.

---

### Problem Definition

Given:

- initial state:  
  $v(0) = v_0,\; a(0) = a_0$

- final state:  
  $v(T) = v_f,\; a(T) = a_f$

- fixed duration:  
  $T$ is known

Unknowns:

- jerk magnitude $j_0$
- switching time $t_1$

The jerk profile is identical to the minimum-time case:

- $j(t) = +j_0 \quad 0 \le t \le t_1$
- $j(t) = -j_0 \quad t_1 < t \le T$

---

### Acceleration Constraint

From the acceleration continuity:

$a_f = a_0 + j_0 (2 t_1 - T)$

When $2 t_1 \ne T$, the jerk magnitude is:

$j_0 = \frac{a_f - a_0}{2 t_1 - T}$

---

### Velocity Constraint and Switching Time

Substituting the kinematic expressions into $v(T)=v_f$ and eliminating $j_0$
yields a **quadratic equation** in $t_1$:

$A t_1^2 + B t_1 + C = 0$

with:

$A = a_0 - a_f$

$B = 2 T a_f + 2 (v_0 - v_f)$

$C = -\frac{1}{2} T^2 (a_0 + a_f) - T (v_0 - v_f)$

A valid solution must satisfy $0 \le t_1 \le T$.

---

### Degenerate and Symmetric Cases

If $2 t_1 \approx T$, the acceleration constraint becomes numerically
ill-conditioned (typically when $a_0 \approx a_f$).

In this case, the jerk magnitude is computed from the velocity constraint:

$j_0 = \frac{v_f - v_0 - a_0 T}
           {-\frac{1}{2} T^2 + 2 T t_1 - t_1^2}$

For symmetric profiles ($t_1 = T/2$), this reduces to the familiar
minimum-jerk solution for fixed-time motion.

---

### Notes

- All equations are solved analytically.
- No numerical optimization or integration is used.
- Zero-motion and near-degenerate cases reduce naturally to $j_0 = 0$.

## Acceleration-Limited Extension (`BangBangLimitedAccel`)

The acceleration-limited profile is derived directly from the unconstrained
bang-bang jerk solution.

The unconstrained solution is always computed first, and only modified
if it violates the acceleration limit.

---

### Step 1: Solve the Unconstrained Problem

An unconstrained bang-bang jerk trajectory is computed using the same model
described above, producing:

- jerk magnitude $j_0$
- switching time $t_1$
- total duration $T$

From this solution, the acceleration profile $a(t)$ is evaluated analytically.

---

### Step 2: Detect Acceleration Clipping

The peak acceleration of the unconstrained trajectory is evaluated:

- at $t = 0$
- at the switching time $t = t_1$
- at $t = T$

If:

$|a(t)| \le a_{\max}$ for all $t \in [0, T]$

then the solution is already valid, and no modification is required.

In this case, the acceleration-limited solver reduces exactly to the
unconstrained bang-bang solution.

---

### Step 3: Construct the Clipped S-Curve

If the unconstrained solution exceeds the acceleration limit, the trajectory
is modified by inserting a **flat acceleration segment**.

The resulting profile consists of three segments:

1. **Acceleration ramp-up**  
   Constant jerk $j = \pm j_{\max}$  
   Acceleration increases from $a_0$ to $\pm a_{\max}$

2. **Constant acceleration segment**  
   Jerk $j = 0$  
   Acceleration is held at $\pm a_{\max}$

3. **Acceleration ramp-down**  
   Constant jerk $j = \mp j_{\max}$  
   Acceleration decreases to the final value $a_f$

This produces the classic **S-curve** acceleration profile.

---

### Segment Durations

The segment durations are computed analytically:

- Ramp-up time:  
  $t_1 = \frac{a_{\max} - a_0}{j_{\max}}$

- Ramp-down time:  
  $t_3 = \frac{a_f - a_{\max}}{-j_{\max}}$

The required velocity change $\Delta v = v_f - v_0$ is decomposed as:

$\Delta v = \Delta v_1 + \Delta v_2 + \Delta v_3$

where:
- $\Delta v_1$ is contributed by the ramp-up segment
- $\Delta v_3$ is contributed by the ramp-down segment
- $\Delta v_2 = a_{\max} \cdot t_{\text{flat}}$

Solving for the flat segment duration:

$t_{\text{flat}} = \frac{\Delta v - (\Delta v_1 + \Delta v_3)}{a_{\max}}$

If $t_{\text{flat}} < 0$, it is clamped to zero.

The total duration becomes:

$T = t_1 + t_{\text{flat}} + t_3$

---

### Switching Times

The resulting trajectory has either:

- **one switching time** (no acceleration clipping), or
- **two switching times** (acceleration limit reached)

These switching times are returned by the solver and can be used directly
as spline knot boundaries in downstream motion controllers.

---

### Notes

- The acceleration-limited profile is a minimal modification of the
  unconstrained solution.
- Jerk limits are always respected.
- The solution remains fully analytic.
- No re-optimization or iterative search is required.

