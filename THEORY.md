# Theory

The mathematics behind the simulator: the model, how it is fitted, how it is simulated, how the bond is priced and its risk measured, and how the code is verified. Numbers in square brackets refer to [REFERENCES.md](REFERENCES.md).

Results quoted here come from the notebook run with seed 42 on data to 9 July 2026.

## Contents

1. [Notation](#1-notation)
2. [The Vasicek model](#2-the-vasicek-model)
3. [Calibration](#3-calibration)
4. [Simulation: the Euler–Maruyama scheme](#4-simulation-the-eulermaruyama-scheme)
5. [Monte Carlo bond pricing](#5-monte-carlo-bond-pricing)
6. [One-year Value-at-Risk and Expected Shortfall](#6-one-year-value-at-risk-and-expected-shortfall)
7. [Exact results used for checking](#7-exact-results-used-for-checking)
8. [Verification](#8-verification)
9. [Extension: the Cox–Ingersoll–Ross model](#9-extension-the-coxingersollross-model)
10. [Assumptions](#10-assumptions)

## 1. Notation

| Symbol | Meaning | Unit |
|---|---|---|
| $`r_t`$ | short-term interest rate at time $`t`$, continuously compounded | 1/year |
| $`\kappa`$ | speed of mean reversion | 1/year |
| $`\theta`$ | long-run mean rate | 1/year |
| $`\sigma`$ | volatility of the rate | 1/year$`^{3/2}`$ |
| $`W_t`$ | standard Brownian motion | $`\sqrt{\text{year}}`$ |
| $`\Delta`$ | spacing of the observed data, 1/52 | year |
| $`\Delta t`$ | simulation time step, 1/52 by default | year |
| $`N`$ | number of simulated paths | count |
| $`T`$ | bond maturity, 5 | year |
| $`F`$ | face value, 100 | LKR |
| $`c`$ | annual coupon rate, 0.10 | 1/year |
| $`P(r,\tau)`$ | price of a zero-coupon bond paying 1 after $`\tau`$ years, given rate $`r`$ | dimensionless |

## 2. The Vasicek model

The short rate follows the stochastic differential equation [1]

$$dr_t = \kappa(\theta - r_t)\,dt + \sigma\,dW_t .$$

The first term pulls the rate towards $`\theta`$: when $`r_t > \theta`$ the drift is negative, and when $`r_t < \theta`$ it is positive. The second term adds random shocks whose size does not depend on the level of the rate.

Without the noise the equation is the linear ordinary differential equation $`dr/dt = \kappa(\theta - r)`$, with solution

$$r(t) = \theta + (r_0 - \theta)e^{-\kappa t},$$

so the gap to $`\theta`$ halves every $`\ln 2/\kappa`$ years.

With the noise, the solution is an Ornstein–Uhlenbeck process [3]. Given $`r_0`$, the rate at time $`t`$ is normally distributed with

$$\mathbb{E}[r_t] = \theta + (r_0-\theta)e^{-\kappa t}, \qquad \mathrm{Var}[r_t] = \frac{\sigma^2}{2\kappa}\left(1 - e^{-2\kappa t}\right).$$

## 3. Calibration

### 3.1 From published yields to rates

The published 91-day Treasury bill yield $`y`$ is a simple annualised yield. It is converted to a continuously compounded rate over the bill's life $`\tau_b = 91/365`$:

$$r = \frac{\ln(1 + y\,\tau_b)}{\tau_b}.$$

### 3.2 Exact discretisation

Because $`r_{t+\Delta}`$ given $`r_t`$ is normal with the mean and variance above, the model sampled every $`\Delta`$ years is exactly a first-order autoregression [4]:

$$r_{t+\Delta} = a + b\,r_t + \varepsilon_t, \qquad a = \theta(1-e^{-\kappa\Delta}), \quad b = e^{-\kappa\Delta}, \quad \varepsilon_t \sim \mathcal{N}\left(0,\ \sigma^2\frac{1-e^{-2\kappa\Delta}}{2\kappa}\right).$$

### 3.3 Estimators

Fit the line $`r_{t+\Delta} = a + b\,r_t`$ by ordinary least squares, let $`\hat{s}^2`$ be the residual variance, and invert:

$$\hat\kappa = -\frac{\ln\hat b}{\Delta}, \qquad \hat\theta = \frac{\hat a}{1-\hat b}, \qquad \hat\sigma = \sqrt{\frac{2\hat\kappa\,\hat{s}^2}{1-\hat{b}^2}} .$$

This requires $`0 < \hat b < 1`$. If $`\hat b \ge 1`$ the data show no mean reversion, the logarithm gives $`\hat\kappa \le 0`$, and the model is rejected for that window. The code raises an error in that case.

### 3.4 Result and its reliability

| Window | Weeks | $`\hat\kappa`$ | $`\hat\theta`$ | $`\hat\sigma`$ |
|---|---|---|---|---|
| 15 Feb 2024 to 9 Jul 2026 (baseline) | 126 | 0.766 | 8.62% | 0.0141 |
| 1 Jan 2024 to 9 Jul 2026 | 132 | 4.224 | 8.20% | 0.0181 |
| 1 Jan 2020 to 9 Jul 2026 | 332 | 0.214 | 13.94% | 0.0509 |
| 1 Jan 2025 to 9 Jul 2026 | 80 | rejected ($`\hat b > 1`$) | | |

The baseline starts in the first week in which the 91-day yield was back in single digits after the 2022 to 2023 crisis.

The estimate of $`\kappa`$ is the weakest of the three. Its accuracy depends on the length of the sample in years, not on the number of observations, and in short samples it is biased upwards [11]. A parameter-recovery experiment with 1,000 synthetic samples of 126 weeks, generated with $`\kappa = 0.766`$, gave a median estimate of 2.30 and a 95% range of 0.28 to 7.62. The same experiment recovered $`\theta`$ and $`\sigma`$ closely.

## 4. Simulation: the Euler–Maruyama scheme

Time is divided into steps of length $`\Delta t`$. With independent $`Z_n \sim \mathcal{N}(0,1)`$ [7]:

$$r_{n+1} = r_n + \kappa(\theta - r_n)\,\Delta t + \sigma\sqrt{\Delta t}\,Z_n .$$

This is Euler's method applied to the drift, plus a random increment with the variance of a Brownian step, $`\Delta t`$.

**Stability.** The deterministic part multiplies the gap $`r_n - \theta`$ by $`(1 - \kappa\Delta t)`$ each step. It decays without oscillating when $`\kappa\Delta t < 1`$ and is stable when $`\kappa\Delta t < 2`$. With $`\kappa = 0.766`$ and $`\Delta t = 1/52`$, $`\kappa\Delta t = 0.015`$.

**Accuracy.** Euler–Maruyama has strong order 1/2 in general and weak order 1 [7]. Prices are expectations, so the weak order is the relevant one: the bias of a price falls in proportion to $`\Delta t`$. For the Vasicek model the noise does not depend on the rate, and the scheme's only error is in the drift.

**Random numbers.** The notebook draws $`Z_n`$ with NumPy's `Generator.standard_normal`, seeded through `SeedSequence` so that each replication has an independent stream. The web demo uses a seeded uniform generator and the Box–Muller transform [9].

## 5. Monte Carlo bond pricing

### 5.1 Zero-coupon bond

A payment of 1 at time $`T`$ is worth the expected discount factor along the rate path:

$$P(0,T) = \mathbb{E}\left[\exp\left(-\int_0^T r_s\,ds\right)\right].$$

On each simulated path the integral is computed with the trapezoidal rule, and the expectation is replaced by an average over $`N`$ paths [5], [6]:

$$\hat{P} = \frac{1}{N}\sum_{i=1}^{N} D^{(i)}, \qquad D^{(i)} = \exp\left(-\sum_{n=0}^{M-1}\tfrac{1}{2}\left(r^{(i)}_n + r^{(i)}_{n+1}\right)\Delta t\right).$$

### 5.2 Coupon bond

A bond with face value $`F`$, annual coupon rate $`c`$ and half-yearly payments at $`t_1, \dots, t_m = T`$ is a set of zero-coupon payments:

$$V_0 = \sum_{j=1}^{m}\frac{cF}{2}\,P(0,t_j) + F\,P(0,T).$$

In the code the discount factors at all payment dates are recorded on the same paths, so one simulation prices the whole bond.

### 5.3 Error of the estimate

By the law of large numbers $`\hat P \to P`$ as $`N \to \infty`$. By the central limit theorem the error is approximately normal with standard error $`s/\sqrt{N}`$, where $`s`$ is the sample standard deviation of the per-path values, giving the 95% confidence interval

$$\hat{P} \pm 1.96\,\frac{s}{\sqrt{N}} .$$

The total error has two parts: this sampling error, of order $`N^{-1/2}`$, and the discretisation bias, of order $`\Delta t`$. Results are also reported across 10 independent replications, with a Student-t interval.

## 6. One-year Value-at-Risk and Expected Shortfall

Let $`V_1`$ be the value of the position after one year: the value of the payments still to come, given the rate $`r_1`$ at that time, plus the coupons received during the year (held as cash). The loss is

$$L = V_0 - V_1 ,$$

so a negative $`L`$ is a gain. At confidence level $`\alpha = 0.95`$:

$$\mathrm{VaR}_\alpha = \text{the } \alpha\text{-quantile of } L, \qquad \mathrm{ES}_\alpha = \mathbb{E}\left[L \mid L \ge \mathrm{VaR}_\alpha\right].$$

Value-at-Risk is the loss exceeded in only 5% of scenarios. Expected Shortfall is the average loss within that worst 5%, and unlike Value-at-Risk it is a coherent risk measure [12].

**Revaluation grid.** Finding $`V_1`$ for every scenario would need a second simulation inside each path. Instead:

1. Simulate $`r_1`$ on 50,000 paths.
2. Take 41 equally spaced rates covering the range of simulated $`r_1`$, and price the remaining payments by Monte Carlo from each of them, using the same random numbers at every grid point so the resulting curve is smooth.
3. Read each scenario's value from that curve by linear interpolation.

**Measure.** Pricing uses risk-neutral dynamics and the one-year scenarios use real-world dynamics. With the market price of interest-rate risk set to zero (assumption 4), the two coincide, and the same parameters are used for both.

## 7. Exact results used for checking

The Vasicek zero-coupon price has a closed form [1], [4]:

$$P(r,\tau) = A(\tau)\,e^{-B(\tau)\,r}, \qquad B(\tau) = \frac{1-e^{-\kappa\tau}}{\kappa},$$

$$\ln A(\tau) = \left(\theta - \frac{\sigma^2}{2\kappa^2}\right)\left(B(\tau)-\tau\right) - \frac{\sigma^2 B(\tau)^2}{4\kappa}.$$

Because $`r_1 \sim \mathcal{N}(m, s^2)`$, with $`m`$ and $`s^2`$ from Section 2 at $`t = 1`$, and the bond value falls as the rate rises, the risk measures are also exact. With $`z = \Phi^{-1}(0.95)`$, remaining payments $`c_j`$ at remaining times $`\tau_j`$, $`A_j = A(\tau_j)`$, $`B_j = B(\tau_j)`$, and coupons received $`C`$:

$$\mathrm{VaR}_{0.95} = V_0 - C - \sum_j c_j\,P(m + z s,\ \tau_j),$$

$$\mathrm{ES}_{0.95} = V_0 - C - \sum_j c_j\,A_j\,e^{-B_j m + B_j^2 s^2/2}\,\frac{\Phi(-z - B_j s)}{1-0.95}.$$

The second formula follows from $`\mathbb{E}\left[e^{-bZ}\,;\,Z \ge z\right] = e^{b^2/2}\,\Phi(-z-b)`$ for a standard normal $`Z`$.

These formulas are used only to check the simulation. Every reported result comes from the simulation itself.

## 8. Verification

Each check is an `assert` in the notebook, so the notebook stops if one fails.

| Check | What is compared | Pass rule | Result |
|---|---|---|---|
| 1a | Monte Carlo zero-coupon price, $`\Delta t = 1/364`$, against the exact price | within 3 standard errors | 63.8215 against 63.8277; 3 SE = 0.0196 |
| 1b | The same with the weekly step used for all results | within 0.05 LKR | difference +0.0151 |
| 2 | Mean and variance of the simulated rate at 1 and 5 years against Section 2 | mean within 3 SE; variance within 2% | variance errors +0.46% and +0.52% |
| 3 | $`\sigma = 0`$ against the exact ODE solution | identical paths; error below 1 basis point; error halves when the step halves | error $`4.0\times10^{-5}`$; ratio 2.006 |
| 4 | Parameters recovered from synthetic data (200 paths of 1,000 weeks) | $`\kappa`$ within 10%; $`\theta`$, $`\sigma`$ within 2% | +4.6%, −0.4%, +0.2% |
| 5 | Monte Carlo VaR and ES against Section 7 | within 0.05 LKR | differences −0.025 and −0.029 |
| 6a | Sampling error against $`N`$, log–log slope | between −0.6 and −0.4 | −0.503 (theory −0.5) |
| 6b | Time-step bias against $`\Delta t`$, log–log slope | between 0.8 and 1.2 | 0.981 (theory 1) |
| 7 | CIR Monte Carlo price against the exact CIR price | within 3 SE (fine step); within 0.05 LKR (weekly) | 63.6436 against 63.6359 |

**Why check 1 has two parts.** With a fine step the discretisation bias is negligible, so the estimate should agree with the exact price to within sampling error. With the weekly step the estimate carries a small, systematic bias of about 0.015 LKR, which a three-standard-error test would eventually detect as the number of paths grows. It is therefore tested against a stated tolerance.

**How the scaling is measured.** For 6a, one run of 400,000 paths is cut into blocks of $`N`$ paths, and the root-mean-square deviation of the block averages is measured for $`N`$ from 100 to 40,000. For 6b, the price is computed with 2, 4, 12, 26 and 52 steps per year, using 400,000 paths each so that the bias stands clear of the sampling error.

## 9. Extension: the Cox–Ingersoll–Ross model

The CIR model [2] scales the noise by the square root of the rate:

$$dr_t = \kappa(\theta - r_t)\,dt + \sigma\sqrt{r_t}\,dW_t .$$

Shocks shrink as the rate approaches zero. If the Feller condition $`2\kappa\theta \ge \sigma^2`$ holds, the rate stays strictly positive.

**Simulation.** A plain Euler step can produce a negative rate, under which the square root is undefined. The full-truncation scheme [8] uses $`r^+ = \max(r, 0)`$ inside the drift and the square root:

$$r_{n+1} = r_n + \kappa(\theta - r_n^+)\,\Delta t + \sigma\sqrt{r_n^+\,\Delta t}\;Z_n .$$

**Calibration.** Dividing the Euler form by $`\sqrt{r_t}`$ gives a regression with constant error variance:

$$\frac{r_{t+\Delta}-r_t}{\sqrt{r_t}} = \kappa\theta\,\frac{\Delta}{\sqrt{r_t}} - \kappa\,\Delta\sqrt{r_t} + \sigma\sqrt{\Delta}\;\varepsilon_t .$$

Least squares on the two regressors gives $`\kappa\theta`$ and $`\kappa`$, and the residual standard deviation divided by $`\sqrt{\Delta}`$ gives $`\sigma`$.

**Exact price.** With $`\gamma = \sqrt{\kappa^2 + 2\sigma^2}`$ [2]:

$$P(r,\tau) = A(\tau)e^{-B(\tau)r}, \qquad B(\tau) = \frac{2(e^{\gamma\tau}-1)}{(\gamma+\kappa)(e^{\gamma\tau}-1)+2\gamma},$$

$$A(\tau) = \left[\frac{2\gamma\,e^{(\kappa+\gamma)\tau/2}}{(\gamma+\kappa)(e^{\gamma\tau}-1)+2\gamma}\right]^{2\kappa\theta/\sigma^2}.$$

**Result.** Fitted to the baseline window: $`\kappa = 0.678`$, $`\theta = 8.64\%`$, $`\sigma = 0.0470`$, with the Feller condition satisfied (0.117 against 0.002).

| Model | Coupon bond price | 95% VaR | 95% ES |
|---|---|---|---|
| Vasicek | 103.03 ± 0.02 | −8.22 ± 0.02 | −7.73 ± 0.03 |
| CIR | 102.77 ± 0.02 | −7.84 ± 0.02 | −7.20 ± 0.03 |

At current rate levels the two models give similar answers. The choice between them matters far less than the choice of calibration window in Section 3.4.

## 10. Assumptions

1. The 91-day Treasury bill yield is a proxy for the instantaneous short rate.
2. Observations are treated as equally spaced, $`\Delta = 1/52`$ year, although the published week-end dates are 6 to 10 days apart.
3. $`\kappa`$, $`\theta`$ and $`\sigma`$ are constant over the calibration window and over the life of the bond.
4. The market price of interest-rate risk is zero, so the same parameters serve for pricing and for the one-year horizon. Observed 5-year auction yields exceed the model's 5-year par yield by 2.15 ± 0.56 percentage points on average, which measures what this assumption leaves out.
5. Every cash flow is paid in full and on time.
6. Coupons received during the year are held as cash and not reinvested.