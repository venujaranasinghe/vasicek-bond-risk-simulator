# References

All sources behind the model, the numerical methods and the data. Numbers in square brackets are the ones used in [THEORY.md](THEORY.md).

## Interest-rate models

[1] O. Vasicek, "An equilibrium characterization of the term structure," *J. Financ. Econ.*, vol. 5, no. 2, pp. 177–188, 1977, doi: [10.1016/0304-405X(77)90016-2](https://doi.org/10.1016/0304-405X(77)90016-2).
The short-rate model used here and its closed-form bond price.

[2] J. C. Cox, J. E. Ingersoll and S. A. Ross, "A theory of the term structure of interest rates," *Econometrica*, vol. 53, no. 2, pp. 385–407, 1985, doi: [10.2307/1911242](https://doi.org/10.2307/1911242).
The square-root (CIR) model used in the extension, and its closed-form bond price.

[3] G. E. Uhlenbeck and L. S. Ornstein, "On the theory of the Brownian motion," *Phys. Rev.*, vol. 36, no. 5, pp. 823–841, 1930, doi: [10.1103/PhysRev.36.823](https://doi.org/10.1103/PhysRev.36.823).
The mean-reverting process underlying the Vasicek model, with its mean and variance.

[4] D. Brigo and F. Mercurio, *Interest Rate Models – Theory and Practice*, 2nd ed. Berlin, Germany: Springer, 2006, doi: [10.1007/978-3-540-34604-3](https://doi.org/10.1007/978-3-540-34604-3).
Reference text for short-rate models, the exact discretisation and bond-price formulas.

## Monte Carlo and numerical methods

[5] P. P. Boyle, "Options: A Monte Carlo approach," *J. Financ. Econ.*, vol. 4, no. 3, pp. 323–338, 1977, doi: [10.1016/0304-405X(77)90005-8](https://doi.org/10.1016/0304-405X(77)90005-8).
The original use of Monte Carlo simulation for pricing financial contracts.

[6] P. Glasserman, *Monte Carlo Methods in Financial Engineering*. New York, NY, USA: Springer, 2003, doi: [10.1007/978-0-387-21617-1](https://doi.org/10.1007/978-0-387-21617-1).
Monte Carlo error and confidence intervals, discretisation bias, simulation of the Vasicek and CIR models.

[7] P. E. Kloeden and E. Platen, *Numerical Solution of Stochastic Differential Equations*. Berlin, Germany: Springer, 1992, doi: [10.1007/978-3-662-12616-5](https://doi.org/10.1007/978-3-662-12616-5).
The Euler–Maruyama scheme and its strong and weak orders of convergence.

[8] R. Lord, R. Koekkoek and D. van Dijk, "A comparison of biased simulation schemes for stochastic volatility models," *Quant. Finance*, vol. 10, no. 2, pp. 177–194, 2010, doi: [10.1080/14697680802392496](https://doi.org/10.1080/14697680802392496).
The full-truncation Euler scheme used to simulate the CIR model.

[9] G. E. P. Box and M. E. Muller, "A note on the generation of random normal deviates," *Ann. Math. Statist.*, vol. 29, no. 2, pp. 610–611, 1958, doi: [10.1214/aoms/1177706645](https://doi.org/10.1214/aoms/1177706645).
The method used in the web demo to turn uniform random numbers into normal ones.

[10] A. M. Law, *Simulation Modeling and Analysis*, 5th ed. New York, NY, USA: McGraw-Hill, 2015.
Input modelling, random variate generation, replications and confidence intervals.

## Estimation and risk measures

[11] C. Y. Tang and S. X. Chen, "Parameter estimation and bias correction for diffusion processes," *J. Econometrics*, vol. 149, no. 1, pp. 65–81, 2009, doi: [10.1016/j.jeconom.2008.11.001](https://doi.org/10.1016/j.jeconom.2008.11.001).
Why the speed of mean reversion is hard to estimate and biased in short samples, for the Vasicek and CIR models.

[12] P. Artzner, F. Delbaen, J.-M. Eber and D. Heath, "Coherent measures of risk," *Math. Finance*, vol. 9, no. 3, pp. 203–228, 1999, doi: [10.1111/1467-9965.00068](https://doi.org/10.1111/1467-9965.00068).
Properties of risk measures; the basis for reporting Expected Shortfall alongside Value-at-Risk.

## Data

[13] Central Bank of Sri Lanka, "Rates on Government Securities (per cent per annum)," Data Library. [Online]. Available: https://www.cbsl.lk/eResearch/MoneyMarketRatesDefault.aspx?ReportId=6169
Weekly Treasury bill and Treasury bond yields, January 2020 to July 2026. See [data/README.md](data/README.md).