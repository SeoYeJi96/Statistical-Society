# Statistical-Society

Abstract : We propose a sequential importance sampling for estimating the multi-period market risk. 1
In the Monte Carlo method for estimating the market risk when holding a portfolio of assets over 2
multiple time periods, the commonly applied method is to fit a GARCH-type model to the observed 3
log returns over a single time-period, and simulate the process of future log returns by applying the 4
fitted model. We propose to apply the importance sampling in generation of the white noises of the 5
fitted GARCH model instead of sampling from the nominal distribution. We first approximate the 6
distribution of the white noises of the fitted GARCH model as a finite mixture of normal distributions, 7
and then adopt the exponentially twisted distribution of it as the importance sampling distribution 8
of white noises for the simulation of future log return processes. We estimate the risk measures such 9
as the value-at-risk and the expected short-fall from the simulated log returns over the multiple 10
time periods of interest. We prove that the optimal twisting parameter is unique, and show that an 11
approximate value to it can be found by applying the stochastic approximation. Numerical results are 12
given to compare the performance of proposed method to that of the crude Monte Carlo simulation 13
in terms of the accuracy of the estimated risk measures.
