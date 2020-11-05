# fastlogitME
Fast but Basic Marginal Effects for Logit Models in R.

Calculates marginal effects based on logistic model objects such as "glm" or "speed.glm" at the average (default) or at given values indicated by at. It also returns confidence intervals for said marginal effects and the p-values, which can easily be inputed in stargazer. The function only returns the essentials and is therefore much faster but not as detailed as other functions available to calculate marginal effects. As a result, it is highly suitable for large datasets for which other packages may require too much time or calculating power.

Background:
I wrote this function to work with large datasets in R. Estimating logit models and the respective marginal effects can in these cases take a long time and a lot of CPU. Therefore I make use of speed.glm to estimate the logit model and I use this function to calculate the marginal effects, as well as the confidence intervals and p-values of said marginal effects. In the help file of the function I included details on how to export these results with stargazer.
