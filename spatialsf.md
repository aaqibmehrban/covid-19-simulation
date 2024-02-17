Moran Scatter Plots: Each subplot shows the Moran scatter plot for a particular variable or PCA component. The Moran scatter plot is a way of visualizing spatial autocorrelation, where the plot's x-axis represents the values of the variable or PCA component, and the y-axis represents the spatially lagged values. Points clustered around the diagonal line indicate a positive autocorrelation, while dispersion from the line indicates negative autocorrelation.
Violin Plot: The violin plot compares the distribution of R-squared values across spatial folds for models with and without the inclusion of interactions, providing a visual assessment of the improvement in model performance due to the interactions.
Statistical Metrics: Above each Moran plot, there is a summary of the variable’s contribution to the model:
R2 (R-squared): Indicates how much of the variance in the response variable is explained by the model. A higher R-squared value indicates a better fit.
Imp. (%): Reflects the relative importance of the variable or PCA component in the model.
Max cor: Represents the maximum Pearson correlation between the variable and the predictors, suggesting the degree of linear association.
