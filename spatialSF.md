## Intercity Travel Resilience Evaluation Indicators

The evaluation indicators for intercity travel resilience primarily encompass fluctuation ratio, recovery ratio, recovery elasticity, and recovery index. These concepts and their corresponding formulas are standardized and organized below.

### 1. Fluctuation Ratio (FR) and Recovery Ratio (RR)

#### **Fluctuation Ratio (FR)**
The Fluctuation Ratio (FR) is defined as the proportion of the difference between the baseline value of the population migration index before the pandemic outbreak and its minimum value during the pandemic (i.e., fluctuation amount) to the baseline value. This ratio reflects the city's level of impact due to the pandemic.

#### **Recovery Ratio (RR)**
The Recovery Ratio (RR) refers to the proportion of the difference between the migration index from its minimum value during the pandemic to its value at a specific post-pandemic time (i.e., recovery amount) to the baseline value. This ratio measures the extent of the city's post-pandemic recovery.

**Formulas:**

$$FR_i = \frac{Trav_{i,t_1} - Trav_{i,t_m}}{Trav_{i,t_1}} \times 100\%$$

$$RR_i = \frac{Trav_{i,t_2} - Trav_{i,t_m}}{Trav_{i,t_1}} \times 100\%$$

Where:
- \(FR_i\) represents the fluctuation ratio of the ith city.
- \(RR_i\) represents the recovery ratio of the ith city.
- \(Trav_{i,t}\) denotes the population migration index of the ith city at time \(t\).
- \(t_1\), \(t_2\), and \(t_m\) respectively signify the times corresponding to seven days before high-risk area designation, seven days after the lifting of high-risk areas, and the lowest point of the migration index during the pandemic.

### 2. City Recovery Elasticity (CRE)

#### **City Recovery Elasticity (CRE)**
CRE is used to describe the ratio of the city's recovery amount relative to the fluctuation amount, equivalent to the ratio of the recovery ratio to the fluctuation ratio. A higher CRE indicates that the city can recover more quickly and to a greater extent under the same fluctuation ratio, demonstrating stronger resilience.

**Formula:**

$$CRE_{i,t} = \frac{Trav_{i,t_2} - Trav_{i,t_m}}{Trav_{i,t_1} - Trav_{i,t_m}} \times 100\% = \frac{RR_i}{FR_i} \times 100\%$$

### 3. City Recovery Index (CRI)

#### **City Recovery Index (CRI)**
The CRI reflects the recovery level of the city at a specific post-pandemic time relative to the baseline value, characterized by the change in the population migration index. Values closer to 100% indicate a better recovery effect of the city.

**Formula:**

$$CRI_{i,t} = \frac{Trav_{i,t_2}}{Trav_{i,t_1}} \times 100\%$$

Here, \(CRI_{i,t}\) represents the recovery index of the ith city at time \(t\). Other variables are defined as before.

---
## Moran Scatter Plots Analysis

### Overview
The Moran scatter plots generated using `spatialRF::plot_training_df_moran` function display the relationship between selected variables or interactions and the response variable (CRE). These plots are instrumental in visualizing spatial autocorrelation and understanding the significance of variable interactions in explaining the variation in the response variable.

### Moran Scatter Plots
- **Subplots**: Each subplot represents a Moran scatter plot for a particular variable or PCA component, with the x-axis representing the variable or PCA component values, and the y-axis representing the spatially lagged values.
- **Positive Autocorrelation**: Clustering of points around the diagonal line indicates a positive spatial autocorrelation.
- **Negative Autocorrelation**: Dispersion of points from the line indicates negative spatial autocorrelation.

### Violin Plot
- **Model Comparison**: The violin plot compares the distribution of R-squared values for models with and without the inclusion of interactions over spatial folds, providing insight into the improvement in model performance due to interactions.

### Statistical Metrics
- **R-squared (+R2)**: Indicates the explanatory power of the model, with higher values suggesting a better model fit.
- **Importance (Imp. (%))**: Reflects the relative importance of the variable or PCA component within the model.
- **Maximum Correlation (Max cor)**: Shows the maximum Pearson correlation between the variable and the predictors, suggesting the degree of linear association.

### Interaction Identification with `the_feature_engineer()`
The `the_feature_engineer()` function identifies beneficial interactions between variables by evaluating potential interactions and recommending those that are not highly correlated with other predictors and improve the model's R-squared value on independent data.

### Results Application
Identified interactions and their importance are added to the dataset, allowing for the extended dataset to be used for further model training and analysis.

### Code and Plot Generation
The code attempts to explore and visualize the impact of these interactions. The Moran scatter plots show the relationships between predictor variables and the response variable, while the violin plot contrasts the performance of models with and without these interactions.

### Final Notes
The decision to use specific thresholds for interactions should be guided by domain knowledge and specific analysis requirements. It is crucial to interpret these interactions in the context of the study area and ensure they are sensible from a biological, ecological, or geographical standpoint.

### Image Analysis
![alt text](<figure/Comparing models with and without interactions via spatial cross-validation..png>)
*The image shows a series of Moran scatter plots indicating the relationship between different predictors and the response variable, along with a violin plot comparing models with and without interactions.*

The Moran scatter plots indicate varying degrees of spatial autocorrelation for different variables. The violin plot demonstrates that certain interactions can significantly impact the predictive power of the spatial model, as evidenced by changes in the R-squared values across spatial folds.

