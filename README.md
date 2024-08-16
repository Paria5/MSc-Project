# MSc-Project
This is the edition where we separate the snr plots from the data obtained in another function called plotResults for more modularity and a neater code.
In this version. the results will also be stored at the end of the code in a .mat file.
ber and throughput graphs need refining.

MCS indexes without RIS has been obtained
MCS indexing is corrected.

Actual vs. Predicted SNR Plot:
Purpose: This plot helps you assess how well your model's predictions match the actual SNR values.
Interpretation: The scatter plot shows the relationship between the actual SNR values (y_test) and the predicted SNR values (y_pred). The red dashed line represents the ideal scenario where predictions perfectly match actual values. If most points lie close to this line, it indicates good model performance. Significant deviations from this line suggest areas where the model's predictions are off.

Residuals Plot:
Purpose: This plot helps diagnose the model's prediction errors by examining residuals (the difference between actual and predicted values).
Interpretation: The residuals plot shows residuals on the y-axis and predicted values on the x-axis. Ideally, residuals should be randomly scattered around zero. If there's a pattern or systematic structure in the residuals, it may indicate that the model isn't capturing some underlying trend or that there's a problem with the data or model. The red horizontal line at zero helps visualize where the residuals are centered.

Distribution of Residuals Plot:
Purpose: This plot provides insight into the distribution and magnitude of the residuals.
Interpretation: The histogram shows the frequency of different residual values. A roughly normal distribution of residuals (bell-shaped curve) is often a sign of a well-behaved model. If the residuals are skewed or have heavy tails, it may suggest issues with model assumptions or the presence of outlier
