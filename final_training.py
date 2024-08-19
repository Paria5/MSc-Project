import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt

# Load the dataset
df = pd.read_csv('SNR_results.csv')

# Reshape the data
snr_columns = df.columns[:15]  # Assuming first 15 columns are SNR values for different MCS indexes
data = []

for _, row in df.iterrows():
    for i, snr in enumerate(row[:15]):
        data.append([i + 1, snr])  # i + 1 represents the MCS index (assuming 1-based index)

# Convert to DataFrame
df_long = pd.DataFrame(data, columns=['MCS_index', 'SNR'])

X = df_long[['MCS_index']]  # MCS index as feature
y = df_long['SNR']  # SNR as target

# Split the data using the best configuration
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train the Gradient Boosting Regressor with the best configuration
model = GradientBoostingRegressor(n_estimators=200, learning_rate=0.2, random_state=42)
model.fit(X_train, y_train)

# Predict on the test set
y_pred = model.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
r2 = r2_score(y_test, y_pred)

# Print evaluation metrics
print(f'MSE: {mse}')
print(f'RMSE: {rmse}')
print(f'R^2 Score: {r2}')

# Plotting actual vs predicted SNR values
plt.figure(figsize=(10, 6))
plt.scatter(y_test, y_pred, color='blue', alpha=0.5)
plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], color='red', linewidth=2)
plt.xlabel('Actual SNR')
plt.ylabel('Predicted SNR')
plt.title('Actual vs Predicted SNR')
plt.grid(True)
plt.show()

# Define a function to calculate throughput based on SNR and bandwidth
def calculate_throughput(snr, bandwidth_mhz):
    return bandwidth_mhz * np.log2(1 + snr)

# Generate SNR values for plotting throughput
snr_values = np.linspace(y.min(), y.max(), 100)

# Calculate throughput for different bandwidths
bandwidths = [10, 15, 20]
throughput_dict = {bw: calculate_throughput(snr_values, bw) for bw in bandwidths}

# Plot throughput vs SNR for different bandwidths
plt.figure(figsize=(10, 6))
for bw in bandwidths:
    plt.plot(snr_values, throughput_dict[bw], label=f'{bw} MHz Bandwidth')

plt.xlabel('SNR')
plt.ylabel('Throughput (bps)')
plt.title('Throughput vs SNR for Different Bandwidths')
plt.legend()
plt.grid(True)
plt.show()
