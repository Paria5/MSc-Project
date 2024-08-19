import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error
import numpy as np

# Load the dataset
df = pd.read_csv('SNR_results.csv')

# Reshape the data
snr_columns = df.columns[:15]  # SNR columns
data = []

for _, row in df.iterrows():
    for i, snr in enumerate(row[:15]):
        data.append([i + 1, snr])  # i + 1 represents the MCS index (assuming 1-based index)

# Convert to DataFrame
df_long = pd.DataFrame(data, columns=['MCS_index', 'SNR'])

X = df_long[['MCS_index']]  # MCS index as feature
y = df_long['SNR']  # SNR as target

# Define the range of the SNR values
snr_min = y.min()
snr_max = y.max()
snr_range = snr_max - snr_min

# Store results
results = []

# Gradient Boosting Regressor configurations
test_sizes = [0.1, 0.2, 0.3, 0.4, 0.5]
n_estimators_list = [50, 100, 200]
learning_rates = [0.01, 0.1, 0.2]

for test_size in test_sizes:
    for n_estimators in n_estimators_list:
        for learning_rate in learning_rates:
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
            model = GradientBoostingRegressor(n_estimators=n_estimators, learning_rate=learning_rate, random_state=42)
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            mse = mean_squared_error(y_test, y_pred)
            rmse = np.sqrt(mse)
            relative_error = (rmse / snr_range) * 100
            estimated_accuracy = 100 - relative_error
            results.append({
                'model': 'GBM',
                'test_size': test_size,
                'n_estimators': n_estimators,
                'learning_rate': learning_rate,
                'mse': mse,
                'rmse': rmse,
                'relative_error (%)': relative_error,
                'estimated_accuracy (%)': estimated_accuracy
            })

# MLP Regressor configurations
hidden_layer_sizes_list = [(10,), (50,), (100,), (50, 50), (100, 50)]
activation_functions = ['identity', 'logistic', 'tanh', 'relu']

for test_size in test_sizes:
    for hidden_layer_sizes in hidden_layer_sizes_list:
        for activation in activation_functions:
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
            model = MLPRegressor(hidden_layer_sizes=hidden_layer_sizes, activation=activation, max_iter=500, random_state=42)
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            mse = mean_squared_error(y_test, y_pred)
            rmse = np.sqrt(mse)
            relative_error = (rmse / snr_range) * 100
            estimated_accuracy = 100 - relative_error
            results.append({
                'model': 'MLP',
                'test_size': test_size,
                'hidden_layer_sizes': hidden_layer_sizes,
                'activation': activation,
                'mse': mse,
                'rmse': rmse,
                'relative_error (%)': relative_error,
                'estimated_accuracy (%)': estimated_accuracy
            })

# Support Vector Regressor configurations
kernels = ['linear', 'rbf', 'poly']
C_values = [0.1, 1, 10]

for test_size in test_sizes:
    for kernel in kernels:
        for C in C_values:
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
            model = SVR(kernel=kernel, C=C)
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            mse = mean_squared_error(y_test, y_pred)
            rmse = np.sqrt(mse)
            relative_error = (rmse / snr_range) * 100
            estimated_accuracy = 100 - relative_error
            results.append({
                'model': 'SVR',
                'test_size': test_size,
                'kernel': kernel,
                'C': C,
                'mse': mse,
                'rmse': rmse,
                'relative_error (%)': relative_error,
                'estimated_accuracy (%)': estimated_accuracy
            })

# Random Forest Regressor configurations
for test_size in test_sizes:
    for n_estimators in n_estimators_list:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
        model = RandomForestRegressor(n_estimators=n_estimators, random_state=42)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        mse = mean_squared_error(y_test, y_pred)
        rmse = np.sqrt(mse)
        relative_error = (rmse / snr_range) * 100
        estimated_accuracy = 100 - relative_error
        results.append({
            'model': 'RF',
            'test_size': test_size,
            'n_estimators': n_estimators,
            'mse': mse,
            'rmse': rmse,
            'relative_error (%)': relative_error,
            'estimated_accuracy (%)': estimated_accuracy
        })

# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Sort by estimated accuracy and select top 10 configurations
top_10_results = results_df.sort_values(by='estimated_accuracy (%)', ascending=False).head(10)

# Print the top 10 configurations
print("Top 10 Configurations Based on Estimated Accuracy:")
print(top_10_results)
