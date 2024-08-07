import numpy as np
import torch
import os
import scipy
import matplotlib.pyplot as plt

from uq360.algorithms.actively_learned_model import ActivelyLearnedModel
from uq360.algorithms.ensemble_heteroscedastic_regression import EnsembleHeteroscedasticRegression

# Samples new x points, array of shape (n_points, 9)
def sample_(n_points):

    print(f'sample n_points: {n_points}')

    dimension = np.linspace(0, 1, 200)

    X_data = np.random.choice(dimension, size=(n_points, 9))

    return X_data

# Samples new y points
def querry_(X_data):

    # Save designs to a .mat file
    designs = (X_data * 0.2e-6) + 0.05e-6
    scipy.io.savemat('designs.mat', {'designs': designs})
    print(f'designs: {designs}')

    # RCWA on designs, saved in z_real.mat and z_imag.mat
    #os.system(f'matlab -batch "run(\'matlab/get_TiO2_9_pillars.m\')"')
    os.system('srun matlab -nosplash -nodesktop  -r "run(\'matlab/get_TiO2_9_pillars\'); exit"')

    # Load the results
    z_real = scipy.io.loadmat('z_real.mat')['z_real']
    z_imag = scipy.io.loadmat('z_imag.mat')['z_imag']

    return z_real

# Set device
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"device: {device}")

lamb = 0.488

# define config for model
config_model = {"num_features": 9, 
                "num_hidden_1": 256, "num_hidden_2": 256, "num_hidden_3": 256,
                "num_outputs": 1, 
                "batch_size": 16, "num_epochs": 10, "lr": 0.001}

# define config for ensemble
config_ensemble = {"num_models": 5, 
                   "batch_size": 16,
                   "model_kwargs":{"model_type":'mlp',
                                   "config": config_model,
                                   "device": device}, }

# define config for active learning object
config_AL = {"num_init": 8, 
 "T": 1, 
 "K": 8, 
 "M": 4, 
 "sampling_function": sample_, 
 "querry_function" : querry_,
 "model_function": EnsembleHeteroscedasticRegression,
 "model_kwargs": {"model_type":'ensembleheteroscedasticregression', 
                  "config":config_ensemble, 
                  "device":device}, }

# Instantiate the class object and train the model
uq_model = ActivelyLearnedModel(config=config_AL, device=device, online=True)
print('ActivelyLearnedModel created!')
#uq_model = uq_model.fit()
#print('ActivelyLearnedModel fitted!')
uq_model.save(f'neural_networks/9_pillars/lamb_{lamb}um.pth')
print('ActivelyLearnedModel saved!')

# Create a test dataset
X_test = sample_(4)
y_test = querry_(X_test)

# Predict on the test dataset
res = uq_model.predict(X_test)
y_test_pred = np.squeeze(res.y_mean, axis=1)

# Convert to phase and transmission
#phase_test = np.mod(np.angle(y_test), 2*np.pi)
#phase_test_pred = np.mod(np.angle(y_test_pred), 2*np.pi)
#T_test = np.abs(y_test)**2
#T_test_pred = np.abs(y_test_pred)**2

# Calculate the error
mae = np.mean(np.abs(y_test - y_test_pred))
#mae_phase = np.mean(np.abs(phase_test - phase_test_pred))
#mae_T = np.mean(np.abs(T_test - T_test_pred))
#frac_err_AL_ens = np.sqrt(np.sum(np.square(y_test - y_test_pred)))/np.sqrt(np.sum(np.square(y_test)))

# Print the results
print(f'y_test: {y_test}')
print(f'y_test_pred: {y_test_pred}')
print(f'mae z_real: {mae:.03f}') 
#print(f'mae phase: {mae_phase:.03f}')
#print(f'mae T: {mae_T:.03f}')
#print(f'frac_err_AL_ens: {frac_err_AL_ens}')"""

# Create a test in 1 dimension
X_test_1 = sample_(1)
X_test = np.tile(X_test_1, (200, 1))
X_test[:, 0] = np.linspace(0, 1, 200)
y_test = querry_(X_test)

# Predict on the test dataset
res = uq_model.predict(X_test)
y_test_pred = np.squeeze(res.y_mean, axis=1)

print('making figure')
plt.plot((X_test[:, 0] * 200) + 50 , y_test, 'RCWA')
plt.plot((X_test[:, 0] * 200) + 50, y_test_pred, label='pred')
plt.xlabel('Diameter 1 (nm)')
plt.ylabel('z_real')
plt.savefig('figs/active_learning_z_real.png')
print('figure saved')
