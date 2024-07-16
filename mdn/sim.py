#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 19:44:21 2024

@author: yidong
"""

import os
import numpy as np
import torch
import torch.optim as optim
from mdn_models import MDNModel
from trainers import MDNTrainer
from my_data import get_mimic_dataloader  # Assuming the provided data loader code is saved in 'my_data.py'
import sys
import time

device = 'cuda' if torch.cuda.is_available() else 'cpu'

def generate_and_save_data(n_samples, file_path):
    # Generate binary covariates x from a Bernoulli distribution with p = 0.5
    x = (torch.rand(n_samples, 1) > 0.5).long().numpy()
    
    # Generate event times t1 from an exponential distribution with rate parameter 2 (scale = 0.5)
    t1 = np.random.exponential(scale=0.5, size=n_samples).reshape(-1, 1)
    
    # Generate event times t2 from a Rayleigh distribution with scale parameter 0.5
    t2 = np.random.rayleigh(scale=0.5, size=n_samples).reshape(-1, 1)
    
    # Assign t1 or t2 based on the value of x
    t = np.where(x == 0, t1, t2).reshape(-1)
    
    # Generate censoring times uniformly distributed over [0, 2]
    censoring_time = np.random.uniform(0, 2, size=n_samples)
    
    # Determine the observed time as the minimum of event time and censoring time
    observed_time = np.minimum(t, censoring_time)
    
    # Determine the event indicator: 1 if event is observed, 0 if censored
    delta = (t <= censoring_time).astype(float)
    
    # Save the data to a .npz file with the covariates (x) and a stack of observed times and event indicators
    np.savez(file_path, arr_0=x, arr_1=np.stack([observed_time, delta], axis=1))

def ground_truth_quantile_0(p):
    return -np.log(1 - p) / 2

def ground_truth_quantile_1(p):
    return np.sqrt(-np.log(1 - p) / 2)

def predicted_quantile(trainer, x_val, quantiles):
    t_vals = np.linspace(0, 5, 100)
    inputs = {
        "features": torch.tensor([[x_val]], dtype=torch.float32).to(device),
        "t": torch.tensor(t_vals, dtype=torch.float32).to(device)
    }
    outputs = trainer.model(inputs)
    survival_func = outputs["survival_func"].cpu().detach().numpy()
    cdf = 1 - survival_func
    pred_quantiles = np.interp(quantiles, cdf, t_vals)
    return pred_quantiles

def compute_2_wasserstein_distance(trainer):
    quantiles = np.linspace(1/5000, 4999/5000, 4999)
    
    # Predicted quantiles
    pred_quantiles_0 = predicted_quantile(trainer, 0, quantiles)
    pred_quantiles_1 = predicted_quantile(trainer, 1, quantiles)
    
    # True quantiles
    true_quantiles_0 = ground_truth_quantile_0(quantiles)
    true_quantiles_1 = ground_truth_quantile_1(quantiles)
    
    # Compute 2-Wasserstein distance
    wasserstein_0 = np.sqrt(np.mean((pred_quantiles_0 - true_quantiles_0) ** 2))
    wasserstein_1 = np.sqrt(np.mean((pred_quantiles_1 - true_quantiles_1) ** 2))
    
    # Average Wasserstein distance
    avg_wasserstein = (wasserstein_0 + wasserstein_1) / 2
    return avg_wasserstein

def survival_loss(outputs, labels):
    if torch.isnan(torch.abs(outputs["lambda"]).max()):
        sys.exit()
    batch_loss = -labels * torch.log(
        outputs["lambda"].clamp(min=1e-8)) + outputs["Lambda"]
    return torch.mean(batch_loss)

class SurvivalLossMeter(object):
    def __init__(self):
        super(SurvivalLossMeter, self).__init__()
        self.reset()

    def add(self, outputs, labels):
        self.values.append(survival_loss(outputs, labels).item())

    def value(self):
        return [np.mean(self.values)]

    def reset(self):
        self.values = []

# Simulation parameters
sample_sizes = [5000, 20000, 125000]
n_runs = 100
results = {}
data_dir = "./simulation_data"
os.makedirs(data_dir, exist_ok=True)

# Start timer
start_time = time.time()

for sample_size in sample_sizes:
    wasserstein_distances = []
    
    for run in range(n_runs):
        # Determine the number of samples for training, validation, and testing
        n_train = int(sample_size * 0.7)  # 70% for training
        n_valid = int(sample_size * 0.2)  # 20% for validation
        n_test = sample_size - n_train - n_valid  # 10% for testing
        
        # Generate and save training, validation, and test data
        train_file_path = os.path.join(data_dir, f"train_{sample_size}_{run}.npz")
        valid_file_path = os.path.join(data_dir, f"valid_{sample_size}_{run}.npz")
        test_file_path = os.path.join(data_dir, f"test_{sample_size}_{run}.npz")
        generate_and_save_data(n_train, train_file_path)
        generate_and_save_data(n_valid, valid_file_path)
        generate_and_save_data(n_test, test_file_path)
        
        # Get the dataloaders
        train_dataloader, feature_size = get_mimic_dataloader(train_file_path, batch_size=n_train, random_state=np.random.RandomState(seed=0), is_eval=False)
        valid_dataloader, _ = get_mimic_dataloader(valid_file_path, batch_size=n_valid, random_state=np.random.RandomState(seed=0), is_eval=True)
        test_dataloader, _ = get_mimic_dataloader(test_file_path, batch_size=n_test, random_state=np.random.RandomState(seed=0), is_eval=True)
        
        # Define the model
        model_config = {
            "hidden_size": 32,
            "init_type": "default",
            "num_components": 10
        }
        model = MDNModel(model_config=model_config, feature_size=1).to(device)
        
        # Define the trainer
        optimizer = optim.RMSprop(
            model.parameters(),
            lr=0.001,
            weight_decay=1e-5)
        
        criterions = {"survival_loss": survival_loss}
        metrics = {"survival_loss": SurvivalLossMeter()}
        
        trainer = MDNTrainer(
            model=model,
            device=device,
            criterions=criterions,
            optimizer=optimizer,
            dataloaders={"train": train_dataloader, "valid": valid_dataloader, "test": test_dataloader},
            metrics=metrics,
            earlystop_metric_name="survival_loss",
            batch_size=n_train,
            num_epochs=1000,
            patience=10,
            grad_clip=5,
            result_path=None,
            model_path=None,
            log_path=None,
            log_step=1000,
            exp_name=f"survival_mdn_simulation_{sample_size}_{run}",
            verbose=0,
            fine_tune=False,
            debug=False)
        
        # Train the model
        trainer.train()
        
        # Compute 2-Wasserstein distance
        if hasattr(trainer.model, "set_last_eval"):
            trainer.model.set_last_eval(True)
        
        avg_wasserstein_distance = compute_2_wasserstein_distance(trainer)
        wasserstein_distances.append(avg_wasserstein_distance)
        print(f"Run {run + 1}/{n_runs} for sample size {sample_size}: 2-Wasserstein Distance = {avg_wasserstein_distance}")
    
    # Compute average and standard error
    avg_distance = np.mean(wasserstein_distances)
    std_error = np.std(wasserstein_distances) / np.sqrt(n_runs)
    
    results[sample_size] = {
        "average_distance": avg_distance,
        "standard_error": std_error
    }
    
# End timer
end_time = time.time()
elapsed_time = (end_time - start_time) / 60

# Print final results
for sample_size, result in results.items():
    print(f"Sample Size: {sample_size}")
    print(f"Average 2-Wasserstein Distance: {result['average_distance']}")
    print(f"Standard Error: {result['standard_error']}")
    print()
    
# Print elapsed time
print(f"Total elapsed time: {elapsed_time:.2f} minutes")