import subprocess
import optuna
import numpy as np
from optuna.integration import BoTorchSampler
import optuna.visualization as vis
import matplotlib.pyplot as plt
import gc
import psutil
import os
import time

count = 0

def monitor_memory():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024  # MB

def cleanup_processes():
    current_user = os.getenv('USER') or os.getenv('USERNAME')
    
    try:
        # Kill ALL MPI processes from this user (more aggressive)
        subprocess.run(["pkill", "-9", "-u", current_user, "mpirun"], check=False)
        subprocess.run(["pkill", "-9", "-u", current_user, "harps.exe"], check=False)
        
        # Clean up shared memory segments
        subprocess.run(["ipcrm", "--all=shm"], check=False)
        
        time.sleep(2)
    except:
        pass

def test_function(yr):
    # maximum at around M(0.166, 4.102) and ml(0.031, 1.161)

    global count
    count = (count + 1)

    x = yr*100 - 5
    
    return -(x + 0.1*x*x - 0.05*x*x*x + 0.0027*x*x*x*x)

def evaluate_absorbed_power(yr):
    global count
    count = (count + 1)
    #print(f"Full Evaluation with yr={yr}. Count: {count}")

    cleanup_processes()
    
    try:
        subprocess.run(["python3", "./analysis/write_input_files.py", "6", str(yr)], check=True, timeout=100 )

        subprocess.run(["mpirun", "--bind-to", "socket", "-np", "1", "harps.exe", "input/example6_symm.in"], check=True, timeout=1000)

        with open("Outputs/p_abs.txt", "r") as f:
            result = float(f.read().strip())

        cleanup_processes()
        gc.collect()

        return result

    except subprocess.TimeoutExpired:
        print("Subprocess timed out - cleaning up")
        cleanup_processes()
        raise
    except Exception as e:
        print(f"Error in evaluation: {e}")
        cleanup_processes()
        raise
    

# Objective for Optuna
def objective(trial):
    mem_start = monitor_memory()

    yr = trial.suggest_float("yr", 0.06, 0.14)

    try:
        result = evaluate_absorbed_power(yr)
        #result = test_function(yr)

        with open("Outputs/results_yr_pabs.txt", "a") as f:
            f.write(f"{yr:.6f}\t{result:.6e}\n")

        # Monitor memory after evaluation
        mem_end = monitor_memory()
        print(f"Trial {trial.number}: Memory usage: {mem_start:.1f} -> {mem_end:.1f} MB (Δ{mem_end-mem_start:.1f})")
        
        # Force garbage collection if memory usage is high
        if mem_end > 100*1024:  # more than 100GB
            gc.collect()
    except subprocess.CalledProcessError as e:
        print("External program failed:", e)
        return 1e-25 # Large penalty when the simulation fails

    return result


if __name__ == "__main__":
    print(f"Initial memory usage: {monitor_memory():.1f} MB")

    #sampler = BoTorchSampler(n_startup_trials=10)
    sampler=optuna.samplers.TPESampler( n_startup_trials=10, n_ei_candidates=24)

    study = optuna.create_study(direction="maximize", sampler=sampler)
    study.optimize(objective, n_trials=35)

    print("Best params:", study.best_params)
    print("Best value:", study.best_value)
    print(f"Final memory usage: {monitor_memory():.1f} MB")

    #cleanup_processes()

    #fig1 = vis.plot_optimization_history(study)
    #fig1.show()

    fig2 = vis.plot_slice(study, params=['yr'])
    fig2.show()

