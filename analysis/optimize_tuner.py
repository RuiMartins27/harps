import subprocess
import optuna
import numpy as np
from optuna.integration import BoTorchSampler
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

def test_function(z1, z2, z3):
    # maximum at around M(0.2, 0.3, 0.15) and ml(0.4, 0.1, 0.35)

    global count
    count = (count + 1)
    #print(f"Full Evaluation with z1={z1}, z2={z2}, z3={z3}. Count: {count}")

    G1 = np.exp(-200*((z1-0.2)**2 + (z2-0.3)**2 + (z3-0.15)**2))
    G2 = 0.8 * np.exp(-100*((z1-0.4)**2 + (z2-0.1)**2 + (z3-0.35)**2))
    osc = 0.5 * np.sin(40*(z1 + z2 + z3)) * np.exp(-50*((z1-0.25)**2 + (z2-0.25)**2 + (z3-0.25)**2))
    quad = -0.5 * ((z1-0.25)**2 + (z2-0.25)**2 + (z3-0.25)**2)
    
    return G1 + G2 + osc + quad

def evaluate_absorbed_power(z1, z2, z3):
    global count
    count = (count + 1)
    #print(f"Full Evaluation with z1={z1}, z2={z2}, z3={z3}. Count: {count}")

    cleanup_processes()
    
    try:
        subprocess.run(["python3", "./analysis/write_input_files.py", "8", str(z1), str(z2), str(z3)], check=True, timeout=100 )

        subprocess.run(["mpirun", "--bind-to", "socket", "-np", "64", "harps.exe", "input/example8.in"], check=True, timeout=1000)

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

    if(trial.number == 0):
        z1 = trial.suggest_float("z1", 0.1, 0.12); z2 = trial.suggest_float("z2", 0.2, 0.4); z3 = trial.suggest_float("z3", 0.35, 0.355)
    else:
        z1 = trial.suggest_float("z1", 0.0, 0.5)
        z2 = trial.suggest_float("z2", 0.0, 0.5)
        z3 = trial.suggest_float("z3", 0.0, 0.5)

    try:
        result = evaluate_absorbed_power(z1, z2, z3)
        #result = test_function(z1, z2, z3)

        # Monitor memory after evaluation
        mem_end = monitor_memory()
        print(f"Trial {trial.number}: Memory usage: {mem_start:.1f} -> {mem_end:.1f} MB (Δ{mem_end-mem_start:.1f})")
        
        # Force garbage collection if memory usage is high
        if mem_end > 250*1000:  # If using more than 500GB
            gc.collect()
    except subprocess.CalledProcessError as e:
        print("External program failed:", e)
        return 1e-25 # Large penalty when the simulation fails

    return result


if __name__ == "__main__":
    print(f"Initial memory usage: {monitor_memory():.1f} MB")

    #sampler = BoTorchSampler(n_startup_trials=10)
    sampler=optuna.samplers.TPESampler( n_startup_trials=10, n_ei_candidates=24, multivariate=True)

    study = optuna.create_study(direction="maximize", sampler=sampler)
    study.optimize(objective, n_trials=100)

    print("Best params:", study.best_params)
    print("Best value:", study.best_value)
    print(f"Final memory usage: {monitor_memory():.1f} MB")
    cleanup_processes()
