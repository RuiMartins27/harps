from pathlib import Path
import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

N_POINTS = 30

Y_REFLECTOR_MIN = 0.08
Y_REFLECTOR_MAX = 0.16

MPI_PROCESSES = 8

INPUT_WRITER = Path("analysis/write_input_files.py")
HARPS_EXECUTABLE = Path("harps.exe")
HARPS_INPUT = Path("input/example17.in")
RZ_SCRIPT = Path("analysis/rz_graph.py")

OUTPUT_FILE = Path("analysis/analysis_output/y_reflector_power.txt")
PLOT_FILE = Path("analysis/analysis_output/y_reflector_power.png")


# ----------------------------------------------------------------------
# Run command
# ----------------------------------------------------------------------

def run_command(command, description):
    print()
    print("=" * 70)
    print(description)
    print("Command:", " ".join(map(str, command)))
    print("=" * 70)

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if result.stdout:
        print(result.stdout)

    if result.stderr:
        print(result.stderr)

    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with return code {result.returncode}: "
            f"{' '.join(map(str, command))}"
        )

    return result.stdout


# ----------------------------------------------------------------------
# Extract integrated power
# ----------------------------------------------------------------------

def extract_total_power(output):
    """
    Extract:

        Total volume-integrated absorbed power: XXX W

    from rz_graph.py output.
    """

    pattern = (
        r"Total volume-integrated absorbed power:\s*"
        r"([-+0-9.eE]+)\s*W"
    )

    match = re.search(pattern, output)

    if match is None:
        raise RuntimeError(
            "Could not find total integrated power in rz_graph.py output."
        )

    return float(match.group(1))


# ----------------------------------------------------------------------
# Main sweep
# ----------------------------------------------------------------------

def main():

    y_values = np.linspace(
        Y_REFLECTOR_MIN,
        Y_REFLECTOR_MAX,
        N_POINTS,
    )

    results = []

    for i, y_reflector in enumerate(y_values, start=1):

        print()
        print()
        print("#" * 70)
        print(
            f"POINT {i}/{N_POINTS}: "
            f"y_reflector = {y_reflector:.8f} m"
        )
        print("#" * 70)

        # --------------------------------------------------------------
        # 1. Generate HARPS input
        # --------------------------------------------------------------

        run_command(
            [
                "python3",
                str(INPUT_WRITER),
                "-17",
                f"{y_reflector:.8f}",
            ],
            f"Generating input file for y_reflector = {y_reflector:.8f} m",
        )

        # --------------------------------------------------------------
        # 2. Run HARPS
        # --------------------------------------------------------------

        run_command(
            [
                "mpirun",
                "--bind-to",
                "socket",
                "-np",
                str(MPI_PROCESSES),
                str(HARPS_EXECUTABLE),
                str(HARPS_INPUT),
            ],
            f"Running HARPS for y_reflector = {y_reflector:.8f} m",
        )

        # --------------------------------------------------------------
        # 3. Run RZ analysis
        # --------------------------------------------------------------

        rz_output = run_command(
            [
                "python3",
                str(RZ_SCRIPT),
            ],
            "Calculating integrated absorbed power",
        )

        # --------------------------------------------------------------
        # 4. Extract power
        # --------------------------------------------------------------

        total_power = extract_total_power(rz_output)

        results.append((y_reflector, total_power))

        print()
        print(
            f">>> y_reflector = {y_reflector:.8f} m"
        )
        print(
            f">>> P_abs = {total_power:.6f} W"
        )

    # ------------------------------------------------------------------
    # Save numerical results
    # ------------------------------------------------------------------

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    np.savetxt(
        OUTPUT_FILE,
        np.array(results),
        header="y_reflector_m total_absorbed_power_W",
    )

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------

    results = np.array(results)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=300)

    ax.plot(
        results[:, 0] * 100,
        results[:, 1],
        marker="o",
        linewidth=1.5,
        markersize=4,
    )

    ax.set_xlabel(r"Reflector position $y_{\mathrm{reflector}}$ [cm]")
    ax.set_ylabel(r"Total absorbed power [W]")

    ax.grid(True, alpha=0.25)

    fig.tight_layout()

    PLOT_FILE.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        PLOT_FILE,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print()
    print("=" * 70)
    print("SWEEP COMPLETE")
    print("=" * 70)

    print(f"Results saved to: {OUTPUT_FILE}")
    print(f"Plot saved to:    {PLOT_FILE}")

    print()
    print("Results:")
    for y_reflector, power in results:
        print(
            f"  {y_reflector:.8f} m  ->  "
            f"{power:.6f} W"
        )


if __name__ == "__main__":
    main()