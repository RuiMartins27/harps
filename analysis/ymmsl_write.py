def generate_ymmsl(n_macros=40, filename="../model_polar.ymmsl"):
    with open(filename, "w") as f:
        f.write("ymmsl_version: v0.1\n\n")
        f.write("model:\n")
        f.write(f"  name: muscle3_harps_polar{n_macros}\n\n")
        f.write("  components:\n")

        # Macros
        for i in range(n_macros):
            f.write(f"    macro_{i}:\n")
            f.write("      implementation: macro\n")
            f.write("      ports:\n")
            f.write("        o_i:\n")
            f.write("          - \"run_harps\"\n")
            f.write("        s:\n")
            f.write("          - \"read_harps\"\n\n")

        # Micro
        f.write("    micro:\n")
        f.write("      implementation: micro\n")
        f.write("      ports:\n")
        f.write("        f_init:\n")
        for i in range(n_macros):
            f.write(f"          - \"begin_harps_{i}\"\n")
        f.write("        o_f:\n")
        for i in range(n_macros):
            f.write(f"          - \"end_harps_{i}\"\n")
        f.write("\n")

        # Conduits
        f.write("  conduits:\n")
        for i in range(n_macros):
            f.write(f"    macro_{i}.run_harps: micro.begin_harps_{i}\n")
            f.write(f"    micro.end_harps_{i}: macro_{i}.read_harps\n")

        # Implementations
        f.write("\nimplementations:\n")
        f.write("  macro:\n")
        f.write("    env:\n")
        f.write("      +LD_LIBRARY_PATH: :/home/martins/muscle3/lib\n")
        f.write("    executable: ~/power_coupling/1d_codes/main_polar\n\n")
        f.write("  micro:\n")
        f.write("    env:\n")
        f.write("      +LD_LIBRARY_PATH: :/home/martins/muscle3/lib\n")
        f.write("    executable: mpirun\n")
        f.write("    args: --bind-to socket -np 1 ~/power_coupling/harps/harps_polar.exe input/coupled_with_1D_radial_muscle_symm.in\n\n")

        # Settings
        f.write("settings:\n")
        f.write("  s_1: 1e-6\n\n")

        # Resources
        f.write("resources:\n")
        for i in range(n_macros):
            f.write(f"  macro_{i}:\n")
            f.write("    threads: 1\n")
        f.write("  micro:\n")
        f.write("    threads: 1\n")


if __name__ == "__main__":
    generate_ymmsl(40)
