import os
import csv
import matplotlib.pyplot as plt

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# =======================
# USER SETTINGS
# =======================
directoryNames = ["from_Ref1_C3_0.3_SST", "from_Ref1_C3_0.5_SST", "from_Ref1_C3_0.7_SST"]       # simulation directories
PlotLabelNames = ["C3=0.3", "C3=0.5", "C3=0.7"]         # labels for legend
SORT_FIELD = "R"

# =======================
# PATHS
# =======================
try:
    base_path = os.path.dirname(os.path.abspath(__file__))
except NameError:
    base_path = os.getcwd()

csv_file = os.path.join(base_path, "experiments.csv")

# =======================
# HELPERS
# =======================
def is_time_dir(name):
    try:
        return float(name) > 0.0
    except ValueError:
        return False

# =======================
# READ EXPERIMENT CSV
# =======================
experiment_data = {}
with open(csv_file, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        experiment_data[row["Run"].strip()] = row

# =======================
# EXPERIMENTAL DATA
# =======================
# We'll assume all directories share the same experimental lengths
exp_lengths_dict = {}  # key: run_id, value: experimental length

# =======================
# PROCESS SIMULATION DIRECTORIES
# =======================
sim_data = {}  # key: label, value: dict with 'labels' and 'xmin_values'

for dir_name, plot_label in zip(directoryNames, PlotLabelNames):
    dir_path = os.path.join(base_path, dir_name)
    if not os.path.isdir(dir_path):
        print(f"{dir_name} not found, skipping...")
        continue

    # Collect valid cases
    case_folders = [
        f for f in os.listdir(dir_path)
        if os.path.isdir(os.path.join(dir_path, f))
    ]

    valid_cases = []
    for folder in case_folders:
        full_path = os.path.join(dir_path, folder)
        entries = os.listdir(full_path)

        has_processor = any(d.startswith("processor") for d in entries)
        has_time = any(is_time_dir(d) for d in entries)

        if not (has_processor or has_time):
            continue

        run_id = "_".join(folder.split("_")[-2:])
        if run_id not in experiment_data:
            continue

        row = experiment_data[run_id]

        valid_cases.append({
            "folder": folder,
            "run": run_id,
            "R": float(row["R"]),
            "Fo": float(row["Fo"]),
            "K": float(row["K"])
        })

    # Sort
    valid_cases.sort(key=lambda c: c[SORT_FIELD])

    # Data collection
    labels = []
    xmin_values = []

    for case in valid_cases:
        folder = case["folder"]
        caseFolder = os.path.join(dir_path, folder)

        print(f"\nProcessing {dir_name}/{folder}")

        # --- OpenFOAM reader ---
        reader = OpenFOAMReader(FileName=os.path.join(caseFolder, "case.foam"))
        reader.MeshRegions = ['internalMesh']
        reader.CellArrays = ['alpha.sludge']
        reader.CaseType = 'Decomposed Case'

        reader.UpdatePipelineInformation()
        timesteps = reader.TimestepValues
        if not timesteps:
            continue

        latest_time = timesteps[-1]

        clip = Clip(Input=reader)
        clip.ClipType = 'Scalar'
        clip.Scalars = ['POINTS', 'alpha.sludge']
        clip.Value = 0.5
        clip.Invert = 0

        reader.UpdatePipeline(time=latest_time)
        clip.UpdatePipeline(time=latest_time)

        bounds = clip.GetDataInformation().GetBounds()
        xmin = bounds[0]

        if xmin > 1e290:
            continue

        # --- Experimental length ---
        Lfile = os.path.join(caseFolder, "experimentalLength.txt")
        if os.path.exists(Lfile):
            with open(Lfile, "r") as f:
                L_value = float(f.read().strip())
                exp_lengths_dict[case["run"]] = L_value

        labels.append(str(int(case[SORT_FIELD])))
        xmin_values.append(7.34 - xmin)

    # Store data for plotting
    sim_data[plot_label] = {
        "labels": labels,
        "xmin_values": xmin_values
    }

# =======================
# SAVE DATA TO CSV
# =======================
out_csv = os.path.join(base_path, "interface_positions.csv")

with open(out_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        SORT_FIELD,
        "Simulation",
        "InterfacePosition",
        "Experimental"
    ])

    # Loop over x-labels (R values)
    for idx, r_label in enumerate(labels):
        # Experimental value for this R
        exp_val = None
        for run_id, row in experiment_data.items():
            if str(int(row[SORT_FIELD])) == r_label:
                exp_val = exp_lengths_dict.get(run_id, None)
                break

        # Each simulation contributes one row
        for sim_label, data in sim_data.items():
            sim_val = data["xmin_values"][idx]
            writer.writerow([
                r_label,
                sim_label,
                sim_val,
                exp_val
            ])

print(f"Saved plotted data → {out_csv}")

# =======================
# PLOT ALL SIMULATIONS WITH DISTINCT COLORS
# =======================
plt.figure(figsize=(10, 5))

x = range(len(labels))  # assume all dirs share same labels

# Define distinct colors for each simulation
sim_colors = ["tab:blue", "tab:green", "tab:orange"]  # one per simulation
width = 0.25

for i, ((label, data), color) in enumerate(zip(sim_data.items(), sim_colors)):
    offset = (i - len(sim_data)/2) * width + width/2
    plt.bar(
        [xi + offset for xi in x],
        data["xmin_values"],
        width=width,
        color=color,
        label=label
    )

# Experimental points only (black)
exp_y = [
    exp_lengths_dict.get(run_id, None)
    for run_id in experiment_data.keys()
    if str(int(experiment_data[run_id][SORT_FIELD])) in labels
]

plt.scatter(
    x,          # x positions
    exp_y,      # y positions
    color="black",
    marker="o",
    s=60,
    zorder=10,          # <-- key line
    edgecolors="white",
    linewidths=0.8,
    label="Experiment"
)

plt.xticks(x, labels)
plt.xlabel(SORT_FIELD)
plt.ylabel("Interface position")
plt.title("Interface position at latest time (multiple simulations)")
plt.grid(axis="y", alpha=0.3)
plt.legend()
plt.tight_layout()

out_png = os.path.join(base_path, "xmin_bar_chart_multi.png")
plt.savefig(out_png, dpi=300)
plt.close()

print(f"\nSaved multi-simulation bar chart → {out_png}")
