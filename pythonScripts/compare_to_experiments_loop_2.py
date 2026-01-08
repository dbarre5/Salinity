import os
import csv

# =======================
# COLOR SETTINGS
# =======================
# 5 Viridis colors (RGB 0–1), evenly spaced
VIRIDIS_COLORS = [
    [0.267, 0.005, 0.329],  # dark purple
    [0.230, 0.322, 0.545],  # blue
    [0.128, 0.566, 0.551],  # teal
    [0.369, 0.788, 0.383],  # green
    [0.993, 0.906, 0.144],  # yellow
]

colorSelection = 0

# =======================
# USER SETTINGS
# =======================
try:
    base_path = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # Fallback for ParaView GUI
    base_path = os.getcwd()

csv_file = os.path.join(base_path, "experiments.csv")

# =======================
# HELPERS
# =======================
def is_time_dir(name):
    """Return True if name is numeric AND > 0."""
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
        run = row["Run"].strip()
        experiment_data[run] = row

# =======================
# SCAN CASE FOLDERS
# =======================
case_folders = [
    f for f in os.listdir(base_path)
    if os.path.isdir(os.path.join(base_path, f))
]

valid_cases = []

for folder in case_folders:
    full_path = os.path.join(base_path, folder)
    entries = os.listdir(full_path)

    # Rule 1: decomposed case
    has_processor = any(
        d.startswith("processor") and os.path.isdir(os.path.join(full_path, d))
        for d in entries
    )

    # Rule 2: serial case with time > 0
    has_time = any(
        is_time_dir(d) and os.path.isdir(os.path.join(full_path, d))
        for d in entries
    )

    if not (has_processor or has_time):
        print(f"Skipping {folder}: no processor dirs and no time dirs > 0.")
        continue

    # Extract experiment ID (last two underscore tokens)
    parts = folder.split("_")
    run_id = "_".join(parts[-2:])   # e.g. 4_10a

    if run_id not in experiment_data:
        print(f"Skipping {folder}: run {run_id} not found in experiments.csv")
        continue

    row = experiment_data[run_id]

    valid_cases.append({
        "folder": folder,
        "run": run_id,
        "Re": float(row["R"]),
        "Fo": float(row["Fo"]),
        "K": float(row["K"])
    })

# =======================
# SORT CASES BY PHYSICS
# =======================
SORT_FIELD = "Re"   # "Re", "Fo", or "K"

valid_cases.sort(key=lambda c: c[SORT_FIELD])

print(f"Valid cases (sorted by {SORT_FIELD}):")
for c in valid_cases:
    value = c[SORT_FIELD]
    print(f"{c['folder']}  {SORT_FIELD}={value}")


# =======================
# PARAVIEW PIPELINE
# =======================
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.InteractionMode = '2D'

translationIndex = 0

for case in valid_cases:

    folder = case["folder"]
    run_id = case["run"]

    sort_value = case[SORT_FIELD]

    caseFolder = os.path.join(base_path, folder)

    # --- Read experimental length ---
    Lfile = os.path.join(caseFolder, "experimentalLength.txt")
    if not os.path.exists(Lfile):
        print(f"{folder}: experimentalLength.txt not found, skipping.")
        continue

    with open(Lfile, "r") as f:
        L_value = float(f.read().strip())



    # --- Translation offset ---
    translation = 0.2 * translationIndex
    translationIndex += 1

    # --- Load OpenFOAM case ---
    fileName = os.path.join(caseFolder, "case.foam")
    casefoam = OpenFOAMReader(FileName=fileName)
    casefoam.MeshRegions = ['internalMesh']
    casefoam.CellArrays = [
        'cellRegion', 'U', 'alpha.sludge', 'ddt0(alpha.sludge)',
        'ddt0(epsilon)', 'ddt0(k)', 'ddt0(rho,U)', 'ddtCorrDdt0(U)',
        'epsilon', 'k', 'multU', 'nut', 'p', 'p_rgh'
    ]
    casefoam.CaseType = 'Reconstructed Case'

    # --- Transform ---
    transform1 = Transform(Input=casefoam)
    transform1.Transform.Translate = [0.0, 0.0, -translation]

    casefoamDisplay = Show(transform1, renderView1)
    casefoamDisplay.SetRepresentationType('Feature Edges')
    ColorBy(casefoamDisplay, None)
    casefoamDisplay.DiffuseColor = [0.49, 0.49, 0.49]
    casefoamDisplay.AmbientColor = [0.49, 0.49, 0.49]

    # --- Clip alpha ---
    clip1 = Clip(Input=transform1)
    clip1.ClipType = 'Scalar'
    clip1.Scalars = ['POINTS', 'alpha.sludge']
    clip1.Value = 0.5
    clip1.Invert = 0

    clip1Display = Show(clip1, renderView1)
    clip1Display.Representation = 'Feature Edges'
    clip1Display.LineWidth = 2.0 
    ColorBy(clip1Display, None)
    clip1Display.DiffuseColor = VIRIDIS_COLORS[colorSelection]
    clip1Display.AmbientColor = VIRIDIS_COLORS[colorSelection]
    clip1Display.SetScalarBarVisibility(renderView1, False)

    # --- Clip alpha ---
    clip2 = Clip(Input=transform1)
    clip2.ClipType = 'Scalar'
    clip2.Scalars = ['POINTS', 'alpha.sludge']
    clip2.Value = 0.5
    clip2.Invert = 0

    clip2Display = Show(clip2, renderView1)
    clip2Display.Representation = 'Surface'
    ColorBy(clip2Display, None)
    clip2Display.DiffuseColor = VIRIDIS_COLORS[colorSelection]
    clip2Display.AmbientColor = VIRIDIS_COLORS[colorSelection]
    clip2Display.Opacity = 0.3
    clip2Display.SetScalarBarVisibility(renderView1, False)


    # --- Experimental line ---
    line1 = Line(
        Point1=[7.33, 0.0, 0.1 - translation],
        Point2=[7.33 - L_value, 0.0, 0.1 - translation]
    )
    line1Display = Show(line1, renderView1)
    line1Display.LineWidth = 3.0
    line1Display.DiffuseColor = [0.666667, 0.0, 0.0]
    line1Display.AmbientColor = [0.666667, 0.0, 0.0]

    # --- Text annotation ---
    value = case[SORT_FIELD]

    value = case[SORT_FIELD]

    # Decide formatting based on field type
    if SORT_FIELD in ("Re", "K"):
        value_str = f"{int(value)}"       # integer for large numbers
    else:
        value_str = f"{value:.3f}"       # 3 decimal places for small floats

    text1 = Text(
        Text=f"{SORT_FIELD} = {value_str}"
    )
    
    text1Display = Show(text1, renderView1)
    text1Display.TextPropMode = 'Billboard 3D Text'
    text1Display.BillboardPosition = [7.35, 0.0, 0.055 - translation]

# --- Camera ---
renderView1.CameraPosition = [6.686899, -3.804107, -0.377216]
renderView1.CameraFocalPoint = [6.686899, -1.5, -0.377216]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 1.019187
