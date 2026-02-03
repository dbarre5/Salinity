# =======================
# IMPORTS
# =======================
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # pvpython-safe
import matplotlib.pyplot as plt
import seaborn as sns

# =======================
# PORTABLE BASE PATH
# =======================
def get_base_path():
    try:
        # Works in scripts
        return os.path.dirname(os.path.abspath(__file__))
    except NameError:
        # Works in notebook / pvpython
        return os.getcwd()

base_path = get_base_path()
csv_file = os.path.join(base_path, "interface_positions.csv")

# =======================
# READ CSV
# =======================
df = pd.read_csv(csv_file)
df["R"] = df["R"].astype(int)

# Average duplicates to prevent seaborn from computing error bars
df = df.groupby(["R", "Simulation"], as_index=False).mean()

# Extract experimental values (assume same R has same Experimental)
exp_df = df.drop_duplicates("R")[["R", "Experimental"]].sort_values("R")

# =======================
# PLOT (seaborn style, no error bars)
# =======================
plt.figure(figsize=(10, 5))

# Seaborn barplot with edge outlines
bars = sns.barplot(
    x="R",
    y="InterfacePosition",
    hue="Simulation",
    data=df,
    palette="Blues",
    errorbar=None  # seaborn >=0.12
    # ci=None      # older seaborn
)

# Add dark blue outline to each bar
for bar in bars.patches:
    bar.set_edgecolor("darkblue")  # outline color
    bar.set_linewidth(0.5)         # outline thickness

# =======================
# EXPERIMENTAL POINTS
# =======================
x_positions = range(len(exp_df))
plt.scatter(
    x=x_positions,
    y=exp_df["Experimental"],
    color="black",
    s=60,
    zorder=10,         
    edgecolors="white",
    linewidths=0.8,
    label="Experiment"
)

# =======================
# FORMATTING
# =======================
plt.grid(False)
plt.xlabel("R")
plt.ylabel("Wedge Length (m)")
plt.title("Wedge Length with Ascending Reynolds Number Across Various C3 Values")
plt.legend(title="")  # remove legend title

# =======================
# SAVE PLOT
# =======================
out_png = os.path.join(base_path, "xmin_bar_chart_multi_seaborn.png")
plt.savefig(out_png, dpi=300)
plt.close()

print(f"Saved plot → {out_png}")
# =======================
# COMPUTE % ERROR
# =======================
df_error = df.copy()

# Merge experimental values into df_error
exp_map = exp_df.set_index("R")["Experimental"].to_dict()
df_error["Experimental"] = df_error["R"].map(exp_map)

# Compute percent error
df_error["PercentError"] = (df_error["InterfacePosition"] - df_error["Experimental"]) / df_error["Experimental"] * 100

# =======================
# PLOT % ERROR
# =======================
plt.figure(figsize=(10, 5))

bars = sns.barplot(
    x="R",
    y="PercentError",
    hue="Simulation",
    data=df_error,
    palette="Blues",
    errorbar=None  # seaborn >=0.12
)

# Dark blue outlines
for bar in bars.patches:
    bar.set_edgecolor("darkblue")
    bar.set_linewidth(0.5)

# Optional: add horizontal line at 0%
plt.axhline(0, color="black", linestyle="--", linewidth=1)

# Formatting
plt.grid(False)
plt.xlabel("R")
plt.ylabel("Percent Error (%)")
plt.title("Simulation Percent Error with Ascending Reynolds Number Across Various C3 Values")
plt.legend(title="")

# Save
out_png_error = os.path.join(base_path, "xmin_bar_chart_percent_error.png")
plt.savefig(out_png_error, dpi=300)
plt.close()

print(f"Saved percent error plot → {out_png_error}")
