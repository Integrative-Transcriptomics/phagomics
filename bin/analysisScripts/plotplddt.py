import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def collectPlddtValuesByCategory(folderPath):
    knownValues = []
    unknownValues = []

    for filename in os.listdir(folderPath):
        if filename.endswith(".json"):
            filepath = os.path.join(folderPath, filename)
            try:
                with open(filepath, 'r') as f:
                    data = json.load(f)
                    plddtValues = data.get("plddt", [])

                    if "_known" in filename:
                        knownValues.extend(plddtValues)
                    elif "_unknown" in filename:
                        unknownValues.extend(plddtValues)
            except Exception as e:
                print(f"Error reading {filename}: {e}")

    return knownValues, unknownValues

def plotAreaByCategory(ax, known, unknown, title=""):
    print(np.average(unknown + known))

    if known:
        kdeKnown = gaussian_kde(known)
        xKnown = np.linspace(min(known), max(known), 500)
        yKnown = kdeKnown(xKnown)
        ax.fill_between(xKnown, yKnown, alpha=0.5, label="Known", color="#24b064")
        ax.plot(xKnown, yKnown, color="#24b064")

    if unknown:
        kdeUnknown = gaussian_kde(unknown)
        xUnknown = np.linspace(min(unknown), max(unknown), 500)
        yUnknown = kdeUnknown(xUnknown)
        ax.fill_between(xUnknown, yUnknown, alpha=0.5, label="Unknown", color="#F18046")
        ax.plot(xUnknown, yUnknown, color="#F18046")

    ax.set_xlabel("pLDDT Score", fontsize=14)
    ax.set_ylabel("Density", fontsize=14)
    ax.tick_params(axis='y', labelsize=12)
    ax.tick_params(axis='x', labelsize=12)
    # ax.set_title(title)
    # ax.legend()
    ax.grid(True)

fig, axes = plt.subplots(1, 2, figsize=(12, 6))

knownData1, unknownData1 = collectPlddtValuesByCategory("resultsT4/colabfold")
plotAreaByCategory(axes[0], knownData1, unknownData1, title="resultsT4")

knownData2, unknownData2 = collectPlddtValuesByCategory("resultsFull/colabfold")
plotAreaByCategory(axes[1], knownData2, unknownData2, title="resultsFull")

# Add panel labels
axes[0].text(-0.1, 1.05, "A", transform=axes[0].transAxes, fontweight="bold",
             fontsize=16, va='top', ha='left')
axes[0].legend(fontsize=16)
axes[0].set_title("T4", fontsize=16)

axes[1].text(-0.1, 1.05, "B", transform=axes[1].transAxes, fontweight="bold",
             fontsize=16, va='top', ha='left')
axes[1].set_title("BASEL collection", fontsize=16)

plt.tight_layout()
plt.savefig("figures/plddtPlots.pdf", bbox_inches='tight')
plt.show()

# Single plot example (commented)
# fig, axes = plt.subplots(1, 1, figsize=(6, 6))

# knownData1, unknownData1 = collectPlddtValuesByCategory("resultsT4/foldseek/validate/notHit")
# plotAreaByCategory(axes, knownData1, unknownData1, title="resultsT4")

# axes.text(-0.1, 1.05, "A", transform=axes.transAxes,
#           fontsize=16, fontweight='bold', va='top', ha='left')

# plt.tight_layout()
# plt.show()
