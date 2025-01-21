import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.size"] = 16
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["lines.linewidth"] = 1.5
plt.rcParams["lines.markersize"] = 1
plt.rcParams["axes.grid"] = True
plt.rcParams["axes.grid.which"] = "major"
plt.rcParams["grid.linestyle"] = "--"
plt.rcParams["grid.linewidth"] = 1
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.major.size"] = 6.0
plt.rcParams["ytick.major.size"] = 6.0
plt.rcParams["xtick.minor.size"] = 3.0
plt.rcParams["ytick.minor.size"] = 3.0
plt.rcParams["figure.autolayout"] = True
plt.rcParams['savefig.dpi'] = 600


groups = ("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6")
data = {
'IAEA-ORIGEN Yields': [1.84, 1.49, 8.72, 2.85, 2.26, 0.51],
'IAEA-ORIGEN Half-lives':	[1.52, 1.07, 7.32, 3.42, 2.03, 2.53],
'Pure ORIGEN Yields':     [1.11, 1.76, 2.49, 0.14, 3.41, 5.92],
'Pure ORIGEN Half-lives':	[0.88, 0.42, 0.09, 0.16, 0.16, 4.36]
}

x = np.arange(len(groups))  # the label locations
width = 1/5  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')
patterns = ["-", "o"]
for attribute, measurement in data.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute, 
                   hatch=patterns[multiplier%len(patterns)])
    #ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Absolute Percent Difference [%]')
ax.set_xticks(x + 1.5*width, groups, fontsize=12, weight='bold')
ax.legend(fontsize=12)
plt.tight_layout()

plt.savefig('group_bar.png')
plt.close()