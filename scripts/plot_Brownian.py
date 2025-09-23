import numpy as np
import matplotlib.pyplot as plt


def simulate_brownian_motion(n_steps, sigma, n_replicates=10):
    """
    Simulate Brownian motion with given parameters
    """
    time = np.arange(n_steps + 1)
    brownian_paths = []

    for _ in range(n_replicates):
        # Generate random increments
        increments = np.random.normal(0, sigma, n_steps)
        # Start at 0 and cumulative sum
        path = np.concatenate([[0], np.cumsum(increments)])
        brownian_paths.append(path)

    return time, brownian_paths


def simulate_substitutions(time, rate=0.05):
    """
    Simulate substitution events as Poisson process
    """
    substitution_times = []
    current_time = 0

    while current_time < time[-1]:
        # Exponential waiting time between events
        wait_time = np.random.exponential(1 / rate)
        current_time += wait_time
        if current_time < time[-1]:
            substitution_times.append(current_time)

    return substitution_times


# Set random seed for reproducibility
np.random.seed(42)

# Parameters
n_steps = 100
n_replicates = 10

# Long generation time (smaller sigma) - 10x smaller than short generation time
sigma_long = 0.15
time_long, paths_long = simulate_brownian_motion(n_steps, sigma_long, n_replicates)

# Short generation time (bigger sigma) - 10x bigger than long generation time
sigma_short = 0.5
time_short, paths_short = simulate_brownian_motion(n_steps, sigma_short, n_replicates)

# Create figure with 4 subplots stacked vertically with custom height ratios
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(16 * .6, 9 * .6),
                                         gridspec_kw={'height_ratios': [4, 1, 4, 1]})

# Plot 1: Long generation time Brownian motion (top)
ax1.set_title('Long Generation Time', fontsize=14)
color_long = 'blue'

for i, path in enumerate(paths_long):
    ax1.plot(time_long, path, color=color_long, alpha=0.7, linewidth=1.5)

ax1.set_xlabel('Time')
ax1.set_ylabel('Trait')
ax1.grid(False)
for y in np.linspace(-5, 15, 5):
    ax1.axhline(y, color='black', linestyle='-', linewidth=0.5, alpha=0.2)
for x in np.linspace(0, 100, 30):
    ax1.axvline(x, color='black', linestyle='-', linewidth=0.5, alpha=0.2)
ax1.axhline(y=0, color='black', linestyle='--', alpha=0.5)
ax1.set_xlim(0, 100)
# Plot 2: Substitutions for long generation time
# Generate fewer substitutions with lower rate
substitutions_long = simulate_substitutions(time_long, rate=0.02)  # Very low rate

# Draw horizontal line only - no frame, labels, or ticks
ax2.axhline(y=0, color='black', linewidth=2, alpha=0.8)

if substitutions_long:
    # Plot substitutions as dots on the horizontal line
    ax2.scatter(substitutions_long, [0] * len(substitutions_long),
                color='darkred', s=50, zorder=5)

# Remove all frame elements
ax2.set_xlim(ax1.get_xlim())  # Match x-axis with Brownian plot above
ax2.set_ylim(-0.5, 0.5)
ax2.axis('off')  # Remove all axes, ticks, and labels

# Plot 3: Short generation time Brownian motion
ax3.set_title('Short Generation Time', fontsize=14)
color_short = 'red'

for i, path in enumerate(paths_short):
    ax3.plot(time_short, path, color=color_short, alpha=0.7, linewidth=1.5)

ax3.set_xlabel('Time')
ax3.set_ylabel('Trait')
ax3.grid(False)
for y in np.linspace(-5, 15, 5):
    ax3.axhline(y, color='black', linestyle='-', linewidth=0.5, alpha=0.2)
for x in np.linspace(0, 100, 90):
    ax3.axvline(x, color='black', linestyle='-', linewidth=0.5, alpha=0.2)
ax3.axhline(y=0, color='black', linestyle='--', alpha=0.5)
ax3.set_xlim(0, 100)

# Make both Brownian plots use the same y-axis limits
y_min = min(min(min(path) for path in paths_long), min(min(path) for path in paths_short))
y_max = max(max(max(path) for path in paths_long), max(max(path) for path in paths_short))
y_range = y_max - y_min
y_margin = y_range * 0.1  # Add 10% margin
ax1.set_ylim(y_min - y_margin, y_max + y_margin)
ax3.set_ylim(y_min - y_margin, y_max + y_margin)

# Plot 4: Substitutions for short generation time
# Generate fewer substitutions with slightly higher rate
substitutions_short = simulate_substitutions(time_short, rate=0.2)  # Low rate

# Draw horizontal line only - no frame, labels, or ticks
ax4.axhline(y=0, color='black', linewidth=2, alpha=0.8)

if substitutions_short:
    # Plot substitutions as dots on the horizontal line
    ax4.scatter(substitutions_short, [0] * len(substitutions_short),
                color='darkred', s=50, zorder=5)

# Remove all frame elements
ax4.set_xlim(ax3.get_xlim())  # Match x-axis with Brownian plot above
ax4.set_ylim(-0.5, 0.5)
ax4.axis('off')  # Remove all axes, ticks, and labels

plt.tight_layout()

# Save as PDF
plt.savefig('../manuscript/figures/brownian_motion_comparison.pdf',
            format='pdf', dpi=300, bbox_inches='tight')
