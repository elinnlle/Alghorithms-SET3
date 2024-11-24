import numpy as np
import matplotlib.pyplot as plt

# Function to check if a point is in a circle
def is_in_circle(cx, cy, r, x, y):
    return (x - cx) ** 2 + (y - cy) ** 2 <= r ** 2

# Monte Carlo simulation function
def monte_carlo_area(a, b, c, narrow_bounds, n_points):
    # Define bounds for random point generation
    if narrow_bounds:
        x_min = max(a[0] - a[2], b[0] - b[2], c[0] - c[2])
        x_max = min(a[0] + a[2], b[0] + b[2], c[0] + c[2])
        y_min = max(a[1] - a[2], b[1] - b[2], c[1] - c[2])
        y_max = min(a[1] + a[2], b[1] + b[2], c[1] + c[2])
    else:
        x_min = min(a[0] - a[2], b[0] - b[2], c[0] - c[2])
        x_max = max(a[0] + a[2], b[0] + b[2], c[0] + c[2])
        y_min = min(a[1] - a[2], b[1] - b[2], c[1] - c[2])
        y_max = max(a[1] + a[2], b[1] + b[2], c[1] + c[2])

    # Random points generation
    x_rand = np.random.uniform(x_min, x_max, n_points)
    y_rand = np.random.uniform(y_min, y_max, n_points)

    # Count points inside all three circles
    inside_count = np.sum(
        is_in_circle(a[0], a[1], a[2], x_rand, y_rand) &
        is_in_circle(b[0], b[1], b[2], x_rand, y_rand) &
        is_in_circle(c[0], c[1], c[2], x_rand, y_rand)
    )

    # Calculate estimated area
    rect_area = (x_max - x_min) * (y_max - y_min)
    estimated_area = (inside_count / n_points) * rect_area

    return estimated_area

# True area for comparison (from provided analytical value)
true_area = 0.25 * np.pi + 1.25 * np.arcsin(0.8) - 1.0

# Circles' parameters: (x, y, r)
circle_a = (1.0, 1.0, 1.0)
circle_b = (1.5, 2.0, np.sqrt(5) / 2)
circle_c = (2.0, 1.5, np.sqrt(5) / 2)

# Number of points for Monte Carlo and ranges
n_points_range = range(100, 100001, 500)

# Results storage
results_narrow = []
results_wide = []

# Run experiments
for n_points in n_points_range:
    area_narrow = monte_carlo_area(circle_a, circle_b, circle_c, True, n_points)
    area_wide = monte_carlo_area(circle_a, circle_b, circle_c, False, n_points)
    results_narrow.append((n_points, area_narrow, abs(area_narrow - true_area) / true_area))
    results_wide.append((n_points, area_wide, abs(area_wide - true_area) / true_area))

# Convert results to numpy arrays for easy handling
results_narrow = np.array(results_narrow)
results_wide = np.array(results_wide)

# Plotting
plt.figure(figsize=(14, 6))

# Graph 1: Estimated area vs N
plt.subplot(1, 2, 1)
plt.plot(results_narrow[:, 0], results_narrow[:, 1], label="Narrow Bounds", color='blue')
plt.plot(results_wide[:, 0], results_wide[:, 1], label="Wide Bounds", color='red')
plt.axhline(y=true_area, color='green', linestyle='--', label="True Area")
plt.xlabel("Number of Points (N)")
plt.ylabel("Estimated Area")
plt.title("Estimated Area vs N")
plt.legend()
plt.grid()

# Graph 2: Relative error vs N
plt.subplot(1, 2, 2)
plt.plot(results_narrow[:, 0], results_narrow[:, 2], label="Narrow Bounds", color='blue')
plt.plot(results_wide[:, 0], results_wide[:, 2], label="Wide Bounds", color='red')
plt.xlabel("Number of Points (N)")
plt.ylabel("Relative Error")
plt.title("Relative Error vs N")
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()
