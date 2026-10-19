from collections import defaultdict
from statistics import mean, stdev

INPUT_FILE = "NERO_HISTchi2.txt"
OUTPUT_FILE = "NERO_HISTchi2_outliers.txt"
SIGMA_THRESHOLD = 5.0

# Store values for each channel
channel_values = defaultdict(list)

# Store all measurements
measurements = []

# Read file
with open(INPUT_FILE) as f:
    for line in f:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        run, channel, value = line.split()

        run = int(run)
        channel = int(channel)
        value = float(value)

        measurements.append((run, channel, value))
        channel_values[channel].append(value)

# Calculate mean and sigma for each channel
channel_mean = {}
channel_sigma = {}

for ch, values in channel_values.items():
    channel_mean[ch] = mean(values)
    if len(values) > 1:
        channel_sigma[ch] = stdev(values)
    else:
        channel_sigma[ch] = 0.0

# Find outliers
outliers = []

for run, channel, value in measurements:
    sigma = channel_sigma[channel]

    if sigma == 0:
        continue

    deviation = value - channel_mean[channel]

    if deviation > SIGMA_THRESHOLD * sigma:
        z = deviation / sigma
        outliers.append((run, channel, value, channel_mean[channel], sigma, z))

# Save results
with open(OUTPUT_FILE, "w") as f:
    f.write("# Run Channel Value Mean Sigma Zscore\n")

    for o in outliers:
        f.write(
            f"{o[0]} {o[1]} {o[2]:.6f} {o[3]:.6f} {o[4]:.6f} {o[5]:.2f}\n"
        )

print(f"Found {len(outliers)} outliers.")
print(f"Results written to {OUTPUT_FILE}")
