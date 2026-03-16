import matplotlib.pyplot as plt
import itertools

# Define the output file paths for the two files
output_file_path1 = "Benchmark/Benchmark_logs/Functional/B3LYP/homos_and_lumos_01.txt"
output_file_path2 = "Benchmark/Benchmark_logs/Functional/BP86/homos_and_lumos_01.txt"

# Create empty lists to store data
file_names = []
homos_values = []

# Read data from the first output file
with open(output_file_path1, 'r') as output_file:
    lines = output_file.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("File: "):
            file_name = line[len("File: "):]
            file_name = file_name.replace(".log", "")  # Remove ".log" extension
            file_names.append(file_name)
            i += 1
            homo_line = lines[i + 1].strip()
            homos_values.append(list(map(float, homo_line.split(','))))
        i += 1

# Create a list of unique colors for plotting
colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

# Plot the data from the first file
fig, ax = plt.subplots(figsize=(10, 6))

for i, file_name in enumerate(file_names):
    color = next(colors)
    ax.plot([i] * len(homos_values[i]), homos_values[i], 'o', label=file_name, color=color)

# Read data from the second output file
file_names = []  # Reset file_names list
homos_values = []  # Reset homos_values list

with open(output_file_path2, 'r') as output_file:
    lines = output_file.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("File: "):
            file_name = line[len("File: "):]
            file_name = file_name.replace("01.log", "")  # Remove ".log" extension
            file_names.append(file_name)
            i += 1
            homo_line = lines[i + 1].strip()
            homos_values.append(list(map(float, homo_line.split(','))))
        i += 1

# Plot the data from the second file with different colors
for i, file_name in enumerate(file_names):
    color = next(colors)
    ax.plot([i] * len(homos_values[i]), homos_values[i], 'o', label="linestyle='--' belongs to BP86", color=color, linestyle='--')

ax.set_xlabel("File Name")
ax.set_ylabel("HOMO Values")
ax.set_title("HOMO Values for Each File")
ax.set_xticks(range(len(file_names)))
ax.set_xticklabels(file_names, rotation=45, fontsize=8)
ax.legend()

plt.tight_layout()
plt.show()
