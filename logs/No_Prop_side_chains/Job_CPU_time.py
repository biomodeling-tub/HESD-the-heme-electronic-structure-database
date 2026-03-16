import os

# Get the current directory
current_directory = os.getcwd()

# Name of the output file
output_file_name = "Job_cpu_time_01.txt"

# Open the output file in write mode
with open(output_file_name, "w") as output_file:
    # Iterate through all files in the current directory
    for filename in os.listdir(current_directory):
        if os.path.isfile(filename):
            try:
                # Open the current file
                with open(filename, "r") as current_file:
                    # Loop through each line in the current file
                    for line in current_file:
                        # Check if the line contains "Job cpu time: "
                        if "Job cpu time: " in line:
                            # Write the line to the output file
                            output_file.write(f"{filename}: {line}")
            except FileNotFoundError:
                print(f"File {filename} not found.")

print(f"Lines containing 'Job cpu time:' from all files in the directory have been copied to {output_file_name}.")
