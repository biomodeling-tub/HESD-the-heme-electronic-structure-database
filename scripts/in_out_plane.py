import os
import numpy as np

def read_xyz(file_path):
    atoms = []
    coords = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
    for line in lines[2:]:  # Skip the first two lines which contain header information
        parts = line.split()
        if len(parts) >= 4:
            atom, x, y, z = parts[:4]
            atoms.append(atom)
            coords.append([float(x), float(y), float(z)])

    return atoms, np.array(coords)


def save_xyz(file_path, atoms, coords):
    with open(file_path, 'w') as f:
        f.write(f"{len(atoms)}\n\n")
        for atom, coord in zip(atoms, coords):
            f.write(f"{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")

def rotate_out_plane(coords, fe_coords, angle_deg):
    angle_rad = np.deg2rad(angle_deg)
    rotation_matrix = np.array([[np.cos(angle_rad), -np.sin(angle_rad), 0],
                                [np.sin(angle_rad), np.cos(angle_rad), 0],
                                [0, 0, 1]])
    # Subtract Fe-Atom coordinates to perform rotation around it
    coords_shifted = coords - fe_coords
    # Rotate the coordinates
    rotated_coords = np.dot(coords_shifted, rotation_matrix.T) + fe_coords
    return rotated_coords


def find_rotation_axis(coords, fe_index):
    """Finds the rotation axis along the line between atoms 2 and 4 and perpendicular to the plane formed by this line and the Fe atom."""
    atom_2 = coords[1]
    atom_4 = coords[3]
    fe_atom = coords[fe_index]
    
    # Define vectors representing the line between atoms 2 and 4 and from Fe to the midpoint of 2 and 4
    vector_2_to_4 = atom_4 - atom_2
    midpoint_2_4 = (atom_2 + atom_4) / 2
    vector_fe_to_midpoint = midpoint_2_4 - fe_atom
    
    # Calculate cross product to find a vector perpendicular to the plane formed by the line and the Fe atom
    rotation_axis = np.cross(vector_2_to_4, vector_fe_to_midpoint)
    return rotation_axis


def rotate_atoms(coords, angle_degrees, axis, fe_coords):
    """Rotates atoms around the Fe atom along the specified axis."""
    rotation_matrix = rotation_matrix_from_axis_and_angle(axis, np.radians(angle_degrees))
    rotated_coords = []
    for coord in coords:
        rotated_coord = np.dot(rotation_matrix, coord - fe_coords) + fe_coords
        rotated_coords.append(rotated_coord)
    return np.array(rotated_coords)


def rotation_matrix_from_axis_and_angle(axis, angle):
    """Compute rotation matrix from axis and angle."""
    axis = axis / np.linalg.norm(axis)
    a = np.cos(angle / 2.0)
    b, c, d = -axis * np.sin(angle / 2.0)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])



def outplane(directory):
    file_pathes = [file for file in os.listdir(directory) if file.endswith('_g16.xyz')]
    for file_path in file_pathes:
        try:
            atoms, coords = read_xyz(os.path.join(directory, file_path))
        except (IOError, EOFError) as e:
            print(f"Error reading the XYZ file {file_path} in directory {directory}: {e}")
            continue

        # Define the normal vector of the plane using the first two nitrogen atoms and the Fe-Atom
        normal_vector = np.cross(coords[1] - coords[0], coords[2] - coords[0])

        # Find Fe-Atom's coordinates
        fe_index = atoms.index("Fe")
        fe_coords = coords[fe_index]

        # Define rotation angles
        rotation_angles = np.random.permutation(np.arange(360))[:10]
        #rotation_angles = range(0, 361, 90)

        for rotation_angle in rotation_angles:
            # Rotate all atoms in the plane
            rotated_coords = rotate_out_plane(coords, fe_coords, rotation_angle)
            # Save rotated coordinates to a new XYZ file
            rotated_file_name = f"{file_path[:4]}_{rotation_angle}_out_plane.xyz"
            rotated_file_path = os.path.join(directory, rotated_file_name)
            save_xyz(rotated_file_path, atoms, rotated_coords)
            print(f"Rotation out of plane completed and saved to: {rotated_file_name}")

def inplane(directory):
    input_files = [file for file in os.listdir(directory) if file.endswith('_out_plane.xyz')]
    for input_file in input_files:
        try:
            # Read XYZ file
            atoms, coords = read_xyz(os.path.join(directory, input_file))
        except (IOError, EOFError) as e:
            print(f"Error reading the XYZ file {input_file} in directory {directory}: {e}")
            continue

        # Define x-axis (between atoms 2 and 4) and y-axis (between atoms 3 and 5)
        x_axis = coords[3] - coords[1]
        y_axis = coords[4] - coords[2]

        # Define z-axis as the cross product of x-axis and y-axis
        z_axis = np.cross(x_axis, y_axis)

        # Find the index of Fe atom
        fe_index = atoms.index("Fe")

        # Get Fe atom coordinates
        fe_coords = coords[fe_index]

        # Generate random rotation angles for each rotation
        rotation_angles_in = np.random.permutation(np.arange(360))[:10]

        # Perform rotations
        for angle in rotation_angles_in:
            rotated_coords = rotate_atoms(coords, angle, z_axis, fe_coords)
            coords = rotated_coords
            output_file = f"{input_file[:-4]}_{angle}_in_plane.xyz"
            output_file_path = os.path.join(directory, output_file)
            # Write rotated coordinates to output file
            save_xyz(output_file_path, atoms, coords)
            print(f"Rotation in plane completed and saved to: {output_file}")


        
for root, dirs, files in os.walk("../../PDB/"):
    if root != "../../PDB/":
        outplane(root)
        inplane(root)