from pathlib import Path
import numpy as np

# Configuration
N_JITTER = 1000     # how many jittered copies per file
STD_DEV   = 0.02    # standard deviation of the Gaussian noise
SKIP_ATOMS = 5      # number of leading atoms to leave untouched

def read_xyz(path):
    """Read atoms and coordinates from an .xyz file."""
    lines = path.read_text().splitlines()
    atoms = []
    coords = []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) >= 4:
            atoms.append(parts[0])
            coords.append([float(c) for c in parts[1:4]])
    return atoms, np.array(coords)

def write_xyz(path, atoms, coords):
    """Write atoms and coords to an .xyz file."""
    with path.open('w') as f:
        f.write(f"{len(atoms)}\n\n")
        for atom, (x, y, z) in zip(atoms, coords):
            f.write(f"{atom} {x:.6f} {y:.6f} {z:.6f}\n")

def jitter_directory(directory):
    """For each *_plane.xyz file under `directory`, produce N_JITTER jittered copies."""
    directory = Path(directory)
    for infile in directory.rglob("*_g16.xyz"):
        atoms, coords = read_xyz(infile)
        # Prepare slice of coords that will be noised
        fixed, to_jitter = coords[:SKIP_ATOMS], coords[SKIP_ATOMS:]
        for i in range(1, N_JITTER + 1):
            noise = np.random.normal(0, STD_DEV, size=to_jitter.shape)
            jittered = np.vstack([fixed, to_jitter + noise])
            outfile = infile.with_name(
                infile.stem + f"_jittered_{i}.xyz"
            )
            write_xyz(outfile, atoms, jittered)

if __name__ == "__main__":
    base = Path("PDB")
    for subdir in base.iterdir():
        if subdir.is_dir():
            jitter_directory(subdir)
