def is_comment_line(line):
    return len(line.split()) != 11

def parse_ice_structure(file_path):
    frames = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    frame_lines = []
    for line in lines:
        if is_comment_line(line):
            if frame_lines:
                frames.append(frame_lines)
                frame_lines = []
            frame_lines.append(line.strip())
        else:
            frame_lines.append(line.strip())
    
    if frame_lines:
        frames.append(frame_lines)
    
    return frames

import re
def process_frame(frame):
    header = frame[0].strip().replace("= ", "=").replace("= ", "=")
    metadata = header.split()
    
    print(header, metadata)
    # Extract the dimensions and shifts
    dimensions = metadata[0]
    xl = float(metadata[1].split('=')[1])
    yl = float(metadata[2].split('=')[1])
    zl = float(metadata[3].split('=')[1])
    x_shift = float(metadata[4].split('=')[1])
    
    atoms = []
    
    for line in frame[2:]:  # Skip the header and the column names line
        parts = line.split()
        if len(parts) != 11:
            continue
        pattern = int(parts[0])
        atom_type = int(parts[1])
        ox, oy, oz = map(float, parts[2:5])
        h1x, h1y, h1z = map(float, parts[5:8])
        h2x, h2y, h2z = map(float, parts[8:11])
        
        # Add Oxygen atom
        atoms.append(f"O {ox:.5f} {oy:.5f} {oz:.5f} {pattern} {atom_type} 0")
        
        # Add first Hydrogen atom
        atoms.append(f"H {h1x:.5f} {h1y:.5f} {h1z:.5f} {pattern} {atom_type} 1")
        
        # Add second Hydrogen atom
        atoms.append(f"H {h2x:.5f} {h2y:.5f} {h2z:.5f} {pattern} {atom_type} 1")
    
    return atoms, xl, yl, zl, x_shift

def write_xyz(file_path, frames_data):
    with open(file_path, 'w') as file:
        for atoms, xl, yl, zl, x_shift in frames_data:
            num_atoms = len(atoms)
            # Write the number of atoms and a comment line
            file.write(f"{num_atoms}\n")
            file.write(f"Lattice=\"{xl} 0.0 0.0 0.0 {yl} 0.0 0.0 0.0 {zl}\" Properties=species:S:1:pos:R:3:pattern:I:1:atom_type:I:1:H_or_O:I:1\n")
            
            for atom in atoms:
                file.write(f"{atom}\n")

input_file_path = 'hayw-reim97jcp-cells.dat'  # Replace with your input file path
output_file_path = 'hayw-reim97jcp-cells.xyz'  # Replace with your desired output file path


frames = parse_ice_structure(input_file_path)
frames_data = [process_frame(frame) for frame in frames]
write_xyz(output_file_path, frames_data)

print(f"Conversion complete! The data has been written to {output_file_path}")

