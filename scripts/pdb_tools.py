import itertools

class file_operations():
    def __init__(self):
        return self

    def read_file(pdb_file):
        """
        file_operations.readfile() reads a .pdb file using the readlines method. It also filters out all lines containing the keywords
        REMARK, TER, TITLE, CRYST1, SCALE
        and returns the list containing the lines (strings)
        """
        f=open(pdb_file, 'r')
        lines = f.readlines()
        f.close()
        lines = [k for k in lines if "HEADER" not in k]
        lines = [k for k in lines if "TITLE" not in k]
        lines = [k for k in lines if "CRYST1" not in k]
        lines = [k for k in lines if "REMARK" not in k]
        lines = [k for k in lines if "SCALE" not in k]
        lines = [k for k in lines if "MODEL" not in k]
        lines = [k for k in lines if "ENDMDL" not in k]
        lines = [k for k in lines if "TER" not in k]
        lines = [k for k in lines if "END" not in k]
        lines = [k for k in lines if "ORIGX1" not in k]
        lines = [k for k in lines if "ORIGX2" not in k]
        lines = [k for k in lines if "ORIGX3" not in k]
        #lines = [k.replace("\n", '') for k in lines]
        return lines

    def write_file(file, lines):
        """
        file_operations.writefile takes a list of strings, e.g. from a read file like with the file_operations.readfile function and opens a writer.
        Using this write method, file_operations.writefile writes the lines into a file that is saved, and closes the write-method.
        It returns nothing.
        """
        #lines_ = line_operations.add_terminus(lines=lines)
        #if "END" not in lines_[-1]:
        #    lines.append("END")
        lines_ = line_operations.add_ending(lines=lines)
        f = open(file, mode="w", encoding="utf-8")
        f.writelines(lines_)
        f.close()
        return

class line_operations():
    """
    The class lines contains all the functions that operate on a line, which are iterated over in the operations functions.
    It relies on being handed a single line or line_dict object as input.
    """
    def __init__(self):
        return None

    def read_pdb_line(line):
        """
        line_operations.read_pdb_line() creates a dictionary (line_dict) which is filled with the content of the line it is given based on
        string indexing. Based on the typical .pdb format, line_dict then knows:
        atom, serial_no, atom_name, resname, chainID, resi_no, x_coord, y_coord, z_coord, occupancy,
        temp_fac, segment, element_symbol.
        line_operations.read_pdb_line() returns the dictionary line_dict.
        """
        line_dict = {
            "atom": line[0:6],
            "serial_no": line[6:11],
            "atom_name": line[12:16],
            "resname": line[17:21],
            "chainID": line[21],
            "resi_no": line[22:26],
            "ins_code": line[26],
            "x_coord": line[31:38],
            "y_coord": line[39:46],
            "z_coord": line[47:54],
            "occupancy": line[55:60],
            "temp_fac": line[60:66],
            "segment": line[72:76],
            "elem_symb": line[77:79],
            "charge": line[79:81]
        }
        if "\n" in line_dict["elem_symb"]:
            line_dict["elem_symb"] = line_dict["elem_symb"].replace("\n", "")
        line_dict["elem_symb"] = line_dict["elem_symb"].strip()
        if line_dict["elem_symb"] == "":
            line_dict["elem_symb"] = "  "
        if "\n" in line_dict["charge"]:
            line_dict["charge"] = line_dict["elem_symb"].replace("\n", "")
        if line_dict["charge"]=="":
            line_dict["charge"]="  "
        return line_dict

    def create_line(line_dict):
        """
        line_operations.create_line() takes a line_dict and creates the PDB-style line with the information contained in the dictionary.
        It returns "line", an object containing the string that was produced.
        """
        if len(line_dict["elem_symb"]) == 1:
            line = f'{line_dict["atom"]}{line_dict["serial_no"]} {line_dict["atom_name"]} {line_dict["resname"]}{line_dict["chainID"]}{line_dict["resi_no"]}{line_dict["ins_code"]}    {line_dict["x_coord"]} {line_dict["y_coord"]} {line_dict["z_coord"]} {line_dict["occupancy"]}{line_dict["temp_fac"]}      {line_dict["segment"]} {line_dict["elem_symb"]}      '
        if len(line_dict["elem_symb"]) ==2:
            line = f'{line_dict["atom"]}{line_dict["serial_no"]} {line_dict["atom_name"]} {line_dict["resname"]}{line_dict["chainID"]}{line_dict["resi_no"]}{line_dict["ins_code"]}    {line_dict["x_coord"]} {line_dict["y_coord"]} {line_dict["z_coord"]} {line_dict["occupancy"]}{line_dict["temp_fac"]}      {line_dict["segment"]}{line_dict["elem_symb"]}      '
        #line = f'{line_dict["atom"]}{line_dict["serial_no"]} {line_dict["atom_name"]} {line_dict["resname"]}{line_dict["chainID"]}{line_dict["resi_no"]}{line_dict["ins_code"]}   {line_dict["x_coord"]} {line_dict["y_coord"]} {line_dict["z_coord"]} {line_dict["occupancy"]} {line_dict["temp_fac"]}       {line_dict["segment"]} {line_dict["elem_symb"]}{line_dict["charge"]}\n'
        if len(line) > 82:
            line=line.strip()
        if len(line) < 82:
            line=f"{line: <82}"
        line = line + "\n"
        return line

    def fill_serial(serial_no: int, line_dict: dict):
        """
        line_operations.fill_serial() takes a serial number (serial_no) and a line_dict and creates line_dict["serial_no"] objects with
        the appropriate number of spaces inserted in front, so that the serial number is inserted at the place the .pdb-format
        dictates. It returns the line_dict with the appropriate serial_no.
        """
        if isinstance(serial_no, int):
            if serial_no >= 1e5:
                raise ValueError("Only serial numbers until 99.999 allowed. ")
            line_dict["serial_no"] = f"{serial_no: >5}"
        else:
            raise TypeError(f"The serial number {serial_no} is not an integer !!")
        return line_dict

    def fill_resi_no(resi_no, line_dict):
        """
        line_operations.fill_resi_no() takes a residue number (resi_no) and a line_dict and creates a line_dict with the serial
        number and the appropriate number of spaces inserted into line_dict["resi_no"]. It returns the line_dict.
        """
        if isinstance(resi_no, int):
            if resi_no >= 1e4:
                raise ValueError("Only residue numbers until 9.999 allowed. ")
            line_dict["resi_no"] = f"{resi_no: >4}"
        else:
            raise TypeError(f"The residue number {resi_no} is not an integer !!")
        return line_dict

    def add_terminus(lines):
        """
        line_operations.add_terminus() takes a list of strings (lines) and adds the string "TER" as the very last string in the list.
        """
        if lines[-1] != "TER":
            lines.append("TER")
        return lines

    def add_ending(lines):
        """
        line_operations.add_terminus() takes a list of strings (lines) and adds the string "TER" as the very last string in the list.
        """
        #print("add the ending:")
        #print(lines)
        if lines[-1] != "END":
            lines.append("END")
        return lines

    def exchange_segment(line_dict, segment):
        """
        line_operations.exchange_segment() takes a line_dict and a segment name; then exchanges the previous segment name with the new
        name and returns the line_dict with updated segment name. If the given segment name is too long, it will raise a
        ValueError. If the given segment name is too short, it will start filling them up with whitespaces from the left
        until the desired length of 4 characters is reached.
        """
        if len(segment) > 4:
            raise ValueError("segment value given is longer than 4. ")
        if len(segment) < 4:
            whitespaces=" "*(4-len(segment))
            segment=f'{whitespaces}{segment}'
        line_dict["segment"] = segment
        return line_dict

    def exchange_chainID(line_dict, chainID):
        """
        line_operations.exchange_chainID() takes a line_dict and a name for a chainID (maximum of 1 character) and exchanges the previous
        chainID with the new chainID. It then returns the line_dict with the updated chainID.
        """
        if len(chainID) > 1:
            raise ValueError("chainID value given is longer than 1. ")
        line_dict["chainID"] = chainID
        return line_dict

    def exchange_HETATM_ATOM(line_dict):
        line_dict["atom"] = "ATOM  "
        return line_dict

class operations():
    def __init__(self):
        return

    def _filter_segment(lines, segname):
        lines=[k for k in lines if segname in k]
        return lines

    def _split_segment(pdb_file, segname, pdb_id):
        """
        operations.split_segment() takes a pdb_file, a segname, and a pdb_id; then reads the file and drops all instances that
        do not have the segname within them. It writes the file as "coords/{pdb_id}_{segname}.pdb. Not having a folder "coords"
        will produce an error until the writefile function checks if the folder exists.
        """
        lines=file_operations.read_file(pdb_file=pdb_file)
        lines=operations._filter_segment(lines=lines, segname=segname)
        file_operations.write_file(file=f'PDB/{pdb_id}_{segname}.pdb', lines=lines)
        return

    def split_segments(pdb_file, segnames, pdb_id):
        """
        operations.split_segments() takes a pdb_file and a list of segname strings (segnames) and a pdb_id; then for each
        instance in segnames calls the operations.split_segment() function, producing a single file that contains only the
        lines in the pdb that had the segname in them.
        """
        for segname in segnames:
            operations._split_segment(pdb_file=pdb_file, segname=segname, pdb_id=pdb_id)
        return

    def split_waterchains(pdb_file, output_name):
        """
        operations.split_waterchains() takes a pdb_file and an output_name; first it reads the file (which should only contain
        H2O molecules in the pdb format and TIP3 water model/other water model that has EXACTLY three atoms in the water) and
        splits them into chains of 10.000 water molecules each. It writes the files based on the operations.renumber_tip3 method.
        """
        lines=file_operations.read_file(pdb_file=pdb_file)
        length, counter, filenames=len(lines), 0, []
        while counter < length:
            lines_=lines[0:29997]
            filename="coords/"+output_name+f"{counter//29997}.pdb"
            filenames.append(filename)
            file_operations.write_file(file=filename, lines=lines_)
            lines=lines[29997:]
            counter +=29997
        for filename in filenames:
            operations.renumber_tip3(pdb_file=filename, pdb_file_output=filename, segment=filename[12:16])
        return

    def fuse_segments(pdb_files, pdb_output):
        """
        operations.fuse_segments() takes a list of strings that point to .pdb files and a name for a pdb_output; then it
        appends all the files into one, removing potential "TER" lines, and writes a single fused file.
        """
        lines_=[]
        for pdb_file in pdb_files:
            lines=file_operations.read_file(pdb_file=pdb_file)
            lines_.append(lines)
            lines=[]
        lines_ = list(itertools.chain(*lines_))
        file_operations.write_file(file=pdb_output, lines=lines_)
        return

    def add_segment(pdb_file, pdb_file_output, segment):
        """
        operations.add_segment() takes a pdb_file, the name of a pdb_file_output, and a segment name (segment); then it
        calls the line_operations.exchange_segment() function to exchange the segment identifier in the line_dict. In the end it
        writes the file based on pdb_file_output.
        """
        lines=file_operations.read_file(pdb_file=pdb_file)
        lines_ = []
        serial_no=1
        for line in lines:
            line_dict=line_operations.read_pdb_line(line=line)
            line_dict=line_operations.exchange_segment(line_dict=line_dict, segment=segment)
            line_=line_operations.create_line(line_dict=line_dict)
            lines_.append(line_)
            serial_no+=1
        file_operations.write_file(file=pdb_file_output, lines=lines_)
        return

    def add_chainID(pdb_file, pdb_file_output, chainID):
        """
        operations.add_chainID() takes a pdb_file and a output name pdb_file_output and a segment name; then iterates over
        all lines in the PDB file changing the chainID. In the end it saves the new file according to pdb_file_output.
        """
        lines=file_operations.read_file(pdb_file=pdb_file)
        lines_ = []
        for line in lines:
            line_dict=line_operations.read_pdb_line(line=line)
            line_dict=line_operations.exchange_chainID(line_dict=line_dict, chainID=chainID)
            line_=line_operations.create_line(line_dict=line_dict)
            lines_.append(line_)
        file_operations.write_file(file=pdb_file_output, lines=lines_)
        return

    def change_temp_factors(pdb_file, restraints_file):
        """
        operations.change_temp_factors() takes a pdb_file and a name for a restraints_file; then it
        creates a restraints file based on the specifications of CHARMM-GUIs NAMD constraint files.
        H-atoms will recieve the temp_fac 0.00
        All C-atoms but CA (C-alphas) will recieve temp_fac 0.50
        All other atoms will recieve temp_fac 1.00.
        """
        lines=file_operations.read_file(pdb_file=pdb_file)
        lines_ = []
        for line in lines:
            line_dict=line_operations.read_pdb_line(line)
            if line_dict["atom_name"].startswith("H"):
                line_dict["temp_fac"] = "  0.00"
            else:
                if line_dict["atom_name"].startswith("C") and not line_dict["atom_name"].startswith("CA"):
                    line_dict["temp_fac"] = "  0.50"
                else:
                    line_dict["temp_fac"] = "  1.00"
            line_ = line_operations.create_line(line_dict=line_dict)
            lines_.append(line_)
            if line.startswith("TER"):
                line_ = line
                lines_.append("line_")
        file_operations.write_file(file=restraints_file, lines=lines_)
        lines=lines_
        return

    def renumber(pdb_file, pdb_file_output):
        """
        operations.renumber() takes a pdb_file and a pdb_file_output name; then it checks if there are more than 99.999 atoms. If
        so it will raise a ValueError, if not it will use the line_operations.fill_serial() method to renumber the atoms starting at 1. In
        the end it will write a file based on pdb_file_output.
        """
        lines=file_operations.read_file(pdb_file=pdb_file)
        if len(lines) > 99999:
            raise ValueError("len(lines)>99999. Try again with less atoms.")
        lines_ = []
        serial_no=1
        for line in lines:
            line_dict=line_operations.read_pdb_line(line=line)
            line_dict=line_operations.fill_serial(serial_no=serial_no, line_dict=line_dict)
            line_=line_operations.create_line(line_dict=line_dict)
            lines_.append(line_)
            serial_no+=1
        file_operations.write_file(file=pdb_file_output, lines=lines_)
        return

    def renumber_tip3(pdb_file, pdb_file_output, segment):
        """
        operations.renumber_tip3() takes a pdb_file, a pdb_file_output and a segment name; then it will read the file and
        check if there are more than 99.999 atoms in the current file. If so it will raise a ValueError, if not it will
        not only renumber the serial numbers of all atoms but also the residue numbers. This only works with a water model
        that has exactly three atoms per water molecule. In the end it writes a file based on the pdb_file_output specifications.
        """
        lines=file_operations.read_file(pdb_file=pdb_file)
        if len(lines) > 99999:
            raise ValueError("len(lines)>99999. Try again with less atoms.")
        lines_ = []
        serial_no=1
        for line in lines:
            line_dict=line_operations.read_pdb_line(line=line)
            line_dict=line_operations.fill_serial(serial_no=serial_no, line_dict=line_dict)
            resi_no=((serial_no-1)//3)+1
            line_operations.fill_resi_no(resi_no=resi_no, line_dict=line_dict)
            line_dict["segment"] = segment
            line_=line_operations.create_line(line_dict=line_dict)
            lines_.append(line_)
            serial_no+=1
        file_operations.write_file(file=pdb_file_output, lines=lines_)
        return