# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #


# TODO: test with K2
def create_dat_file(line, path, names=None):

    # helper function
    def _get_spaces(str1,str2, total_spaces):
        length1 = len(str(str1))
        length2 = len(str(str2))
        middle_spaces = total_spaces - length2 - length1
        return middle_spaces
 
    with open(f'{path}.dat', 'w') as file:
        onesided = []
        crystal  = []
        for name in names:
            if line[name].side != 'both':
                onesided.append(name)
                continue
            # if line[name].__class__ == EverestCrystal:
            #     crystal.append(name)
            #     continue
            gap = line[name].gap
            if gap == None:
                gap = "null"
            mat = line[name].material
            length = line[name].length
            angle = line[name].angle
            offset = ( (line[name]._jaw_LU + line[name]._jaw_LD)/2 + 
                       (line[name]._jaw_RU + line[name]._jaw_RD)/2 )/2
            file.write(f"{name}" + " " * _get_spaces(name,gap, 65) + f"{gap}" + f" {mat}" \
                       + " " * _get_spaces(mat, "1.111", 15) + f"{length}" \
                       + " " * _get_spaces(length, angle, 14) + f"{angle}" \
                       + " " * 8 + f"{offset} \n")
        if len(onesided) > 0:
            file.write("SETTINGS \n")
            for name in onesided:
                file.write(f"ONESIDED {name}" + " " * _get_spaces(name,"1", 50) \
                           + f"{line[name]._side} \n")
        if len(crystal) > 0:
            if len(onesided) == 0:
                file.write("SETTINGS \n")
            for name in crystal:
                file.write(f"CRYSTAL {name}" + f" {line[name].bending_radius}"\
                           + f" {line[name].width}" + f" {line[name].height}" \
                           + f" {line[name].miscut}" + f" {line[name]._orient}")
    print("File created.")
    file.close()
    return file  