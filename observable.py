import math
import matplotlib.pyplot as plt
import numpy as np
import sys
import spglib as spg

def read_poscar(filename):
    cell = []
    atomic_positions = []
    direct_var = "Direct"
    poscar_file = open(str(filename),"r")
    poscar_file.readline()
    poscar_file.readline()

    line = poscar_file.readline()
    cell.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    line = poscar_file.readline()
    cell.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    line = poscar_file.readline()
    cell.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    atom_types = poscar_file.readline().split()
    atoms_type_numbers_degen =  poscar_file.readline().split()
    total_n_atoms = 0
    for i in range(0,len(atoms_type_numbers_degen)):
        atoms_type_numbers_degen[i] = int(atoms_type_numbers_degen[i])
        total_n_atoms = total_n_atoms + atoms_type_numbers_degen[i]
    direct_var = poscar_file.readline()
    for i in range(0,total_n_atoms):
        line = poscar_file.readline()
        atomic_positions.append([ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ])

    poscar_file.close()

    return cell, atomic_positions, atoms_type_numbers_degen, atom_types, direct_var

def read_ICOHP(filename,TRESHOLD,DIST_TRESHOLD):
    ICOHP_values_list = []
    NEW_ICOHP_FLAG = False

    ICOHP_file = open(str(filename),"r")
    lines = ICOHP_file.readlines()

    for i in range(1,len(lines)):
        temp_values = [float(lines[i].split()[3]),  float(lines[i].split()[7]),lines[i].split()[1],lines[i].split()[2]]
        if abs(temp_values[1]) > TRESHOLD: 
            if DIST_TRESHOLD >= temp_values[0]:
                ICOHP_values_list.append([float(lines[i].split()[3]),  float(lines[i].split()[7]),lines[i].split()[1],lines[i].split()[2]])

    return ICOHP_values_list






if len(sys.argv) <4:
    print("Ths code requires the following inputs in order: [POSCAR] [iCOBI_lobster_file] [iCOBI_cutoff (0.018 recommended)] [Interatomic Distance Cutoff]")
    exit()


input_poscar  = sys.argv[1]
input_ICOHP   = sys.argv[2]
ICOHP_CUT     = float(sys.argv[3])
DISTANCE_CUT  = float(sys.argv[4])


cell = []
atoms = []
pos = []
numbers = []
types =[]



cell, atoms, tmp_num, types, input_type = read_poscar(input_poscar)
variab =                      read_ICOHP(input_ICOHP,ICOHP_CUT,DISTANCE_CUT)






for i in range(0,len(atoms)):
    pos.append([0,0,0])
    if input_type == "Direct\n":
        pos[i][0] = atoms[i][0]*cell[0][0] + atoms[i][1]*cell[1][0] +atoms[i][2]*cell[2][0] 
        pos[i][1] = atoms[i][0]*cell[0][1] + atoms[i][1]*cell[1][1] +atoms[i][2]*cell[2][1]
        pos[i][2] = atoms[i][0]*cell[0][2] + atoms[i][1]*cell[1][2] +atoms[i][2]*cell[2][2]
    else:
        pos[i][0] = atoms[i][0]
        pos[i][1] = atoms[i][1]
        pos[i][2] = atoms[i][2]

count = 0
for i in range(0,len(tmp_num)):
    count = count + 1
    for j in range(0,tmp_num[i]):
        numbers.append(count)
    



cell_SPG_LIB = (cell,atoms,numbers)
equivalent_atoms = spg.get_symmetry(cell_SPG_LIB, symprec=5e-02, angle_tolerance=-1.0, mag_symprec=-1.0, is_magnetic=True)['equivalent_atoms']



data = []
stored_int =  []
count = 0

print("A1  Id1 A2  Id2 Dist    iCOBI   D_iCOBI  Repeated Cell")
for i in range(0,len(pos)):
    count = 0
    data.append([0,0,0])#,0,0,0])
    for j in range(0,len(pos)):
        for a in range(-3,4):
            for b in range(-3,4):
                for c in range(-3,4):
                    distance = [ pos[j][0]+cell[0][0]*a+cell[1][0]*b+cell[2][0]*c-pos[i][0],
                           pos[j][1]+cell[0][1]*a+cell[1][1]*b+cell[2][1]*c-pos[i][1],
                           pos[j][2]+cell[0][2]*a+cell[1][2]*b+cell[2][2]*c-pos[i][2]]
                    ddd = math.sqrt(distance[0]*distance[0]+ distance[1]*distance[1] + distance[2]*distance[2])
                    for q in range(0,len(variab)):
                        if ddd > variab[q][0] - 0.0001 and ddd < variab[q][0] + 0.0001:
                            if variab[q][2] == str(types[numbers[i]-1]) + str(i+1) or variab[q][3] == str(types[numbers[i]-1]) + str(i+1):
                                if variab[q][2] == str(types[numbers[j]-1]) + str(j+1) or variab[q][3] == str(types[numbers[j]-1]) + str(j+1):
                                    count = count + 1 
                                    print( f"{types[numbers[i]-1]:<3} " f"{i+1:<3} " f"{types[numbers[j]-1]:<3} " f"{j+1:<3} "     f"{ddd:<2.5f} " f"{variab[q][0]:<2.5f} " f"{variab[q][1]:<2.6f} " f"{a:<1} " f"{b:<1} " f"{c:<1}"
)
                                    data[i][0] = data[i][0] - distance[0]/ddd*variab[q][1] 
                                    data[i][1] = data[i][1] - distance[1]/ddd*variab[q][1] 
                                    data[i][2] = data[i][2] - distance[2]/ddd*variab[q][1]

    stored_int.append([i,count])
    if count !=0 :
        data[i][0] = data[i][0]/count
        data[i][1] = data[i][1]/count
        data[i][2] = data[i][2]/count

    else:
        data[i][0] = 0
        data[i][1] = 0
        data[i][2] = 0

print("")
print("---------Interactions---------")
print("")
for i in stored_int: print('Atom: ' + f"{types[numbers[i[0]]-1]:<3}" +  f"({i[0]:<3}) " + '  Interactions: ' f"({i[1]:<3})" +   '  Equivalent to ' + str(equivalent_atoms[i[0]]) )
print("")
print("---------Symmetry Indexes---------")
print("")


count = 0
S_VAR = 0


for i in range(0,len(equivalent_atoms)):
    tmp_atm = equivalent_atoms[i]
    for j in range(0,len(equivalent_atoms)):

        if tmp_atm != -1 and  tmp_atm == equivalent_atoms[j]:
            S_VAR = S_VAR + math.sqrt(data[j][0]*data[j][0]+ data[j][1]*data[j][1] + data[j][2]*data[j][2])
            equivalent_atoms[j] = -1
            count = count +1
    if tmp_atm != -1:    print("Atom [Multiplicity]: " + f"{types[numbers[i]-1]:<2}" + " [" + str(count) + "] " + " - SI: " + f"{S_VAR/count:<1.5f}")
    count = 0
    S_VAR = 0
