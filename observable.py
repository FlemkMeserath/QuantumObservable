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
        if abs(temp_values[1]) > TRESHOLD: # and len(ICOHP_values_list) == 0
            if DIST_TRESHOLD >= temp_values[0]:
                ICOHP_values_list.append([float(lines[i].split()[3]),  float(lines[i].split()[7]),lines[i].split()[1],lines[i].split()[2]])
#        NEW_ICOHP_FLAG = True
#        for j in range(0,len(ICOHP_values_list)):
#            if temp_values[0] == ICOHP_values_list[j][0] or abs(temp_values[1]) < TRESHOLD:
#                NEW_ICOHP_FLAG = False
#        if NEW_ICOHP_FLAG == True and len(ICOHP_values_list) > 0: 
#            if DIST_TRESHOLD >= temp_values[0]:
#                ICOHP_values_list.append([float(lines[i].split()[3]), - float(lines[i].split()[7])])
    return ICOHP_values_list



cell = []
atoms = []
pos = []
numbers = []
types =[]


input_poscar  = sys.argv[1]
input_ICOHP   = sys.argv[2]
ICOHP_CUT     = float(sys.argv[3])
DISTANCE_CUT  = float(sys.argv[4])

cell, atoms, tmp_num, types, input_type = read_poscar(input_poscar)
variab =                      read_ICOHP(input_ICOHP,ICOHP_CUT,DISTANCE_CUT)

print("variab",variab)

print(variab)


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
    
print(numbers,types)


cell_SPG_LIB = (cell,atoms,numbers)
equivalent_atoms = spg.get_symmetry(cell_SPG_LIB, symprec=5e-03, angle_tolerance=-1.0, mag_symprec=-1.0, is_magnetic=True)['equivalent_atoms']
print(equivalent_atoms)

#print(pos,atoms,cell)
data = []
count = 0
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
                                    count = count + 1 #*variab[q][1]
                                    print(types[numbers[i]-1],i+1,types[numbers[j]-1],j+1,variab[q][0],variab[q][1],ddd,pos[j][0],pos[j][1],pos[j][2],a,b,c)
                                    data[i][0] = data[i][0] - distance[0]/ddd*variab[q][1] 
                                    data[i][1] = data[i][1] - distance[1]/ddd*variab[q][1] 
                                    data[i][2] = data[i][2] - distance[2]/ddd*variab[q][1]
#                                    data[i][3] = data[i][3] - distance[0]
#                                    data[i][4] = data[i][4] - distance[1]
#                                    data[i][5] = data[i][5] - distance[2]
    print(i,count)
    if count !=0 :
        data[i][0] = data[i][0]/count
        data[i][1] = data[i][1]/count
        data[i][2] = data[i][2]/count
#        data[i][3] = data[i][3]/count
#        data[i][4] = data[i][4]/count
#        data[i][5] = data[i][5]/count
    else:
        data[i][0] = 0
        data[i][1] = 0
        data[i][2] = 0
#        data[i][3] = 0
#        data[i][4] = 0
#        data[i][5] = 0

for i in data: print(i)

count = 0
eq_count = equivalent_atoms[0]
S_VAR = 0
#S_VAR2 = 0
for i in range(0,len(equivalent_atoms)):
    if eq_count == equivalent_atoms[i]:
        S_VAR = S_VAR + math.sqrt(data[i][0]*data[i][0]+ data[i][1]*data[i][1] + data[i][2]*data[i][2])
#        S_VAR2 = S_VAR2 + math.sqrt(data[i][3]*data[i][3]+ data[i][4]*data[i][4] + data[i][5]*data[i][5])
        eq_count = equivalent_atoms[i]
        count = count +1
    else:
        print(S_VAR/count,count)
#        print(S_VAR2/count,count,"d")
        count = 1
        eq_count = equivalent_atoms[i]
        S_VAR = math.sqrt(data[i][0]*data[i][0]+ data[i][1]*data[i][1] + data[i][2]*data[i][2])
#        S_VAR2 = math.sqrt(data[i][3]*data[i][3]+ data[i][4]*data[i][4] + data[i][5]*data[i][5])
print(S_VAR/count,count)
#print(S_VAR2/count,count,"d")
