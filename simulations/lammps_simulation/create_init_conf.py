class bead:
    ''' class bead contains all info about bead in chain: global number, 
    the remaining valence, type of the bead and coordinates '''
    def __init__(self):
        self.numbd=0
        self.valence=0
        self.typep=0
        self.x=0.0
        self.y=0.0
        self.z=0.0

    def __lt__(self, other):
        return self.numbd < other.numbd
    
class chain:
    '''class chain has two lists of beads and bonds and general info about system such as total number of particles, density and box size along each axis'''
    def __init__(self):
        self.bd=[]
        self.bnd=[]
        self.number_of_beads=float()
        self.number_of_bonds=float()
        self.density=float()
        self.xbox=float()
        self.ybox=float()
        self.zbox=float()

def lst_bd3(num):
    lst=[]
    for i in range(num):
        a = bead()
        a.typep = 3
        lst.append(a)
    return lst

def count_c(part, s, f):
    cb = 0
    for i in range(s,f):
        if part.bd[i].typep == 3:
            cb += 1
    return cb

def create_system_drosophila_with_elongated_centromer_telomer(system, path_of_marking_ab, path_of_cent_telo):
    import numpy as np
    '''parse the file with A B compartments and another file with centromeric/telomeric regions.
    A - 1 - active
    B - 2 - inactive
    C - 3 - centromeric-telomeric regions
    6 chromosomes
    2 3 2 3 X X
    '''
    c_frac = 0.3
    resolution = 10000
    f=open(path_of_marking_ab)
    system.extend([chain(), chain(), chain(), chain(), chain(), chain()])
    for line in f:
        to_add = bead()
        name, start, end, val = line.split()
        if name=='chr2L' or name=='chr2R':
            if float(val) > 0:
                to_add.typep = 1
                system[0].bd.append(to_add)
            else:
                to_add.typep = 2
                system[0].bd.append(to_add)
        elif name=='chr3L' or name=='chr3R':
            if float(val) > 0:
                to_add.typep = 1
                system[1].bd.append(to_add)
            else:
                to_add.typep = 2
                system[1].bd.append(to_add)
        elif name=='chrX':
            if float(val) > 0:
                to_add.typep = 1
                system[2].bd.append(to_add)
            else:
                to_add.typep = 2
                system[2].bd.append(to_add)
    f.close()

    system[4] = system[2]
    system[2] = system[0]
    system[3] = system[1]
    system[5] = system[4] #new Xchr

    f = open(path_of_cent_telo)
    cent_tel = []
    for line in f:
        name, start, end, val = line.split()
        if name == 'chr2L':
            cent_b = int(int(start) / resolution)
            cent_e = int(np.ceil(int(end) / resolution))
        if name == 'chr2R':
            cent_e += int(np.ceil(int(end) / resolution))
            for i in range(cent_b, cent_e):
                system[0].bd[i].typep = 3
                system[2].bd[i].typep = 3
        if name == 'chr3L':
            cent_b = int(int(start) / resolution)
            cent_e = int(np.ceil(int(end) / resolution))
        if name == 'chr3R':
            cent_e += int(np.ceil(int(end) / resolution))
            for i in range(cent_b, cent_e):
                system[1].bd[i].typep = 3
                system[3].bd[i].typep = 3
    for i in range(6):
        system[i].bd[0].typep = 3
        system[i].bd[-1].typep = 3
    f.close()

    #lengths chr2L = 2351, chr2R = 2528, chr3L = 2811, chr3R = 3207, chrX = 2354
    n_l = 2351
    nc = count_c(system[0], 0, n_l)
    ins_tel_l = int((c_frac*n_l-nc)/(1-c_frac)*0.1)
    ins_cent_l = int((c_frac*n_l-nc)/(1-c_frac)*0.9)
    system[0].bd[0:0] = lst_bd3(ins_tel_l)
    system[0].bd[ins_tel_l+n_l-1:ins_tel_l+n_l-1] = lst_bd3(ins_cent_l)
    
    chr2_mid = ins_tel_l+ins_cent_l+n_l
        
    n_r = 2528
    nc = count_c(system[0], ins_tel_l+ins_cent_l+n_l, len(system[0].bd))
    ins_tel_r = int((c_frac*n_r-nc)/(1-c_frac)*0.9)
    ins_cent_r = int((c_frac*n_r-nc)/(1-c_frac)*0.1)
    system[0].bd[ins_tel_l+ins_cent_l+n_l:ins_tel_l+ins_cent_l+n_l] = lst_bd3(ins_tel_r)
    system[0].bd.extend(lst_bd3(ins_cent_r))
    
    n_l = 2811
    nc = count_c(system[1], 0, n_l)
    ins_tel_l = int((c_frac*n_l-nc)/(1-c_frac)*0.1)
    ins_cent_l = int((c_frac*n_l-nc)/(1-c_frac)*0.9)
    system[1].bd[0:0] = lst_bd3(ins_tel_l)
    system[1].bd[ins_tel_l+n_l-1:ins_tel_l+n_l-1] = lst_bd3(ins_cent_l)
    
    chr3_mid = ins_tel_l+ins_cent_l+n_l
    
    n_r = 3207
    nc = count_c(system[1], ins_tel_l+ins_cent_l+n_l, len(system[1].bd))
    ins_tel_r = int((c_frac*n_r-nc)/(1-c_frac)*0.9)
    ins_cent_r = int((c_frac*n_r-nc)/(1-c_frac)*0.1)
    system[1].bd[ins_tel_l+ins_cent_l+n_l:ins_tel_l+ins_cent_l+n_l] = lst_bd3(ins_tel_r)
    system[1].bd.extend(lst_bd3(ins_cent_r))
    
    n_l = 2354
    nc = count_c(system[4], 0, len(system[4].bd))
    ins_tel_l = int((c_frac*n_l-nc)/(1-c_frac)*0.1)
    ins_cent_l = int((c_frac*n_l-nc)/(1-c_frac)*0.9)
    system[4].bd[0:0] = lst_bd3(ins_tel_l)
    system[4].bd.extend(lst_bd3(ins_cent_l))
    
    cent_tel.append([chr2_mid])
    cent_tel.append([chr3_mid])
    cent_tel.append([chr2_mid])
    cent_tel.append([chr3_mid])

    system[2] = system[0]
    system[3] = system[1]
    system[5] = system[4] #new Xchr
                                
    f.close()
    return cent_tel

def rw_in_cyllinder(beads_number=1000, cell_size_sf=10, xchr=False, mid=0):
    import random
    import numpy as np
    r = 1 # chain step
    d = 10 # radius of cyllinder
    if xchr:
        step = (cell_size_sf - 4) / beads_number
    else:
        step1 = (cell_size_sf - 4) / mid
        step2 = (cell_size_sf - 4) / (beads_number - mid)
    x1 = np.zeros(beads_number)
    y1 = np.zeros(beads_number)
    z1 = np.zeros(beads_number)
    z1[0] = abs(cell_size_sf - 4)
    if xchr:
        for i in range(1, beads_number):
            phi = random.uniform(0, 2 * np.pi)
            x1[i] = x1[i - 1] + r * np.cos(phi)
            y1[i] = y1[i - 1] + r * np.sin(phi)
            while (x1[i])**2 + (y1[i])**2 > d:
                phi = random.uniform(0, 2 * np.pi)
                x1[i] = x1[i - 1] + r * np.cos(phi)
                y1[i] = y1[i - 1] + r * np.sin(phi)
            z1[i] = abs(cell_size_sf - 4 - step * i)
    else:
        for i in range(1, mid):
            phi = random.uniform(0, 2 * np.pi)
            x1[i] = x1[i - 1] + r * np.cos(phi)
            y1[i] = y1[i - 1] + r * np.sin(phi)
            while (x1[i])**2 + (y1[i])**2 > d:
                phi = random.uniform(0, 2 * np.pi)
                x1[i] = x1[i - 1] + r * np.cos(phi)
                y1[i] = y1[i - 1] + r * np.sin(phi)
            z1[i] = abs(cell_size_sf - 4 - step1 * i)
        for i in range(mid, beads_number):
            phi = random.uniform(0, 2 * np.pi)
            x1[i] = x1[i - 1] + r * np.cos(phi)
            y1[i] = y1[i - 1] + r * np.sin(phi)
            while (x1[i])**2 + (y1[i])**2 > d:
                phi = random.uniform(0, 2 * np.pi)
                x1[i] = x1[i - 1] + r * np.cos(phi)
                y1[i] = y1[i - 1] + r * np.sin(phi)
            z1[i] = z1[i-1] - step2

    z1 = abs(z1)
    z1=z1+2

    return x1, y1, z1

def write_lmpdat_drosophila_elongated_pericentromeris_regions(system, path, cent_tel):
    import numpy as np
    '''I assume that volume concentration of chromatin is [concentration.
    cell_size has +2 to expand a bit boarders, the simulation will be in a sphere'''
    concentration = 0.2
    h_cell_size = 4 * np.power((len(system[0].bd)+len(system[1].bd)+len(system[2].bd)+len(system[3].bd)+len(system[4].bd)+len(system[5].bd))/concentration/np.pi/4, 1/3)
    print ("cell size is %f"%h_cell_size)
    r_cell_size = h_cell_size/2
    shift_x=[r_cell_size/2-7, r_cell_size/2-7, r_cell_size/2, r_cell_size/2+7, r_cell_size/2+7, r_cell_size/2-3]
    shift_y=[r_cell_size/2-7, r_cell_size/2+7, r_cell_size/2, r_cell_size/2-7, r_cell_size/2+7, r_cell_size/2-3]
    shift_z=[0, 0, 0, 0, 0, 0, 0, 0, 0,0]
    num_beads=len(system[0].bd)+len(system[1].bd)+len(system[2].bd)+len(system[3].bd)+len(system[4].bd)+len(system[5].bd)
    #fix telomeric regions
    for i in range(6):
        if i < 4:
            system[i].bd[0].x = shift_x[i] - 1
            system[i].bd[0].y = shift_y[i] - 1
            system[i].bd[0].z = h_cell_size - shift_z[i]
            system[i].bd[0].typep = 4
            system[i].bd[-1].x = shift_x[i] + 1
            system[i].bd[-1].y = shift_y[i] + 1
            system[i].bd[-1].z = h_cell_size - shift_z[i]
            system[i].bd[-1].typep = 4
        else:
            system[i].bd[0].x = shift_x[i] - 1
            system[i].bd[0].y = shift_y[i] - 1
            system[i].bd[0].z = h_cell_size - shift_z[i]
            system[i].bd[0].typep = 4
            system[i].bd[-1].x = shift_x[i] + 1
            system[i].bd[-1].y = shift_y[i] + 1
            system[i].bd[-1].z = shift_z[i]
            system[i].bd[-1].typep = 4
    #write lmp file        
    with open(path,'w') as thefile:
        thefile.write("#drosophila melanogaster 10kb\n\n")
        thefile.write("%8d atoms\n"%num_beads)
        thefile.write("%8d bonds\n\n\n"%(num_beads-6))
        thefile.write("  4 atom types\n")
        thefile.write("  1 bond types\n\n")
        thefile.write("0.000 %5.3f xlo xhi\n"%r_cell_size)
        thefile.write("0.000 %5.3f ylo yhi\n"%r_cell_size)
        thefile.write("0.000 %5.3f zlo zhi\n\n"%h_cell_size)
        thefile.write("Masses\n\n1	1.00\n2	1.00\n3	1.00\n4	1.00\n\n Atoms\n\n")
    ttl_nm=0
    with open(path,'a') as f:
        for chr in range(6):
            if chr < 4:
                x,y,z = rw_in_cyllinder(beads_number=len(system[chr].bd), cell_size_sf=h_cell_size, xchr=False, mid=cent_tel[chr][0])
                z[cent_tel[chr][0]] = 2
                system[chr].bd[cent_tel[chr][0]].typep = 4
            else:
                x,y,z = rw_in_cyllinder(beads_number=len(system[chr].bd), cell_size_sf=h_cell_size, xchr=True, mid = 0)
            for i in range(len(x)):
                x[i] += shift_x[chr] + i * 10 / float(len(x)) - 5
                y[i] += shift_y[chr] + i * 10 / float(len(x)) - 5
                z[i] += shift_z[chr]
            for i in range(0, len(system[chr].bd)):
                f.write("%8d 0 %3d %5.4f %5.4f %5.4f\n"%(ttl_nm+i+1, system[chr].bd[i].typep, x[i], y[i], z[i]))
            ttl_nm=ttl_nm+len(system[chr].bd)
    with open(path,'a') as thefile:
        thefile.write("\n Bonds\n\n")
    f=open(path,'a')
    bnd_nm=1
    bnd=1
    for chr in range(6):
        for i in range(1,len(system[chr].bd)):
            f.write("%8d 1 %5d %5d\n"%(bnd_nm,bnd,bnd+1))
            bnd_nm+=1
            bnd+=1
        bnd+=1
    f.close()
    print ('Done!')

def main():
    system=[]
    cent_tel=create_system_drosophila_with_elongated_centromer_telomer(system, 
                  './comp_HiTC_sub_10000.bed', 
                  './dm6_centromeres_telomeres')
    write_lmpdat_drosophila_elongated_pericentromeris_regions(system, './rwincyl', cent_tel)
main()
