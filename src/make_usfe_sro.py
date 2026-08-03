import numpy as np
import numpy.random as rand
import time
from scipy.interpolate import griddata
import numpy.linalg as la
import sys

args=sys.argv
material=args[1] #monbti_300, monbti_1673, tanbti_300, tanbti_1673
N1=int(args[2])
N2=int(args[3])
N3=int(args[4])
NN=int(args[5]) #number of nearest neighbor shells
num_surf=int(args[6]) #id number of surface

print("N1 = %d, N2 = %d, N3 = %d, NN Shells = %d"%(N1*32,N2,N3,NN))

#DEFINE SRO PARAMETERS
#Mo or Ta = 0
#Nb = 1
#Ti = 2
prob_monbti_300 = np.zeros((3,3))
prob_monbti_1673 = np.zeros((3,3))
prob_tanbti_300 = np.zeros((3,3))
prob_tanbti_1673 = np.zeros((3,3))

prob_rand = np.zeros((3,3))
prob_rand[0,0] = 0.333333333333 #Mo-Mo
prob_rand[1,1] = 0.333333333333 #Nb-Nb
prob_rand[2,2] = 0.333333333333 #Ti-Ti
prob_rand[0,1] = prob_rand[1,0] = 0.333333333333 #Mo-Nb
prob_rand[0,2] = prob_rand[2,0] = 0.333333333333 #Mo-Ti
prob_rand[1,2] = prob_rand[2,1] = 0.333333333333 #Nb-Ti

prob_monbti_300[0,0] = 0.10768229 #Mo-Mo
prob_monbti_300[1,1] = 0.29674479 #Nb-Nb
prob_monbti_300[2,2] = 0.25247396 #Ti-Ti
prob_monbti_300[0,1] = prob_monbti_300[1,0] = 0.42402344 #Mo-Nb
prob_monbti_300[0,2] = prob_monbti_300[2,0] = 0.46829427 #Mo-Ti
prob_monbti_300[1,2] = prob_monbti_300[2,1] = 0.27923177 #Nb-Ti

prob_monbti_1673[0,0] = 0.28072917 #Mo-Mo
prob_monbti_1673[1,1] = 0.32552083 #Nb-Nb
prob_monbti_1673[2,2] = 0.30585937 #Ti-Ti
prob_monbti_1673[0,1] = prob_monbti_1673[1,0] = 0.34980469 #Mo-Nb
prob_monbti_1673[0,2] = prob_monbti_1673[2,0] = 0.36946615 #Mo-Ti
prob_monbti_1673[1,2] = prob_monbti_1673[2,1] = 0.32467448 #Nb-Ti

prob_tanbti_300[0,0] = 0.40598958 #Ta-Ta
prob_tanbti_300[1,1] = 0.30846354 #Nb-Nb
prob_tanbti_300[2,2] = 0.38046875 #Ti-Ti
prob_tanbti_300[0,1] = prob_tanbti_300[1,0] = 0.33300781 #Ta-Nb
prob_tanbti_300[0,2] = prob_tanbti_300[2,0] = 0.2610026 #Ta-Ti
prob_tanbti_300[1,2] = prob_tanbti_300[2,1] = 0.35852865 #Nb-Ti

prob_tanbti_1673[0,0] = 0.34973958 #Ta-Ta
prob_tanbti_1673[1,1] = 0.32539062 #Nb-Nb
prob_tanbti_1673[2,2] = 0.3265625 #Ti-Ti
prob_tanbti_1673[0,1] = prob_tanbti_1673[1,0] = 0.32571615 #Ta-Nb
prob_tanbti_1673[0,2] = prob_tanbti_1673[2,0] = 0.32454427 #Ta-Ti
prob_tanbti_1673[1,2] = prob_tanbti_1673[2,1] = 0.34889323 #Nb-Ti

#READ IN USFE VALUES
filename='local_composition_vs_SFE.csv'
f=open(filename,'r')
usfe_data = np.zeros((580,6)) #Mo Nb Ta Ti Temp USFE
l = f.readline()
for i in range(580):
    l = f.readline().split(',')
    if l[2] != "":
        usfe_data[i,2] = float(l[2]) #Ta
    usfe_data[i,3] = float(l[3]) #Ti
    usfe_data[i,1] = float(l[4]) #Nb
    if l[5] != "":
        usfe_data[i,0] = float(l[5]) #Mo
    usfe_data[i,4] = float(l[6][:-2])
    usfe_data[i,5] = float(l[1])
usfe_monbti_300K = usfe_data[(usfe_data[:,0] != 0) & (usfe_data[:,-2] == 300)]
usfe_tanbti_300K = usfe_data[(usfe_data[:,2] != 0) & (usfe_data[:,-2] == 300)]
usfe_monbti_1673K = usfe_data[(usfe_data[:,0] != 0) & (usfe_data[:,-2] == 1673)]
usfe_tanbti_1673K = usfe_data[(usfe_data[:,2] != 0) & (usfe_data[:,-2] == 1673)]
usfe_monbti_rand = usfe_data[(usfe_data[:,0] != 0) & (usfe_data[:,-2] == -1)]
usfe_tanbti_rand = usfe_data[(usfe_data[:,2] != 0) & (usfe_data[:,-2] == -1)]

usfe_mo = 1232.10
usfe_nb = 777.74
usfe_ta = 763.08
usfe_ti = 439.22


#interpolate USFE values
grid_x, grid_y = np.mgrid[0:1:2001j, 0:np.sqrt(3)/2:1001j]

#MoNbTi 1673
xy_monbti_1673 = np.zeros((len(usfe_monbti_1673K)+3,2))
z_monbti_1673 = np.zeros(len(usfe_monbti_1673K)+3)
xy_monbti_1673[:-3,0] = 0.5*(2*usfe_monbti_1673K[:,1]+usfe_monbti_1673K[:,3])
xy_monbti_1673[:-3,1] = np.sqrt(3)/2*usfe_monbti_1673K[:,3]
z_monbti_1673[:-3] = usfe_monbti_1673K[:,-1]
xy_monbti_1673[-3] = [0,0] #Mo
z_monbti_1673[-3] = usfe_mo
xy_monbti_1673[-2] = [1,0] #Nb
z_monbti_1673[-2] = usfe_nb
xy_monbti_1673[-1] = [0.5,np.sqrt(3)/2]#Ti
z_monbti_1673[-1] = usfe_ti
monbti_1673 = griddata(xy_monbti_1673, z_monbti_1673, (grid_x, grid_y), method='linear')

#MoNbTi 300
xy_monbti_300 = np.zeros((len(usfe_monbti_300K)+3,2))
z_monbti_300 = np.zeros(len(usfe_monbti_300K)+3)
xy_monbti_300[:-3,0] = 0.5*(2*usfe_monbti_300K[:,1]+usfe_monbti_300K[:,3])
xy_monbti_300[:-3,1] = np.sqrt(3)/2*usfe_monbti_300K[:,3]
z_monbti_300[:-3] = usfe_monbti_300K[:,-1]
xy_monbti_300[-3] = [0,0] #Mo
z_monbti_300[-3] = usfe_mo
xy_monbti_300[-2] = [1,0] #Nb
z_monbti_300[-2] = usfe_nb
xy_monbti_300[-1] = [0.5,np.sqrt(3)/2]#Ti
z_monbti_300[-1] = usfe_ti
monbti_300 = griddata(xy_monbti_300, z_monbti_300, (grid_x, grid_y), method='linear')

#TaNbTi 1673
xy_tanbti_1673 = np.zeros((len(usfe_tanbti_1673K)+3,2))
z_tanbti_1673 = np.zeros(len(usfe_tanbti_1673K)+3)
xy_tanbti_1673[:-3,0] = 0.5*(2*usfe_tanbti_1673K[:,1]+usfe_tanbti_1673K[:,3])
xy_tanbti_1673[:-3,1] = np.sqrt(3)/2*usfe_tanbti_1673K[:,3]
z_tanbti_1673[:-3] = usfe_tanbti_1673K[:,-1]
xy_tanbti_1673[-3] = [0,0] #Ta
z_tanbti_1673[-3] = usfe_ta
xy_tanbti_1673[-2] = [1,0] #Nb
z_tanbti_1673[-2] = usfe_nb
xy_tanbti_1673[-1] = [0.5,np.sqrt(3)/2]#Ti
z_tanbti_1673[-1] = usfe_ti
tanbti_1673 = griddata(xy_tanbti_1673, z_tanbti_1673, (grid_x, grid_y), method='linear')

#TaNbTi 300
xy_tanbti_300 = np.zeros((len(usfe_tanbti_300K)+3,2))
z_tanbti_300 = np.zeros(len(usfe_tanbti_300K)+3)
xy_tanbti_300[:-3,0] = 0.5*(2*usfe_tanbti_300K[:,1]+usfe_tanbti_300K[:,3])
xy_tanbti_300[:-3,1] = np.sqrt(3)/2*usfe_tanbti_300K[:,3]
z_tanbti_300[:-3] = usfe_tanbti_300K[:,-1]
xy_tanbti_300[-3] = [0,0] #Ta
z_tanbti_300[-3] = 763.08
xy_tanbti_300[-2] = [1,0] #Nb
z_tanbti_300[-2] = 777.74
xy_tanbti_300[-1] = [0.5,np.sqrt(3)/2]#Ti
z_tanbti_300[-1] = 439.22
tanbti_300 = griddata(xy_tanbti_300, z_tanbti_300, (grid_x, grid_y), method='linear')

#MoNbTi random
xy_monbti_rand = np.zeros((len(usfe_monbti_rand)+3,2))
z_monbti_rand = np.zeros(len(usfe_monbti_rand)+3)
xy_monbti_rand[:-3,0] = 0.5*(2*usfe_monbti_rand[:,1]+usfe_monbti_rand[:,3])
xy_monbti_rand[:-3,1] = np.sqrt(3)/2*usfe_monbti_rand[:,3]
z_monbti_rand[:-3] = usfe_monbti_rand[:,-1]
xy_monbti_rand[-3] = [0,0] #Mo
z_monbti_rand[-3] = usfe_mo
xy_monbti_rand[-2] = [1,0] #Nb
z_monbti_rand[-2] = usfe_nb
xy_monbti_rand[-1] = [0.5,np.sqrt(3)/2]#Ti
z_monbti_rand[-1] = usfe_ti
monbti_rand = griddata(xy_monbti_rand, z_monbti_rand, (grid_x, grid_y), method='linear')

#TaNbTi rand
xy_tanbti_rand = np.zeros((len(usfe_tanbti_rand)+3,2))
z_tanbti_rand = np.zeros(len(usfe_tanbti_rand)+3)
xy_tanbti_rand[:-3,0] = 0.5*(2*usfe_tanbti_rand[:,1]+usfe_tanbti_rand[:,3])
xy_tanbti_rand[:-3,1] = np.sqrt(3)/2*usfe_tanbti_rand[:,3]
z_tanbti_rand[:-3] = usfe_tanbti_rand[:,-1]
xy_tanbti_rand[-3] = [0,0] #Ta
z_tanbti_rand[-3] = 763.08
xy_tanbti_rand[-2] = [1,0] #Nb
z_tanbti_rand[-2] = 777.74
xy_tanbti_rand[-1] = [0.5,np.sqrt(3)/2]#Ti
z_tanbti_rand[-1] = 439.22
tanbti_rand = griddata(xy_tanbti_rand, z_tanbti_rand, (grid_x, grid_y), method='linear')

if material == "monbti_300":
    prob=prob_monbti_300
    usfe_interp=monbti_300
    xy_test=xy_monbti_300
elif material == "monbti_1673":
    prob=prob_monbti_1673
    usfe_interp=monbti_1673
    xy_test=xy_monbti_1673
elif material == "tanbti_300":
    prob=prob_tanbti_300
    usfe_interp=tanbti_300
    xy_test=xy_tanbti_300
elif material == "tanbti_1673":
    prob=prob_tanbti_1673
    usfe_interp=tanbti_1673
    xy_test=xy_tanbti_1673
elif material == "monbti_rand":
    prob=prob_rand
    usfe_interp=monbti_rand
    xy_test=xy_monbti_rand
elif material == "tanbti_rand":
    prob=prob_rand
    usfe_interp=tanbti_rand
    xy_test=xy_tanbti_rand
else:
    print("ERROR material not available")


#GENERATE LATTICE WITH SRO -- (110) PLANE
#110 plane -- nonorthogonal
a1 = np.array([1,1,-1])/np.sqrt(3)
a2 = np.array([-1,1,1])/np.sqrt(3)
a3 = np.array([1,-1,1])/np.sqrt(3)

coords = np.zeros((N1*N2*N3,3))
pos = np.zeros((N1*N2*N3,3))
neighbors = np.zeros((N1*N2*N3,8),dtype=np.int)
index = np.zeros((N1*N2*N3),dtype=np.int)
for i in range(N1):
    for j in range(N2):
        for k in range(N3):
            idx = N1*N2*k + N1*j + i
            coords[idx] = i,j,k
            pos[idx] = i*a1 +  j*a2 + k*a3
            neighbors[idx,0] = N1*N2*k + N1*j + (i+1)%N1
            neighbors[idx,1] = N1*N2*k + N1*j + (i-1)%N1
            neighbors[idx,2] = N1*N2*k + N1*((j+1)%N2) + i
            neighbors[idx,3] = N1*N2*k + N1*((j-1)%N2) + i
            neighbors[idx,4] = N1*N2*((k+1)%N3) + N1*j + i
            neighbors[idx,5] = N1*N2*((k-1)%N3) + N1*j + i
            neighbors[idx,6] = N1*N2*((k+1)%N3) + N1*((j+1)%N2) + (i+1)%N1
            neighbors[idx,7] = N1*N2*((k-1)%N3) + N1*((j-1)%N2) + (i-1)%N1

print('lattice and neighbor lists built')

goal_bonds = np.zeros((3,3))
curr_bonds = np.zeros((3,3))

expected_upper = np.zeros((3,3))
expected_lower = np.zeros((3,3))

for i in range(3):
    for j in range(3):
        goal_bonds[i,j] = N1*N2*N3*8*prob[i,j]/3

        expected_upper[i,j] = np.ceil(prob[i,j]*8)
        expected_lower[i,j] = np.floor(prob[i,j]*8)

        if expected_lower[i,j] == 0:
            expected_lower[i,j] = 1
        elif expected_upper[i,j] == 8:
            expected_upper[i,j] = 7
print("goal bonds set")

neigh_types = np.zeros((len(coords),3),dtype=np.int)
lattice = np.zeros(len(coords),dtype=np.int)
lattice[0:len(coords)//3] = 1
lattice[len(coords)//3:2*len(coords)//3] = 2
rand.shuffle(lattice)
print('assigned initial atom types')

lattice_random = np.copy(lattice)

for i in range(len(coords)):
    for j in neighbors[i]:
        curr_bonds[lattice[i],lattice[j]] += 1
        neigh_types[i,lattice[j]] += 1
print('calculated current bond numbers')

pairs = [(0,1),(2,0),(1,2)]
swaps = 0
rejected = 0
tol = 0.0001*N1*N2*N3*8/3
current = [[],[],[]]
init_bonds = curr_bonds.copy()
init_time = time.time()
print('bond number tolerance = %f'%tol)

print('start swapping...')
#print(curr_bonds)
#print(goal_bonds)
#print(np.sum(lattice == 0),np.sum(lattice == 1),np.sum(lattice == 2))
while abs(goal_bonds[1,0] - curr_bonds[1,0]) > tol or abs(goal_bonds[2,0] - curr_bonds[2,0]) > tol or abs(goal_bonds[2,1] - curr_bonds[2,1]) > tol:
    for ip,p in enumerate(pairs):
        if curr_bonds[p] > goal_bonds[p]: #need to DECREASE unlike bonds. find 0 with more 1s than average
            first = np.argwhere((neigh_types[:,p[1]] > expected_upper[p]) & (lattice==p[0]))
            num = rand.randint(0,len(first))
            idx_0 = first[num][0]

            second = np.argwhere((lattice==p[1]) & (neigh_types[:,p[0]] > expected_upper[p]))
            num = rand.randint(0,len(second))
            idx_1 = second[num][0]
        else:
            first = np.argwhere((neigh_types[:,p[1]] < expected_lower[p]) & (lattice==p[0]))
            num = rand.randint(0,len(first))
            idx_0 = first[num][0]

            second = np.argwhere((lattice==p[1]) & (neigh_types[:,p[0]] < expected_lower[p]))
            num = rand.randint(0,len(second))
            idx_1 = second[num][0]

        m = np.zeros((2,2)) #neighbors around central atoms

        m[0,0] = neigh_types[idx_0,p[0]]
        m[0,1] = neigh_types[idx_0,p[1]]
        m[1,0] = neigh_types[idx_1,p[0]]
        m[1,1] = neigh_types[idx_1,p[1]]

        new = curr_bonds[p] + (m[1,1]+m[0,0]-m[0,1]-m[1,0])
        if abs(new - goal_bonds[p]) < abs(curr_bonds[p]-goal_bonds[p]):
            swaps += 1
            neighbors_0 = neighbors[idx_0]
            neighbors_1 = neighbors[idx_1]

            for n in neighbors_0:
                neigh_types[n,p[0]] -= 1
                neigh_types[n,p[1]] += 1
                curr_bonds[p[0],lattice[n]] -= 1
                curr_bonds[p[1],lattice[n]] += 1
                curr_bonds[lattice[n],p[0]] -= 1
                curr_bonds[lattice[n],p[1]] += 1
            for n in neighbors_1:
                neigh_types[n,p[1]] -= 1
                neigh_types[n,p[0]] += 1
                curr_bonds[p[1],lattice[n]] -= 1
                curr_bonds[p[0],lattice[n]] += 1
                curr_bonds[lattice[n],p[1]] -= 1
                curr_bonds[lattice[n],p[0]] += 1

            if idx_0 in neighbors_1:
                curr_bonds[p] += 2
                curr_bonds[p[1],p[0]] +=2
                curr_bonds[p[0],p[0]] -= 2
                curr_bonds[p[1],p[1]] -= 2

            lattice[idx_0] = p[1]
            lattice[idx_1] = p[0]

        else:
            rejected +=1
        current[ip].append(curr_bonds[p])
        if swaps % 10000 == 0:
            print("%d swaps done"%swaps)

print("Time to generate SRO lattice (seconds): %f" % (time.time()-init_time))
#print(curr_bonds)
#print(goal_bonds)

#DETERMINE LOCAL COMPOSITION

n_NN = [0,6,10,16,26,30,38,44,52,62,64,72,82]
neigh_comp = np.zeros((N1*N2*N3,n_NN[NN]),dtype=np.int)

for i in range(N1):
    for j in range(N2):
        for k in range(N3):
            idx = N1*N2*k + N1*j + i

            #1st NN
            neigh_comp[idx,0] = N1*N2*k + N1*j + (i+1)%N1
            neigh_comp[idx,1] = N1*N2*k + N1*((j+1)%N2) + i
            neigh_comp[idx,2] = N1*N2*k + N1*((j-1)%N2) + i
            neigh_comp[idx,3] = N1*N2*((k+1)%N3) + N1*j + i
            neigh_comp[idx,4] = N1*N2*((k-1)%N3) + N1*j + i
            neigh_comp[idx,5] = N1*N2*((k+1)%N3) + N1*((j+1)%N2) + (i+1)%N1

            if NN > 1:#2nd NN
                neigh_comp[idx,6] = N1*N2*((k+1)%N3) + N1*((j+1)%N2) + i
                neigh_comp[idx,7] = N1*N2*((k-1)%N3) + N1*((j-1)%N2) + i
                neigh_comp[idx,8] = N1*N2*k + N1*((j+1)%N2) + (i+1)%N1
                neigh_comp[idx,9] = N1*N2*((k+1)%N3) + N1*j + (i+1)%N1

            if NN > 2:#3rd NN
                neigh_comp[idx,10] = N1*N2*((k-1)%N3) + N1*((j+1)%N2) + i
                neigh_comp[idx,11] = N1*N2*((k+1)%N3) + N1*((j-1)%N2) + i
                neigh_comp[idx,12] = N1*N2*((k+1)%N3) + N1*((j+2)%N2) + (i+1)%N1
                neigh_comp[idx,13] = N1*N2*((k+2)%N3) + N1*((j+1)%N2) + (i+1)%N1
                neigh_comp[idx,14] = N1*N2*k + N1*((j-1)%N2) + (i+1)%N1
                neigh_comp[idx,15] = N1*N2*((k-1)%N3) + N1*j + (i+1)%N1

            if NN > 3:#4th NN
                neigh_comp[idx,16] = N1*N2*((k+1)%N3) + N1*((j+2)%N2) + i
                neigh_comp[idx,17] = N1*N2*((k+2)%N3) + N1*((j+1)%N2) + i
                neigh_comp[idx,18] = N1*N2*((k-1)%N3) + N1*((j-2)%N2) + i
                neigh_comp[idx,19] = N1*N2*((k-2)%N3) + N1*((j-1)%N2) + i
                neigh_comp[idx,20] = N1*N2*k + N1*((j+2)%N2) + (i+1)%N1
                neigh_comp[idx,21] = N1*N2*((k+2)%N3) + N1*j + (i+1)%N1
                neigh_comp[idx,22] = N1*N2*((k-1)%N3) + N1*((j+1)%N2) + (i+1)%N1
                neigh_comp[idx,23] = N1*N2*((k+1)%N3) + N1*((j-1)%N2) + (i+1)%N1
                neigh_comp[idx,24] = N1*N2*((k+2)%N3) + N1*((j+2)%N2) + (i+1)%N1
                neigh_comp[idx,25] = N1*N2*((k-1)%N3) + N1*((j-1)%N2) + (i+1)%N1

            if NN > 4:#5th NN
                neigh_comp[idx,26] = N1*N2*((k+2)%N3) + N1*j + i
                neigh_comp[idx,27] = N1*N2*((k-2)%N3) + N1*j + i
                neigh_comp[idx,28] = N1*N2*k + N1*((j+2)%N2) + i
                neigh_comp[idx,29] = N1*N2*k + N1*((j-2)%N2) + i

            if NN > 5:#6th NN
                neigh_comp[idx,30] = N1*N2*((k+1)%N3) + N1*((j+3)%N2) + (i+1)%N1
                neigh_comp[idx,31] = N1*N2*((k+3)%N3) + N1*((j+1)%N2) + (i+1)%N1
                neigh_comp[idx,32] = N1*N2*((k-1)%N3) + N1*((j+2)%N2) + i
                neigh_comp[idx,33] = N1*N2*((k+2)%N3) + N1*((j-1)%N2) + i
                neigh_comp[idx,34] = N1*N2*((k-2)%N3) + N1*j + (i+1)%N1
                neigh_comp[idx,35] = N1*N2*k + N1*((j-2)%N2) + (i+1)%N1
                neigh_comp[idx,36] = N1*N2*((k-2)%N3) + N1*((j+1)%N2) + i
                neigh_comp[idx,37] = N1*N2*((k+1)%N3) + N1*((j-2)%N2) + i

            if NN > 6:#7th NN
                neigh_comp[idx,38] = N1*N2*((k+2)%N3) + N1*((j+3)%N2) + (i+1)%N1
                neigh_comp[idx,39] = N1*N2*((k+3)%N3) + N1*((j+2)%N2) + (i+1)%N1
                neigh_comp[idx,40] = N1*N2*((k-1)%N3) + N1*((j+2)%N2) + (i+1)%N1
                neigh_comp[idx,41] = N1*N2*((k+2)%N3) + N1*((j-1)%N2) + (i+1)%N1
                neigh_comp[idx,42] = N1*N2*((k-1)%N3) + N1*((j-2)%N2) + (i+1)%N1
                neigh_comp[idx,43] = N1*N2*((k-2)%N3) + N1*((j-1)%N2) + (i+1)%N1

            if NN > 7:#8th NN
                neigh_comp[idx,44] = N1*N2*((k+0)%N3) + N1*((j+3)%N2) + (i+1)%N1
                neigh_comp[idx,45] = N1*N2*((k+3)%N3) + N1*((j+0)%N2) + (i+1)%N1
                neigh_comp[idx,46] = N1*N2*((k+1)%N3) + N1*((j+3)%N2) + (i+0)%N1
                neigh_comp[idx,47] = N1*N2*((k+3)%N3) + N1*((j+1)%N2) + (i+0)%N1
                neigh_comp[idx,48] = N1*N2*((k-3)%N3) + N1*((j-1)%N2) + (i+0)%N1
                neigh_comp[idx,49] = N1*N2*((k-1)%N3) + N1*((j-3)%N2) + (i+0)%N1
                neigh_comp[idx,50] = N1*N2*((k-2)%N3) + N1*((j+1)%N2) + (i+1)%N1
                neigh_comp[idx,51] = N1*N2*((k+1)%N3) + N1*((j-1)%N2) + (i+1)%N1

            if NN > 8:#9th NN
                neigh_comp[idx,52] = N1*N2*((k+0)%N3) + N1*((j+3)%N2) + (i+0)%N1
                neigh_comp[idx,53] = N1*N2*((k+3)%N3) + N1*((j+0)%N2) + (i+0)%N1
                neigh_comp[idx,54] = N1*N2*((k+0)%N3) + N1*((j-3)%N2) + (i+0)%N1
                neigh_comp[idx,55] = N1*N2*((k-3)%N3) + N1*((j+0)%N2) + (i+0)%N1
                neigh_comp[idx,56] = N1*N2*((k+3)%N3) + N1*((j+2)%N2) + (i+0)%N1
                neigh_comp[idx,57] = N1*N2*((k+2)%N3) + N1*((j+3)%N2) + (i+0)%N1
                neigh_comp[idx,58] = N1*N2*((k+3)%N3) + N1*((j+3)%N2) + (i+1)%N1
                neigh_comp[idx,59] = N1*N2*((k-3)%N3) + N1*((j-2)%N2) + (i+0)%N1
                neigh_comp[idx,60] = N1*N2*((k-2)%N3) + N1*((j-3)%N2) + (i+0)%N1
                neigh_comp[idx,61] = N1*N2*((k-2)%N3) + N1*((j-2)%N2) + (i+1)%N1

            if NN > 9:#10th NN
                neigh_comp[idx,62] = N1*N2*((k+2)%N3) + N1*((j-2)%N2) + (i+0)%N1
                neigh_comp[idx,63] = N1*N2*((k-2)%N3) + N1*((j+2)%N2) + (i+0)%N1

            if NN > 10:#11th NN
                neigh_comp[idx,64] = N1*N2*((k-1)%N3) + N1*((j+3)%N2) + (i+1)%N1
                neigh_comp[idx,65] = N1*N2*((k+3)%N3) + N1*((j-1)%N2) + (i+1)%N1
                neigh_comp[idx,66] = N1*N2*((k+4)%N3) + N1*((j+2)%N2) + (i+1)%N1
                neigh_comp[idx,67] = N1*N2*((k+2)%N3) + N1*((j+4)%N2) + (i+1)%N1
                neigh_comp[idx,68] = N1*N2*((k+2)%N3) + N1*((j-2)%N2) + (i+1)%N1
                neigh_comp[idx,69] = N1*N2*((k-2)%N3) + N1*((j+2)%N2) + (i+1)%N1
                neigh_comp[idx,70] = N1*N2*((k-3)%N3) + N1*((j-1)%N2) + (i+1)%N1
                neigh_comp[idx,71] = N1*N2*((k-1)%N3) + N1*((j-3)%N2) + (i+1)%N1

            if NN > 12:#12th NN
                neigh_comp[idx,72] = N1*N2*((k+4)%N3) + N1*((j+1)%N2) + (i+1)%N1
                neigh_comp[idx,73] = N1*N2*((k+1)%N3) + N1*((j+4)%N2) + (i+1)%N1
                neigh_comp[idx,74] = N1*N2*((k+3)%N3) + N1*((j+3)%N2) + (i+0)%N1
                neigh_comp[idx,75] = N1*N2*((k+3)%N3) + N1*((j-1)%N2) + (i+0)%N1
                neigh_comp[idx,76] = N1*N2*((k-1)%N3) + N1*((j+3)%N2) + (i+0)%N1
                neigh_comp[idx,77] = N1*N2*((k+1)%N3) + N1*((j-3)%N2) + (i+0)%N1
                neigh_comp[idx,78] = N1*N2*((k-3)%N3) + N1*((j+1)%N2) + (i+0)%N1
                neigh_comp[idx,79] = N1*N2*((k-3)%N3) + N1*((j-3)%N2) + (i+0)%N1
                neigh_comp[idx,80] = N1*N2*((k+0)%N3) + N1*((j-3)%N2) + (i+1)%N1
                neigh_comp[idx,81] = N1*N2*((k-3)%N3) + N1*((j+0)%N2) + (i+1)%N1

#translate atomic coords to compositions
comp = np.zeros((coords.shape))
for i in range(N1):
    for j in range(N2):
        for k in range(N3):
            idx = N1*N2*k + N1*j + i
            for n in neigh_comp[idx]:
                comp[idx,int(lattice[n])] += 1
            comp[idx,int(lattice[idx])] += 1
            comp[idx] /= len(neigh_comp[0])+1
print("Calculated local compositions")

#CONVERT COMPOSITION TO USFE
usfe = np.zeros((N1*N2*N3))
x = np.round(0.5*(2*comp[:,1]+comp[:,2])*2000)
y = np.round(comp[:,2]*1000)

z_test = np.zeros(len(xy_test))
z_test[-3] = 1
z_test[-2] = 1
z_test[-1] = 1
interp_test = griddata(xy_test, z_test, (grid_x, grid_y), method='linear')
sampled = np.zeros(3)

for i in range(len(usfe)):
    usfe[i] = usfe_interp[int(x[i]),int(y[i])]
    if not usfe[i] > 0: #checks if not nan
        print("error interpolating usfe")
    usfe_test = interp_test[int(x[i]),int(y[i])]
    if usfe_test == 0:
        sampled[0] += 1 #interpolated from inner data
    elif usfe_test > 0:
        sampled[1] += 1 #interpolated from corners
    else:
        sampled[2] += 1 #nan
        print(x[i],y[i])
print("Calculated local USFE values")
print("Percentage of USFE values covered by alloy USFE: %.3f"%(sampled[0]/np.sum(sampled)*100))

#WRITE INPUT FILE AND LATTICE FILE
filename="INPUT/%s_%d_%d_%d_NN%d_%d.txt" % (material,N1*32,N2,N3,NN,num_surf)
f = open(filename,'w')
#f.write("6291456\n")
f.write("%d %d %d %d\n"%(N1*32,N2,N3,3))
for li in range(32):
    for i in range(N1):
        for j in range(N2):
            for k in range(N3):
                for l in range(3):
                    f.write("%d %d %d %d %f\n"%(li*N1+i,j,k,l,usfe[N1*N2*k + N1*j + i]))
f.close()

filename2="INPUT/lattice_%s_%d_%d_%d_NN%d_%d.txt" % (material,N1*32,N2,N3,NN,num_surf)
f2 = open(filename2,'w')
f2.write("%d\n\n"%(len(coords)))
for i in range(len(coords)):
    f2.write("%f %f %f %d %f %f %f %f\n"%(pos[i,0],pos[i,1],pos[i,2],lattice[i],comp[i,0],comp[i,1],comp[i,2],usfe[i]))
f2.close()

print("created files %s and %s"%(filename,filename2))
