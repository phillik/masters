#!/usr/bin/env python3.7
import pymatgen_diffusion
import pymatgen
import io
import sys
import os

from pymatgen_diffusion import neb
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifWriter
from pymatgen_diffusion.neb import (full_path_mapper, pathfinder, io)
from pymatgen_diffusion.neb.full_path_mapper import FullPathMapper
from pymatgen_diffusion.neb.pathfinder import (IDPPSolver, MigrationPath, DistinctPathFinder)
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import (Specie, Element)
from pymatgen.core.sites import PeriodicSite

#Description: This script is a hodgepodge of different functions from the pymatgen_diffusion.neb module.
#Includes functionality for manual Path creation and file creation for NEB
#images using the IDPP algorithm. It also allows you to visualize all symmetrically possible migration
#pathways for a given species in a given structure. Instructions are given with each section.

#A quick tester I incorporated if you need to make sure you're actually loading the desired module correctly
#modulename = 'pymatgen_diffusion.neb.full_path_mapper'
#if modulename not in sys.modules:
#	print('You have not imported the module')

#mig_species is used as the input for full_path_mapper and DistinctPathFinder 
mig_species = sys.argv[1]

print('Modules Loaded')
#POSCAR should be CONTCAR from a relaxation of the target cell
poscar = Poscar.from_file("POSCAR")
LLZO = poscar.structure
print("Structure Loaded")
#I haven't quite figured out a full functionality for the FullPathMapper, but I used it here to map all sites in my structure 


#If I can find it, add a function to obtain all geometrically possible interstitial sites. Maybe call Voronoi_interstitial Finder

#The following 23 lines down to num+=1 are used to manually find the target sites you're looking for and create a MigrationPath class object
#This is necessary if you want to manually force a migration between two vacancies. 
#The MigrationPath container object requires two PeriodicSite objects for the initial and end sites as well as a symmetrized_structure 
txt1= input("Create a MigrationPath Object Manually?(y/n): ")
if txt1=="y":

	#I haven't quite figured out a full functionality for the FullPathMapper, but I used it here to map all sites in my structure
	LLZO_full_paths = FullPathMapper(LLZO, mig_species)
	all_sites = LLZO_full_paths.get_only_sites() #stores a structure type object of all possible mig_species vacancy sites
	print("Full Paths Loaded")
	
	t=all_sites.sites #Stores structure class attribute PeriodicSite in a list of Tuples
	num=1 
	for z in t:
	 print(f'{num}: {z}') #Use to find placeholder for target sites
	 num+=1

	print("The placeholder number for each atom is found at the far left next to each item in the list above.")
	s1= int(input("Enter Initial Site number: "))
	s2= int(input("Enter End Site number: "))
	check1 = isinstance(s1, int)
	check2 = isinstance(s2, int)
	if not check1 or check2:
	 print("Input is not an Integer")
	 sys.exit()
	num=1
	for z in t: 
	 if num==s1: #Initial Site
	  print(num,z)
	  site1=z #Stores a PeriodicSite object. This object stores cartesian coordinates, Lattice params, and atom type. Found in pymatgen.core.sites
	  num+=1
	  continue
	 if num==s2: #end site
	  print(num,z)
	  site2=z
	  num+=1
	 else:
	  num+=1
	LLZO_temp = SpacegroupAnalyzer(LLZO) #SpacegroupAnalyzer Class includes many really cool functionalities. Pymatgen.symmetry.analyzer module
	LLZO_SS = LLZO_temp.get_symmetrized_structure() #Here I use it to obtain a symmetrized structure to create my MP object
	LLZO_MP = MigrationPath(site1,site2,LLZO_SS)  #Create MP object
	n_images = int(input("Number of images (Total structures will be n+2): "))
	IDPP_check = input("Use IDPP algorithm?(y/n): ")
	if IDPP_check=="y":
	 c = True
	else:
	 c= False
	u = LLZO_MP.get_structures(n_images,True,c)  #This is the great part. (nimages=5,vac_mode=True,idpp=True) This gives an improved guess
	#for the initial path using the Independent Pair Potential algorithm. This should drastically cut down on convergence times. This also
	#stores the migrating atom in the zero position of the POSCAR file for simplified analysis. 
	num=0
	LLZO_MP.write_path('test_path.cif') #So that you can visualize your full path 
	for e in u: #loops through the list of structure objects and stores each into your home directory with their correct number. 0 is initial.
	 e.to(filename=f"pos{num}.vasp",fmt='poscar')
	 num+=1
	#sites = CifWriter(all_sites) #if you prefer a cif file
	#sites.write_file('all_sites.cif')

Dpaths = DistinctPathFinder(LLZO, mig_species) #Finds all symmetrically inequivalent migration pathways

txt2= input("Write all possible Migration Paths?(y/n): ")
if txt2=="y":
	
 print("1: All symmetrically inequivalent paths")
 print("2: All possible symmetrical paths")
 option = str(input("Enter 1 or 2: "))
 if option=="1":
  Dpaths.write_all_paths('all_paths_found.cif') #uses hydrogens as placeholders to show all distinct migration pathways possible in the given structure

 if option=="2":
  LLZO_temp = SpacegroupAnalyzer(LLZO) #SpacegroupAnalyzer Class includes many really cool functionalities. Pymatgen.symmetry.analyzer module
  LLZO_SS = LLZO_temp.get_symmetrized_structure()
  
  paths = []
  ms = mig_species + '1' #get_el_sp(mig_species)
  for s in LLZO_SS.equivalent_sites:
   for s2 in s:
    if str(s2.species) == ms:
     for nn in LLZO_SS.get_neighbors(s2, r=round(3, 3) + 0.01):
      if str(nn.species) == ms:
       path = MigrationPath(s2, nn, LLZO_SS)
       paths.append(path)
  n=1
  sites = []
  for p in paths:
   structures = p.get_structures(5,ms)
   sites.append(structures[0][0])
   sites.append(structures[-1][0])
   for s in structures[1:-1]:
    sites.append(PeriodicSite("H", s[0].frac_coords, s.lattice))
   if n == (len(paths)):
    sites.extend(structures[0].sites[1:])
   n+=1
  Structure.from_sites(sites).to(filename="equiv_paths_found.cif")
  LLZO2 = Structure.from_sites(sites)

 else:
  print("okay, nm")	

#Use the list to determine the placeholder number for the desired MigrationPath. Enter these numbers in the prompt to extract and write
#your NEB images. It once again places the migrating ion first in the POSCAR file for easy analysis. Also useful as it creates a vacancy at the 
#end site.
txt3= input("Obtain NEB images from generated MigrationPath object?(y/n): ")
if txt3=="y":
	   
	MPaths = Dpaths.get_paths() 
	n=1
	if txt2 == 'y':
	 if option == "2":
	  MPaths = paths
	for MP in MPaths:
	 print(n, MP.isite, MP.esite) #Lists each placeholder n alongside the MP initial and final site. Use this list to find target MP
	 MP.write_path(f'{n}path.cif') #Creates a small file to visualize the path and position of target species migration.
	 n+=1
	
	MP_num = input("Enter MigrationPath number (press d when finished: ")
	MP_num_list = []
	while MP_num != "d":
		MP_num_list.append(int(MP_num))
		MP_num = input("Enter another MigrationPath number? (press d when finished: ")	

	n=1
	n_images = int(input("Number of images (Total structures will n+2): "))
	IDPP_check = input("Use IDPP algorithm?(y/n): ")
	if IDPP_check=="y":
	 c = True
	else:
	 c= False
	for MP in MPaths:	
	 if n in MP_num_list: #Insert Placeholder number here
	  x = MP.get_structures(n_images,True,c)
	  print("BOOM BABY!!!")
	  count=0 #for nomenclature
	  for y in x:
	   y.to(filename=f"path{n}_pos{count}.vasp",fmt='poscar')
	   count+=1
	   print("Ba-da Bing")
	 n+=1

# Once in bash,create a new directory for your chosen path and use < for i in <MP_nums>;do mkdir path$i; for x in {0..n_images}; do  mkdir path$i/0$x;t=path$i;v=_pos$x; mv $t$v.vasp path$i/0$x/POSCAR; done; done > to organize your files into their subdirectories
