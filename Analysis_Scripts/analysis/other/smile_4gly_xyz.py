#!/usr/bin/env python
# coding: utf-8

# In[300]:


import numpy as np
import sys
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import os


# In[301]:


path0='/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/020_1e-6_add_water/008_gly128H2O_sd_B35_compress_f'
path1='/0-30ns/gaussian/xyz'
path2='/0-30ns/gaussian/gas/'
path3='/0-30ns/gaussian/smd/'
path4=['anion','cation','neutral','zwitterion']
data1=['animate_anion_10-20ns.xyz','animate_cation_10-20ns.xyz',
      'animate_neutral_10-20ns.xyz','animate_zwitterion_10-20ns.xyz']
data2=['smile_anion.xyz','smile_cation.xyz',
      'smile_neutral.xyz','smile_zwitterion.xyz']


# In[302]:


choose_geo=3
path4=['anion','cation','neutral','zwitterion'][choose_geo]
myinput=os.path.join(path0+path1, data1[choose_geo])
myoutput=os.path.join(path0+path1,data2[choose_geo])
choose_q=[-1,1,0,0][choose_geo]


# In[303]:


with open(myinput, "r") as f:
    lines = f.readlines()
    latom=[]
    lidex=[]
    numfr=0
    j=0
    for line in lines:
        if 'Frame' in line:
            numfr+=1
            latom.append(int(lines[j-1]))
            lidex.append(j-1)
        j+=1


# In[304]:


def read_multi_xyz_file(xyz_filename,atoms,idx):
    with open(xyz_filename, 'r') as f:
        lines = f.readlines()
        num_atoms = atoms
        choose_line = lines[idx:idx+atoms+2]
        #print(choose_line)
        comment = choose_line[1]
        atomic_symbols = []
        atomic_positions = []
        for line in choose_line[2:]:
            fields = line.split()
            #print(fields)
            atomic_symbols.append(fields[0])
            atomic_positions.append([float(x) for x in fields[1:4]])
        assert len(atomic_symbols) == num_atoms
        assert len(atomic_positions) == num_atoms
    return atomic_symbols, atomic_positions, comment


# In[305]:


def count_short_bonds(coords_by_element, O_symbol, H_symbol, max_bond_length):
    count = 0
    for O_coord in coords_by_element.get(O_symbol, []):
        for H_coord in coords_by_element.get(H_symbol, []):
            distance = np.sqrt((H_coord[0] - O_coord[0])**2 + (H_coord[1] - O_coord[1])**2 + (H_coord[2] - O_coord[2])**2)
            if distance < max_bond_length:
                count += 1
    return count


# In[306]:


def write_gaussian_input(coords_by_element, file_path, title="", functional="m062x", 
                         basis_set="def2tzvp", solvent=None, charge=0, multiplicity=1, 
                         nprocshared=20, mem="30GB"):
    """
    This function takes in atomic coordinates and additional 
    Gaussian calculation options and outputs a Gaussian input file.
    """
    # Open the output file and write the start of the Gaussian input file
    with open(file_path, 'w') as f:
        f.write("%"+"nprocshared=%d\n"%nprocshared)
        f.write("%"+"mem=%s\n"%mem)
        if solvent == 'smd':
            f.write("#p %s/%s int=ultrafine freq scrf=%s\n\n"%(functional,basis_set,solvent))
        else:
            f.write("#p %s/%s int=ultrafine freq\n\n"%(functional,basis_set))
        f.write("Frame %s\n\n"%title)
        f.write("%d %d\n"%(charge,multiplicity))
        # Add atomic coordinates to the Gaussian input file
        for element, coords in coords_by_element.items():
            for coord in coords:
                f.write(f" {element}  {'  '.join([f'%16.8f' % x for x in coord])}\n")
        f.write("\n\n")


# In[307]:


def get_mkdir(path):
    dir_path = os.path.dirname(path)
    if not os.path.exists(path):
        os.makedirs(path)


# In[308]:


# Initialize a dictionary to hold the SMILES strings and their frequencies
smiles_dict = {}
for i in range(numfr):#numfr
    #print(latom[i],lidex[i])
    atomic_symbols, atomic_positions, comment = read_multi_xyz_file(myinput,latom[i],lidex[i])
    rd_mol = Chem.rdchem.RWMol()
    for symbol, xyz in zip(atomic_symbols, atomic_positions):
        atom = Chem.rdchem.Atom(symbol)
        atom.SetDoubleProp('x', xyz[0])
        atom.SetDoubleProp('y', xyz[1])
        atom.SetDoubleProp('z', xyz[2])
        rd_mol.AddAtom(atom)
    mol = rd_mol.GetMol()

    # Generate the canonical SMILES representation of the molecule
    smiles = Chem.MolToSmiles(mol)
    
    # Print the SMILES string
    #print(smiles)
    if smiles in smiles_dict:
        smiles_dict[smiles] += 1
    else:
        smiles_dict[smiles] = 1


# In[309]:


# Find the SMILES string that appears the most in the input file
most_common_smiles = max(smiles_dict, key=lambda x: smiles_dict[x])

# Print out all the SMILES strings and their cumulative number of occurrences
print("SMILES string   |   Number of occurrences")
print("-----------------------------------------")
for smiles, count in smiles_dict.items():
    print(smiles + "   |   " + str(count))

# Print out the most common SMILES string
print("\nThe most common SMILES string is: " + most_common_smiles)


# In[310]:


g=0
for i in range(numfr):#numfr
    #print(latom[i],lidex[i])
    atomic_symbols, atomic_positions, comment = read_multi_xyz_file(myinput,latom[i],lidex[i])
    rd_mol = Chem.rdchem.RWMol()
    for symbol, xyz in zip(atomic_symbols, atomic_positions):
        atom = Chem.rdchem.Atom(symbol)
        atom.SetDoubleProp('x', xyz[0])
        atom.SetDoubleProp('y', xyz[1])
        atom.SetDoubleProp('z', xyz[2])
        rd_mol.AddAtom(atom)
    mol = rd_mol.GetMol()
    
    if Chem.MolToSmiles(mol)==most_common_smiles:
        coords_by_element = {}
        for atom in mol.GetAtoms():
            element = atom.GetSymbol()
            x, y, z = float(atom.GetProp('x')), float(atom.GetProp('y')), float(atom.GetProp('z'))
            #print(f"{element}: {x}, {y}, {z}")
            if element not in coords_by_element:
                coords_by_element[str(element)] = [[x, y, z]]
            else:
                coords_by_element[str(element)].append([x, y, z]) 

        symbol1    = ['N','O','C','C','C','C']
        symbol2    = ['H','H','C','O','H','N']
        max_length = [1.3,1.3,1.9,1.8,1.4,1.8]    #Angstroms
        num_shortb = [[2  ,0  ,4  ,2  ,2  ,1  ],  #anion
                      [3  ,1  ,4  ,2  ,2  ,1  ],  #cation
                      [2  ,1  ,4  ,2  ,2  ,1  ],  #neutral
                      [3  ,0  ,4  ,2  ,2  ,1  ]   #zwitterion
                     ][choose_geo]  
        flag=0
        for ii in range(len(symbol1)):
            num_short_bonds = count_short_bonds(coords_by_element, symbol1[ii], symbol2[ii], max_length[ii])
            if num_short_bonds==num_shortb[ii]:
                flag+=1
            #print(f"Number of {symbol1[ii]}-{symbol2[ii]} bonds less than {max_length[ii]} Angstroms: {num_short_bonds}")
        #print(flag)
        
        if flag==6:
            # Create directory if it doesn't exist
            g+=1
            gaspath=os.path.join(path0+path2+path4)
            smdpath=os.path.join(path0+path3+path4)
            get_mkdir(os.path.join(gaspath,'%d'%(10000+g)))
            get_mkdir(os.path.join(smdpath,'%d'%(10000+g)))
            write_gaussian_input(coords_by_element, file_path=os.path.join(gaspath,'%d/gas.gjf'%(10000+g)), 
                                 title="%d"%(i+1),solvent=None, charge=choose_q, multiplicity=1)
            write_gaussian_input(coords_by_element, file_path=os.path.join(smdpath,'%d/smd.gjf'%(10000+g)), 
                                 title="%d"%(i+1),solvent='smd', charge=choose_q, multiplicity=1)            


# In[ ]:




