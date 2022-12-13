def scm(pdb_file,top_file,R):
    '''
    in house code
    author:zhanpeng@CADD
    
    pdb_file:pdb
    top_file:prmtop
    R:neighbor distance,stored in nm (10 angstrom in paper,here should be 1)
    R=1 or 0.5 
    
    this function calculate the SCM_score of Fv model,as defined in paper:
    http://dx.doi.org/10.1080/19420862.2015.1099773
    
    mdtraj default distance stored in nm(nanometers)
    
    '''
    ###################################################
    ###   step 1 compute arr_charge_exposed_maped   ###
    ###################################################
    
    ##### use pytraj to get charge #####
    traj = pt.load(pdb_file,top_file)
    
    # slice Traj obj: remove solvent and return a Traj obj only with protein(Fv)
    traj_prot = pt.strip(traj,":WAT,Cl-")
    
    # get the array,containg the charge of each atoms in Fv,default indexed as the atoms index
    arr_charge = traj_prot.topology.charge
    
    
    ###### use mdtraj to get neighbors,sasa #######
    
    # mdtraj load file,return a Trajectory object
    md_traj = md.load(pdb_file,top=top_file)  
    
    # get atom index of protein sidechain
    prot_sc_idx = md_traj.top.select('protein and sidechain')  
    
    # slice Traj obj,return a Traj only with protein sidechain atoms,with the idx defined in the pre step
    md_traj_prot_sc = md_traj.atom_slice(prot_sc_idx) 
    
    # calculating the SASA of each residue (only sidechain here)
    sasa = md.shrake_rupley(md_traj_prot_sc,mode='residue') 
    
    # identify the exposed residues, and mask each residue as 1 if SASA > sasa_cutoff,otherwise 0,return a mask array
    sasa_cutoff = 0.1
    exposed_residues_map = np.heaviside(sasa[0] - sasa_cutoff,sasa[0] - sasa_cutoff) 
    '''
    ref:mdtraj
    "All of the distances in the Trajectory are stored in nanometers" 
    "https://mdtraj.org/1.9.4/examples/introduction.html?highlight=nano#"
    
    in paper sasa_cutoff=10 square angstrom
    0.1 square nm = 10 square angstrom,so here sasa_cutoff=0.1
    '''
    
    # convert exposed residue map to exposed atom map。here map means mask
    exposed_atoms_map=[]
    for i in range(len(exposed_residues_map)):
        # traverse residus，and get the atom_idx of each residue
        atoms_idx_of_resid=md_traj.top.select('protein and resid '+str(i))
        for j in range(len(atoms_idx_of_resid)):
            # 遍历每个氨基酸的每个原子，并将每个原子的exposed map追加到列表exposed_atoms_map
            exposed_atoms_map.append(exposed_residues_map[i])
    
    # convert list to np array
    arr_exposed_atoms_map = np.asarray(exposed_atoms_map)
    # charge map: the charge of atoms in no exposed residues was masked to 0
    arr_charge_exposed_maped = arr_charge*arr_exposed_atoms_map


    ######################################
    ###   step 2 compute SCM_score     ###
    ######################################
    # build a new enpty list to contain the SCM of each atoms
    l_SCM_atom_i=[]
    
    # get the protein atoms idx
    prot_idx = md_traj.top.select('protein')
    
    # slice Traj obj,return an obj only with protein
    md_traj_prot = md_traj.atom_slice(prot_idx)
    '''print(len(prot_idx))'''

    for i in prot_idx:
        
        #get neighbor index of each atoms
        neighbor_idx_of_atom_i = md.compute_neeighbors(md_traj_prot, R, np.asarray([i]))

        # sum the mapped neighbor charge
        SCM_atom_i=np.sum(arr_charge_exposed_maped[neighbor_idx_of_atom_i])
        
        # append SCM_atom_i to empyty list
        l_SCM_atom_i.append(SCM_atom_i)

    # convert list to np array
    arr_l_SCM_atom_i=np.asanyarray(l_SCM_atom_i)
    
    #paper公式第2步，heaviside，sum，abs
    SCM_score = np.absolute(np.sum(arr_l_SCM_atom_i * np.heaviside(-1 * arr_l_SCM_atom_i,-1 * arr_l_SCM_atom_i)))
    return SCM_score
