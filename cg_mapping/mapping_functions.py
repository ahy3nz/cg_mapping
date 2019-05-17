import time
import os
import pdb
import warnings
import itertools
from multiprocessing import Pool
from collections import OrderedDict
from optparse import OptionParser

from lxml import etree
import numpy as np
import mdtraj
import mbuild as mb

import cg_mapping
from cg_mapping.CG_bead import CG_bead

def load_mapping(mapfile=None):
    """ Load a forward mapping
    
    Parameters
    ----------
    mapfile : str
        Path to mapping file

    Returns
    -------
    mapping_dict : OrderedDict()
        OrderedDict (CG bead index : [beadtype, list of atom indices])
    
    Notes
    -----
    mapping files are ":" delimited
    col0: bead index
    col1: bead type
    col2: atom indices

    xml files are arranged by mapping > molecule > beads + bonds
        """
    extension = mapfile[-4:]
    if 'map' in extension:
        return _load_map_file(mapfile=mapfile)
    elif 'xml' in extension:
        return _load_xml_file(mapfile=mapfile)
    else:
        sys.exit("Invalid mappingfile")

def _load_map_file(mapfile=None):
    """ Load from a .map file """
    mapping_dict = OrderedDict()
    bonding_info = []
    with open(mapfile,'r') as f:
        for line in f:
            if line.rstrip():
                if 'bond' in line.split(":")[0]:
                    atom_i = line.split(":")[1].split()[0].rstrip()
                    atom_j = line.split(":")[1].split()[1].rstrip()
                    bonding_info.append((atom_i, atom_j))

                else:
                    atom_indices = [int(i) for i in line.split(":")[2].strip().split()]
                    mapping_dict.update({line.split(":")[0].rstrip():
                            [line.split(":")[1].rstrip(), atom_indices]})

    return mapping_dict, bonding_info

def _load_xml_file(mapfile=None):
    """ load from a .xml file """
    mapping_dict = OrderedDict()
    bonding_info = []
    tree = etree.parse(mapfile)
    root = tree.getroot()
    molecule = root.find("Molecule")

    beads = molecule.find("Beads")
    for bead in beads.iterchildren():
        if 'Bead' in bead.tag:
            mapping_dict.update({bead.attrib['index'].strip():
                                [bead.attrib['beadtype'], 
                               [int(i) for i in bead.attrib['map'].strip().split()]]})
        else:
            sys.exit("Unidentified tag {}".format(bead.tag))
    
    bonds = molecule.find("Bonds")
    for bond in bonds.iterchildren():
        if "Bond" in bond.tag:
            bonding_info.append((bond.attrib['bead1'], bond.attrib['bead2']))
        else:
            sys.exit("Unidentified tag {}".format(bond.tag))

    return mapping_dict, bonding_info

   
def create_CG_topology(topol=None, all_CG_mappings=None, water_bead_mapping=4,
        all_bonding_info=None):
    """ Create CG topology from given topology and mapping

    Parameters
    ---------
    topol : mdtraj Topology
    all_CG_mappings : dict
        maps residue names to respective CG 
        mapping dictionaries(CG index, [beadtype, atom indices])
    water_bead_mapping : int
        specifies how many water molecules get mapped to a water CG bead
    all_bonding_info : dict
        maps residue names to bonding info arrays 
        np.ndarray (n, 2)


    Returns
    -------
    CG_topology_map : tuple of CG_bead()
    CG_topology : mdtraj topology
        Need to fill in more details for topology creation

    """
    CG_topology_map = []
    CG_topology = mdtraj.Topology()
    water_counter = 0
    # Loop over all residues
    for residue in topol.residues:
        if not residue.is_water:
            # Obtain the correct molecule mapping based on the residue
            molecule_mapping = all_CG_mappings[residue.name]
            temp_residue = CG_topology.add_residue(residue.name, CG_topology.add_chain())
            temp_CG_indices = []
            temp_CG_atoms = []
            temp_CG_beads = [None]*len(molecule_mapping.keys())
            CG_atoms = []

            # For each CG bead, construct CG bead object, relating it to the
            # global atom indices
            # This first requires getting the local atom indices for that residue
            for key in molecule_mapping.keys():
                temp_CG_atoms = [residue.atom(index) for index in molecule_mapping[key][1]]
                new_bead = CG_bead(beadindex=0,
                                   beadtype=molecule_mapping[key][0],
                                   resname=residue.name,
                                   atom_indices=[atom.index for atom in temp_CG_atoms])
                temp_CG_beads[int(key)] = new_bead

            #for index, atom in enumerate(residue.atoms):
            #    temp_CG_indices.append(str(index))
            #    temp_CG_atoms.append(atom)
            #    for key in molecule_mapping.keys():
            #        if set(molecule_mapping[key][1]) == set(temp_CG_indices):
            #            new_bead = CG_bead(beadindex=0, 
            #                               beadtype=molecule_mapping[key][0],
            #                               resname=residue.name,
            #                               atom_indices=[atom.index for atom in temp_CG_atoms])
            #            temp_CG_indices = []
            #            temp_CG_atoms = []
            #            temp_CG_beads[int(key)] = new_bead

            #        else:
            #            pass


            # Add beads to topology by adding atoms
            for index, bead in enumerate(temp_CG_beads):
                bead.beadindex = index
                CG_topology_map.append(bead)
                new_CG_atom = CG_topology.add_atom(bead.beadtype, None, temp_residue)
                CG_atoms.append(new_CG_atom)

            # Add bonds to topolgoy 
            bonding_info = all_bonding_info[residue.name]
            for (index_i, index_j) in bonding_info:
                CG_topology.add_bond(CG_atoms[int(index_i)], CG_atoms[int(index_j)])

        else:
            water_counter +=1
            if water_counter % water_bead_mapping == 0:
                temp_residue = CG_topology.add_residue("HOH", CG_topology.add_chain())
                new_bead = CG_bead(beadindex=0, beadtype="W",
                        resname='HOH')
                CG_topology_map.append(new_bead)
                CG_topology.add_atom("W", None, temp_residue)





    return CG_topology_map, CG_topology

def _map_waters(traj, frame_index):
    """ Worker function to parallelize mapping waters via kmeans

    Parameters
    ----------
    traj : mdtraj trajectory
        full atomstic trajectory
    frame index : int
        parallelizing calculation frame by frame
    
    """
    from sklearn import cluster
    frame = traj[frame_index]
    # Get atom indices of all water oxygens
    waters = frame.topology.select('water and name O')

    # Get coordinates and number of all water oxygens
    n_aa_water = len(waters)
    aa_water_xyz = frame.atom_slice(waters).xyz[0,:,:]

    # Number of CG water molecules based on mapping scheme
    water_bead_mapping = 4
    n_cg_water = int(np.floor(n_aa_water /  water_bead_mapping))
    # Water clusters are a list (n_cg_water) of empty lists
    water_clusters = [[] for i in range(n_cg_water)]

    # Perform the k-means clustering based on the AA water xyz
    k_means = cluster.KMeans(n_clusters=n_cg_water)
    k_means.fit(aa_water_xyz)

    # Each cluster index says which cluster an atom belongs to
    for atom_index, cluster_index in enumerate(k_means.labels_):
        # Sort each water atom into the corresponding cluster
        # The item being added should be an atom index
        
        water_clusters[cluster_index].append(waters[atom_index])


    single_frame_coms = []
    # For each cluster, compute enter of mass
    for cg_index, water_cluster in enumerate(water_clusters):
        com = mdtraj.compute_center_of_mass(frame.atom_slice(water_cluster))
        single_frame_coms.append((frame_index, com))
        #CG_xyz[frame_index, cg_index + water_start,:] = com
    return single_frame_coms
 
def convert_xyz(traj=None, CG_topology_map=None, water_bead_mapping=4,parallel=True):
    """Take atomistic trajectory and convert to CG trajectory

    Parameters
    ---------
    traj : mdtraj Trajectory
        Atomistic trajectory
    CG_topology : list
        list of CGbead()
    parallel : boolean
        True if using parallelized, false if using serial

    Returns
    ------
    CG_xyz : np.ndarray(n_frame, n_CG_beads, 3)

    
    """
    # Iterate through the CG_topology
    # For each bead, get the atom indices
    # Then slice the trajectory and compute hte center of mass for that particular bead
    entire_start = time.time()
    CG_xyz = np.ndarray(shape=(traj.n_frames, len(CG_topology_map),3))
    water_indices = []
    print("Converting non-water atoms into beads over all frames")
    start = time.time()
    for index, bead in enumerate(CG_topology_map):
        if 'HOH' not in bead.resname:
            # Handle non-water residuse with center of mass calculation over all frames
            atom_indices = bead.atom_indices 
            # Two ways to compute center of mass, both are pretty fast
            bead_coordinates = mdtraj.compute_center_of_mass(traj.atom_slice(atom_indices))
            CG_xyz[:, index, :] = bead_coordinates
        else:
            # Handle waters by initially setting the bead coordinates to zero
            # Remember which coarse grain indices correspond to water
            water_indices.append(index)

    end = time.time()
    print("Converting took: {}".format(end-start))

    # Figure out at which coarse grain index the waters start
    print("Converting water beads via k-means")
    start = time.time()
    if len(water_indices)>0:
        # Perform kmeans, frame-by-frame, over all water residues
        # Workers will return centers of masses of clusters, frame index, and cg index
        # Master will assign to CG_xyz
        if parallel:
            all_frame_coms = []
            with Pool() as p:
                all_frame_coms = p.starmap(_map_waters, 
                                        zip(itertools.repeat(traj), 
                                        range(traj.n_frames),
                                        itertools.repeat(water_bead_mapping)))

            end = time.time()
            print("K-means and converting took: {}".format(end-start))

            print("Writing to CG-xyz")
            start = time.time()
            for snapshot in all_frame_coms:
                for element, water_index in zip(snapshot, water_indices):
                    CG_xyz[element[0], water_index , : ] = element[1]
            end =  time.time()
            print("Writing took: {}".format(end-start))

            return CG_xyz
        
        else:
            from sklearn import cluster
            for frame_index, frame in enumerate(traj):
                # Get atom indices of all water oxygens
                waters = traj.topology.select('water and name O')

                # Get coordinates and number of all water oxygens
                n_aa_water = len(waters)
                aa_water_xyz = traj.atom_slice(waters).xyz[frame_index,:,:]

                # Number of CG water molecules based on mapping scheme
                n_cg_water = int(n_aa_water /  water_bead_mapping)
                # Water clusters are a list (n_cg_water) of empty lists
                water_clusters = [[] for i in range(n_cg_water)]

                # Perform the k-means clustering based on the AA water xyz
                start = time.time()
                print("Clustering for frame {}".format(frame_index))
                k_means = cluster.KMeans(n_clusters=n_cg_water)
                k_means.fit(aa_water_xyz)
                end = time.time()
                print("Clustering one frame took: {}".format(end-start))

                # Each cluster index says which cluster an atom belongs to
                print("Assigning water atoms to water clusters for frame {}".format(frame_index))
                start = time.time()
                for atom_index, cluster_index in enumerate(k_means.labels_):
                    # Sort each water atom into the corresponding cluster
                    # The item being added should be an atom index
                    water_clusters[cluster_index].append(waters[atom_index])
                end = time.time()
                print("Assigning took: {}".format(end-start))


                # For each cluster, compute enter of mass
                print("Computing cluster centers for frame {}".format(frame_index))
                start = time.time()
                for cg_index, water_cluster in enumerate(water_clusters):
                    com = mdtraj.compute_center_of_mass(traj.atom_slice(water_cluster)[frame_index])
                    CG_xyz[frame_index, cg_index + water_start,:] = com
                end = time.time()
                print("Computing took: {}".format(end-start))

    entire_end = time.time()
    print("XYZ conversion took: {}".format(entire_end - entire_start))
    return CG_xyz

def compute_avg_box(traj):
    """ Compute average box lengths"""

    return np.mean(traj.unitcell_lengths, axis=0)


