import time
import os
import pdb
import warnings
from collections import OrderedDict
from optparse import OptionParser

import mdtraj
import mbuild as mb

import cg_mapping.mapping_functions as mapping_functions

PATH_TO_MAPPINGS='/raid6/homes/ahy3nz/Programs/cg_mapping/cg_mapping/charmm_mappings/'
HOOMD_FF="/raid6/homes/ahy3nz/Programs/setup/FF/CG/msibi_ff.xml"

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", action="store", type="string", dest = "trajfile", default='last20.xtc')
    parser.add_option("-c", action="store", type="string", dest = "topfile", default='md_pureDSPC.pdb')
    parser.add_option("-o", action="store", type="string", dest = "output", default='cg-traj')
    (options, args) = parser.parse_args()
    
    
    traj = mdtraj.load(options.trajfile, top=options.topfile)
    topol = traj.topology
    start=time.time()
    # Read in the mapping files, could be made more pythonic
    DSPCmapfile = os.path.join(PATH_TO_MAPPINGS,'DSPC.xml')#'mappings/DSPC.map'
    #watermapfile = os.path.join(PATH_TO_MAPPINGS,'water.map')
    #alc16mapfile = os.path.join(PATH_TO_MAPPINGS,'C16OH.map')
    #acd16mapfile = os.path.join(PATH_TO_MAPPINGS,'C16FFA.map')
    # Huge dictionary of dictionaries, keys are molecule names
    # Values are the molecule's mapping dictionary
    # could be made more pythonic
    all_CG_mappings = OrderedDict()
    all_bonding_info = OrderedDict()
    
    molecule_mapping, molecule_bonding = mapping_functions.load_mapping(mapfile=DSPCmapfile)
    all_CG_mappings.update({'DPPC': molecule_mapping})
    all_bonding_info.update({'DPPC': molecule_bonding})
    
    CG_topology_map, CG_topology = mapping_functions.create_CG_topology(topol=topol, 
                            all_CG_mappings=all_CG_mappings, 
                            all_bonding_info=all_bonding_info)
    CG_xyz = mapping_functions.convert_xyz(traj=traj, CG_topology_map=CG_topology_map)
    
    CG_traj = mdtraj.Trajectory(CG_xyz, CG_topology, time=traj.time, 
            unitcell_lengths=traj.unitcell_lengths, unitcell_angles = traj.unitcell_angles)
    
    
    avg_box_lengths = mapping_functions.compute_avg_box(traj)
    print("Avg box length: {}".format(avg_box_lengths))
    
    CG_traj.save('{}.xtc'.format(options.output))
    CG_traj[-1].save('{}.gro'.format(options.output))
    CG_traj[-1].save('{}.h5'.format(options.output))
    CG_traj[-1].save('{}.xyz'.format(options.output))
    CG_traj[-1].save('{}.pdb'.format(options.output))
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mb_compound = mb.Compound()
        mb_compound.from_trajectory(CG_traj, frame=-1, coords_only=False)
        original_box = mb.Box(lengths=[length for length in mb_compound.periodicity])
    
        # Because we are resizing the box based on the average box over the trajectory,
        # we need to scale the coordinates appropriately, since they are taken only
        # from the final frame.
        scaling_ratio = [new/old for old, new in zip(original_box.lengths, avg_box_lengths)]
        print("scaling_ratio: {}".format(scaling_ratio))
        for particle in mb_compound.particles():
            for i, frame in enumerate(particle.xyz):
                particle.xyz[i] = [factor*coord for factor, coord in zip(scaling_ratio, particle.xyz[i])]
    
        # Tile this using grid 3d pattern
        cube = mb.Grid3DPattern(2,2,2)
        # Define a new box based on average box lengths from AA traj, scaled by 8
        new_box = mb.Box(lengths=[2*length for length in avg_box_lengths])
        # Scale pattern lengths based on the new box lengths
        cube.scale([length for length in new_box.lengths])
        replicated = cube.apply(mb_compound)
        mirrored_image = mb.Compound()
        for item in replicated:
            mirrored_image.add(item)
    
        # Particle renaming due to mbuild coarsegrained format
        for particle in mb_compound.particles():
            particle.name = "_"+ particle.name.strip()
        for particle in mirrored_image.particles():
            particle.name = "_"+ particle.name.strip()
    
        #from foyer import Forcefield
        #from mbuild.formats.hoomdxml import write_hoomdxml
        #ff = Forcefield(forcefield_files=HOOMD_FF)
        #kwargs  = {}
        #kwargs['rigid_bodies'] = [p.rigid_id for p in mb_compound.particles()]
        #structure = mb_compound.to_parmed(box=original_box)
        #structure = ff.apply(structure, assert_dihedral_params=False)
        #write_hoomdxml(structure, '{}.hoomdxml'.format(options.output), 
        #        ref_energy = 0.239, ref_distance = 1, **kwargs)

        #ff = Forcefield(forcefield_files=HOOMD_FF)
        #kwargs  = {}
        #kwargs['rigid_bodies'] = [p.rigid_id for p in mirrored_image.particles()]
        #structure = mirrored_image.to_parmed(box=new_box)
        #structure = ff.apply(structure, assert_dihedral_params=False)
        #write_hoomdxml(structure, '{}_2x2x2.hoomdxml'.format(options.output), 
        #        ref_energy = 0.239, ref_distance = 1, **kwargs)
        

        foyer_kwargs = {"assert_angle_params": False,
                "assert_dihedral_params": False}
        mb_compound.save('{}.hoomdxml'.format(options.output), ref_energy = 0.239, ref_distance = 1, forcefield_files=HOOMD_FF, overwrite=True, box=original_box,
                foyerkwargs=foyer_kwargs)
        mirrored_image.save('{}_2x2x2.hoomdxml'.format(options.output), ref_energy = 0.239, ref_distance = 1, forcefield_files=HOOMD_FF, overwrite=True, box=new_box,
                foyerkwargs=foyer_kwargs)
    
        
    end=time.time()
    print(end-start)
    
    
