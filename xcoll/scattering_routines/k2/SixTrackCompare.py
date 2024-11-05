import numpy as np
from pathlib import Path
import pandas as pd
import sys
import json

import xobjects as xo
import xtrack as xt
import xcoll as xc
import xpart as xp

from xcoll.scattering_routines.k2 import K2Engine
from xcoll.beam_elements.k2 import _K2Collimator


num_turns = 5
# num_particles = 100 000
path = Path("/home/ssolstra/Documents/my_scripts/SixTrackComp/")

# Not all dumps have data in collgaps. Had to remove tcsg.b4l7.b1 tcsg.b4r7.b2 tcsg.e5l7.b2 tcsg.e5r7.b1 tcspm.d4r7.b2 
coll_array_b1 = np.array([ 'tcdqa.a4r6.b1', 'tcdqa.b4r6.b1', 'tcdqa.c4r6.b1','tcl.5r1.b1','tcl.5r5.b1', 'tcl.6r1.b1', 'tcl.6r5.b1', 'tcla.6r3.b1', 'tcla.7r3.b1',
                          'tcla.a5r3.b1', 'tcla.a6r7.b1', 'tcla.a7r7.b1','tcla.b5r3.b1','tcla.b6r7.b1','tcla.c6r7.b1','tcla.d6r7.b1','tclpx.4r1.b1','tclpx.4r5.b1','tcp.6l3.b1','tcp.b6l7.b1',
                          'tcp.c6l7.b1','tcp.d6l7.b1','tcsg.4r3.b1','tcsg.5l3.b1','tcsg.a4l7.b1','tcsg.a4r7.b1','tcsg.a5l7.b1','tcsg.a5r3.b1','tcsg.a6l7.b1','tcsg.b5l7.b1','tcsg.b5r3.b1',
                          'tcsg.b5r7.b1','tcsg.d4l7.b1','tcsg.d5r7.b1','tcsp.a4r6.b1','tcspm.6r7.b1','tcspm.b4l7.b1','tcspm.e5r7.b1','tctph.4l2.b1','tctph.4l8.b1','tctph.6l1.b1','tctph.6l5.b1',
                          'tctpv.4l2.b1','tctpv.4l8.b1','tctpv.6l1.b1','tctpv.6l5.b1','tctpxh.4l1.b1', 'tctpxh.4l5.b1', 'tctpxv.4l5.b1'])
                        
coll_array_b2 = np.array([ 'tcdqa.a4l6.b2', 'tcdqa.b4l6.b2', 'tcdqa.c4l6.b2','tcl.5l1.b2','tcl.5l5.b2', 'tcl.6l1.b2', 'tcl.6l5.b2', 'tcla.6l3.b2', 'tcla.7l3.b2',
                          'tcla.a5l3.b2', 'tcla.a6l7.b2', 'tcla.a7l7.b2','tcla.b5l3.b2','tcla.b6l7.b2','tcla.c6l7.b2','tcla.d6l7.b2','tclpx.4l1.b2','tclpx.4l5.b2',
                          'tcp.6r3.b2','tcp.b6r7.b2','tcp.c6r7.b2','tcp.d6r7.b2','tcsg.4l3.b2','tcsg.5r3.b2','tcsg.a4l7.b2','tcsg.a4r7.b2','tcsg.a5l3.b2','tcsg.a5r7.b2',
                          'tcsg.a6r7.b2','tcsg.b5l3.b2','tcsg.b5l7.b2','tcsg.b5r7.b2','tcsg.d5l7.b2','tcsp.a4l6.b2',
                          'tcspm.6l7.b2','tcspm.b4r7.b2','tcspm.e5l7.b2','tctph.4r2.b2','tctph.4r8.b2','tctph.6r1.b2','tctph.6r5.b2',
                          'tctpv.4r2.b2','tctpv.4r8.b2','tctpv.6r1.b2','tctpv.6r5.b2','tctpxh.4r1.b2', 'tctpxh.4r5.b2', 'tctpxv.4r1.b2', 'tctpxv.4r5.b2'])

coll_dup = [['tcdqa.a4r6.b1', 'tcdqa.c4r6.b1', 'tcl.5r1.b1', 'tcla.7r3.b1', 'tcla.a7r7.b1', 'tclpx.4r1.b1', 
             'tcsg.5l3.b1', 'tcsg.a4l7.b1', 'tcsg.a5r3.b1', 'tcsg.b5r3.b1', 'tcsg.b5r7.b1', 'tcsg.d5r7.b1', 
             'tcsp.a4r6.b1', 'tcspm.e5r7.b1', 'tctph.4l8.b1', 'tctph.6l1.b1', 'tctpv.4l2.b1', 'tctpv.6l1.b1', 
             'tctpv.6l5.b1'], 
            ['tcdqa.a4l6.b2', 'tcl.6l5.b2', 'tcla.6l3.b2', 'tcla.7l3.b2', 'tcla.a6l7.b2', 'tcla.a7l7.b2', 
             'tcp.6r3.b2', 'tcp.c6r7.b2', 'tcp.d6r7.b2', 'tcsg.4l3.b2', 'tcsg.a5r7.b2', 'tcsg.a6r7.b2', 
             'tcsg.b5l7.b2', 'tcsg.b5r7.b2', 'tcsg.d5l7.b2', 'tcspm.e5l7.b2', 'tctph.4r2.b2', 'tctph.6r5.b2', 
             'tctpv.4r8.b2', 'tctpxh.4r5.b2', 'tctpxv.4r1.b2', 'tctpxv.4r5.b2']]

for beam in [1,4]:
    coll_array = coll_array_b1 if beam == 1 else coll_array_b2

    # get line
    line = xt.Line.from_json(f"/home/ssolstra/Documents/scripts/generate_sixtrack/job_sample_thin_b{beam}.json")
    file_coll = Path(f"/eos/home-b/bjlindst/share/forSimone/b{beam}_250muradXing15cmBetaNoTCLD_default+OFFSETAPER/runSingle/collgaps.dat")
    collgaps = pd.read_csv(file_coll, delim_whitespace=True)
    collgaps.columns = ['ID', 'name', 'angle[rad]', 'betax[m]', 'betay[m]', 'halfgap[m]',
       'mat.', 'length[m]', 'sigx[m]', 'sigy[m]', 'tilt1[rad]', 'tilt2[rad]',
       'nsig', "#"]
    collgaps = collgaps.drop(columns=['#'])

    # Open colldb 
    colldb_path = Path(path / f"CollDB_HL_relaxed_b{beam}.data")
    colldb = xc.CollimatorDatabase.from_SixTrack(file=colldb_path,nemitt_x=3.5e-6, nemitt_y=3.5e-6)
    colldb._install_k2_collimators(line=line, verbose=True)
    tw = line.twiss()

    # Aperture model check
    print('\nAperture model check after introducing collimators:')
    df_with_coll = line.check_aperture()
    assert not np.any(df_with_coll.has_aperture_problem)

    xc.assign_optics_to_collimators(line=line, twiss=tw)

    # Add values to collimator jaws
    with open(path / 'jaws' /f"jaws_b{beam}.json", "r") as f:
        jaw_values = json.load(f)
    for idx, name in enumerate(coll_array):
        # line[name].jaw_L = collgaps[collgaps['name'] == name]['halfgap[m]'].values[0]

        line[name].jaw_L = jaw_values[idx][0]
        line[name].jaw_R = jaw_values[idx][1]
    print("Values added to jaws")
    
    # Build the tracker
    line.build_tracker()
    K2Engine.start(line=line, cwd='run_1', _capacity = 200000, seed=1336)
    part_ref = xp.reference_from_pdg_id('proton', p0c=7000e9) # alternatively energy=7e12
    for name in coll_array:
        drift_length = line[name].length/2.   
        for i in range(2,6):
            print(f"Getting particles for {name} and turn {i}...") 
            part_data = np.load(path/ f"before/part_{name}_turn_{i}_before.npz")

            part = xp.Particles(x=part_data['x'], px=part_data['px'], y=part_data['y'], py=part_data['py'], delta=part_data['delta'], particle_id=part_data['particle_id'], state=part_data['state'],s=part_data['s'])
            part.particle_ref = part_ref
            # with open(Path(path / f"before/part_{name}_turn_{i}_before.json"), 'r') as fid:
            #     part = xp.Particles.from_dict(json.load(fid))

            xc.enable_scattering(line)
            line[name].track(part)
            xc.disable_scattering(line)

            part.sort(interleave_lost_particles=True)
            mask_alive = (abs(part.x) > 0.) & (abs(part.px) > 0) & (abs(part.y) > 0) & (abs(part.py) > 0) 
            part.state[~mask_alive] = -400
            drift = xt.Drift(length=drift_length)
            drift.track(part)
            np.savez(path/f"after_K2/part_{name}_turn_{i}_after.npz", x=part.x, px=part.px, y=part.y, py=part.py, delta=part.delta, particle_id=part.particle_id, state=part.state,s=part.s)
            # with open(Path(path/ f"after_K2/part_{name}_turn_{i}_after.json"), 'w') as fid:
            #     json.dump(part.to_dict(), fid, cls=xo.JEncoder)

            # delete part for memory
            del part
        print(f"Finished with {name}.")
    # Stop engine
    print("ALL DONE!")
    K2Engine.stop()


















