# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
from pathlib import Path
import pandas as pd
import sys
import json

import xtrack as xt
import xpart as xp
import xcoll as xc
import xobjects as xo

# This script prepares SixTrack data to be compared to K2 coupling results. 
# Not all dumps have data in collgaps. Had to remove tcsg.b4l7.b1 tcsg.b4r7.b2 tcsg.e5l7.b2 tcsg.e5r7.b1 tcspm.d4r7.b2 
coll_array_b1 = np.array([ 'tcdqa.a4r6.b1', 'tcdqa.b4r6.b1', 'tcdqa.c4r6.b1','tcl.5r1.b1','tcl.5r5.b1', 'tcl.6r1.b1', 'tcl.6r5.b1', 'tcla.6r3.b1', 'tcla.7r3.b1',
                          'tcla.a5r3.b1', 'tcla.a6r7.b1', 'tcla.a7r7.b1','tcla.b5r3.b1','tcla.b6r7.b1','tcla.c6r7.b1','tcla.d6r7.b1','tclpx.4r1.b1','tclpx.4r5.b1','tcp.6l3.b1','tcp.b6l7.b1',
                          'tcp.c6l7.b1','tcp.d6l7.b1','tcsg.4r3.b1','tcsg.5l3.b1','tcsg.a4l7.b1','tcsg.a4r7.b1','tcsg.a5l7.b1','tcsg.a5r3.b1','tcsg.a6l7.b1','tcsg.b5l7.b1','tcsg.b5r3.b1',
                          'tcsg.b5r7.b1','tcsg.d4l7.b1','tcsg.d5r7.b1','tcsp.a4r6.b1','tcspm.6r7.b1','tcspm.b4l7.b1','tcspm.e5r7.b1','tctph.4l2.b1','tctph.4l8.b1','tctph.6l1.b1','tctph.6l5.b1',
                          'tctpv.4l2.b1','tctpv.4l8.b1','tctpv.6l1.b1','tctpv.6l5.b1','tctpxh.4l1.b1', 'tctpxh.4l5.b1', 'tctpxv.4l1.b1', 'tctpxv.4l5.b1'])
                        
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

path_before = Path("SixTrackComp/before")
path_after  = Path("SixTrackComp/after")
path_jaw    = Path("SixTrackComp/jaws")

for beam in [4,1]:
    coll_array = coll_array_b1 if beam == 1 else coll_array_b2
    colldb_path = Path(f"/home/ssolstra/Documents/my_scripts/SixTrackComp/CollDB_HL_relaxed_b{beam}.data")
    colldb = xc.CollimatorDatabase.from_SixTrack(file=colldb_path,nemitt_x=3.5e-6, nemitt_y=3.5e-6)
    line = xt.Line.from_json(f"/home/ssolstra/Documents/scripts/generate_sixtrack/job_sample_thin_b{beam}.json")
    colldb._install_k2_collimators(line=line, verbose=True)

    # Jaws 
    file_coll = Path(f"/eos/home-b/bjlindst/share/forSimone/b{beam}_250muradXing15cmBetaNoTCLD_default+OFFSETAPER/runSingle/collgaps.dat")
    collgaps = pd.read_csv(file_coll, delim_whitespace=True)
    file_sigma = Path(f"/eos/home-b/bjlindst/share/forSimone/b{beam}_250muradXing15cmBetaNoTCLD_default+OFFSETAPER/runSingle/sigmasettings.out")
    sigmasettings = pd.read_csv(file_sigma, delim_whitespace=True)

    sigmasettings.columns = ['collimator', 'gap_h1', 'gap_h2', 'gap_h3', 'gap_h4', 'sig_offset',
        'coll_offset', 'nsig', 'gap_rms_error', 'beta_x', 'beta_y', 'orb_x',
        'orb_y', "#"]
    sigmasettings = sigmasettings.drop(columns=['#'])
    collgaps.columns = ['ID', 'name', 'angle[rad]', 'betax[m]', 'betay[m]', 'halfgap[m]',
       'mat.', 'length[m]', 'sigx[m]', 'sigy[m]', 'tilt1[rad]', 'tilt2[rad]',
       'nsig', "#"]
    collgaps = collgaps.drop(columns=['#'])
    jaw_values = []

    for name in coll_array:
        halfgap = collgaps[collgaps['name'] == name]['halfgap[m]'].values[0]                     
        angle   = collgaps[collgaps['name'] == name]['angle[rad]'].values[0]                     
        x       = sigmasettings[sigmasettings['collimator'] == name]['orb_x'].values[0] / 1e4  
        y       = sigmasettings[sigmasettings['collimator'] == name]['orb_y'].values[0] / 1e4   

        jaw_L = x * np.cos(angle) + y * np.sin(angle) + halfgap                                  
        jaw_R = x * np.cos(angle) + y * np.sin(angle) - halfgap                                  
        jaw_values.append([jaw_L, jaw_R])
        print(f"Calculated jaw values for {name}.")

    # with open(path_jaw/f"jaws_b{beam}.json", "w") as f:
    #     json.dump(jaw_values, f, cls=xo.JEncoder)

    # PARTICLES ########################################################
    directory = Path(f"/eos/home-b/bjlindst/share/forSimone/b{beam}_250muradXing15cmBetaNoTCLD_default+OFFSETAPER/runSingle")

    for part_file in directory.glob("dump_m*.dat"):
        if not part_file.exists():
            print(f"Error: The file {part_file} does not exist.")
            sys.exit(1)
        coll = part_file.name.split('_')[2].rsplit('.', 1)[0]
        duplicate = False
        dup = 2 if beam == 4 else 1
        if coll in coll_dup[dup-1]:
            duplicate = True
        # Read the file into a DataFrame
        df = pd.read_csv(part_file, delim_whitespace=True, header=0)
        df.columns = ['particleID', 'turn', 's[m]', 'x[mm]', 'xp[mrad]', 'y[mm]',
       'yp[mrad]', '(E-E0)/E0[1]', 'ktrack',"#"]
        df.drop(columns=["#"], inplace=True)
        
        # then split the turns
        for i in range(2,6):
            print(f"Getting particles for {coll} at turn {i}...")
            part_ref = xp.reference_from_pdg_id('proton', energy=7e12) # alternatively energy=7e12 p0c=7000e9

            # Extract particle data
            part_x  = df[df['turn']  == i]["x[mm]"].values  / 1000                        # [m] 
            part_y  = df[df['turn']  == i]["y[mm]"].values  / 1000                        # [m]
            ptau    = df[df['turn'] == i]["(E-E0)/E0[1]"].values
            delta   = np.sqrt(ptau*ptau + 2*ptau/part_ref.beta0+1) - 1
            part_px = (df[df['turn'] == i]["xp[mrad]"].values / 1000) * (1 + delta)             # [rad]
            part_py = (df[df['turn'] == i]["yp[mrad]"].values / 1000) * (1 + delta)             # [rad]
            part_id = df[df['turn'] == i]["particleID"].values
    
            # Create particle object
            part = xt.Particles(x=part_x, px=part_px, y=part_y, py=part_py, delta=delta, particle_id=part_id, particle_ref=part_ref) # bytta ut ptau med delta
            drift_length = line[coll].length/2.                                          # [m]

            if duplicate and part_file.name.split('_')[1][1] == 'b':    
                part_x     = part.x[:int(len(part.x)/2)]
                part_y     = part.y[:int(len(part.y)/2)]
                part_px    = part.px[:int(len(part.px)/2)]
                part_py    = part.py[:int(len(part.py)/2)]
                part_d     = part.delta[:int(len(part.delta)/2)]
                part_id    = part.particle_id[:int(len(part.particle_id)/2)]
                del part 
                part = xt.Particles(x=part_x, px=part_px, y=part_y, py=part_py, delta=part_d, particle_id=part_id, particle_ref=part_ref)

            mask_alive = (abs(part.x) > 0.) & (abs(part.px) > 0) & (abs(part.y) > 0) & (abs(part.py) > 0) 
            part.state[~mask_alive] = -400
            # Save and drift particle object change back to json and part object
            if part_file.name.split('_')[1][1] == 'b':
                drift = xt.Drift(length= -drift_length)
                drift.track(part)
                # with open(Path(path_before/ f"part_{coll}_turn_{i}_before.json"), 'w') as fid:
                #     json.dump(part.to_dict(), fid, cls=xo.JEncoder)
                np.savez(path_before/f"part_{coll}_turn_{i}_before.npz", x=part.x, px=part.px, y=part.y, py=part.py, delta=part.delta, particle_id=part.particle_id, state=part.state, s=part.s)
                print(f"    Saved particles for {coll} at turn {i}!")
            else:
                drift = xt.Drift(length=drift_length)
                drift.track(part)
                # with open(Path(path_after/ f"part_{coll}_turn_{i}_after.json"), 'w') as fid:
                #     json.dump(part.to_dict(), fid, cls=xo.JEncoder)
                print(f"    Saved particles for {coll} at turn {i}!")
                np.savez(path_after/f"part_{coll}_turn_{i}_after.npz", x=part.x, px=part.px, y=part.y, py=part.py, delta=part.delta, particle_id=part.particle_id,state=part.state, s=part.s)
        print()
        print(f"Finished processing particles for {coll}.")
        print()
    print(f"Finished processing beam {beam}")
print("All done!")
