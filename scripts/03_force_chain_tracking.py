#!/usr/bin/env yade
# -*- coding: utf-8 -*-


# ============================================================================
# Force-chain tracking from saved YADE DEM states
# ============================================================================
# Purpose
# -------
# This script performs force-chain detection and tracking from saved YADE state
# files (.yade). It is designed as a post-processing tool: it does not run a new
# mechanical simulation, but reloads stored states and analyzes the particle
# stress network at each saved step.
#
# Typical use
# -----------
# In this repository, the script is used to analyze states saved during cyclic
# fluid injection after direct shear. However, the same workflow can be applied
# to any saved YADE states produced by other simulations, provided that the state
# files contain the information required for stress and contact-network analysis:
# particles, contacts/interactions, body stress tensors, particle positions, and
# the relevant boundary bodies.
#
# Main stages
# -----------
# 1. Load each .yade state file in the working directory.
# 2. Compute per-particle stress tensors and principal stresses.
# 3. Identify highly stressed grains using |sigma3| > mean(|sigma3|).
# 4. Select force-chain candidates based on the contact topology.
# 5. Build force chains using an alpha-angle alignment criterion.
# 6. Export one tracking text file and one VTU visualization file per state.
#
# Outputs
# -------
# - <state>_forceChainTracking_iterXXXX.txt
#   Compact text file containing stress values, force-chain flags, chain IDs,
#   and chain statistics needed for tracking, plotting, and counting.
#
# - <state>.00000001.vtu
#   VTU file for visualization in ParaView, with detected chain particles colored
#   in red and other particles in blue.
#
# Notes
# -----
# - Units follow the original YADE simulation state.
# - Wall-attached particles can be excluded from the force-chain analysis.
# ============================================================================

from __future__ import print_function

import os
import glob
import math
import numpy as np

from yade import utils, wrapper

Sphere = wrapper.Sphere
Box = wrapper.Box
VTKRecorder = wrapper.VTKRecorder


# =============================================================================
# Force-chain detection algorithm
# =============================================================================

class ForceChainDetector(object):
    """
    Detect force chains from saved YADE states.

    The workflow is:
    1. Compute per-particle stress tensors using utils.bodyStressTensors().
    2. Extract principal stresses and the direction associated with sigma3.
    3. Select highly stressed particles using |sigma3| > mean(|sigma3|).
    4. Keep chain candidates with at least two highly stressed neighbors.
    5. Build chains using an alpha-angle alignment criterion.
    """

    def __init__(self, exclude_ids=None, alpha_deg=45.0):
        self.exclude_ids = set(exclude_ids) if exclude_ids is not None else set()
        self.alpha_deg = float(alpha_deg)
        self.alpha_rad = math.radians(self.alpha_deg)
        self.cos_alpha = math.cos(self.alpha_rad)

    def compute_stresses(self):
        body_tensors = utils.bodyStressTensors()  # index corresponds to body id

        particle_stresses = {}
        sigma3_abs_list = []

        for b in O.bodies:
            if not isinstance(b.shape, Sphere):
                continue
            if b.id in self.exclude_ids:
                continue

            m = body_tensors[b.id]
            sigma_ij = np.array([
                [m[0, 0], m[0, 1], m[0, 2]],
                [m[1, 0], m[1, 1], m[1, 2]],
                [m[2, 0], m[2, 1], m[2, 2]],
            ], dtype=float)

            sigma_ij = 0.5 * (sigma_ij + sigma_ij.T)

            eigvals, eigvecs = np.linalg.eigh(sigma_ij)

            order = np.argsort(eigvals)[::-1]  # sigma1 >= sigma2 >= sigma3
            eigvals = eigvals[order]
            eigvecs = eigvecs[:, order]

            sigma1, sigma2, sigma3 = eigvals

            v3 = eigvecs[:, 2]
            nv = np.linalg.norm(v3)
            if nv > 0.0:
                v3 = v3 / nv
            else:
                v3 = np.array([0.0, 0.0, 0.0], dtype=float)

            v3_xy = np.array([v3[0], v3[1], 0.0])
            norm_xy = np.linalg.norm(v3_xy)
            if norm_xy > 0.0:
                v3_xy /= norm_xy
                theta = math.degrees(math.atan2(v3_xy[0], v3_xy[1]))
            else:
                theta = float("nan")

            pos = b.state.pos

            particle_stresses[b.id] = {
                "sigma1": float(sigma1),
                "sigma2": float(sigma2),
                "sigma3": float(sigma3),
                "theta_deg": float(theta),
                "dir3": [float(v3[0]), float(v3[1]), float(v3[2])],
                "pos": [float(pos[0]), float(pos[1]), float(pos[2])],
            }
            sigma3_abs_list.append(abs(float(sigma3)))

        if len(particle_stresses) == 0:
            print("\n[ForceChainDetector] WARNING: no sphere stresses found.")
            mean_abs_sigma3 = 0.0
        else:
            mean_abs_sigma3 = float(np.mean(sigma3_abs_list))

        print("\n=== PER-PARTICLE STRESS ANALYSIS ===")
        print("  Number of spheres analysed : %d" % len(particle_stresses))
        if sigma3_abs_list:
            print("  |sigma3| range : %.3e .. %.3e" % (min(sigma3_abs_list), max(sigma3_abs_list)))
            print("  Mean |sigma3|  : %.3e" % mean_abs_sigma3)

        return particle_stresses, mean_abs_sigma3

    def Chains_Candidate(self, particle_stresses, mean_abs_sigma3):
        if not particle_stresses:
            print("\n[ForceChainDetector] WARNING: no stress data available.")
            return {}, set(), set()

        neighbors = {}
        for pid in particle_stresses.keys():
            neighbors[pid] = set()

        for intr in O.interactions:
            if not intr.isReal:
                continue
            id1, id2 = intr.id1, intr.id2
            if (id1 in neighbors) and (id2 in neighbors):
                neighbors[id1].add(id2)
                neighbors[id2].add(id1)

        highly_stressed = set()
        for pid, vals in particle_stresses.items():
            if abs(vals["sigma3"]) > mean_abs_sigma3:
                highly_stressed.add(pid)

        print("\n=== HIGHLY STRESSED PARTICLES ===")
        print("  Threshold |sigma3| > mean(|sigma3|) = %.3e" % mean_abs_sigma3)
        print("  # highly stressed grains : %d" % len(highly_stressed))

        chain_candidates = set()
        for pid in highly_stressed:
            n_high = 0
            for nb in neighbors[pid]:
                if nb in highly_stressed:
                    n_high += 1
            if n_high >= 2:
                chain_candidates.add(pid)

        print("\n=== CHAIN CANDIDATE FILTER ===")
        print("  Criterion: particle must have at least two highly stressed neighbors.")
        print("  # chain candidate grains : %d" % len(chain_candidates))

        return neighbors, highly_stressed, chain_candidates

    def Build_Chains(self, particle_stresses, neighbors, highly_stressed, chain_candidates):
        chains = []
        if not chain_candidates:
            print("\n[ForceChainDetector] No chain candidates found.")
            return chains

        cos_alpha = self.cos_alpha
        assigned = set()

        def get_dir3(pid):
            v = np.array(particle_stresses[pid]["dir3"], dtype=float)
            n = np.linalg.norm(v)
            if n == 0.0:
                return None
            return v / n

        def get_pos(pid):
            return np.array(particle_stresses[pid]["pos"], dtype=float)

        def next_in_direction(current_id, direction_sign, local_used):
            u = get_dir3(current_id)
            if u is None:
                return None
            u = direction_sign * u

            pA = get_pos(current_id)
            best_nb = None
            best_cos = cos_alpha

            for nb in neighbors[current_id]:
                if nb not in highly_stressed:
                    continue
                if (nb in assigned) or (nb in local_used):
                    continue

                pB = get_pos(nb)
                l = pB - pA
                norm_l = np.linalg.norm(l)
                if norm_l == 0.0:
                    continue

                cos1 = np.dot(l, u) / (norm_l * np.linalg.norm(u))
                if cos1 < cos_alpha:
                    continue

                vB = get_dir3(nb)
                if vB is None:
                    continue
                l_back = pA - pB
                norm_b = np.linalg.norm(l_back)
                if norm_b == 0.0:
                    continue
                cos2 = np.dot(l_back, vB) / (norm_b * np.linalg.norm(vB))
                if cos2 < cos_alpha:
                    continue

                if cos1 > best_cos:
                    best_cos = cos1
                    best_nb = nb

            return best_nb

        def trace_side(seed_id, direction_sign, local_used):
            side = []
            current = seed_id
            while True:
                nxt = next_in_direction(current, direction_sign, local_used)
                if nxt is None:
                    break
                side.append(nxt)
                local_used.add(nxt)
                current = nxt
            return side

        print("\n=== BUILDING FORCE CHAINS (alpha = %.1f deg) ===" % self.alpha_deg)

        for seed in chain_candidates:
            if seed in assigned:
                continue

            local_used = set([seed])

            forward = trace_side(seed, +1, local_used)
            backward = trace_side(seed, -1, local_used)

            chain_ids = list(reversed(backward)) + [seed] + forward

            if len(chain_ids) >= 3:
                chains.append(chain_ids)
                for pid in chain_ids:
                    assigned.add(pid)

        print("\n=== FORCE CHAIN SUMMARY ===")
        print("  # chains found : %d" % len(chains))
        if chains:
            lengths = [len(c) for c in chains]
            print("  chain length range : %d .. %d" % (min(lengths), max(lengths)))
            print("  mean chain length  : %.2f" % float(np.mean(lengths)))

        return chains


# =============================================================================
# User parameters
# =============================================================================

ALPHA_DEG = 45.0
BASE_COLOR = (0.0, 0.0, 1.0)
CHAIN_COLOR = (1.0, 0.0, 0.0)
EXCLUDE_WALL_ATTACHED = True


# =============================================================================
# Utility functions
# =============================================================================

def detect_wall_attached_spheres():
    """Return sphere ids directly attached to wall bodies."""
    wall_like = set()
    for b in O.bodies:
        if isinstance(b.shape, Box):
            wall_like.add(b.id)

    attached = set()
    for intr in O.interactions:
        if not intr.isReal:
            continue
        i1, i2 = intr.id1, intr.id2
        if (i1 in wall_like) and isinstance(O.bodies[i2].shape, Sphere):
            attached.add(i2)
        elif (i2 in wall_like) and isinstance(O.bodies[i1].shape, Sphere):
            attached.add(i1)

    return attached


def reset_sphere_colors(color_rgb):
    """Assign one color to all spheres."""
    for b in O.bodies:
        if isinstance(b.shape, Sphere):
            b.shape.color = color_rgb


def color_chains_one_color(chains, color_rgb):
    """Color all detected force-chain particles with one color."""
    if not chains:
        return

    for ch in chains:
        for pid in ch:
            if pid < 0 or pid >= len(O.bodies):
                continue
            b = O.bodies[pid]
            if isinstance(b.shape, Sphere):
                b.shape.color = color_rgb


def write_force_chain_tracking_file(prefix, state_name, particle_stresses, mean_abs_sigma3,
                                    highly_stressed, chain_candidates, chains):
    """
    Write one compact text file containing force-chain tracking information.

    The file contains:
    - one summary header for the processed state;
    - one row per analysed particle;
    - particle stress state and chain-membership flags;
    - force-chain id for particles belonging to a detected chain.
    """
    it = int(O.iter)
    output_file = prefix + "_forceChainTracking_iter%d.txt" % it

    chain_id_by_particle = {}
    chain_length_by_particle = {}

    for chain_id, ch in enumerate(chains, start=1):
        chain_length = len(ch)
        for pid in ch:
            chain_id_by_particle[pid] = chain_id
            chain_length_by_particle[pid] = chain_length

    lengths = [len(c) for c in chains]
    chain_len_min = min(lengths) if lengths else 0
    chain_len_max = max(lengths) if lengths else 0
    chain_len_mean = float(np.mean(lengths)) if lengths else 0.0

    with open(output_file, "w") as f:
        f.write("# Force-chain tracking output\n")
        f.write("# state_file\t%s\n" % state_name)
        f.write("# iter\t%d\n" % it)
        f.write("# alpha_deg\t%.6f\n" % ALPHA_DEG)
        f.write("# mean_abs_sigma3\t%.6e\n" % mean_abs_sigma3)
        f.write("# num_spheres_with_stress\t%d\n" % len(particle_stresses))
        f.write("# num_highly_stressed\t%d\n" % len(highly_stressed))
        f.write("# num_chain_candidates\t%d\n" % len(chain_candidates))
        f.write("# num_chains\t%d\n" % len(chains))
        f.write("# chain_len_min\t%d\n" % chain_len_min)
        f.write("# chain_len_max\t%d\n" % chain_len_max)
        f.write("# chain_len_mean\t%.6f\n" % chain_len_mean)
        f.write(
            "# columns: "
            "particle_id\tx\ty\tz\tsigma1\tsigma2\tsigma3\tabsSigma3\t"
            "theta_deg\tdir3_x\tdir3_y\tdir3_z\tmeanAbsSigma3\t"
            "isHighlyStressed\tisChainCandidate\tisInChain\tchainID\tchainLength\n"
        )

        for pid in sorted(particle_stresses.keys()):
            vals = particle_stresses[pid]
            x, y, z = vals["pos"]
            d0, d1, d2 = vals["dir3"]
            is_high = 1 if pid in highly_stressed else 0
            is_candidate = 1 if pid in chain_candidates else 0
            is_in_chain = 1 if pid in chain_id_by_particle else 0
            chain_id = chain_id_by_particle.get(pid, 0)
            chain_length = chain_length_by_particle.get(pid, 0)

            f.write(
                "%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%d\t%d\t%d\t%d\t%d\n" % (
                    pid,
                    x,
                    y,
                    z,
                    vals["sigma1"],
                    vals["sigma2"],
                    vals["sigma3"],
                    abs(vals["sigma3"]),
                    vals["theta_deg"],
                    d0,
                    d1,
                    d2,
                    mean_abs_sigma3,
                    is_high,
                    is_candidate,
                    is_in_chain,
                    chain_id,
                    chain_length,
                )
            )

    print("Tracking file written: %s" % output_file)


def export_vtu(prefix):
    """Export one VTU file for visualization in ParaView."""
    saveSolid = VTKRecorder(
        fileName=prefix + ".",
        iterPeriod=1,
        skipNondynamic=0,
        recorders=["spheres", "colors", "velocity", "bstresses", "intr"],
    )
    saveSolid.dead = True

    # Run only the VTK recorder to avoid modifying the loaded mechanical state.
    O.engines = [saveSolid]

    saveSolid.dead = False
    O.run(1, 1)
    saveSolid.dead = True

    print("VTU exported: %s00000001.vtu" % (prefix + "."))


# =============================================================================
# Main batch processing
# =============================================================================

work_dir = os.getcwd()
yade_files = sorted(glob.glob(os.path.join(work_dir, "*.yade")), key=os.path.getmtime)

if not yade_files:
    raise RuntimeError("No .yade files found in current directory: %s" % work_dir)

print("\nFound %d .yade files in %s" % (len(yade_files), work_dir))

for k, path in enumerate(yade_files):
    fn = os.path.basename(path)
    base = os.path.splitext(fn)[0]

    print("\n" + "=" * 80)
    print("Processing file %d/%d : %s" % (k + 1, len(yade_files), fn))
    print("=" * 80)

    O.load(path)

    # Use the state filename in the output prefix to avoid overwriting outputs
    # when several .yade files are processed in the same directory.
    out_prefix = "FV_cycle_" + base

    reset_sphere_colors(BASE_COLOR)

    exclude_ids = set()
    if EXCLUDE_WALL_ATTACHED:
        exclude_ids = detect_wall_attached_spheres()
    print("Excluded wall-attached spheres : %d" % len(exclude_ids))

    detector = ForceChainDetector(exclude_ids=exclude_ids, alpha_deg=ALPHA_DEG)
    particle_stresses, mean_abs_sigma3 = detector.compute_stresses()
    neighbors, highly_stressed, chain_candidates = detector.Chains_Candidate(
        particle_stresses,
        mean_abs_sigma3,
    )
    chains = detector.Build_Chains(
        particle_stresses,
        neighbors,
        highly_stressed,
        chain_candidates,
    )

    print("\nSummary:")
    print("  iter (state)               : %d" % int(O.iter))
    print("  # spheres with stress data : %d" % len(particle_stresses))
    print("  # highly stressed grains   : %d" % len(highly_stressed))
    print("  # chain candidate grains   : %d" % len(chain_candidates))
    print("  # force chains             : %d" % len(chains))
    print("  Mean |sigma3|              : %.3e" % mean_abs_sigma3)

    write_force_chain_tracking_file(
        out_prefix,
        fn,
        particle_stresses,
        mean_abs_sigma3,
        highly_stressed,
        chain_candidates,
        chains,
    )

    color_chains_one_color(chains, CHAIN_COLOR)
    export_vtu(out_prefix)

print("\nALL DONE.")
