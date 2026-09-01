import numpy as np
from itertools import product
from pyscf.pbc.lib import kpts_helper
import sys
sys.path.insert(0, "/Users/qclove00/block2-preview/build-k")


def _block2_k_labels_1d(cell, kpts, kmesh):
    kmesh = np.asarray(kmesh, dtype=int)
    axes = np.where(kmesh > 1)[0]
    if len(axes) != 1:
        raise ValueError("block2 built-in k_symmetry/k_mod is a single cyclic group; use a 1D kmesh here.")

    axis = int(axes[0])
    nk = int(np.prod(kmesh))
    scaled = cell.get_scaled_kpts(kpts)
    labels = np.rint(scaled[:, axis] * nk).astype(int) % nk

    kconserv = kpts_helper.get_kconserv(cell, kpts)
    check = (labels[:, None, None] - labels[None, :, None]
             + labels[None, None, :] - labels[kconserv]) % nk
    if np.any(check):
        raise RuntimeError("PySCF kconserv and scalar block2 k labels are inconsistent.")

    return labels.tolist(), kconserv


def krhf_to_block2_k_fcidump(
    mf,
    kmesh,
    filename="FCIDUMP.K",
    mo_slice=None,
    n_elec=None,
    twos=0,
    target_k=0,
    total_energy=True,
    tol=1e-10,
):
    """
    Convert PySCF PBC KRHF orbitals/integrals to a block2 K-symmetry FCIDUMP.

    total_energy=True:
        Builds the Born-von Karman supercell Hamiltonian.
        DMRG energy / nkpts is the per-cell energy.
    total_energy=False:
        Divides the whole Hamiltonian by nkpts, so DMRG energy is per-cell.
    """
    from block2.cpx import FCIDUMP
    from block2 import VectorInt, VectorUInt8

    cell = mf.cell
    kpts = np.asarray(mf.kpts)
    nk = len(kpts)

    if int(np.prod(kmesh)) != nk:
        raise ValueError(f"kmesh={kmesh} is inconsistent with nkpts={nk}")

    k_labels, kconserv = _block2_k_labels_1d(cell, kpts, kmesh)

    mo = [np.asarray(c) for c in mf.mo_coeff]
    if mo_slice is not None:
        mo = [c[:, mo_slice] for c in mo]

    nmo_k = [c.shape[1] for c in mo]
    offsets = np.zeros(nk + 1, dtype=int)
    offsets[1:] = np.cumsum(nmo_k)
    n_sites = int(offsets[-1])

    def ks(k):
        return slice(offsets[k], offsets[k + 1])

    hcore = np.asarray(mf.get_hcore(cell, kpts))
    h1e = np.zeros((n_sites, n_sites), dtype=np.complex128)
    g2e = np.zeros((n_sites, n_sites, n_sites, n_sites), dtype=np.complex128)

    for k in range(nk):
        h1e[ks(k), ks(k)] = mo[k].conj().T @ hcore[k] @ mo[k]

    for ki, kj, kk in product(range(nk), repeat=3):
        kl = int(kconserv[ki, kj, kk])
        eri = mf.with_df.ao2mo(
            (mo[ki], mo[kj], mo[kk], mo[kl]),
            (kpts[ki], kpts[kj], kpts[kk], kpts[kl]),
            compact=False,
        )
        eri = np.asarray(eri).reshape(nmo_k[ki], nmo_k[kj], nmo_k[kk], nmo_k[kl])
        g2e[ks(ki), ks(kj), ks(kk), ks(kl)] = eri / nk

    ecore = mf.energy_nuc() * nk

    if not total_energy:
        h1e /= nk
        g2e /= nk
        ecore /= nk

    if n_elec is None:
        n_elec = cell.tot_electrons(nk)

    fd = FCIDUMP()
    fd.initialize_su2(
        n_sites,
        int(n_elec),
        int(twos),
        1,  # C1 trivial irrep in MOLPRO convention
        float(np.real(ecore)),
        np.ascontiguousarray(h1e).ravel(),
        np.ascontiguousarray(g2e).ravel(),
    )

    fd.orb_sym = VectorUInt8([1] * n_sites)
    fd.k_sym = VectorInt([k_labels[k] for k in range(nk) for _ in range(nmo_k[k])])
    fd.k_mod = nk
    fd.k_isym = int(target_k) % nk

    err = fd.symmetrize(fd.k_sym, fd.k_mod)
    if abs(err) > tol:
        raise RuntimeError(f"K-symmetry error is too large: {err}")

    fd.write(filename)
    return fd

def run_block2_su2k_dmrg(fd, scratch="./node0", n_threads=4, memory=2e9):
    import os
    from block2 import (
        SU2K, Global, Threading, ThreadingTypes, DoubleFPCodec,
        init_memory, release_memory,
        VectorUInt8, VectorUBond, VectorDouble, PointGroup,
        Random, QCTypes, SeqTypes, OpNamesSet, OpNames,
        NoiseTypes, DecompositionTypes,
    )
    from block2.su2k import MPSInfo
    from block2.cpx.su2k import (
        HamiltonianQC, MPS, MPOQC, SimplifiedMPO, RuleQC,
        MovingEnvironment, DMRG,
    )

    os.makedirs(scratch, exist_ok=True)
    Random.rand_seed(1234)

    init_memory(isize=int(1e8), dsize=int(memory), save_dir=scratch)
    Global.threading = Threading(
        ThreadingTypes.OperatorBatchedGEMM | ThreadingTypes.Global,
        n_threads, n_threads, 1,
    )
    Global.threading.seq_type = SeqTypes.Nothing
    Global.frame.fp_codec = DoubleFPCodec(1e-16, 1024)
    Global.frame.use_main_stack = False
    # set_mkl_num_threads(1)

    try:
        n_sites = fd.n_sites
        vacuum = SU2K(0)

        pg = PointGroup.swap_c1(fd.isym)  # C1 前提。D2h 等なら swap_d2h に変える
        target_pg = SU2K.pg_combine(pg, fd.k_isym, fd.k_mod)
        target = SU2K(fd.n_elec, fd.twos, target_pg)

        orb_pg = VectorUInt8(map(PointGroup.swap_c1, fd.orb_sym))
        orb_sym = HamiltonianQC.combine_orb_sym(orb_pg, fd.k_sym, fd.k_mod)

        hamil = HamiltonianQC(vacuum, n_sites, orb_sym, fd)

        mpo = MPOQC(hamil, QCTypes.Conventional)
        mpo = SimplifiedMPO(
            mpo, RuleQC(), True, True,
            OpNamesSet((OpNames.R, OpNames.RD)),
        )

        bond_dims = [100] * 4 + [200] * 4
        noises = [1e-5] * 4 + [1e-6] * 2 + [0.0]

        mps_info = MPSInfo(n_sites, vacuum, target, hamil.basis)
        mps_info.tag = "KET"
        mps_info.set_bond_dimension(bond_dims[0])

        mps = MPS(n_sites, 0, 2)
        mps.initialize(mps_info)
        mps.random_canonicalize()
        mps.save_mutable()
        mps.deallocate()
        mps_info.save_mutable()
        mps_info.deallocate_mutable()

        me = MovingEnvironment(mpo, mps, mps, "DMRG")
        me.delayed_contraction = OpNamesSet.normal_ops()
        me.cached_contraction = True
        me.init_environments(True)

        dmrg = DMRG(me, VectorUBond(bond_dims), VectorDouble(noises))
        dmrg.noise_type = NoiseTypes.ReducedPerturbativeCollected
        dmrg.decomp_type = DecompositionTypes.DensityMatrix

        energy = dmrg.solve(len(bond_dims), True, 1e-8)
        return energy

    finally:
        release_memory()

from pyscf.pbc import gto, scf

cell = gto.M(
    a=[[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]],
    atom="He 0 0 0",
    basis="gth-szv",
    pseudo="gth-pade",
    verbose=4,
)

kmesh = [4, 1, 1]
kpts = cell.make_kpts(kmesh)

mf = scf.KRHF(cell, kpts=kpts).density_fit()
mf.kernel()

fd = krhf_to_block2_k_fcidump(mf, kmesh, filename="FCIDUMP.K", target_k=0)
e = run_block2_su2k_dmrg(fd)
print("DMRG energy =", e)
print("per-cell energy =", e / len(mf.kpts))