import numpy as np
from scipy.sparse import dok_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
from itertools import combinations
from scipy.linalg import eigh

def count_bits(x):
    return bin(x).count('1')

def fermionic_sign(state, pos):
    # sign for creating/annihilating at 'pos' with ordering 0..L-1
    mask = (1 << pos) - 1
    return (-1) ** count_bits(state & mask)

def create_state(state, pos):
    if (state >> pos) & 1:
        return None, 0
    sign = fermionic_sign(state, pos)
    new_state = state | (1 << pos)
    return new_state, sign

def annihilate_state(state, pos):
    if ((state >> pos) & 1) == 0:
        return None, 0
    sign = fermionic_sign(state, pos)
    new_state = state & ~(1 << pos)
    return new_state, sign

def build_basis(L, N_up, N_down):
    basis = []
    sites = range(L)
    up_states = list(combinations(sites, N_up))
    down_states = list(combinations(sites, N_down))
    for up_occ in up_states:
        up_bits = 0
        for pos in up_occ:
            up_bits |= (1 << pos)
        for down_occ in down_states:
            down_bits = 0
            for pos in down_occ:
                down_bits |= (1 << pos)
            basis.append((up_bits, down_bits))
    return basis

def state_index_map(basis):
    return {state: i for i, state in enumerate(basis)}

def get_neighbors_1d(L, PBC):
    # list of  nearest-neighbor pairs 
    pairs = []
    for i in range(L - 1):
        pairs.append((i, i+1))
        pairs.append((i+1, i))
    if PBC and L>1:
        pairs.append((L-1, 0))
        pairs.append((0, L-1))
    return pairs

def get_nnn_1d(L, PBC):
    pairs = []
    for i in range(L - 2):
        pairs.append((i, i+2))
        pairs.append((i+2, i))
    if PBC and L>2:
        pairs.append((L-2, 0))
        pairs.append((0, L-2))
        pairs.append((L-1, 1))
        pairs.append((1, L-1))
    return pairs

def Hubbard_ED_Hamiltonian_1D(t, t_prime, U, mu, L, N_up, N_down, PBC=False):
    basis = build_basis(L, N_up, N_down)
    dim = len(basis)
    print(f"Basis dimension: {dim}")
    H = dok_matrix((dim, dim), dtype=np.float64)
    basis_map = state_index_map(basis)

    nn = get_neighbors_1d(L, PBC)
    nnn = get_nnn_1d(L, PBC)

    # Build Hamiltonian
    for i, (up_state, down_state) in enumerate(basis):
        # On-site terms
        for site in range(L):
            n_up = (up_state >> site) & 1
            n_down = (down_state >> site) & 1
            H[i, i] += -mu * (n_up + n_down) + U * n_up * n_down

        # Hopping spin-up: sum_{ij} -t_{ij} c†_i c_j
        for (s, d) in nn:
            amp = -t
            # annihilate at j = s (we treat pair as (s->d) meaning destroy at s, create at d)
            new_up1, sign1 = annihilate_state(up_state, s)
            if new_up1 is None: 
                continue
            new_up2, sign2 = create_state(new_up1, d)
            if new_up2 is None:
                continue
            new_state = (new_up2, down_state)
            j = basis_map.get(new_state, None)
            if j is not None:
                H[i, j] += amp * sign1 * sign2

        for (s, d) in nnn:
            amp = -t_prime
            new_up1, sign1 = annihilate_state(up_state, s)
            if new_up1 is None: 
                continue
            new_up2, sign2 = create_state(new_up1, d)
            if new_up2 is None:
                continue
            new_state = (new_up2, down_state)
            j = basis_map.get(new_state, None)
            if j is not None:
                H[i, j] += amp * sign1 * sign2

        # Hopping spin-down
        for (s, d) in nn:
            amp = -t
            new_dn1, sign1 = annihilate_state(down_state, s)
            if new_dn1 is None:
                continue
            new_dn2, sign2 = create_state(new_dn1, d)
            if new_dn2 is None:
                continue
            new_state = (up_state, new_dn2)
            j = basis_map.get(new_state, None)
            if j is not None:
                H[i, j] += amp * sign1 * sign2

        for (s, d) in nnn:
            amp = -t_prime
            new_dn1, sign1 = annihilate_state(down_state, s)
            if new_dn1 is None:
                continue
            new_dn2, sign2 = create_state(new_dn1, d)
            if new_dn2 is None:
                continue
            new_state = (up_state, new_dn2)
            j = basis_map.get(new_state, None)
            if j is not None:
                H[i, j] += amp * sign1 * sign2

    return H.tocsr(), basis

def double_occupancy(basis, vec, L):
    doubl_occ = 0.0
    for coeff, (up_state, down_state) in zip(vec, basis):
        for site in range(L):
            n_up = (up_state >> site) & 1
            n_down = (down_state >> site) & 1
            doubl_occ += abs(coeff)**2 * (n_up * n_down)
    return doubl_occ.real

def total_density(basis, vec, L):
    dens = 0.0
    for coeff, (up_state, down_state) in zip(vec, basis):
        for site in range(L):
            n_up = (up_state >> site) & 1
            n_down = (down_state >> site) & 1
            dens += abs(coeff)**2 * (n_up + n_down)
    return dens.real

# compute free-fermion ground energy by filling single-particle levels
def free_fermion_energy(L, t, t_prime, N_particles, PBC):
    # construct single-particle Hamiltonian (L x L)
    H1 = np.zeros((L, L), dtype=float)
    for i in range(L):
        # nearest neighbor hopping
        j = i+1
        if j < L:
            H1[i, j] = -t
            H1[j, i] = -t
        elif PBC and j==L:
            H1[i, 0] = -t
            H1[0, i] = -t
    # next-nearest
    for i in range(L):
        j = i+2
        if j < L:
            H1[i, j] = -t_prime
            H1[j, i] = -t_prime
        else:
            if PBC:
                H1[i, j % L] = -t_prime
                H1[j % L, i] = -t_prime
    eigs = np.linalg.eigvalsh(H1)
    eigs.sort()
    # fill lowest N_particles (spinless)
    return eigs[:N_particles].sum(), eigs  # return sum and spectrum

def test():
    L = 4
    N_up = int(L/2)
    N_down = L - N_up
    t = 1.0
    t_prime = 0.
    U = 8.0  
    mu = 0.0

    for PBC in (True, False):
        print("\n==== Test:  PBC =", PBC, "====")
        H, basis = Hubbard_ED_Hamiltonian_1D(t, t_prime, U, mu, L, N_up, N_down, PBC)
        # ground state energy with ED
        vals, vecs = eigsh(H, k=1, which='SA')
        E_mb = vals[0]
        print("Many-body ED ground energy (U=0):", E_mb)
        
        #full diagonalization
        eig1, vecc =  eigh(H.toarray()) 
        print("Full energy spectrum",  eig1)

        # free-fermion results: single-particle fill for each spin
        E1sum, spectrum = free_fermion_energy(L, t, t_prime, N_up, PBC)
        E_free_total = 2 * E1sum  # two spin species
        print("Free fermion (filled) energy  :", E_free_total)

        print("Difference (ED - free):", E_mb - E_free_total)


test()

