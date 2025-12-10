import numpy as np
from scipy.sparse import dok_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
from itertools import combinations
from scipy.linalg import eigh

def mod(a, b):
    return (a + b) % b
def site_index(m, n, Lx, L2):
    return m * L2 + n
def get_neighbors(Lx, L2, PBC_x, PBC_y):
    neighbors = []
    for m in range(Lx):
        for n in range(L2):
            i = site_index(m, n, Lx, L2)
            print(i,m,n)

            # nearest neighbors x
            if PBC_x or m + 1 < Lx:
                m_x = mod(m+1, Lx) if PBC_x else m+1
                if m_x < Lx:
                    j = site_index(m_x, n, Lx, L2)
                    neighbors.append((i, j))

            # nearest neighbors y
            if PBC_y or n + 1 < L2:
                n_y = mod(n+1, L2) if PBC_y else n+1
                if n_y < L2:
                    j = site_index(m, n_y, Lx, L2)
                    neighbors.append((i, j))

            # next-nearest diagonal +1, +1
            if (PBC_x or m + 1 < Lx) and (PBC_y or n + 1 < L2):
                m_d = mod(m+1, Lx) if PBC_x else m+1
                n_d = mod(n+1, L2) if PBC_y else n+1
                if m_d < Lx and n_d < L2:
                    j = site_index(m_d, n_d, Lx, L2)
                    neighbors.append((i, j))

            # next-nearest diagonal -1, -1
            if (PBC_x or m - 1 >= 0) and (PBC_y or n - 1 >= 0):
                m_d = mod(m-1, Lx) if PBC_x else m-1
                n_d = mod(n-1, L2) if PBC_y else n-1
                if m_d >= 0 and n_d >= 0:
                    j = site_index(m_d, n_d, Lx, L2)
                    neighbors.append((i, j))
    return neighbors

def count_bits(x):
    return bin(x).count('1')

def fermionic_sign(state, pos):
    mask = (1 << pos) - 1
    return (-1)**count_bits(state & mask)

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

def build_basis(Lx, L2, N_up, N_down):
    N = Lx * L2
    basis = []

    up_states = list(combinations(range(N), N_up))
    down_states = list(combinations(range(N), N_down))

    for up_occ in up_states:
        up_bits = 0
        for pos in up_occ:
            up_bits |= (1 << pos)
        for down_occ in down_states:
            down_bits = 0
            for pos in down_occ:
                down_bits |= (1 << pos)
            basis.append( (up_bits, down_bits) )
    return basis

def state_index_map(basis):
    return {state: i for i, state in enumerate(basis)}

def is_next_nearest(i, j, Lx, L2):
    # Test if sites i and j are next-nearest neighbors along diagonals (+1,+1) or (-1,-1)
    m_i, n_i = i // L2, i % L2
    m_j, n_j = j // L2, j % L2
    return (abs(m_i - m_j) == 1 and abs(n_i - n_j) == 1)

def Hubbard_ED_Hamiltonian(t, t_prime, U, mu, Lx, L2, N_up, N_down, PBC_x=False, PBC_y=False):
    N = Lx * L2
    basis = build_basis(Lx, L2, N_up, N_down)
    dim = len(basis)
    print(f"Basis dimension: {dim}")
    H = dok_matrix((dim, dim), dtype=np.float64)
    neighbors = []
    # Separate nearest neighbors and next nearest neighbors for hopping
    for m in range(Lx):
        for n in range(L2):
            i = site_index(m, n, Lx, L2)
            print(i,m,n)

            # nearest neighbors x
            if PBC_x or m + 1 < Lx:
                m_x = mod(m+1, Lx) if PBC_x else m+1
                if m_x < Lx:
                    j = site_index(m_x, n, Lx, L2)
                    neighbors.append( (i,j,'nn') )

            # nearest neighbors y
            if PBC_y or n + 1 < L2:
                n_y = mod(n+1, L2) if PBC_y else n+1
                if n_y < L2:
                    j = site_index(m, n_y, Lx, L2)
                    neighbors.append( (i,j,'nn') )

            # next-nearest diagonal +1, +1
            if (PBC_x or m + 1 < Lx) and (PBC_y or n + 1 < L2):
                m_d = mod(m+1, Lx) if PBC_x else m+1
                n_d = mod(n+1, L2) if PBC_y else n+1
                if m_d < Lx and n_d < L2:
                    j = site_index(m_d, n_d, Lx, L2)
                    neighbors.append( (i,j,'nnn') )

            # next-nearest diagonal -1, -1
            if (PBC_x or m - 1 >= 0) and (PBC_y or n - 1 >= 0):
                m_d = mod(m-1, Lx) if PBC_x else m-1
                n_d = mod(n-1, L2) if PBC_y else n-1
                if m_d >= 0 and n_d >= 0:
                    j = site_index(m_d, n_d, Lx, L2)
                    neighbors.append( (i,j,'nnn') )

    basis_map = state_index_map(basis)

    for i, (up_state, down_state) in enumerate(basis):
        # chemical potential and Hubbard interaction(onsite term)
        for site in range(N):
            n_up = (up_state >> site) & 1
            n_down = (down_state >> site) & 1
            H[i, i] += -mu * (n_up + n_down) + U * n_up * n_down

        # Hopping terms for spin up electrons
        for (s, d, hop_type) in neighbors:
            new_up_state1, sign1 = annihilate_state(up_state, s)
            if new_up_state1 is None:
                continue
            new_up_state2, sign2 = create_state(new_up_state1, d)
            if new_up_state2 is None:
                continue
            new_state = (new_up_state2, down_state)
            j = basis_map.get(new_state, None)
            if j is not None:
                amp = -t if hop_type == 'nn' else -t_prime
                H[i, j] += amp * sign1 * sign2

        # Hopping terms for spin down electrons
        for (s, d, hop_type) in neighbors:
            new_down_state1, sign1 = annihilate_state(down_state, s)
            if new_down_state1 is None:
                continue
            new_down_state2, sign2 = create_state(new_down_state1, d)
            if new_down_state2 is None:
                continue
            new_state = (up_state, new_down_state2)
            j = basis_map.get(new_state, None)
            if j is not None:
                amp = -t if hop_type == 'nn' else -t_prime
                H[i, j] += amp * sign1 * sign2

    return H.tocsr(), basis

def double_occupancy(basis, vec, Lx, L2):
    N = Lx * L2
    doubl_occ = 0.0
    for coeff, (up_state, down_state) in zip(vec, basis):
        for site in range(N):
            n_up = (up_state >> site) & 1
            n_down = (down_state >> site) & 1
            doubl_occ += coeff.conj() * coeff * (n_up * n_down)
    return doubl_occ.real

def total_density(basis, vec, Lx, L2):
    N = Lx * L2
    dens = 0.0
    for coeff, (up_state, down_state) in zip(vec, basis):
        for site in range(N):
            n_up = (up_state >> site) & 1
            n_down = (down_state >> site) & 1
            dens += coeff.conj() * coeff * (n_up + n_down)
    return dens.real



def spin_spin_correlation(basis, vec, Lx, L2, site_i=None, site_j=None):
    N = Lx * L2
    szsz = np.zeros((N, N), dtype=np.float64) if site_i is None else 0.0

    for coeff, (up, down) in zip(vec, basis):
        for i in range(N):
            szi = ((up >> i) & 1) - ((down >> i) & 1)
            if site_i is not None and i != site_i:
                continue

            for j in range(N):
                if site_j is not None and j != site_j:
                    continue

                szj = ((up >> j) & 1) - ((down >> j) & 1)
                val = (coeff.conj() * coeff * szi * szj).real

                if site_i is not None:
                    szsz += val
                else:
                    szsz[i, j] += val
    return szsz

if __name__ == "__main__":
    Lx, L2 = 4, 2
    N_up = int(Lx/1.)
    N_down = int(L2/1.)
    t = 1
    t_prime = 0.
    U = 4.0
    mu = 0.0
    PBC_x = True
    PBC_y = True

    H, basis = Hubbard_ED_Hamiltonian(t, t_prime, U, mu, Lx, L2, N_up, N_down, PBC_x, PBC_y)
    print("Diagonalizing Hamiltonian sparse from")
    vals, vecs = eigsh(H, k=5, which='SA')
    print("Lowest 5 eigenvalues:", vals)
        
    with open(f'Energies.txt', 'w') as file:
        file.write(f"{vals}\n")
     
    #full digoanlizations    
    #eig1, vecc =  eigh(H.toarray()) 
    #print("Ground state energy",  eig1)

    with open(f'Densities.txt', 'w') as file:
        for i in range(5):
            print(f"Eigenstate {i} double occupancy: {double_occupancy(basis, vecs[:, i], Lx, L2)}")
            print(f"Eigenstate {i} total density: {total_density(basis, vecs[:, i], Lx, L2)}")
            file.write(f"Eigenstate {i} double occupancy: {double_occupancy(basis, vecs[:, i], Lx, L2)}\n")
            file.write(f"Eigenstate {i} total density: {total_density(basis, vecs[:, i], Lx, L2)}\n")
                    

