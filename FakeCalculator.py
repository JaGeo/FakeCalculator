import numpy as np

from ase.calculators.calculator import PropertyNotImplementedError
from ase.calculators.lj import LennardJones
from ase.neighborlist import NeighborList


# This is adapted code from the ase package (LennardJones class)




# This Potential is not intended to give correct energies, forces or stresses, it just delivers something to test the rest of the code.
# So, don't worry, if anything is wrong.
class FakeCalc(LennardJones):
    implemented_properties = ['energy', 'forces', 'stress']

    nolabel = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def get_potential_energy(self, atoms=None, force_consistent=False):

        natoms = len(self)
        sigma = 1.0
        epsilon = 1.0
        rc = None

        if rc is None:
            rc = 3 * sigma

        self.nl = NeighborList([rc / 2] * natoms, self_interaction=False)

        self.nl.update(self)

        positions = self.positions
        cell = self.cell

        e0 = 4 * epsilon * ((sigma / rc) ** 12 - (sigma / rc) ** 6)

        energy = 0.0
        forces = np.zeros((natoms, 3))
        stress = np.zeros((3, 3))

        for a1 in range(natoms):
            neighbors, offsets = self.nl.get_neighbors(a1)
            cells = np.dot(offsets, cell)
            d = positions[neighbors] + cells - positions[a1]
            r2 = (d ** 2).sum(1)
            c6 = (sigma ** 2 / r2) ** 3
            c6[r2 > rc ** 2] = 0.0
            energy -= e0 * (c6 != 0.0).sum()
            c12 = c6 ** 2
            energy += 4 * epsilon * (c12 - c6).sum()
            f = (24 * epsilon * (2 * c12 - c6) / r2)[:, np.newaxis] * d
            forces[a1] -= f.sum(axis=0)
            for a2, f2 in zip(neighbors, f):
                forces[a2] += f2
            stress += np.dot(f.T, d)

        if self.number_of_lattice_vectors == 3:
            stress += stress.T.copy()
            stress *= -0.5 / self.get_volume()
            stress = stress.flat[[0, 4, 8, 5, 2, 1]]
        else:
            raise PropertyNotImplementedError

        return energy

    def get_stress(self, atoms=None):

        natoms = len(self)
        sigma = 1.0
        epsilon = 1.0
        rc = None

        if rc is None:
            rc = 3 * sigma

        self.nl = NeighborList([rc / 2] * natoms, self_interaction=False)

        self.nl.update(self)

        positions = self.positions
        cell = self.cell

        e0 = 4 * epsilon * ((sigma / rc) ** 12 - (sigma / rc) ** 6)

        energy = 0.0
        forces = np.zeros((natoms, 3))
        stress = np.zeros((3, 3))

        for a1 in range(natoms):
            neighbors, offsets = self.nl.get_neighbors(a1)
            cells = np.dot(offsets, cell)
            d = positions[neighbors] + cells - positions[a1]
            r2 = (d ** 2).sum(1)
            c6 = (sigma ** 2 / r2) ** 3
            c6[r2 > rc ** 2] = 0.0
            energy -= e0 * (c6 != 0.0).sum()
            c12 = c6 ** 2
            energy += 4 * epsilon * (c12 - c6).sum()
            f = (24 * epsilon * (2 * c12 - c6) / r2)[:, np.newaxis] * d
            forces[a1] -= f.sum(axis=0)
            for a2, f2 in zip(neighbors, f):
                forces[a2] += f2
            stress += np.dot(f.T, d)

        if self.number_of_lattice_vectors == 3:
            stress += stress.T.copy()
            stress *= -0.5 / self.get_volume()
            stress = stress.flat[[0, 4, 8, 5, 2, 1]]
        else:
            raise PropertyNotImplementedError
        return stress

    def get_forces(self, atoms=None):

        natoms = len(self)
        sigma = 1.0
        epsilon = 1.0
        rc = None

        if rc is None:
            rc = 3 * sigma

        self.nl = NeighborList([rc / 2] * natoms, self_interaction=False)

        self.nl.update(self)

        positions = self.positions
        cell = self.cell

        e0 = 4 * epsilon * ((sigma / rc) ** 12 - (sigma / rc) ** 6)

        energy = 0.0
        forces = np.zeros((natoms, 3))
        stress = np.zeros((3, 3))

        for a1 in range(natoms):
            neighbors, offsets = self.nl.get_neighbors(a1)
            cells = np.dot(offsets, cell)
            d = positions[neighbors] + cells - positions[a1]
            r2 = (d ** 2).sum(1)
            c6 = (sigma ** 2 / r2) ** 3
            c6[r2 > rc ** 2] = 0.0
            energy -= e0 * (c6 != 0.0).sum()
            c12 = c6 ** 2
            energy += 4 * epsilon * (c12 - c6).sum()
            f = (24 * epsilon * (2 * c12 - c6) / r2)[:, np.newaxis] * d
            forces[a1] -= f.sum(axis=0)
            for a2, f2 in zip(neighbors, f):
                forces[a2] += f2
            stress += np.dot(f.T, d)

        if self.number_of_lattice_vectors == 3:
            stress += stress.T.copy()
            stress *= -0.5 / self.get_volume()
            stress = stress.flat[[0, 4, 8, 5, 2, 1]]
        else:
            raise PropertyNotImplementedError
        return forces
