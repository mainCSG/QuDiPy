'''
File used to generate and plot charge stability diagrams from provided potentials using the Hubbrad model
For more information about the method, see the references https://doi.org/10.1103/PhysRevB.83.235314 and https://doi.org/10.1103/PhysRevB.83.161301
'''

import math
import sys
import copy
import numpy as np
from numpy.linalg.linalg import eigvals
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import itertools
from scipy import linalg as la
from scipy.ndimage import gaussian_filter

class HubbardCSD:
    '''
    Initialize the charge stability diagram class which generates charge stability diagrams based on given capacitance parameters.
    Based on section III in https://doi.org/10.1103/PhysRevB.83.235314.
    This class is intended for use with NextNano potentials to simulate charge stability diagrams of various designs, but can also be used with analytic potentials.
    '''
    def __init__(self, n_sites, n_e, h_mu=False, h_t=False, h_u=False):
        '''

        Parameters
        ----------
        n_sites: Number of sites that a hamilt
        n_e: Maximum number of electrons in the system. Must be less than or equal to 2*n_sites

        Keyword Arguments
        -----------------
        h_mu: Whether to include chemical potential term in Hamiltonian when Hamiltonian is created (default False)
        h_t: Whether to include tunnel coupling term in Hamiltonian when Hamiltonian is created (default False)
        h_u: Whether to include coulomb repulsion repulsion in Hamiltoninan when hamiltonian is created (default False)

        Returns
        -------
        None
        '''

        # Save which parts of Hamiltonian to be used for later
        self.h_mu = h_mu
        self.h_t = h_t
        self.h_u = h_u

        if self.h_mu is False:
            raise Exception("Hamiltonian will be independent of gate volatges so no charge stability diagram can be generated")

        # Check that number of electrons doesn't exceed twice the amount of sites
        if n_e >= 2 * n_sites:
            raise Exception(f"Number of electrons ({n_e}) exceeds twice the amount of sites ({n_sites}) allowed")
        else:
            self.n_e = n_e
            self.n_sites = n_sites

        # Generates the basis to be used
        self.basis, self.basis_labels = self._generate_basis()

        # These next steps generate the fixed portion of the Hamiltonian, which is created on initialization

        # First, generate the matrix of the correct size
        self.fixed_hamiltonian = np.zeros((len(self.basis),len(self.basis)))

        # Then add the component to the fixed portion of the Hamiltonian that you want to consider
        if h_t is True:
            self.fixed_hamiltonian += self._generate_h_t()

        if h_u is True:
            self.fixed_hamiltonian += self._generate_h_u()

    def generate_csd(self, v_g1_max, v_g2_max, v_g1_min=0, v_g2_min=0, c_cs_1=None, c_cs_2=None, num=100, plotting=False, blur=False, blur_sigma=0):

        self.num = num
        self.v_g1_min = v_g1_min
        self.v_g1_max = v_g1_max
        self.v_g2_min = v_g2_min
        self.v_g2_max = v_g2_max

        self.v_1_values = [round(self.v_g1_min + i/num * (v_g1_max - self.v_g1_min), 4) for i in range(num)]
        self.v_2_values = [round(self.v_g2_min + j/num * (v_g2_max - self.v_g2_min), 4) for j in range(num)]

        data = []
        for v_1 in self.v_1_values:
            for v_2 in self.v_2_values:
                h_mu = np.zeros(self.fixed_hamiltonian.shape)

                for i in range(self.fixed_hamiltonian.shape[0]):
                    state_1 = self.basis[i]

                    result = 0
                    for k in range(len(self.basis_labels)):
                        if k == 0 or k == 1:
                            result += - v_1 * self._inner_product(state_1, self._number(state_1, k))
                        if k == 2 or k == 3:
                            result += - v_2 * self._inner_product(state_1, self._number(state_1, k))

                    h_mu[i][i] = result

                eigenvals, eigenvects = la.eig(self.fixed_hamiltonian + h_mu)
                eigenvals = np.round(eigenvals, 4) # To avoid numerical imprecision
                eigenvects = np.round(eigenvects, 4)
                lowest_eigenvect = eigenvects[np.where(eigenvals == eigenvals.min())]
                occupation = (lowest_eigenvect * np.array([sum(x) for x in self.basis])).sum(1).mean()
                data.append([v_1, v_2, occupation])

        data = np.real(data)
        df = pd.DataFrame(data, columns=['V_g1', 'V_g2', 'Current'])
        self.csd = df.pivot_table(index='V_g1', columns='V_g2', values='Current')

        # Show plots if flag is set to True
        if plotting is True:

            # Toggles colorbar if charge sensor information is given
            cbar_flag = False
            if (c_cs_1 is not None) and (c_cs_2 is not None):
                cbar_flag = True

            # Plot the chagre stability diagram
            p1 = sb.heatmap(self.csd, cbar=cbar_flag, xticklabels=int(
                num/5), yticklabels=int(num/5), cbar_kws={'label': 'Current (arb.)'})
            p1.axes.invert_yaxis()
            p1.set(xlabel=r'V$_1$', ylabel=r'V$_2$')
            plt.show()

            # # Plot the "derivative" of the charge stability diagram
            # p2 = sb.heatmap(df_der, cbar=cbar_flag, xticklabels=int(
            # num/5), yticklabels=int(num/5))
            # p2.axes.invert_yaxis()
            # p2.set(xlabel=r'V$_1$', ylabel=r'V$_2$')
            # plt.show()

    def _generate_basis(self):

        # Compute all possible occupations with the cartesian product, and then
        # remove all states that exceed the number of electron specified
        # TODO make this more efficient so we only generate the states we want
        all_combos = list(itertools.product(*[[0,1] for i in range(self.n_sites * 2)])) # * 2 is for spin degeneracy (could add another *2 for valleys)

        basis = []
        for combo in all_combos:
            if sum(combo) <= self.n_e:
                basis.append(list(combo))

        # Labels each index in the basis state with site number and spin direction (could add valley states)
        site = [f'site_{n+1}' for n in range(self.n_sites)]
        spin = ['spin_up', 'spin_down']
        basis_labels = list(itertools.product(site, spin))

        return basis, basis_labels

    def _generate_h_t(self):
        h_t = np.zeros(self.fixed_hamiltonian.shape)
        for i in range(self.fixed_hamiltonian.shape[0]):
            for j in range(i):
                state_1 = self.basis[i]
                state_2 = self.basis[j]

                if sum(state_1) != sum(state_2):
                    continue # No tunnel coupling between states with different charge occupation
                
                if sum(state_1[::2]) != sum(state_2[::2]):
                    continue # No tunnel coupling between states with different number of spins in each orientation

                result = 0
                term_1 = 0
                term_2 = 0
                for k in range(len(self.basis_labels)):
                    for l in range(k):
                        term_1 += 0.1 * self._inner_product(state_1, self._create(self._annihilate(state_2, k), l))
                        term_2 += 0.1 * self._inner_product(state_1, self._create(self._annihilate(state_2, l), k))

                result += term_1 + term_2
                h_t[i][j] = -result
                h_t[j][i] = -result #Since matrix is symmetric
        return h_t

    def _generate_h_u(self):
        h_u = np.zeros(self.fixed_hamiltonian.shape)
        for i in range(self.fixed_hamiltonian.shape[0]):
                state_1 = self.basis[i]

                result = 0
                for k in range(len(self.basis_labels)):
                    for l in range(k):
                        if self.basis_labels[k][0] == self.basis_labels[l][0]: # check if electrons are on same site
                            result += 0.5 * self._inner_product(state_1, self._number(self._number(state_1, k), l))
                        else:
                            result += 0.25 * self._inner_product(state_1, self._number(self._number(state_1, k), l))

                h_u[i][i] = result
        return h_u

    def _inner_product(self, state_1, state_2):
        '''
        docstring
        '''
        if state_1==None or state_2==None: # Deals with cases where eigenvalue of number is 0, so the inner product is multiplied by 0
            return 0
        elif state_1 == state_2:
            return 1
        else:
            return 0

    def _create(self, state, position):
        '''
        docstring
        '''
        state = copy.copy(state)
        if state == None:
            pass # keep state as None
        elif state[position] == 0:
            state[position] = 1
        else:
            state = None
        return state

    def _annihilate(self, state, position):
        '''
        docstring
        '''
        state = copy.copy(state)
        if state == None:
            pass
        elif state[position] == 1:
            state[position] = 0
        else:
            state = None
        return state

    def _number(self, state, position):
        '''
        docstring
        '''
        return self._create(self._annihilate(state, position), position)

    def load_potential(self, file_path):
        # TODO
        pass
