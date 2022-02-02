from GaussianInput import *
from SynapticConstants import *
from ConnectivityMatrix import *

from tqdm import tqdm
import numpy as np

class PING:

    def __init__(self, excit_input, simulation_time, dt, nr_excit, nr_inhibit):
        self.nr_excit = nr_excit
        self.nr_inhibit = nr_inhibit
        self.nr_neurons = self.nr_excit + self.nr_inhibit

        # timescale of recover variable u
        self.a = np.array([0.02 for _ in range(self.nr_excit)] + [0.1 for _ in range(self.nr_inhibit)])
        # sensitivity of  u to sub-threshold oscillations
        self.b = np.array([0.2 for _ in range(self.nr_excit)] + [0.2 for _ in range(self.nr_inhibit)])
        # membrane voltage after spike (reset)
        self.c = np.array([-65 for _ in range(self.nr_excit)] + [-65 for _ in range(self.nr_inhibit)])
        # spike reset of recover variable u
        self.d = np.array([8 for _ in range(self.nr_excit)] + [2 for _ in range(self.nr_inhibit)])

        # initial values of v = voltage
        self.v = np.array([-65 for _ in range(self.nr_excit + self.nr_inhibit)])
        # initial values of u = membrane recovery variable
        self.u = np.multiply(self.b, self.v)

        # spike timings
        self.firings = []

        # spike times
        self.simulation_time = simulation_time
        self.dt = dt

        # gaussian input
        self.gaussian_input = GaussianInput(
            excit_input=2,
            mean_excit_input=[excit_input]*nr_excit,
            inhibit_input=1.5,
            mean_inhibit_input=[4]*nr_inhibit
        )

        # synaptic constants
        self.synaptic_constants = SynapticConstants(
            nr_excit=self.nr_excit,
            nr_inhibit=self.nr_inhibit
        )

        # connectivity matrix
        self.connect_matrix = ConnectivityMatrix(
            nr_neurons=self.nr_neurons,
            nr_excit=self.nr_excit,
            nr_inhibit=self.nr_inhibit
        )

    def run(self):
        print("Simulation started")

        for t in tqdm(range(self.simulation_time)):

            # thalamic input
            I = self._get_new_thalamic_input()

            # indices of spikes
            fired = np.argwhere(self.v > 30).flatten()
            self.firings = [self.firings, [t for _ in range(len(fired))] + fired]
            for f in fired:
                self.v[f] = self.c[f]
                self.u[f] += self.d[f]

            # synaptic potentials
            self.synaptic_constants.gampa = self._get_new_gampa()
            self.synaptic_constants.gaba = self._get_new_gaba()
            gsyn = np.append(self.synaptic_constants.gampa, self.synaptic_constants.gaba)

            # defining input to eah neuron as the summation of all synaptic input
            # form all connected neurons
            I = np.add(I, np.matmul(self.connect_matrix.S, gsyn))

            self.v = np.add(self.v, self._addon_to_v(I=I))
            self.v = np.add(self.v, self._addon_to_v(I=I))
            self.u = np.add(self.u, self._addon_to_u())

        print("Simulation ended")

    def _addon_to_u(self):
        return np.multiply(
            self.a,
            np.multiply(self.b, self.v) - self.u
        )

    def _addon_to_v(self, I):
        return 0.5 * (0.04 * self.v**2 + 5 * self.v + 140 - self.u + I)

    def _get_new_gaba(self):
        alpha = self.v[self.nr_excit:] / 10.0 + 2
        z = np.tanh(alpha)
        comp1 = (z + 1) / 2.0
        comp2 = (1 - self.synaptic_constants.gaba) / self.synaptic_constants.rise_gaba
        comp3 = self.synaptic_constants.gaba / self.synaptic_constants.decay_gaba
        new_comp = self.dt * 0.3 * (np.multiply(comp1, comp2) - comp3)
        return np.add(
            self.synaptic_constants.gaba,
            new_comp
        )

    def _get_new_gampa(self):
        alpha = self.v[:self.nr_excit] / 10.0 + 2
        z = np.tanh(alpha)
        comp1 = (z + 1) / 2.0
        comp2 = (1 - self.synaptic_constants.gampa) / self.synaptic_constants.rise_ampa
        comp3 = self.synaptic_constants.gampa / self.synaptic_constants.decay_ampa
        new_comp = self.dt * 0.3 * (np.multiply(comp1, comp2) - comp3)
        return np.add(
            self.synaptic_constants.gampa,
            new_comp
        )

    def _get_new_thalamic_input(self):
        return np.append(
            self.gaussian_input.excit * np.random.randn(self.nr_excit) + self.gaussian_input.mean_excit,
            self.gaussian_input.inhibit * np.random.randn(self.nr_inhibit) + self.gaussian_input.mean_inhibit
        )