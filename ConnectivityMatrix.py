import numpy as np

class ConnectivityMatrix():

    def __init__(self, nr_neurons, nr_excit, nr_inhibit):
        # excitatory  to excitatory
        self.EE = 0.05

        # excitatory  to inhibitory
        self.EI = 0.4

        # inhibitory to excitatory
        self.IE = 0.3

        # inhibitory to inhibitory
        self.II = 0.2

        self.S = self.get_S(
            nr_neurons=nr_neurons,
            nr_excit=nr_excit,
            nr_inhibit=nr_inhibit
        )

    def get_S(self, nr_neurons, nr_excit, nr_inhibit):
        S = np.zeros((nr_neurons, nr_neurons))
        S[:nr_excit, :nr_excit] = self.EE * np.random.standard_normal((nr_excit, nr_excit))
        S[nr_excit:nr_neurons, nr_excit:nr_neurons] = self.II * np.random.standard_normal((nr_inhibit, nr_inhibit))
        S[:nr_excit, nr_excit:nr_neurons] = -self.IE * np.random.standard_normal((nr_excit, nr_inhibit))
        S[nr_excit:nr_neurons, :nr_excit] = -self.EI * np.random.standard_normal((nr_inhibit, nr_excit))
        S = np.nan_to_num(S)
        return S





