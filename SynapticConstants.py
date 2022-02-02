import numpy as np

class SynapticConstants:

    def __init__(self, nr_excit, nr_inhibit):
        self.gampa = np.zeros(nr_excit)
        self.gaba = np.zeros(nr_inhibit)

        self.decay_ampa = 1.0
        self.rise_ampa = 0.15

        self.decay_gaba = 7.0
        self.rise_gaba = 0.2