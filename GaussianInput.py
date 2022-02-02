class GaussianInput:

    def __init__(self, excit_input, mean_excit_input, inhibit_input, mean_inhibit_input):
        # gaussian input to excit. neurons
        self.excit = excit_input
        self.mean_excit = mean_excit_input

        # gaussian input to inhibit. neurons
        self.inhibit = inhibit_input
        self.mean_inhibit = mean_inhibit_input