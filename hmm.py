initial_probs = {"H": 0.1, "S": 0.1, "L": 0.7, "T": 0.1}

emission_probs = {'H': {'A': 0.10936768965100001, 'C': 0.0120779196309, 'E': 0.08765545587009999, 'D': 0.0494356886845,
                        'G': 0.036320624799500004, 'F': 0.0410824185458, 'I': 0.0613439825582, 'H': 0.021650095459,
                        'K': 0.06429801648729999, 'M': 0.027698864187200002, 'L': 0.119967292066,
                        'N': 0.032890607674000004, 'Q': 0.0464852589605, 'P': 0.023702530621900003,
                        'S': 0.0496523059773, 'R': 0.0590054465015, 'T': 0.0430491283732, 'W': 0.015331969423,
                        'V': 0.0642414897759, 'Y': 0.034743214753},
                  'S': {'A': 0.061978103675699996, 'C': 0.020655746687, 'E': 0.0443725923907, 'D': 0.0302529632084,
                        'G': 0.0483672537899, 'F': 0.0566702831483, 'I': 0.0961072519038, 'H': 0.022512165280600002,
                        'K': 0.0468276119082, 'M': 0.022289618619399997, 'L': 0.103110741013, 'N': 0.0258673691175,
                        'Q': 0.030366007780599998, 'P': 0.0188204255483, 'S': 0.0551174159963, 'R': 0.0432645823956,
                        'T': 0.0694092885796, 'W': 0.0185181336562, 'V': 0.133828556728, 'Y': 0.051663888572700004},
                  'L': {'A': 0.0629678572883, 'C': 0.014113526723200001, 'E': 0.0543433607893, 'D': 0.0764986318714,
                        'G': 0.111476427641, 'F': 0.0314745483818, 'I': 0.0353650322419, 'H': 0.0250797215404,
                        'K': 0.057279623429799996, 'M': 0.0177023456215, 'L': 0.0646919014075, 'N': 0.0593511915713,
                        'Q': 0.034062859409700004, 'P': 0.08744285923650001, 'S': 0.0759237901519,
                        'R': 0.046280824851600004, 'T': 0.0608880580334, 'W': 0.010841071203499999,
                        'V': 0.0461452252754, 'Y': 0.0280711433305},
                  'T': {'A': 0.0555519635253, 'C': 0.0127708829702, 'E': 0.059153932649299995, 'D': 0.0771780383593,
                        'G': 0.133118751387, 'F': 0.0347317458288, 'I': 0.0366403292525, 'H': 0.0261878452502,
                        'K': 0.0602976566945, 'M': 0.0158100646744, 'L': 0.061276093331599994, 'N': 0.0591169860522,
                        'Q': 0.0343238748424, 'P': 0.0516644685062, 'S': 0.0771712324072, 'R': 0.0494739242099,
                        'T': 0.0643567589343, 'W': 0.011443884357800001, 'V': 0.0474702195033, 'Y': 0.0322613472641}}

transition_probs = {
    "H": {"H": 0.59, "S": 0.01, "L": 0.2, "T": 0.2},
    "S": {"H": 0.01, "S": 0.59, "L": 0.2, "T": 0.2},
    "L": {"H": 0.2, "S": 0.2, "L": 0.59, "T": 0.01},
    "T": {"H": 0.2, "S": 0.2, "L": 0.01, "T": 0.59}
}


class HiddenMarkovModel(object):

    def __init__(self, init_probs, emission_probs, transition_probs):
        self.init_probs = init_probs
        self.emission_probs = emission_probs
        self.transition_probs = transition_probs
        self.states = self.init_probs.keys()
        self.symbols = self.emission_probs[self.states[0]].keys()

    def get_init_prob(self, state):
        if state in self.init_probs.keys():
            return self.init_probs[state]
        else:
            return 0

    def get_emission_prob(self, state, symbol):
        if state in self.states and symbol in self.symbols:
            return self.emission_probs[state][symbol]
        else:
            return 0

    def get_transition_prob(self, original, dest):
        if original in self.states and dest in self.states:
            return self.transition_probs[original][dest]
        else:
            return 0

    def calculate_joint_prob(self, path, sequence):
        """
        This method will calculate the joint probability for the given sequence and state paths defined by path variable
        :param path: the state paths for this particular sequence
        :param sequence: the symbol sequence ,  Nucleic acid sequence
        :return:
        """
        if len(sequence) < 1:
            return 0
        if len(path) < 1:
            return 0

        if len(sequence) != len(path):
            print("Sequence and state paths don't match each others' length")
            return 0

        prob = self.get_init_prob(path[0]) * self.get_emission_prob(path[0], sequence[0])
        for pos in range(1, len(path)):
            prob = prob * self.get_emission_prob(path[pos], sequence[pos]) * self.get_transition_prob(path[pos - 1],
                                                                                                      path[pos])

        return prob

    def viterbi(self, sequence):

        if len(sequence) == 0:
            return []

        viterbi = {}
        state_path = {}
        for state in self.states:
            viterbi[state] = self.get_init_prob(state) * self.get_emission_prob(state, sequence[0])
            state_path[state] = [state]

        for t in range(1, len(sequence)):
            new_state_path = {}
            new_path = {}
            viterbi_tmp = {}
            for state_dest in self.states:
                intermediate_probs = []
                for state_orig in self.states:
                    prob = viterbi[state_orig] * self.get_transition_prob(state_orig, state_dest)
                    intermediate_probs.append((prob, state_orig))

                max_prob, max_state = max(intermediate_probs)
                prob = self.get_emission_prob(state_dest, sequence[t]) * max_prob
                viterbi_tmp[state_dest] = prob
                new_state_path[state_dest] = max_state
                new_path[state_dest] = state_path[max_state] + [state_dest]
            viterbi = viterbi_tmp
            state_path = new_path

        max_state = None
        max_prob = -1
        for state in self.states:
            if viterbi[state] > max_prob:
                max_prob = viterbi[state]
                max_state = state

        return max_prob, state_path[max_state]


if __name__ == '__main__':
    hmm = HiddenMarkovModel(initial_probs, emission_probs, transition_probs)
    max_prob, paath = hmm.viterbi(
        "MDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGGGDPGETLRFADAIIVGASIYLADNPAAAAAGIIESIKMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAAGGDPGETLRFADAIIVGASIYLADNPAAAAAGIIESIK")
    print(max_prob)
    print("".join(paath))
