from PCFG import PCFG
import numpy as np


def load_sents_to_parse(filename):
    sents = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line:
                sents.append(line)
    return sents


def recursive_derivation_string(pcfg, non_terminals, nt_indices,  rules_bp, s_bp, start_index, end_index, lhs_index):
    rhs = rules_bp[(start_index, end_index, lhs_index)]
    if len(rhs) == 1:
        return non_terminals[lhs_index] + " " + rhs[0]
    lhs_string = non_terminals[lhs_index]
    s = s_bp[start_index][end_index][lhs_index]
    first_string = recursive_derivation_string(pcfg, non_terminals, nt_indices,
                                               rules_bp, s_bp, start_index, s, nt_indices[rhs[0]])
    second_string = recursive_derivation_string(pcfg, non_terminals, nt_indices,
                                                rules_bp, s_bp, s + 1, end_index, nt_indices[rhs[1]])

    return lhs_string + " (" + first_string + ") (" + second_string + ")"


def table_to_string(pcfg, non_terminals, nt_indices, rules_bp, s_bp, n, s_symbol_index, end_symbol):
    inner_str = recursive_derivation_string(pcfg, non_terminals, nt_indices, rules_bp, s_bp, 0, n-1, s_symbol_index)
    return "(ROOT (" + inner_str + ") " + end_symbol + ")"


def q_prob(pcfg, lhs, rhs):
    for rule_rhs, rule_weight in pcfg._rules[lhs]:
        if rule_rhs == rhs or rule_rhs == [rhs]:
            return np.float64(rule_weight / pcfg._sums[lhs])
    return 0.0


def cky(pcfg, sent):
    sent = sent.split()
    end_symbol = sent[-1]
    sent = sent[:-1]

    n = len(sent)

    non_terminals = list(pcfg._rules.keys())
    non_terminals.extend(['.', '!'])
    nt_indices = {nt: i for i, nt in enumerate(non_terminals)}
    N = len(non_terminals)

    table = np.zeros((n, n, N), np.float64)
    rules_bp = dict()
    s_bp = np.zeros((n, n, N), np.int16)
    for i, word in enumerate(sent):
        for x, nt in enumerate(non_terminals):
            table[i][i][x] = q_prob(pcfg, nt, word)
            rules_bp[(i, i, x)] = [word]

    for l in xrange(1, n):
        for i, word in enumerate(sent):
            j = i + l
            if j >= n: break
            for x, nt in enumerate(non_terminals):
                curr_max, arg_max_s, arg_max_rhs = 0.0, 0, 0
                for rhs, w in pcfg._rules[nt]:
                    if not pcfg.is_preterminal(rhs):
                        for s in xrange(i, j):
                            y, z = nt_indices[rhs[0]], nt_indices[rhs[1]]
                            val = q_prob(pcfg, nt, rhs) * table[i][s][y] * table[s + 1][j][z]
                            if val > curr_max:
                                curr_max = val
                                arg_max_s = s
                                arg_max_rhs = rhs
                table[i][j][x] = curr_max
                rules_bp[(i, j, x)] = arg_max_rhs
                s_bp[i][j][x] = arg_max_s

    if table[0][n-1][nt_indices["S"]] == 0:
        return "FAILED TO PARSE!"
    return table_to_string(pcfg, non_terminals, nt_indices, rules_bp, s_bp, n, nt_indices['S'], end_symbol)


if __name__ == '__main__':
    import sys
    pcfg = PCFG.from_file_assert_cnf(sys.argv[1])
    sents_to_parse = load_sents_to_parse(sys.argv[2])
    for sent in sents_to_parse:
        print cky(pcfg, sent)
