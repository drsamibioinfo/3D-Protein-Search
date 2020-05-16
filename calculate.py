from propensities import residues
import sys, os, pickle, argparse


def main():
    pass


if __name__ == '__main__':
    main()

file = open('propensities.db', 'rb')
stats = pickle.load(file)
file.close()
secondary_structures = ['H', 'S', 'L', 'T']
ss_props = {s: [(stats[s][i]) / float(stats["{0}_total".format(s)]) for i in range(len(residues.keys()))] for s in
            secondary_structures}
fi_props = []
for i in range(len(residues.keys())):
    total = 0
    for s in secondary_structures:
        total += stats[s][i]
    fi_props.append(total / float(stats['total_residues']))
props = {}
for ss, pp in ss_props.items():
    props[ss] = [x / y for x, y in zip(pp[:-1], fi_props[:-1])]

############################################### Calculate Propensities######################################################################################
output = open('ss_propensities.csv', 'w')
output.write(",".join(["class"] + residues.keys()[:-1]) + "\r")
for ss, pps in props.items():
    oneLine = ",".join([ss] + [str(x) for x in pps]) + "\r"
    output.write(oneLine)
output.close()
####################################################### Outputs Probabilities ###############################################################################
output = open('ss_propabilities.csv', 'w')
output.write(",".join(["class"] + residues.keys()[:-1]) + "\r")
for ss, pps in ss_props.items():
    oneLine = ",".join([ss] + [str(x * 100) for x in pps[:-1]]) + "\r"
    output.write(oneLine)
output.close()

output = open('BinaryPatterning_polarity.csv', 'w')
output.write(",".join(["class"] + residues.keys()[:-1]) + "\r")
polarity = [x / float(stats['total_residues']) for x in stats['Polar'][:-1]]
hydrophobicity = [x / float(stats['total_residues']) for x in stats['Hydrophobic'][:-1]]
output.write(",".join(["Polar"] + [str(x) for x in polarity]) + "\r")
output.write(",".join(["Hydrophobicity"] + [str(x) for x in hydrophobicity]) + "\r")
output.close()

print("Statistics has been written.")
