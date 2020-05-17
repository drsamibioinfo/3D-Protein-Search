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

############################################################# Calculate Propensities for SASA (Hydrophobicity/Polarity) ###################################

SASA = ['Hydrophobic','Polar']
sasa_types = []
for type in SASA:
    for ss in secondary_structures:
        sasa_types.append((ss,type))

residues_sasa_propensities = {}
residues_sasa_probabilities = {}

for ss,type in sasa_types:
    propss = [0] * (len(residues.keys()) - 1)
    probabilis = [0] * (len(residues.keys()) - 1)
    attr = "{0}_{1}".format(ss,type)
    for i in range(0,len(residues.keys())-1):
        residue_props = (stats[attr][i] / float(stats[type][i])) / (float(sum(stats[attr]))/float(sum(stats[type])))
        residue_probabilities = stats[attr][i] / float(stats[ss][i])
        propss[i] = residue_props
        probabilis[i] = residue_probabilities
    residues_sasa_propensities[attr] = propss
    residues_sasa_probabilities[attr] = probabilis

output = open('ss_sasa_propensities.csv','w')
headers = ['Type'] + residues.keys()
output.write(",".join(headers) + "\r")
for k , v in residues_sasa_propensities.items():
    output.write(",".join([k] + [str(x) for x in v])+"\r")
output.close()
output = open('ss_sasa_probabilities.csv','w')
output.write(",".join(headers) + "\r")
for k , v in residues_sasa_probabilities.items():
    output.write(",".join([k] + [str(x) for x in v])+"\r")
output.close()
print("Secondary Structures propensities in different SASA has been written.")


print("Statistics has been written.")
