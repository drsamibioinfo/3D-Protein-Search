from Bio.PDB import PDBParser , DSSP
import os,urllib2
from Bio.PDB.Polypeptide import PPBuilder


download_url = "https://files.rcsb.org/download/{0}.pdb".format("6w41")
response = urllib2.urlopen(download_url)
output_file = os.path.join(os.curdir,"6w41.pdb")
with open(output_file,mode='w') as output:
    output.write(response.read())
print("6w41 was downloaded successfully")

parser = PDBParser()
structure = parser.get_structure('6w41',output_file)
dssp = DSSP(structure[0],output_file,acc_array="Wilke")
ss = "".join([aa[2] for aa in dssp])
builder = PPBuilder()
seq = ""
for chain in builder.build_peptides(structure):
    seq += chain.get_sequence()

print("""
Primary Sequence: {0}
Primary Sequence Length: {1}
Secondary Structure : {2}
Length of SS : {3}
""".format(seq,len(seq),ss,len(ss)))
