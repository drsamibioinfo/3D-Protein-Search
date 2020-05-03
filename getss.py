#!/usr/bin/python
from Bio.PDB import PDBParser , DSSP
import os,urllib2,sys
from Bio.PDB.Polypeptide import PPBuilder

if len(sys.argv) <= 1:
    print("You have to supply the name of the PDB")
else:
    param = sys.argv[1]
    download_url = "https://files.rcsb.org/download/{0}.pdb".format(param)
    response = urllib2.urlopen(download_url)
    output_file = os.path.join(os.curdir,"{0}.pdb".format(param))
    with open(output_file,mode='w') as output:
        output.write(response.read())
    print("{0} was downloaded successfully".format(param))

    parser = PDBParser()
    structure = parser.get_structure(param,output_file)
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
