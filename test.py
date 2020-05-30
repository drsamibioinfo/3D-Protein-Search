from __future__ import print_function
from Bio.PDB import PDBParser, DSSP, MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder, Polypeptide
import Bio.PDB, os


class EPolyPeptide(Polypeptide):
    def __init__(self, pp):
        for res in pp:
            self.append(res)

    @property
    def start(self):
        start = self[0].get_id()[1]
        return int(start)

    @property
    def end(self):
        end = self[-1].get_id()[1]
        return int(end)

    @property
    def chain_id(self):
        full_id = self[0].full_id
        return full_id[2]

    def __repr__(self):
        sup = super(EPolyPeptide, self).__repr__()
        return "<Chain: {0}, {1}>".format(self.chain_id,sup)


def get_primary_sequence(input_file):
    file_name, _ = os.path.splitext(input_file)
    file_name = file_name.replace('./', '')
    parser = PDBParser()
    structure = parser.get_structure(file_name, input_file)
    builder = PPBuilder()
    seq = ""
    for chain in builder.build_peptides(structure, aa_only=False):
        seq += chain.get_sequence()
    return seq


def get_sequence_position(input_file, chain_id, start_position, end_position):
    file_name, _ = os.path.splitext(input_file)
    file_name = file_name.replace('./', '')
    parser = PDBParser()
    structure = parser.get_structure(file_name, input_file)
    builder = PPBuilder()
    peptides = builder.build_peptides(structure,aa_only=False)
    pps = [EPolyPeptide(pp) for pp in peptides]
    seq_leftover = 0
    start = None
    end = None
    for pp in pps:
        if not pp.chain_id == chain_id:
            seq_leftover += len(pp)
            continue
        start = int(start_position) - pp.start
        end = int(end_position) - pp.start
        break

    if not start and not end:
        return -1
    else:
        return seq_leftover + start


def get_chain_position(input_file, global_index):
    chain = None
    position_in_chain = -1
    file_name, _ = os.path.splitext(input_file)
    file_name = file_name.replace('./', '')
    parser = PDBParser()
    structure = parser.get_structure(file_name, input_file)
    builder = PPBuilder()
    peptides = builder.build_peptides(structure, aa_only=False)
    pps = [EPolyPeptide(pp) for pp in peptides]
    total_length = sum([len(pp) for pp in pps])
    if global_index >= total_length:
        return None, -1
    distance = 0
    offset = 0
    global_index = int(global_index) + 1
    while distance < global_index:
        pp = pps[offset]
        distance += len(pp)
        offset += 1
        if global_index <= distance:
            position_in_chain = global_index - (distance - len(pp))
            chain = pp.chain_id
            break
    return chain , position_in_chain



def main():
    input_file = os.path.join(".", "6w41.pdb")
    print("Primary sequence : ", get_primary_sequence(input_file))
    # get chain sequence
    sequence_index = get_sequence_position(input_file,"C",333,343)
    print("Residue in chain C with positions {0} , {1} has a sequence index of {2}".format("C","333:343",sequence_index))
    chain , chain_position = get_chain_position(input_file,443)
    print("Residue exist in chain {0} and chain position {1} ".format(chain,chain_position))


if __name__ == '__main__':
    main()
