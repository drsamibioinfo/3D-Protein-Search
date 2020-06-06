import os, sys, argparse, traceback, urllib2
from infi.clickhouse_orm.database import Database
from infi.clickhouse_orm.models import Model
from infi.clickhouse_orm.fields import *
from datetime import datetime, date
from terminaltables import AsciiTable
from Bio.PDB import PDBParser, DSSP, MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder, Polypeptide
import Bio.PDB
from Bio.PDB.vectors import calc_dihedral
from collections import OrderedDict
from math import sqrt, degrees, atan2, sin, cos
from fuzzysearch import find_near_matches

# borrowed from SARI Sabban <http://www.github.com/sarisabban>
# These numbers represent different residues molecular weights
residues = OrderedDict({'A': 129, 'P': 159, 'N': 195, 'H': 224,
                        'V': 174, 'Y': 263, 'C': 167, 'K': 236,
                        'I': 197, 'F': 240, 'Q': 225, 'S': 155,
                        'L': 201, 'W': 285, 'E': 223, 'T': 172,
                        'M': 224, 'R': 274, 'G': 104, 'D': 193, 'X': 0})


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
        return "<Chain: {0}, {1}>".format(self.chain_id, sup)


def get_sequence_position(structure, chain_id, start_position, end_position):
    builder = PPBuilder()
    peptides = builder.build_peptides(structure, aa_only=False)
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


def get_chain_position(structure, global_index):
    chain = None
    position_in_chain = -1
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
    return chain, position_in_chain


class ProteinModel(Model):
    protein_id = StringField()
    size = UInt64Field()
    helices = UInt64Field()
    sheets = UInt64Field()
    loops = UInt64Field()
    bends = UInt64Field()
    sequence = StringField()
    ss = StringField()
    pattern = StringField()
    enhanced = StringField()
    processdate = DateField(default=date.today())

    @classmethod
    def table_name(cls):
        return "proteins"


class Patterning(object):
    def __init__(self):
        self.p = argparse.ArgumentParser(
            description="Patterning is a program that allows to perform de-novo protein design " \
                        "by searching for possible proteins or parts of protein which match a given 3D structure or a possible 1D secondary structure with possible binary patterning fingerprint")

    def setup(self):
        self.db = None
        self.common_length = 0.0
        self.p.add_argument("-w", "--workdir", default=".",
                            help="The working directory for this script to save and load temporary. By default it "
                                 "takes the current working directory. "
                                 "memory files required to function properly.")
        self.p.add_argument("-q", "--pdbs", type=str, default=None,
                            help="Local repository containing PDB Files to read from")
        self.p.add_argument("-d", "--database", default="propensities", help="Database Name to use")
        self.p.add_argument("-t", "--host", default="localhost", help="Database Host to use. defaults to localhost")
        self.p.add_argument("-r", "--port", default=8123, help="Database Default port to use")
        self.p.add_argument("-u", "--username", help="Username of the database to use")
        self.p.add_argument("-o", "--password", help="Database password to use to login")
        self.p.add_argument("-s", "--source", type=str, default=None
                            ,
                            help="Source PDB which contains the 3D structure of a protein or part of a protein (Fragment), that you want to search for similars in the database.")
        self.p.add_argument("-p", "--pattern", type=str, default=None, help="Binary Pattern of a fragment")
        self.p.add_argument("-e", "--structure", type=str, default=None,
                            help="Possible Secondary structure of the fragment/protein, multiples of (H,B,E,G,I,T,S,-).")
        self.p.add_argument("-m", "--rmsd", type=float, default=-1,
                            help="Return matched proteins with this RMSD value only and exclude any others.")
        self.p.add_argument("-l", "--limit", type=int, default=10,
                            help="Total Number of results to include. Defaults to 10.")
        self.p.add_argument("-x", "--start", type=int, default=-1,
                            help="When searching by a whole protein containing a "
                                 "specific scaffold. This parameter denotes the "
                                 "location of the first residue of the fragment.")
        self.p.add_argument("-y", "--end", type=int, default=-1, help="When searching by a whole protein containing a "
                                                                      "specific scaffold. This parameter denotes the "
                                                                      "location of the last residue of the fragment.")
        self.p.add_argument("-a", "--chain", type=str, default=None,
                            help="Chain Identifier if your start and end are relative to a particular chain.")
        self.p.add_argument('-c', '--cutoff', default=23.9, type=float, help="Cutoff value to use to mark a residue "
                                                                             "Solvent-accessible surface area as polar or "
                                                                             "hydrophobic, if the SASA of a residue equal "
                                                                             "to or higher than this value it will be "
                                                                             "considered polar otherwise it will be "
                                                                             "marked as hydrophobic")
        self.p.add_argument("-f", "--fuzzy", nargs="?", const=True, default=False, help="Perform Fuzzy Search. Allow "
                                                                                        "matching similar but not "
                                                                                        "identical results.")
        self.p.add_argument("-z", "--fuzzylevel", type=float, default=0.90, help="Include only results with equal or "
                                                                                 "higher value than the following "
                                                                                 "fuzziness level .Defaults to 0.90")
        self.p.add_argument("-v","--distance",type=int,default=1,help="Possible Lovenstein distance for fuzzy string search")
        self.p.add_argument("-j","--deletions",type=int,default=0,help="Number of allowed string deletions when fuzzy string search is enabled. Defaults to Zero.")
        self.p.add_argument("-n","--insertions",type=int,default=0,help="Number of allowed string insertions when fuzzy string search is enabled. Defaults to Zero.")
        self.p.add_argument("-k","--substitutions",type=int,default=0,help="Number of allowed string substitutions when fuzzy string search is enabled. Defaults to Zero.")


        if len(sys.argv) <= 1:
            self.p.print_help()
            return False
        self.args = self.p.parse_args()
        return True

    def connect(self):
        db_url = "http://{0}:{1}".format(self.args.host or "localhost", self.args.port or 8123)
        self.db = Database(db_name=self.args.database, db_url=db_url, username=self.args.username,
                           password=self.args.password)
        return True

    def print_results(self, headers, rows):
        limit = self.args.limit if self.args.limit < len(rows) else len(rows)
        data = [headers] + [[getattr(row, x) for x in headers if hasattr(row, x)] for row in rows]
        table = AsciiTable(data)
        table.title = "Possible Matche(s)"
        table.inner_row_border = True
        print(table.table)
        output = 'Total : {0}'.format(len(data) - 1)
        print(output)

    def search_by_common_pattern(self, pattern):
        if not self.db:
            print("Couldn't connect to the database.")
            self.p.print_help()
            return
        if self.args.fuzzy:
            query = "select p.* , ngramSearch(enhanced,'{0}') as ngram from proteins as p where ngram > {2} order by ngram DESC limit {3}".format(
                pattern, pattern, self.args.fuzzylevel, self.args.limit
            )
        else:
            query = "select p.* , position(enhanced,'{0}') as pos from proteins as p where position(enhanced,'{1}') > 0 limit {2}".format(
                pattern, pattern, self.args.limit)
        rows = []
        headers = ["protein_id"]
        if not self.args.fuzzy:
            headers += ["pos"]
        else:
            headers += ["ngram"]
        for row in self.db.select(query):

            rows.append(row)
        # if len(rows) > 0:
        #     self.print_results(headers, rows)
        return rows , headers

    def get_secondary_structure_details(self, name, pdb_file, aa_only=False):
        parser = PDBParser()
        structure = parser.get_structure(name, pdb_file)
        dssp = DSSP(structure[0], pdb_file, acc_array="Wilke")
        ss = "".join([aa[2] for aa in dssp])
        sasa = [residues[aa[1]] * aa[3] for aa in dssp]
        builder = PPBuilder()
        seq = ""
        for chain in builder.build_peptides(structure, aa_only=aa_only):
            seq += chain.get_sequence()
        return name, seq, ss, sasa, structure

    def get_enhanced(self, ss, pattern):
        sequence = ""
        for index, letter in enumerate(ss):
            sasa = pattern[index]
            if int(sasa) == 0:
                if letter == 'H':
                    sequence += 'I'
                elif letter == 'B':
                    sequence += 'J'
                elif letter == 'E':
                    sequence += 'K'
                elif letter == 'G':
                    sequence += 'L'
                elif letter == 'I':
                    sequence += 'M'
                elif letter == 'T':
                    sequence += 'N'
                elif letter == 'S':
                    sequence += 'O'
                else:
                    sequence += 'P'
            else:
                if letter == 'H':
                    sequence += 'A'
                elif letter == 'B':
                    sequence += 'B'
                elif letter == 'E':
                    sequence += 'C'
                elif letter == 'G':
                    sequence += 'D'
                elif letter == 'I':
                    sequence += 'E'
                elif letter == 'T':
                    sequence += 'F'
                elif letter == 'S':
                    sequence += 'G'
                else:
                    sequence += 'H'
        return sequence

    def start(self):
        if not self.setup():
            return
        self.connect()

        if not self.args.source is None:
            self.load_source_pdb()
        elif self.args.pattern is not None and self.args.structure is not None:
            self.process_patterning()
        else:
            self.p.print_help()
            return

    def load_source_pdb(self):
        source_file = self.args.source
        base_name = os.path.basename(source_file)
        name, _ = os.path.splitext(base_name)
        _, seq, ss, sasa, structure = self.get_secondary_structure_details(name, source_file)
        if self.args.start != -1 and self.args.end != -1 and not self.args.chain:
            seq_start = self.args.start - 1
            fragment_length = self.args.end - self.args.start
            seq = seq[seq_start:self.args.end + 1]
            ss = ss[seq_start:self.args.end + 1]
            sasa = sasa[seq_start:self.args.end + 1]
        else:
            start = self.args.start
            end = self.args.end
            fragment_length = end - start
            chain_id = self.args.chain
            seq_start = get_sequence_position(structure, chain_id, start, end)
            if seq_start == -1:
                self.p.error(
                    "Unable to get the sequence position from the chain identifier , start position and end position")
                self.p.print_help()
                return
            seq = seq[seq_start:seq_start + fragment_length + 1]
            ss = ss[seq_start:seq_start + fragment_length + 1]
            sasa = sasa[seq_start:seq_start + fragment_length + 1]
        asa = [1 if a >= self.args.cutoff else 0 for a in sasa] if self.args.pattern is None else [int(x) for x in
                                                                                                   self.args.pattern]
        common_sequence = self.get_enhanced(ss, asa)
        self.common_length = len(common_sequence)
        found_rows , incoming_headers = self.search_by_common_pattern(common_sequence)
        incoming_headers_set = set(incoming_headers)
        if len(found_rows) <= 0:
            print("No Records found.")
            return
        if self.args.fuzzy and len(found_rows) > 0:
            new_matches = []
            deleting = []
            for row in found_rows:
                matches = find_near_matches(common_sequence,row.enhanced,max_l_dist=self.args.distance,max_deletions=self.args.deletions,max_insertions=self.args.insertions,
                                            max_substitutions=self.args.substitutions)
                if len(matches) > 0:
                    repeats = 1
                    for match in matches:
                        if repeats <= 1:
                            setattr(row, "pos", match.start)
                            setattr(row, "end_pos", match.end)
                        else:
                            keys = row.to_dict()
                            ngram = keys['ngram']
                            del keys['ngram']
                            new_model = ProteinModel(**keys)
                            setattr(new_model, "pos", match.start)
                            setattr(new_model, "end_pos", match.end)
                            setattr(new_model,"ngram",ngram)
                            new_matches.append(new_model)
                        repeats += 1
                else:
                    deleting.append(row)
            found_rows.extend(new_matches)
            for todelete in deleting:
                found_rows.remove(todelete)
        print ("Calculating Elenated Score values. Please Wait....")
        deviated_rows = self.calculate_elenated_topology_score(seq, found_rows, seq_start, fragment_length)
        deviated_rows = sorted([x for x in deviated_rows if x.rmsd > -1], key=lambda x: x.rmsd, reverse=False)
        proper_headers_set = ["protein_id", "pos", "chain", "chain_pos","end_pos", "chain_length", "deviation", "rmsd"]
        if self.args.fuzzy:
            proper_headers_set += ['ngram']
        self.print_results(headers=proper_headers_set,
                           rows=deviated_rows)

    def calculate_elenated_topology_score(self, seq, rows, seq_start, fragment_length):
        try:
            source_structure = self.__get_structure__(self.args.source)
            source_residues = [res for res in source_structure.get_residues()]
            for row in rows:
                position = row.pos
                try:
                    target_file = self.get_target_file(row.protein_id)
                except Exception as e:
                    setattr(row, 'rmsd', -1)
                    setattr(row, 'chain_length', -1)
                    setattr(row, 'chain_pos', -1)
                    setattr(row, 'deviation', -1)
                    setattr(row, 'chain', "N/A")
                    continue
                target_structure = self.__get_structure__(target_file)
                target_residues = [res for res in target_structure.get_residues()]
                start_offset_residue = target_residues[position]
                chain, chain_position = get_chain_position(target_structure, position)
                setattr(row, "chain",
                        "{0}".format(chain))
                chain_length = self.get_chain_length(target_structure, start_offset_residue.full_id[2])
                setattr(row, "chain_length", chain_length)
                setattr(row, "chain_pos", chain_position)
                setattr(row,"end_pos",chain_position+fragment_length)
                current_deviation = self.__get_elenated_topology(
                    source_residues[seq_start: (seq_start + fragment_length) + 1],
                    target_residues[position - 1:position + fragment_length + 1], len(seq))
                setattr(row, "deviation", current_deviation)
                self.calculate_rmsd_deviation(row, seq_start, fragment_length, aa_only=False)

        except Exception as e:
            print(e.message)
            raise e

        finally:
            return rows

    def get_chain_length(self, target_structure, chain_name):
        length = 0
        for chain in target_structure.get_chains():
            if chain.id == chain_name:
                length = len(chain)
                break
        return length

    def get_chain_polypeptide(self, structure, pps, global_index):
        chain, position_in_chain = get_chain_position(structure, global_index)
        polypeptides = [EPolyPeptide(pp) for pp in pps]
        current_chain = None
        for pp in polypeptides:
            if pp.chain_id == chain:
                current_chain = pp
                break
        return chain, current_chain, position_in_chain

    def calc_rmsd(self, source_atoms, target_atoms):
        from math import sqrt
        import numpy as np
        from Bio.PDB.QCPSuperimposer import QCPSuperimposer
        if len(source_atoms) != len(target_atoms):
            return -1
        source_arr = []
        for atom in source_atoms:
            xyz = [atom.coord[0], atom.coord[1], atom.coord[2]]
            source_arr.append(xyz)
        source_arr = np.array(source_arr)
        target_arr = []
        for atom in target_atoms:
            xyz = [atom.coord[0], atom.coord[1], atom.coord[2]]
            target_arr.append(xyz)
        target_arr = np.array(target_arr)
        sup = QCPSuperimposer()
        sup.set(source_arr, target_arr)
        sup.run()
        return sup.get_rms()

    def calculate_rmsd_deviation(self, row, source_position, fragment_length, aa_only=False):
        if self.args.source is None:
            setattr(row, "rmsd", -1)
        target_position = row.pos
        source_structure = self.__get_structure__(self.args.source)
        builder = PPBuilder()
        source_pps = [x for x in builder.build_peptides(source_structure, aa_only=aa_only)]
        source_length = sum([len(x) for x in source_pps])
        source_residues = []
        for pp in source_pps:
            source_residues += [x for x in pp]
        # source_backbones = [atom['CA'] for atom in source_residues[source_position:source_position + fragment_length+1]]
        source_chain_name, source_chain, source_position_in_chain = self.get_chain_polypeptide(source_structure,
                                                                                               source_pps,
                                                                                               source_position)
        source_backbones = [atom['CA'] for atom in
                            source_chain[source_position_in_chain - 1:source_position_in_chain + fragment_length]]
        source_backbone_residues = " ".join(
            [x.resname for x in source_chain[source_position_in_chain - 1:source_position_in_chain + fragment_length]])
        builder = PPBuilder()
        target_file = self.get_target_file(row.protein_id)
        if target_file is None:
            setattr(row, "rmsd", -1)
            return
        target_structure = self.__get_structure__(target_file)
        target_pps = [x for x in builder.build_peptides(target_structure, aa_only=aa_only)]
        target_length = sum([len(x) for x in target_pps])
        target_residues = []
        for pp in target_pps:
            target_residues += [x for x in pp]
        # target_backbone = [atom['CA'] for atom in target_residues[target_position:target_position + fragment_length+1]]
        target_chain_name, target_chain, target_position_in_chain = self.get_chain_polypeptide(target_structure,
                                                                                               target_pps,
                                                                                               target_position)
        target_backbone = [atom['CA'] for atom in
                           target_chain[target_position_in_chain - 1:target_position_in_chain + fragment_length]]
        target_backbone_residues = " ".join(
            [x.resname for x in target_chain[target_position_in_chain - 1:target_position_in_chain + fragment_length]])
        lengths = [source_length, target_length]
        smallest = min(int(item) for item in lengths)
        # find RMSD
        if len(source_backbones) != len(target_backbone):
            setattr(row, 'rmsd', -1)
            return
        # sup = Bio.PDB.Superimposer()
        # sup.set_atoms(source_backbones, target_backbone)
        # sup.apply(source_structure[0].get_atoms())
        # RMSD = round(sup.rms, 4)
        # print("RMSD For : {0} - {1} : {2}".format(source_backbone_residues,target_backbone_residues,RMSD))
        RMSD = self.calc_rmsd(source_backbones, target_backbone)
        RMSD = round(RMSD, 4)
        setattr(row, "rmsd", RMSD)

    def calculate_RMSD(self, row, source_position, fragment_length, aa_only=False):
        if self.args.source is None:
            setattr(row, "rmsd", -1)
        target_position = row.pos
        source_structure = self.__get_structure__(self.args.source)
        builder = PPBuilder()
        type1 = builder.build_peptides(source_structure, aa_only=aa_only)
        length1 = type1[-1][-1].get_full_id()[3][1]
        fixed_residues = []
        for pp in type1:
            fixed_residues += [x for x in pp]
        fixed = [atom['CA'] for atom in fixed_residues][source_position:source_position + fragment_length]
        builder = PPBuilder()
        target_file = self.get_target_file(row.protein_id)
        if target_file is None:
            setattr(row, "rmsd", -1)
            return
        target_structure = self.__get_structure__(target_file)
        type2 = builder.build_peptides(target_structure, aa_only=aa_only)
        length2 = type2[-1][-1].get_full_id()[3][1]
        moving_residues = []
        for pp in type2:
            moving_residues += [x for x in pp]
        moving = [atom['CA'] for atom in moving_residues][target_position:target_position + fragment_length]
        lengths = [length1, length2]
        smallest = min(int(item) for item in lengths)
        # find RMSD
        if len(fixed) != len(moving):
            setattr(row, "rmsd", -1)
            return
        sup = Bio.PDB.Superimposer()
        sup.set_atoms(fixed, moving)
        sup.apply(target_structure[0].get_atoms())
        RMSD = round(sup.rms, 4)
        setattr(row, "rmsd", RMSD)

    def __get_structure__(self, file_path):
        base_name = os.path.basename(file_path)
        name, ext = os.path.splitext(base_name)
        if 'cif' in ext:
            parser = MMCIFParser()
        else:
            parser = PDBParser()
        return parser.get_structure(name, file_path)

    def get_target_file(self, protein_id):
        if not self.args.pdbs is None:
            which_file = os.path.join(self.args.pdbs, protein_id)
            if os.path.exists(which_file + ".pdb"):
                return which_file + ".pdb"
            elif os.path.exists(which_file + ".cif"):
                return which_file + ".cif"
            else:
                return which_file + ".pdb"
        elif os.path.exists(os.path.join(self.args.workdir, "{0}.pdb".format(protein_id))):
            return os.path.join(self.args.workdir, "{0}.pdb".format(protein_id))
        else:
            print("Downloading File : {0}".format(protein_id))
            download_url = "https://files.rcsb.org/download/{0}.pdb".format(protein_id)
            response = urllib2.urlopen(download_url)
            output_file = os.path.join(self.args.workdir, "{0}.pdb".format(protein_id))
            with open(output_file, mode='w') as output:
                output.write(response.read())
            print("Downloaded.")
            return output_file

    def get_phi(self, previous, source_res):
        try:
            C_1 = previous['C'].get_vector()
            N = source_res['N'].get_vector()
            CA = source_res['CA'].get_vector()
            C = source_res['C'].get_vector()
            return degrees(calc_dihedral(C_1, N, CA, C))
        except Exception as e:
            return 0.0

    def get_psi(self, target_res, next_res):
        try:
            N = target_res['N'].get_vector()
            CA = target_res['CA'].get_vector()
            C = target_res['C'].get_vector()
            N1_1 = next_res['N'].get_vector()
            return degrees(calc_dihedral(N, CA, C, N1_1))
        except Exception as e:
            return 0.0

    def __get_elenated_topology(self, source_residues, target_residues, seq_length):
        """
        This method will calculate the elenated topology mean square deviation
        :param source_residues: Fragment source residues
        :param target_residues: Target protein residues
        :return: Elenated T(msd) Value between these two proteins
        """
        # target residues should be longer than the source residues in most of cases
        if len(source_residues) > len(target_residues):
            return 0.0
        deviation = 0.0
        total = 0.0
        for index, res in enumerate(source_residues[:-1]):
            if index == 0:
                continue
            source_res = source_residues[index]
            target_res = target_residues[index]
            # calculate Phi and Psi torsional angles for the current residues
            source_phi = self.get_phi(source_residues[index - 1], source_res)
            target_phi = self.get_phi(target_residues[index - 1], target_res)
            source_psi = self.get_psi(source_res, source_residues[index + 1])
            target_psi = self.get_psi(target_res, target_residues[index + 1])
            deviation += atan2(sin(source_phi - target_phi), cos(source_phi - target_phi)) + atan2(
                sin(source_psi - target_psi), cos(source_psi - target_psi))
            total += 1
        if total == 0:
            return 0.0
        else:
            return (abs(deviation) / float(seq_length)) * 100.0

    def process_patterning(self):
        ss = self.args.structure
        pattern = self.args.pattern
        if len(ss) != len(pattern):
            print("Length of Both Secondary Structure and Binary Patterning should equal.")
            self.p.print_help()
            return
        common_sequence = self.get_enhanced(ss, pattern)
        self.common_length = len(common_sequence)
        self.search_by_common_pattern(common_sequence)


if __name__ == '__main__':
    searcher = Patterning()
    searcher.start()
