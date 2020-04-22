#!/usr/bin/python
import os, sys, argparse, traceback, urllib2
from collections import OrderedDict
import pickle
from multiprocessing import Pool
from infi.clickhouse_orm.database import Database
from infi.clickhouse_orm.models import Model
from infi.clickhouse_orm.fields import *
from datetime import date
from Bio.PDB import PDBParser, DSSP
from Bio.PDB.Polypeptide import PPBuilder

HELICES = ["H", "G", "I"]
BETA_STRANDS = ["E"]
LOOPS = ["-", "T"]
BENDS = ["S", "B"]

LBL_HELIX = "HELIX"
LBL_SHEETS = "SHEET"
LBL_LOOPS = "LOOPS"
LBL_BENDS = "BENDS"

def get_type(l):
    if l in HELICES:
        return "H"
    elif l in BETA_STRANDS:
        return "S"
    elif l in LOOPS:
        return "L"
    else:
        return "T"


class ProteinModel(Model):
    protein_id = StringField()
    size = UInt64Field()
    helices = UInt64Field()
    sheets = UInt64Field()
    loops = UInt64Field()
    bends = UInt64Field()
    sequence = StringField()
    ss = StringField()
    pattern = ArrayField(inner_field=UInt64Field())
    processdate = DateField(default=date.today())

    @classmethod
    def table_name(cls):
        return "Proteins"


class StructureModel(Model):
    protein_id = StringField()
    ss_id = StringField()
    type = StringField()
    pos = UInt64Field()
    size = UInt64Field()
    sequence = StringField()
    residues = ArrayField(inner_field=UInt64Field())
    sasa = ArrayField(inner_field=Float64Field())
    pattern = ArrayField(inner_field=UInt64Field())
    processdate = DateField(default=date.today())

    @classmethod
    def table_name(cls):
        return "structures"


class MemoryDatabase(dict):
    def __init__(self, seq=None, **kwargs):
        if seq:
            super(MemoryDatabase, self).__init__(seq, **kwargs)
        else:
            super(MemoryDatabase, self).__init__()

    def contains(self, k):
        return k in self.keys()

    def size(self):
        return len(self.keys())

    def __getattr__(self, item):
        if isinstance(self[item], dict):
            return MemoryDatabase(self[item])
        else:
            return self[item]

    def __setattr__(self, key, value):
        self[key] = value


# borrowed from SARI Sabban <http://www.github.com/sarisabban>
# These numbers represent different residues molecular weights
residues = OrderedDict({'A': 129, 'P': 159, 'N': 195, 'H': 224,
                        'V': 174, 'Y': 263, 'C': 167, 'K': 236,
                        'I': 197, 'F': 240, 'Q': 225, 'S': 155,
                        'L': 201, 'W': 285, 'E': 223, 'T': 172,
                        'M': 224, 'R': 274, 'G': 104, 'D': 193})


class PropensityManager(object):
    def __init__(self):
        self.args = None
        self.db = MemoryDatabase()  # empty in memory database to use
        self.p = argparse.ArgumentParser(description="This program will calculate different residues propensities in "
                                                     "3D structures of proteins to extract statistical knowledge and "
                                                     "their distribution in different proteins structural motifs")
        self.p.add_argument("-l", "--location", type=str, default=None,
                            help="Directory which contains all PDBs to include in this study")
        self.p.add_argument("-w", "--workdir", default=".",
                            help="The working directory for this script to save and load temporary. By default it "
                                 "takes the current working directory. "
                                 "memory files required to function properly.")
        self.p.add_argument("-d", "--database", default="propensities", help="Database Name to use")
        self.p.add_argument("-s", "--host", default="localhost", help="Database Host to use. defaults to localhost")
        self.p.add_argument("-r", "--port", default=8123, help="Database Default port to use")
        self.p.add_argument('-c', '--cutoff', default=23.9, type=float, help="Cutoff value to use to mark a residue "
                                                                             "Solvent-accessible surface area as polar or "
                                                                             "hydrophobic, if the SASA of a residue equal "
                                                                             "to or higher than this value it will be "
                                                                             "considered polar otherwise it will be "
                                                                             "marked as hydrophobic")
        self.p.add_argument("-u", "--username", help="Username of the database to use")
        self.p.add_argument("-p", "--password", help="Database password to use to login")
        self.p.add_argument("-m", "--multiprocessing", nargs="?", const=False, default=False,
                            help="Turn multiprocessing on")

    def load_db(self):
        print("Loading Existing Database")
        db_path = os.path.join(self.args.workdir, "{0}.pdb".format(self.args.database))
        if not os.path.exists(db_path):
            print("Creating Fresh Database")
            self.db = MemoryDatabase()
        else:
            with open(db_path, mode='rb') as db_file:
                self.db = MemoryDatabase(pickle.load(db_file))
            print("Existing Database has been loaded.")

    def __insert__(self, models=[]):
        try:
            db_url = "http://{0}:{1}".format(self.args.host or "localhost", self.args.port or 8123)
            db = Database(db_name=self.args.database, db_url=db_url, username=self.args.username,
                          password=self.args.password)
            ctime = datetime.datetime.now()
            db.insert(models)
            took = datetime.datetime.now() - ctime
            print("{0} Records have been written into clickhouse , Took : {1}".format(len(models), str(took)))
        except Exception as e:
            print (e)
            print (traceback.format_exc())
            raise e

    def save_db(self):
        db_path = os.path.join(self.args.workdir, "{0}.db".format(self.args.database))
        if os.path.exists(db_path):
            print("Removing Old and Existing Database")
            os.remove(db_path)
        print("Saving Recent Memory Database onto Desk")
        with open(db_path, mode='wb') as output:
            pickle.dump(dict(self.db), output)
        print("Done saving the memory database onto Desk")

    def start(self):
        try:
            if len(sys.argv) <= 1:
                self.p.print_help()
                return
            self.args = self.p.parse_args()
            if self.args.location is None or len(self.args.location) == 0:
                self.p.print_help()
                return
            self.load_db()
            for root, dirs, files in os.walk(self.args.location):
                if len(files) <= 0:
                    continue
                self.parse_files(root, files)

        except Exception as e:
            print (e.message)
            self.save_db()
        finally:
            self.save_db()

    def parse_files(self, root, files):
        # Process multiple PDB files simultaneously on multiple cores
        if self.args.multiprocessing:
            with Pool() as p:
                for file in files:
                    if not str(file).endswith('pdb'):
                        continue
                    pdb_file = os.path.join(root, file)
                    p.map(self.process_single_pdb, pdb_file)
        else:
            for file in files:
                if not str(file).endswith('pdb'):
                    continue
                pdb_file = os.path.join(root, file)
                self.process_single_pdb(pdb_file)

    def process_single_pdb(self, pdb_file):
        print("Processing File : {0}".format(pdb_file))
        name, sequence, ss, sasa = self.get_secondary_structure_details(pdb_file)
        if len(sequence) != len(ss):
            return
        protein = ProteinModel()
        protein.protein_id = name
        protein.size = len(sequence)
        protein.helices = sequence.count('H') + sequence.count('G') + sequence.count('I')
        protein.sheets = sequence.count('E')
        protein.loops = sequence.count('-') + sequence.count('T')
        protein.bends = sequence.count('S') + sequence.count('B')
        protein.sequence = str(sequence)
        protein.ss = ss
        protein.pattern = [1 if a >= self.args.cutoff else 0 for a in sasa]
        self.__insert__(models=[protein])
        self.db[name] = {
            "protein": protein.to_dict(),
            "structures": {}
        }
        index = 1
        # create secondary structures for this protein
        for keys , structure , solvents in self.prepare_structures(name,sequence,ss,sasa):
            newStructure = StructureModel()
            newStructure.protein_id = name
            newStructure.ss_id = "{0}_{1}".format(name,index)
            type = get_type(structure.values()[0])
            newStructure.type = type
            newStructure.pos = keys[0]
            newStructure.size = len(keys)
            seq = sequence[keys[0]:keys[-1]+1]
            newStructure.sequence = str(seq)
            frequencies = [0] * len(residues.keys())
            for res in residues.keys():
                frequencies[residues.keys().index(res)] = str(seq).count(res)
            newStructure.residues = frequencies
            newStructure.sasa = solvents.values()
            newStructure.pattern = [1 if a >= self.args.cutoff else 0 for a in solvents.values()]
            if type in self.db[name]['structures'].keys():
                structures_list = self.db[name]['structures'][type]
                structures_list.append(newStructure.to_dict())
            else:
                self.db[name]['structures'][type] = [newStructure.to_dict()]
            self.__insert__(models=[newStructure])

            index += 1

    def get_secondary_structure_details(self, pdb_file):
        parser = PDBParser()
        base_name = os.path.basename(pdb_file)
        name, _ = os.path.splitext(base_name)
        structure = parser.get_structure(name, pdb_file)
        dssp = DSSP(structure[0], pdb_file, acc_array="Wilke")
        ss = "".join([aa[2] for aa in dssp])
        sasa = [residues[aa[1]] * aa[3] for aa in dssp]
        builder = PPBuilder()
        seq = ""
        for chain in builder.build_peptides(structure):
            seq += chain.get_sequence()
        return name, seq, ss, sasa

    def prepare_structures(self, name, sequence, ss, sasa):
        first = ss[0]
        structure = OrderedDict()
        solvents = OrderedDict()
        structure[0] = first
        solvents[0] = sasa[0]
        for i in range(1,len(ss)):
            next = ss[i]
            if next == first:
                structure[i] = next
                solvents[i] = sasa[i]
                continue
            else:
                sorted_keys = sorted(structure.keys())
                yield sorted_keys , structure , solvents
                structure.clear()
                solvents.clear()
                first = next
                structure[i] = next
                solvents[i] = sasa[i]

        if len(structure.keys()) > 0:
            sorted_keys = sorted(structure.keys())
            yield sorted_keys , structure , solvents
            structure.clear()
            solvents.clear()




if __name__ == '__main__':
    manager = PropensityManager()
    manager.start()
