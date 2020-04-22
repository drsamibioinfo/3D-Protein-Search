#!/usr/bin/python
import os, sys, argparse
from collections import OrderedDict
import pickle


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
        self.p.add_argument("-u", "--username", help="Username of the database to use")
        self.p.add_argument("-p", "--password", help="Database password to use to login")

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

    def save_db(self):
        db_path = os.path.join(self.args.workdir, "{0}.pdb".format(self.args.database))
        if os.path.exists(db_path):
            print("Removing Old and Existing Database")
            os.remove(db_path)
        print("Saving Recent Memory Database onto Desk")
        with open(db_path, mode='wb') as output:
            pickle.dump(self.db, output)
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

    def parse_files(self, root, files):
        pass


if __name__ == '__main__':
    manager = PropensityManager()
    manager.start()
