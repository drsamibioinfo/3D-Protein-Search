#Protein Propensities

##introduction
This project contains two subprojects, 
### The first project
the first project is a script that downloads the entire PDB database and calculates amino acids propensities in all proteins along with their secondary structure prediction and also it calculates the binary patterning based on structural 
hydrophobicity threshold which was set to 23.9 then, the script builds a database that enable fast ad hoc protein searches against structural fragments/motifs based solely on the structure not based on sequence similarity.

### The second project
The second project is a search engine that allows researchers and protein engineers to find proteins 
that harbors exact and/or fuzzy similar structural scaffolds to a query structural fragment.

## How to install

1. Download , install and start clickhouse database server.
2. Create a new python virtualenv `virtualenv --python=python2.7 <Your_virtual_env>` 
3. Activate the virtualenv `source <Your_virtual_env_location>/bin/activate`
4. Install all python dependencies in requirements.txt `cat requirements.txt | xargs -n1 pip install`
5. Perform Exact or fuzzy protein searches


## How to use

6. Run `python patterning.py --help`
```usage: patterning.py [-h] [-w WORKDIR] [-q PDBS] [-d DATABASE] [-t HOST]
                     [-r PORT] [-u USERNAME] [-o PASSWORD] [-s SOURCE]
                     [-p PATTERN] [-e STRUCTURE] [-m RMSD] [-l LIMIT]
                     [-x START] [-y END] [-a CHAIN] [-c CUTOFF] [-f [FUZZY]]
                     [-z FUZZYLEVEL] [-v DISTANCE] [-j DELETIONS]
                     [-n INSERTIONS] [-k SUBSTITUTIONS]

Patterning is a program that allows to perform de-novo protein design by
searching for possible proteins or parts of protein which match a given 3D
structure or a possible 1D secondary structure with possible binary patterning
fingerprint

optional arguments:
  -h, --help            show this help message and exit
  -w WORKDIR, --workdir WORKDIR
                        The working directory for this script to save and load
                        temporary. By default it takes the current working
                        directory. memory files required to function properly.
  -q PDBS, --pdbs PDBS  Local repository containing PDB Files to read from
  -d DATABASE, --database DATABASE
                        Database Name to use
  -t HOST, --host HOST  Database Host to use. defaults to localhost
  -r PORT, --port PORT  Database Default port to use
  -u USERNAME, --username USERNAME
                        Username of the database to use
  -o PASSWORD, --password PASSWORD
                        Database password to use to login
  -s SOURCE, --source SOURCE
                        Source PDB which contains the 3D structure of a
                        protein or part of a protein (Fragment), that you want
                        to search for similars in the database.
  -p PATTERN, --pattern PATTERN
                        Binary Pattern of a fragment
  -e STRUCTURE, --structure STRUCTURE
                        Possible Secondary structure of the fragment/protein,
                        multiples of (H,B,E,G,I,T,S,-).
  -m RMSD, --rmsd RMSD  Return matched proteins with this RMSD value only and
                        exclude any others.
  -l LIMIT, --limit LIMIT
                        Total Number of results to include. Defaults to 10.
  -x START, --start START
                        When searching by a whole protein containing a
                        specific scaffold. This parameter denotes the location
                        of the first residue of the fragment.
  -y END, --end END     When searching by a whole protein containing a
                        specific scaffold. This parameter denotes the location
                        of the last residue of the fragment.
  -a CHAIN, --chain CHAIN
                        Chain Identifier if your start and end are relative to
                        a particular chain.
  -c CUTOFF, --cutoff CUTOFF
                        Cutoff value to use to mark a residue Solvent-
                        accessible surface area as polar or hydrophobic, if
                        the SASA of a residue equal to or higher than this
                        value it will be considered polar otherwise it will be
                        marked as hydrophobic
  -f [FUZZY], --fuzzy [FUZZY]
                        Perform Fuzzy Search. Allow matching similar but not
                        identical results.
  -z FUZZYLEVEL, --fuzzylevel FUZZYLEVEL
                        Include only results with equal or higher value than
                        the following fuzziness level .Defaults to 0.90
  -v DISTANCE, --distance DISTANCE
                        Possible Lovenstein distance for fuzzy string search
  -j DELETIONS, --deletions DELETIONS
                        Number of allowed string deletions when fuzzy string
                        search is enabled. Defaults to Zero.
  -n INSERTIONS, --insertions INSERTIONS
                        Number of allowed string insertions when fuzzy string
                        search is enabled. Defaults to Zero.
  -k SUBSTITUTIONS, --substitutions SUBSTITUTIONS
                        Number of allowed string substitutions when fuzzy
                        string search is enabled. Defaults to Zero.
```

 - Most often, Protein engineers have a specific protein structural motif that represents an antigenic fragment that
would like to find other proteins that harbor the same structural motif irrespective of sequence similarity, because in protein science, 
sequence similarity not always imply structural similarity.

- You have a PDB of a protein of interest that harbors the structural fragment/motif and you want to find other proteins contain similar structural scaffolds to this one, 
Then you can perform searches like so,

`python patterning.py --source=<Where_is your PDB>/yourpdb.pdb --chain=<Chain_ID> --start=<Start position of your fragment> --end=<End position of your fragment>`

- You can limit the search results by appending `--limit=###` to a number.

- Moreover, You can perform Fuzzy searches using the following command

`python patterning.py --source=<Where_is your PDB>/yourpdb.pdb --chain=<Chain_ID> --start=<Start position of your fragment> --end=<End position of your fragment> --fuzzy --fuzzylevel=0.50 --distance=2`

Here `--fuzzylevel` is the ngram string similarity percentage > 50% and `--distance` is the lovenstein distance between two fuzzy strings

- Moreover, you supply `--deletions` `--insertions` and `--substitutions` to the fuzzy string search engine that suit your needs.


 

