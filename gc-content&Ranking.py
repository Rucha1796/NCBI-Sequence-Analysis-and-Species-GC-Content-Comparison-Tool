# -*- coding: utf-8 -*-
# NAME: cs223_prog_exercise_2.py

"""
AUTHOR: <Rucha Deo>

  ============== VARIABLE, FUNCTION, etc. NAMING CONVENTIONS ==================
<ALL CAPITOL LETTERS>:  Indicates a symbol defined by a
        #define statement or a Macro.

   <Capitalized Word>:  Indicates a user defined global var, fun, or typedef.

   <all small letters>:  A variable or built in functions.
"""

# IMPORTS
import os
import re
from copy import copy
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import GC
import xmltodict

# ENTREZ LIST OF DBs
""""
{'DbList': ['pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide', 'structure', 'genome',
            'annotinfo', 'assembly', 'bioproject', 'biosample', 'blastdbinfo', 'books',
            'cdd', 'clinvar', 'gap', 'gapplus', 'grasp', 'dbvar', 'gene', 'gds', 'geoprofiles',
            'homologene', 'medgen', 'mesh', 'ncbisearch', 'nlmcatalog', 'omim', 'orgtrack', 'pmc',
            'popset', 'proteinclusters', 'pcassay', 'protfam', 'pccompound', 'pcsubstance', 
            'seqannot', 'snp', 'sra', 'taxonomy', 'biocollections', 'gtr']}
"""
# HELPER FUNS
def get_valid_ids(list_of_ids, db_name, return_type, query_type):
    """Validate each Entrez ID# that is in the list_of_ids  argument. Invalid
    ids are removed from the list. A valid id# is one that does not generate an error
    for the given DB, return type, and query type."""

    # Return if there are no ids to check
    if len(list_of_ids) <= 0:
        return list_of_ids

    # Iterate through list and check each id
    copy_list_of_ids = copy(list_of_ids)    # Make a copy of the list of ids that will be modified.
    for lid in list_of_ids:
        if query_type == 'efetch':
            try:
                Entrez.efetch(db = db_name, id = int(lid), rettype = return_type)
            except:
                print("get_valid_ids: The Entrez.efetch id# " + str(lid) + " for NCBI " + str(db_name) +
                      " DB, and return type = " + str(return_type) + " is not valid.")
                copy_list_of_ids.remove(lid)

    # When we reach here, copy_list_of_ids will only have valid ids
    return copy_list_of_ids

def rank_organisms_by_gc(gc_content, species_gc_content):
    rankings = []
    for species, gc_range in species_gc_content.items():
        gc_midpoint = sum(gc_range)/2
        closeness = abs (gc_content-gc_midpoint)
        rankings.append((species, closeness))
    rankings.sort(key=lambda x:x[1])
    return [species for species, closeness in rankings]

species_gc_content = {
    'human': (36,60),
    'house mouse': (42,66),
    'Norway rat': (51,58),
    'dog': (17,45),
    'chimpanzee':(35,45)
}

# MAIN FUNCTION
def main():
    """The main program collects """

    # Init Entrez email addr to use
    Entrez.email = "rucha.deo@sjsu.edu"

    accession_numbers = ["NC_000007.14", "NC_051339.1","NC_000072.7", "NC_051820.1","NC_036886.1"]
  # Extract Entrez info for the specified ACC number
    for acc_num in accession_numbers:
        print("\n\n#####")
        #acc_num = "NC_000007.14"
        #acc_num = 'NC_051339'
        print("PROCESSING ACC #: ", acc_num)

    # Entrez.esearch for info related to the current accession number
    #handle = Entrez.esearch(db='gene', term='human[organism] ' + str(acc_num))
        handle = Entrez.esearch(db='gene', term = str(acc_num))

    # Entrez.read the result returned from NCBI search for the current accession number
        handle_read = Entrez.read(handle)

    # Get Idlist from previous read of accession number info
        handle_read_ids = handle_read['IdList']
        print("handle_read_ids: ", handle_read_ids)

    # Make sure IDs in Idlist are valid. Return a list of only valid IDs.
        valid_handle_read_ids = get_valid_ids(handle_read_ids, 'nucleotide', 'fasta', 'efetch')
        if valid_handle_read_ids == []:
            print("   epigen_pipeline_get_entrez_info: No valid IDs found for ACC #" + str(acc_num) + ".")
            return

    # If we get here, valid_handle_read_ids has usable ID nums.
        # for acc_handle_id in valid_handle_read_ids:
        acc_handle_id = valid_handle_read_ids[0]
        print('accsession handle id:', acc_handle_id)
        # Entrez.efetch based on ID#
        handle_read_id_info = Entrez.efetch(db="nucleotide", id=int(acc_handle_id), rettype="fasta")

        # SeqIO.read  fasta info
        record_seq = SeqIO.read(handle_read_id_info, "fasta")
            
        #GC Content
        gc_content = GC(record_seq.seq)

        organism_rankings = rank_organisms_by_gc(gc_content, species_gc_content)

       
            # Print the fasta sequences
        print()
            # Print sequence ACC number and sequence
        print("record_seq.id: ", record_seq.id)
        print("record_seq.seq: ", record_seq.seq)
        print("GC percent%: ", gc_content)
        print("Organism rankings: ", organism_rankings)
            
        print("==========")


    print()
    print ("############### DONE ###############")
    return


#####################################################################################

