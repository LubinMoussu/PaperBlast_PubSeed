# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:52:42 2018

@author: lubin.moussu
"""

from urllib.parse import urlencode
from urllib.request import Request, urlopen
from pprint import pprint
import re 
from os import startfile
from operator import itemgetter
from itertools import groupby
import os
from re import findall, match
import numpy as np, scipy.stats as st, pandas as pd
from collections import Counter
from html.parser import HTMLParser
from urllib.parse import urlparse, parse_qs
from random import seed, randint, SystemRandom, sample
from webbrowser import open as webbrowser_open
import pickle
from multiprocessing.dummy import Pool as ThreadPool
import itertools
import urllib
from subprocess import Popen
from time import sleep
import operator
from http.client import *

from selenium import webdriver
from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
from selenium.webdriver.common.keys import Keys
import selenium.webdriver.support.ui as ui
from selenium.webdriver.firefox.options import Options

from Similarity_measures import *
from selenium_browser import *

import bs4
try: 
    from BeautifulSoup import BeautifulSoup
except ImportError:
    from bs4 import BeautifulSoup

    
class Protein_Evidence():
    def __init__(self, protein_ID,protein_IDs,genome_ID,evidence,pubmed_IDs):
        self.protein_ID = protein_ID
        self.protein_IDs = protein_IDs
        self.genome_ID = genome_ID
        self.evidence = evidence
        self.pubmed_IDs = pubmed_IDs

def PaperBlast_Query(seq):
        url =  'http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?'
        values = {"query": seq,
                  "Search" : 'Search'
                }

        for key in values:
            assert type(values[key]=="str")
        
        # Encode values in url code
        data = urlencode(values)
        # Create the request
        request = "{0}&{1}".format(url, data)
        opener = urllib.request.build_opener()
        # try:
        #     response = opener.open(request,timeout=60)
        #     return response
        # except urllib.error.URLError, timeout:
        #     return None
        try:
            response = opener.open(request, timeout=60)
            sleep(1)
            return response
        except Exception as ex:
            print('PaperBlast_Query', ex)
            print(seq)
            return None
    
def get_entry_evidences(seq, targeted_sites, genome_exceptions):
    '''
    '''
    assert type(seq)==str
    assert type(genome_exceptions)==list
    page = PaperBlast_Query(seq)
    if page:
        soup = BeautifulSoup(page, 'html.parser')

        evidence_sets = {}
        genomes_features = soup.find_all('p', attrs={'style':"margin-top: 1em; margin-bottom: 0em;"})
        for genome in genomes_features:
            genome_ID = genomes_features.index(genome)

            entry_evidence = {'entry_id': '',
                              'database_sources': {},
                              'accession_code': [],
                              'pubmed_IDs': [],
                              'functional_roles': {},
                              'genome_names': []
                              }
            Accession_Code = ''
            functional_role = ''
            pubmed_IDs = []
            database_name = ''

            database_sources = {site: [] for site in targeted_sites}
            databases_pubmed_IDs = {site: [] for site in targeted_sites}
            databases_functional_role = {site: [] for site in targeted_sites}

            # if genome_ID == 0:
            if genome_ID:
                for child in genome.children:
                    if child.name:
                        if child.name == 'br' and Accession_Code and pubmed_IDs:
                            # add information to entry_evidence corresponding to each database if the entry is a reference protein
                            database_sources[database_name].append(Accession_Code)
                            databases_pubmed_IDs[database_name].extend(pubmed_IDs)
                            databases_functional_role[database_name].append(functional_role)

                            Accession_Code = ''
                            functional_role = ''
                            pubmed_IDs = []

                        # database_entries
                        if child.has_attr('title') and child['title'] in targeted_sites:
                            database_name = str(child['title'])
                            # print(child.contents, child['title'])
                            Accession_Code = child.contents[0]
                            Accession_Code = Accession_Code if '/' not in Accession_Code else Accession_Code.split(' / ')[1]
                            Accession_Code = str(Accession_Code)

                        # functional_role
                        if child.name == 'b':
                            functional_role = str(child.contents[0])
                            # functional_role = str(child.contents[0]).encode('utf-8')

                        # Genome_name
                        if child.name == 'i':
                            genome_name = str(child.string)
                            entry_evidence['genome_names'].append(genome_name)

                        # Pubmed entries
                        if child.has_attr('onmousedown'):
                            pubmed_url = child['href']
                            if pubmed_url and Accession_Code:
                                pubmed_IDs = pubmed_url.lstrip('http://www.ncbi.nlm.nih.gov/pubmed/').split(',')
                                pubmed_IDs = [str(elt) for elt in pubmed_IDs if elt.isdigit()]

            protein_IDs = list(set([protein_ID for site in database_sources for protein_ID in database_sources[site] if
                                    protein_ID[0].isupper()]))

            entry_evidence['pubmed_IDs'] = databases_pubmed_IDs
            entry_evidence['database_sources'] = database_sources
            entry_evidence['accession_code'] = protein_IDs
            entry_evidence['functional_roles'] = databases_functional_role
            if protein_IDs:
                for protein_ID in protein_IDs:
                    if protein_ID not in evidence_sets:
                        if not set(entry_evidence['genome_names']).intersection(genome_exceptions):
                            entry_evidence['entry_id'] = protein_ID
                            evidence_sets[protein_ID] = entry_evidence
        return evidence_sets
   
def query_protein_ID_to_PubSeed(protein_ID):
    url = "http://pubseed.theseed.org//seedviewer.cgi?pattern={0}&page=Find&act=check_search".format(protein_ID)
    # Create the
    request = url
    opener = urllib.request.build_opener()
    try:
        response = opener.open(request, timeout=30)
        sleep(1)
    except urllib.error.URLError:
        return None
    return response

def get_orgid_from_fig(fig):
    org_id = findall('\d+\.\d+', fig)[0]
    return org_id

def get_uniprot_fasta(uniprot_ID):
    url = "https://www.uniprot.org/uniprot/{0}.fasta".format(uniprot_ID)
    request = url
    opener = urllib.request.build_opener()
    try:
        page = opener.open(request, timeout=120)
        sleep(1)
    # except (urllib.error.URLError, http.client.RemoteDisconnected) as error:
    #     return None
    except Exception as ex:
        print('uniprot', uniprot_ID, ex)

    response = [elt.decode("utf-8").rstrip("\n") for elt in page if elt.decode("utf-8").rstrip("\n")]
    seq = ''.join(response[1:])
    ID = response[0]

    OS = findall('OS=(.+) OX=', ID)[0]
    taxon_identifier = findall(' OX=(\d*) ', ID)[0]
    if 'strain' in OS:
        taxo = OS.split(' (')
        genus = taxo[0]
        strains = taxo[1].rstrip(')').lstrip('strain ').split(' / ')
        genome_names = []
        for strain in strains:
            genome_name = '{0} {1}'.format(genus, strain)
            genome_names.append(genome_name)
    else:
        genome_names = [OS]

    uniprot_entry = {'ID':uniprot_ID,
                     'seq':seq,
                     'genome_names':genome_names,
                     'taxon_identifier':taxon_identifier
                     }
    return uniprot_entry

def get_PubSeed_Matches_features_from_ProteinID(protein_ID):

    page = query_protein_ID_to_PubSeed(protein_ID)
    soup = BeautifulSoup(page, 'html.parser')

    matches_features = []
    for link in soup.find_all('input', attrs={'id': "table_data_0"}):
        figs = findall('fig\|\d+\.\d+\.peg.\d+', str(link))
        for fig in figs:
            matches_features.append(fig)
    matches_features = list(set(matches_features))
    return matches_features

def get_PubSeed_genomes_with_Taxon_identifier(all_PubSeed_genomes,
                                              Taxon_identifier):

    matches = [elt for elt in all_PubSeed_genomes if Taxon_identifier in elt]
    return matches

def get_PubSeed_Matches_features_from_Sources_prot_sequences(id_seq,
                                                             Blast_parameters,
                                                             fasta_seq,
                                                             Taxon_identifier,all_PubSeed_genomes):
    '''
    Aim: protein sequence is blasted against PubSeed genome to get the right protein
    :param id_seq:
    :param Blast_parameters:
    :param fasta_seq:
    :param genome_ID:
    :return:
    '''

    matches_genomes = get_PubSeed_genomes_with_Taxon_identifier(all_PubSeed_genomes,Taxon_identifier)
    candidate_figs = []
    for genome_ID in matches_genomes:
        id_seq, peg_evalue, genome_ID = Blast_query_PEG_evalues(id_seq, Blast_parameters, fasta_seq, genome_ID)

        for candidate_fig in peg_evalue:
            if peg_evalue[candidate_fig][0] == 0.0:
                candidate_figs.append(candidate_fig)
    return candidate_figs

def get_PubSeed_Matches_features_from_Sources_IDs(evidence_sets,
                                                  Blast_parameters,
                                                  all_PubSeed_genomes):
    '''
    Aim : Find the PubSeed Fig associated with the protein_ID. Only keep protein_IDs for which PubSeed entries have been found
    :param evidence_sets:
    :param Blast_parameters:
    :param all_PubSeed_genomes:
    :return: evidence_sets
    '''

    for protein_ID in list(evidence_sets.keys()):
        entries = evidence_sets[protein_ID]['accession_code']
        prot = False
        evidence_sets[protein_ID]['matches_PubSeed_features'] = []
        for entry in entries:
            matches_features = get_PubSeed_Matches_features_from_ProteinID(entry)
            if not matches_features:
                # If entry_ID of databases are not found in PubSeed, the protein sequence is retrieved from databases and blast against genomes of PubSeed
                uniprot_entry = get_uniprot_fasta(entry)
                if uniprot_entry['taxon_identifier']:
                    # If there is a identified reference genome, find PubSeed candidates
                    matches_features = get_PubSeed_Matches_features_from_Sources_prot_sequences(entry, Blast_parameters, uniprot_entry['seq'], uniprot_entry['taxon_identifier'], all_PubSeed_genomes)
            if matches_features:
                # Check if genome_ID is within expected genome_IDs in all_PubSeed_genomes
                matches_features = [fig for fig in matches_features if get_orgid_from_fig(fig) in all_PubSeed_genomes]
                if matches_features:
                    evidence_sets[protein_ID]['matches_PubSeed_features'] = matches_features
                    prot = True

        if not prot:
            try:
                del evidence_sets[protein_ID]
            except KeyError:
                pass
    return evidence_sets

def EC_number_similarity(EC1, EC2):

    EC1 = EC1.split('.')
    EC2 = EC2.split('.')

    si = [1 if EC1[i]==EC2[i] else 0 for i in range(len(EC1))]
    print(si)

def uniprot_potential_functional_roles(reference_roles):

    all_uniprot_names = []
    for reference_role in reference_roles:
        names = reference_role.split('; ')
        EC_number_pattern = 'EC (\d+.\d+.[\d-]+.[\d-]+)'
        EC_number = [elt for elt in names if findall(EC_number_pattern,elt)]
        if EC_number:
            EC_number = EC_number[0]
            names = ['{0} ({1})'.format(elt, EC_number) for elt in names if elt!=EC_number and len(elt)>5]
            all_uniprot_names.extend(names)
    return list(set(all_uniprot_names))

def update_functional_role_consistency(current_role_EC,
                                       reference_role_EC,
                                       consistency,
                                       current_role,
                                       reference_role,
                                       jaccard_threshold,
                                       EC_number_sufficiency):

    if current_role_EC and reference_role_EC:
        reference_role_EC = reference_role_EC[0] if reference_role_EC else ''
        if reference_role_EC == current_role_EC:
            # If both EC numbers are equal,
            if EC_number_sufficiency:
                consistency = True
            else:
                jaccard_score = jaccard_similarity(current_role, reference_role)
                if jaccard_score >= jaccard_threshold:
                    consistency = True
    return consistency

def functional_role_consistency(fig_information, evidence_set, jaccard_threshold, EC_number_sufficiency):
    '''
    Aim: check that the functional role fits ontology systems role2
    :param role1:
    :param role2:
    :return:
    '''
    EC_number_pattern = 'EC (\d+.\d+.[\d-]+.[\d-]+)'
    current_role = fig_information['Function']
    current_role_EC = findall(EC_number_pattern,current_role)
    current_role_EC = current_role_EC[0] if current_role_EC else ''
    reference_roles = evidence_set['functional_roles']

    consistency = False
    for site in reference_roles:
        if not consistency:
            if site == 'SwissProt':
                list_reference_roles = uniprot_potential_functional_roles(reference_roles[site])
            else:
                list_reference_roles = list(reference_roles[site])

            for reference_role in list_reference_roles:
                if not consistency:
                    reference_role_EC = findall(EC_number_pattern, reference_role)
                    consistency = update_functional_role_consistency(current_role_EC,
                                                                       reference_role_EC,
                                                                       consistency,
                                                                       current_role,
                                                                       reference_role,
                                                                       jaccard_threshold,
                                                                       EC_number_sufficiency)
    return consistency

def curate_literature_to_fig(fig, current_pubmed_IDs, entry_pubmed_IDs, browser):

    # if the PubSeed fig is not associated with any pubmed_ID which however should be associated, add the pubmed_ID to PubSeed_ID
    if not current_pubmed_IDs:
        # Define pubmed_ID to add to the fig
        pubmed_IDs_to_add = [pubmed_ID for pubmed_ID in entry_pubmed_IDs if pubmed_ID not in current_pubmed_IDs]

        for pubmed_ID in pubmed_IDs_to_add:
            add_pubmedID_to_fig_selenium(fig, pubmed_ID, browser)
            sleep(3)
            break

def update_literature_PubSeed(evidence_sets, EC_number_sufficiency):

    '''
    Aim : Automatically update literature in PubSeed for figs from evidence_sets
    :param evidence_sets:
    :return: figs for which current_role is different from reference_role
    '''

    # init_headless_firefox_browser
    browser = init_headless_firefox_browser()

    # init_pubseed_login
    ID, password =  'moussulubin', 'E2UGSat8z<'
    browser = init_pubseed_login(browser, ID, password)
    sleep(3)
    print('Connected to PubSeed')

    features_to_curate = {}
    for protein_ID in list(evidence_sets.keys()):
        evidence_set = evidence_sets[protein_ID]
        entry_pubmed_IDs = list(set([pubmed_ID for site in evidence_set['pubmed_IDs'] for pubmed_ID in evidence_set['pubmed_IDs'][site]]))

        if protein_ID and evidence_set['matches_PubSeed_features'] and entry_pubmed_IDs:
            # if the protein_ID has an entry in PubSeed and an evidence article, check if the article is assigned to PubSeed entry
            for fig in evidence_set['matches_PubSeed_features']:
                all_fig_information = get_pubseed_fig_information(fig)
                # Curate literature of the fig
                curate_literature_to_fig(fig, all_fig_information[fig]['Ev'], entry_pubmed_IDs, browser)

                # check if the current annotation is consistent with reference annotations. If not, trigger manual curation
                consistency = functional_role_consistency(all_fig_information[fig], evidence_sets[protein_ID], 0.6, EC_number_sufficiency)
                if not consistency:
                    if protein_ID not in features_to_curate:
                        features_to_curate[protein_ID] = []
                    if fig not in features_to_curate[protein_ID]:
                        features_to_curate[protein_ID].append(fig)
    browser.quit()
    return features_to_curate

def get_multiple_entry_evidences(ref_fig_list, targeted_sites, genome_exceptions):

    multiple_evidence_sets = {}
    fig_list = []
    for elt in ref_fig_list:
        seq = elt['Seq']
        fig_list.append(elt['Fig'])
        evidence_sets = get_entry_evidences(seq, targeted_sites, genome_exceptions)
        if evidence_sets:
            multiple_evidence_sets.update(evidence_sets)
    return multiple_evidence_sets, fig_list

def get_multithreading_multiple_entry_evidences(ref_fig_list, targeted_sites, genome_exceptions, workers=4):

    # make the Pool of workers
    pool = ThreadPool(workers)
    # and return the results

    fig_list = []
    seq_list = []
    for elt in ref_fig_list:
        seq_list.append(elt['Seq'])
        fig_list.append(elt['Fig'])

    args = zip(seq_list,
               itertools.repeat(targeted_sites),
               itertools.repeat(genome_exceptions))

    results = pool.starmap(get_entry_evidences, args)
    # close the pool and wait for the work to finish
    pool.close()
    pool.join()

    # Process multithreading process
    multiple_evidence_sets = {}
    for evidence_sets in results:
        if evidence_sets:
            multiple_evidence_sets.update(evidence_sets)
    return multiple_evidence_sets, fig_list

def get_arborescence_PubSeed_Matches_features_from_Sources_IDs(ref_fig_list, targeted_sites, genome_exceptions, Blast_parameters, all_PubSeed_genomes, depth, workers):

    iter = 0
    multiple_arborescence_evidence_sets = {}
    already_explored_fig_list = []
    while iter < depth:
        evidence_sets, fig_list = get_multithreading_multiple_entry_evidences(ref_fig_list, targeted_sites, genome_exceptions, workers)
        already_explored_fig_list.extend(fig_list)

        print('{0} evidence_sets obtained'.format(len(evidence_sets)))

        # get PubSeed Matches only for protein_ID not present in multiple_arborescence_evidence_sets
        evidence_sets = {key:value for key,value in evidence_sets.items() if key not in multiple_arborescence_evidence_sets}
        evidence_sets = get_PubSeed_Matches_features_from_Sources_IDs(evidence_sets, Blast_parameters, all_PubSeed_genomes)
        print('{0} reference external links associated with PubSeed features'.format(len(evidence_sets.keys())))

        # update multiple_arborescence_evidence_sets
        multiple_arborescence_evidence_sets.update(evidence_sets)

        # define new fig to find papers about it and its homologs
        new_ref_fig_list = [fig for protein_ID in evidence_sets for fig in evidence_sets[protein_ID]['matches_PubSeed_features']]
        new_ref_fig_list = [fig for fig in new_ref_fig_list if fig not in ref_fig_list and fig not in already_explored_fig_list]
        ref_fig_list = get_multiple_fasta_information(new_ref_fig_list)

        iter+=1
    return multiple_arborescence_evidence_sets

if __name__ == "__main__":
    import time
    from Multi_project_scripts import *
    start = time.time()
    print("Start")

    from sys import getrecursionlimit, setrecursionlimit
    getrecursionlimit()
    max_rec = 0x10000000
    setrecursionlimit(max_rec)

    targeted_sites = ['SwissProt', 'EcoCyc', 'MetaCyc', 'BRENDA', 'MicrobesOnline']
    # targeted_sites = ['SwissProt', 'MetaCyc', 'BRENDA', 'MicrobesOnline']

    all_PubSeed_genomes = get_csv_reader('Inputs/All_PubSeed_organisms.txt', delimiter=";")[0] # Only contains Bacteria + Archaea
    print('all_PubSeed_genomes', len(all_PubSeed_genomes))

    fig_list = ['fig|511145.12.peg.3012', 'fig|511145.12.peg.3013', 'fig|511145.12.peg.3014', 'fig|511145.12.peg.3015',
                'fig|511145.12.peg.1055', 'fig|511145.12.peg.3212', 'fig|511145.12.peg.3211', 'fig|511145.12.peg.3210',
                'fig|511145.12.peg.3209', 'fig|511145.12.peg.3208', 'fig|511145.12.peg.3207', 'fig|511145.12.peg.3206',
                'fig|511145.12.peg.3205', 'fig|511145.12.peg.3204', 'fig|99287.12.peg.2179', 'fig|907.4.peg.1106',
                'fig|907.4.peg.826', 'fig|907.4.peg.1480', 'fig|511145.12.peg.4192', 'fig|511145.12.peg.4191',
                'fig|511145.12.peg.4190', 'fig|511145.12.peg.339', 'fig|511145.12.peg.340', 'fig|511145.12.peg.341',
                'fig|511145.12.peg.342', 'fig|511145.12.peg.343', 'fig|511145.12.peg.344'
                ]
    fig_list = ['fig|511145.12.peg.3012', 'fig|511145.12.peg.3013', 'fig|511145.12.peg.3014', 'fig|511145.12.peg.3015']
    ref_fig_list = get_multiple_fasta_information(fig_list)

    # f_out = 'Inputs/my_seqs_to_reference.fasta'
    # write_fasta_sequences(list_figs, f_out)

    # # Get fasta sequences to ref in PubSeed
    # my_seqs_to_reference = 'Inputs/my_seqs_to_reference.fasta'
    # # fasta_to_ref = read_fasta(my_seqs_to_reference)
    # # a = [elt for elt in fasta_to_ref]
    # # pprint(a)

    # Define parameters of the Blast
    cutoff = 1e-05
    sequence = ''
    Blast_parameters = {"seq_type":"aa",
              "filter" : "F",
              "evalue" : str(cutoff),
              "wsize" : "0",
              "fasta" : sequence,
              "organism" : '',
              "act" : "BLAST"
            }

    genome_exceptions = ['Homo sapiens (Human)', 'Equus caballus (Horse)', 'Mus musculus (Mouse)', 'Bos taurus (Bovine)', 'Oryza sativa subsp. japonica (Rice)', 'Arabidopsis thaliana (Mouse-ear cress)']
    depth = 1
    workers = len(fig_list)
    if False:
        evidence_sets = get_arborescence_PubSeed_Matches_features_from_Sources_IDs(ref_fig_list, targeted_sites, genome_exceptions,
                                                                   Blast_parameters, all_PubSeed_genomes, depth, workers)
        pickle.dump(evidence_sets, open('Inputs/pickle/evidence_sets.pkl', "wb"), protocol=pickle.HIGHEST_PROTOCOL)
        print('{0} dump to {1}'.format('evidence_sets', 'Inputs/pickle/evidence_sets.pkl'))

    else:
        evidence_sets = pickle.load(open('Inputs/pickle/evidence_sets.pkl', "rb"))

    pprint(evidence_sets)
    if False:
        EC_number_sufficiency = True
        curated_features = update_literature_PubSeed(evidence_sets, EC_number_sufficiency)
        pickle.dump(curated_features, open('Inputs/curated_features.pkl', "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    else:
        curated_features = pickle.load(open('Inputs/curated_features.pkl', "rb"))
        pprint(curated_features)

    end = time.time()
    print(end - start)
    print ("Done")
