"""
Author: Timothy Tran
Date: 4/9/2024

This script retrieves PubMed article metadata and abstracts then saves them to a CSV file.

It accepts the following flags:
    --query:      The search query term (i.e. 'mental health').
    --start-year: The start year for the search.
    --end-year:   The end year for the search.
    --testing:    Enable testing mode (limits num of articles per year to 20).

Example usage:
    python3 pubmed_scraper.py --query "mental health" --start-year 2000 --end-year 2020 --testing

Contains the following functions:
    build_PMID_list: returns a list of PubMed articles ID's based on the query and year
    search: helper of build_PMID_list
    fetch_details: fetches metadata for a list of article ID's
    build_dataframe: creates a dataframe of metadata for PubMed articles based on the query and year
    create_csv_year: creates a csv file containing metadata of all PubMed articles based on the query and year

Inspired by:
    https://gist.github.com/bonzanini/5a4c39e4c02502a8451d#file-search_biopython-py 
    https://medium.com/@felipe.odorcyk/scrapping-data-from-pubmed-database-78a9b53de8ca

Requires
    biopython
    pandas
    numpy
"""

"""
TODOS:
    Update date (month and year) to better format, add specific date
    Implement full-text retrieval
"""

import argparse
from Bio import Entrez
import pandas as pd
import numpy as np
import time

def build_PMID_list(query, year, is_testing_mode):
    id_list = []

    if is_testing_mode:
        results = search(query, year, is_testing_mode=is_testing_mode)
        return results['IdList']

    # PubMed limits results to 9999, divide up search into four quarters
    for i in range(1, 11, 3): # in three month increments
        mindate = str(year) + '/' + str(i)
        maxdate = str(year) + '/' + str(i + 2) # inclusive
        results = search(query, year, mindate, maxdate)
        id_list += results['IdList']

    return id_list

def search(query, year, mindate=None, maxdate=None, is_testing_mode=False):
    Entrez.email = 'timttran@uw.edu'
    query = query + ' ' + str(year)
    retmax = 20 if is_testing_mode else 10000

    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax=retmax,
                            mindate=mindate,
                            maxdate=maxdate,
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'timttran@uw.edu'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

def build_dataframe(PMID_list):
    title_list= []
    abstract_list=[]
    journal_list = []
    language_list =[]
    pubdate_year_list = []
    pubdate_month_list = []

    chunk_size = 10000
    for chunk_i in range(0, len(PMID_list), chunk_size):
        chunk = PMID_list[chunk_i:chunk_i + chunk_size]
        papers = fetch_details(chunk)
        for i, paper in enumerate (papers['PubmedArticle']):
            title_list.append(paper['MedlineCitation']['Article']['ArticleTitle'])
            try:
                abstract_list.append(paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0])
            except:
                abstract_list.append('No Abstract')
            journal_list.append(paper['MedlineCitation']['Article']['Journal']['Title'])
            language_list.append(paper['MedlineCitation']['Article']['Language'][0])
            try:
                pubdate_year_list.append(paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'])
            except:
                pubdate_year_list.append('No Data')
            try:
                pubdate_month_list.append(paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month'])
            except:
                pubdate_month_list.append('No Data')

    df = pd.DataFrame(list(zip(title_list, 
                            abstract_list, 
                            journal_list, 
                            language_list, 
                            pubdate_year_list, 
                            pubdate_month_list)),
        columns=['Title', 'Abstract', 'Journal', 'Language', 'Year','Month'])
    # TODO: update publication date to be in required format (check medium article for possible guide)
        # and also include day)
    return df

def create_csv_year(query, year, is_testing_mode):
    id_list = build_PMID_list(query, year, is_testing_mode)
    print(f'{len(id_list)} articles in {year}')

    dataframe = build_dataframe(id_list[:20]) if is_testing_mode else build_dataframe(id_list)

    dataframe.to_csv(f'{year}_{query.replace(" ", "_")}.csv')

def main(query, start_year, end_year, is_testing_mode):
    for year in range(start_year, end_year + 1):
        create_csv_year(query, year, is_testing_mode)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process start and end years.")
    parser.add_argument("--query", type=str, help="Query", required=True)
    parser.add_argument("--start-year", type=int, help="Start year", required=True)
    parser.add_argument("--end-year", type=int, help="End year", required=True)
    parser.add_argument("--testing", action="store_true", help="Enable testing mode")
    args = parser.parse_args()
    
    start_year = args.start_year
    end_year = args.end_year
    query = args.query
    is_testing_mode = args.testing

    start_time = time.time()

    main(query, start_year, end_year, is_testing_mode)

    end_time = time.time()
    execution_time = end_time - start_time

    print(f"Execution time: {execution_time:.2f} seconds")