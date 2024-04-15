"""
Author: Timothy Tran
Date: 4/9/2024

This script retrieves PubMed article metadata and abstracts then saves them to a CSV file.

It accepts the following flags:
    --start-year: The start year for the search.
    --end-year:   The end year for the search.
    --testing:    Enable testing mode (limits num of articles per year to 20).
and a txt file where each line represents the search terms

Example usage:
    python3 pubmed_scraper.py --start-year 2000 --end-year 2020 --testing

Contains the following functions:
    build_PMID_list: returns a list of PubMed articles ID's based on the query and year
    search: helper of build_PMID_list
    fetch_details: fetches metadata for a list of article ID's
    build_dataframe: creates a dataframe of metadata for PubMed articles based on the query and year
    create_csv_year: creates a csv file containing metadata of all PubMed articles based on the query and year
    fomrat_date: given three strings (day, month, year) formats it into MM/DD/YY

Inspired by:
    https://gist.github.com/bonzanini/5a4c39e4c02502a8451d#file-search_biopython-py 
    https://medium.com/@felipe.odorcyk/scrapping-data-from-pubmed-database-78a9b53de8ca

Requires
    biopython
    pandas
"""

"""
TODOS:
    Update date (month and year) to better format, add specific date
    Implement full-text retrieval
"""

import argparse
from Bio import Entrez
import pandas as pd
import time
import os

def build_PMID_list(queries, year, is_testing_mode):
    PMID_query_mapping = {}

    if is_testing_mode:
        for query in queries:
            results = search(query, year, is_testing_mode=is_testing_mode)
            for PMID in results['IdList']:
                if PMID not in PMID_query_mapping:
                    PMID_query_mapping[PMID] = []
                PMID_query_mapping[PMID].append(query)
        return PMID_query_mapping

    # PubMed limits results to 9999, divide up search into four quarters
    for query in queries:
        for i in range(1, 11, 3): # in three month increments
            mindate = str(year) + '/' + str(i)
            maxdate = str(year) + '/' + str(i + 2) # inclusive
            results = search(query, year, mindate, maxdate)
            for PMID in results['IdList']:
                if PMID not in PMID_query_mapping:
                    PMID_query_mapping[PMID] = []
                PMID_query_mapping[PMID].append(query)

    return PMID_query_mapping

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

def build_dataframe(PMID_query_mapping):
    pmid_list = []
    title_list = []
    abstract_list = []
    journal_list = []
    language_list = []
    pubdate_list = []
    search_terms_list = [] # The search terms that resulted in this article

    PMID_list = list(PMID_query_mapping.keys())

    chunk_size = 10000
    for chunk_i in range(0, len(PMID_list), chunk_size):
        chunk = PMID_list[chunk_i:chunk_i + chunk_size]
        papers = fetch_details(chunk)
        for i, paper in enumerate (papers['PubmedArticle']):
            pmid = paper['MedlineCitation']['PMID']
            pmid_list.append(pmid)
            title_list.append(paper['MedlineCitation']['Article']['ArticleTitle'])
            try:
                abstract_list.append(paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0])
            except:
                abstract_list.append('No Abstract')
            journal_list.append(paper['MedlineCitation']['Article']['Journal']['Title'])
            language_list.append(paper['MedlineCitation']['Article']['Language'][0])
            
            # Get date of publication
            day = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Day')
            month = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Month')
            year = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year')
            pubdate_list.append(format_date(day, month, year))

            search_terms_list.append(PMID_query_mapping[pmid])

    df = pd.DataFrame(list(zip(pmid_list,
                               title_list, 
                               abstract_list, 
                               journal_list, 
                               language_list, 
                               pubdate_list,
                               search_terms_list)),
        columns=['PMID', 'Title', 'Abstract', 'Journal', 'Language', 'PubDate', 'SearchTerms'])
    return df

def format_date(day, month, year):
    # Mapping of month abbreviations to numeric representation
    month_map = {
        'jan': '01', 'feb': '02', 'mar': '03', 'apr': '04',
        'may': '05', 'jun': '06', 'jul': '07', 'aug': '08',
        'sep': '09', 'oct': '10', 'nov': '11', 'dec': '12'
    }

    # Convert month abbreviation to numeric representation
    # Note: sometimes month is stored as "mar-apr", in this case we always get the first month
    month_numeric = month_map.get(month.lower()[:3], 'NA') if month is not None else 'NA'

    # Format the date string
    return f"{month_numeric}/NA/{year}" if day is None else f"{month_numeric}/{day}/{year}"

def create_csv_year(query, year, is_testing_mode):
    id_list = build_PMID_list(query, year, is_testing_mode)
    print(f'{len(id_list)} articles in {year}')

    dataframe = build_dataframe(id_list)

    dataframe.to_csv(f'pubmed_{year}.csv')# _{query.replace(" ", "_")}.csv')

def main(queries, start_year, end_year, is_testing_mode):
    for year in range(start_year, end_year + 1):
        create_csv_year(queries, year, is_testing_mode)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process start and end years.")
    parser.add_argument("--start-year", type=int, help="Start year", required=True)
    parser.add_argument("--end-year", type=int, help="End year", required=True)
    parser.add_argument("--testing", action="store_true", help="Enable testing mode")
    args = parser.parse_args()
    
    start_year = args.start_year
    end_year = args.end_year
    is_testing_mode = args.testing

    search_terms = []
    # Get full path to searchterms.txt
    current_directory = os.path.dirname(os.path.realpath(__file__))
    path_to_keywords = os.path.join(current_directory, 'searchterms.txt')

    with open(path_to_keywords) as file:
        keywords = file.readlines()
        [search_terms.append(keyword.strip()) for keyword in keywords]

    start_time = time.time()

    main(search_terms, start_year, end_year, is_testing_mode)

    end_time = time.time()
    execution_time = end_time - start_time

    print(f"Execution time: {execution_time:.2f} seconds")