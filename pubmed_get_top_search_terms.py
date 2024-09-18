"""
Author: Timothy Tran
Date: 9/18/2024

This script outputs a CSV file "occurences_by_popularity" where each column contains the top search terms 
from the PubMed articles of that year. 

Note that due to the implementation of the PubMed scraper, the max articles from a certain search term is 
60,000. This is not an issue since the scraper searches by relevance, so any term outside of 10,000 is likely irrelevant.

It contains the following global variables
    --start-year: The start year for the CSV file.
    --end-year:   The end year for the CSV file (inclusive).
    --top-n:      The top N number of search terms to show per year
"""

import csv
import pandas as pd
import re
import ast
import numpy as np

from collections import Counter

start_year = 1995
end_year = 2024
top_n = 100 # Top n search terms

def main():
    top_search_terms_per_year = {}

    for year in range(start_year, end_year + 1):
        csv_file = f"data/pubmed_{year}.csv"
        df = pd.read_csv(csv_file)

        top_search_terms = get_top_search_terms_year(df)

        top_search_terms_per_year[year] = top_search_terms

        print(f"{len(top_search_terms)} top search terms for {year} recorded")

    print(top_search_terms_per_year)
    print('Writing info to csv file...')
    write_to_csv(top_search_terms_per_year)

def get_top_search_terms_year(df):
    top_search_terms = {}

    for index, row in df.iterrows():

        search_terms_set = ast.literal_eval(row['SearchTerms'])
        for search_term in search_terms_set:
            top_search_terms[search_term] = top_search_terms.setdefault(search_term, 0) + 1

    # Use Counter to count occurrences
    counter = Counter(top_search_terms)

    # Get the top n most common strings
    top_search_terms = counter.most_common(top_n)

    return top_search_terms

    # print(f"Top {top_n} most common strings:")
    # for string, occurrences in top_n:
    #     print(f"{string}: {occurrences}")

def write_to_csv(top_search_terms_per_year):

    row = ['Year'] + [str(i) for i in range(1, top_n + 1)]
    print(row) #TODO: FIX THIS
    # column_names = [str(year) for year in range(start_year, end_year + 1)]

    # Write the results to a CSV file
    with open('occurrences_by_popularity.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # writer.writerow(column_names)

        data = []
        data.append(row)

        # Iterate over the years
        for year, occurrences in top_search_terms_per_year.items():
            column_data = [year]
            print(year)
            for (searchterm, occurence) in occurrences:
                # print(f"{searchterm}, {occurence}")
                column_data.append(f"{searchterm}, {occurence}")

            data.append(column_data)

        print()
        # print(data)
        # Write the row to the CSV file
        writer.writerows(np.array(data).transpose())

    print("CSV file 'occurrences_by_popularity.csv' has been created successfully.")

if __name__ == "__main__":
    main()