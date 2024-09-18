"""
Author: Timothy Tran
Date: 4/21/2024

This script reads CSV containing PubMed data (created using pubmed_scraper.py) and outputs stats
regarding the metadata (author info, publication date, keywords, etc.) that each article contained 

Example usage:
    python3 pubmed_get_stats_on_csv.py

Note: Currently, must hardcode the csv file name. Can definitely clean this code up and modularize it

Sample output:
    Metadata coverage stats for data/pubmed_2010.csv
    DOI: 89.97% articles have this field
    Abstract: 91.66% articles have this field
    AuthorsInfo: 99.33% articles have this field
    KeywordList: 3.07% articles have this field
    15.54% articles have a full date (##/##/####)
    68.94% articles only have month and year (##/NA/####)
    15.13% articles only have year (NA/NA/####)
    0.39% articles have no publication date (NA/NA/None)
"""

import pandas as pd
import re
import matplotlib.pyplot as plt

start_year = 1995
end_year = 2024

total_articles = 0

# Define the empty keywords for each column
empty_keywords = {
    # 'PMID': 'EMPTY_P', 
    'DOI': 'No DOI', 
    # 'Title': 'EMPTY_T', 
    'Abstract': 'No Abstract', 
    # 'PubDate': 'EMPTY_Pub', 
    'AuthorsInfo': 'No author info', 
    # 'Journal': 'EMPTY_J', 
    # 'Language': 'EMPTY_L', 
    'KeywordList': 'No keywords', 
    #'SearchTerms': 'EMPTY_ST'
}

# Get stats on publication date (different logic than other columns)
# Define the function to check the validity of date format
def is_valid_date_format(pub_date):
    pattern = r'^\d{2}/\d{2}/\d{4}$'
    return bool(re.match(pattern, str(pub_date)))

def does_not_contain_day(pub_date): 
    pattern = r'^\d{2}/NA/\d{4}$'
    return bool(re.match(pattern, str(pub_date)))

def does_not_contain_month_day(pub_date):
    pattern = r'^NA/NA/\d{4}$'
    return bool(re.match(pattern, str(pub_date)))

def get_publication_date_granualarity(df):
    full_date_info_count = 0
    month_year_count = 0
    year_count = 0
    no_date_count = 0

    global total_articles

    # Iterate through all rows
    for index, row in df.iterrows():
        total_articles += 1
        pub_date = row['PubDate']
        if is_valid_date_format(pub_date):
            full_date_info_count += 1
        elif does_not_contain_day(pub_date):
            month_year_count += 1
        elif does_not_contain_month_day(pub_date):
            year_count += 1
        else:
            no_date_count += 1

    return full_date_info_count, month_year_count, year_count, no_date_count

def create_plot(coverage_data, date_granularity_data):
    # Create the metadata coverage plot
    plt.figure(figsize=(10, 6))
    for col in empty_keywords.keys():
        plt.plot(coverage_data['Year'], coverage_data[col], marker='o', label=col)
    
    plt.xlabel('Year')
    plt.ylabel('Percentage')
    plt.title('Metadata Coverage Over Years')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Create the publication date granularity plot
    plt.figure(figsize=(10, 6))
    for key in ['Full Date', 'Month and Year', 'Year Only', 'No Date']:
        plt.plot(date_granularity_data['Year'], date_granularity_data[key], marker='o', label=key)
    
    plt.xlabel('Year')
    plt.ylabel('Percentage')
    plt.title('Publication Date Granularity Over Years')
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    # Dictionary to store the coverage data
    coverage_data = {key: [] for key in empty_keywords.keys()}
    coverage_data['Year'] = []

    date_granularity_data = {
        'Full Date': [],
        'Month and Year': [],
        'Year Only': [],
        'No Date': [],
        'Year': []
    }

    for year in range(start_year, end_year + 1):
        csv_file = f"data/pubmed_{year}.csv"
        df = pd.read_csv(csv_file)

        # Initialize a dictionary to store the counts
        column_counts = {}

        # Iterate over each column
        for col in empty_keywords.keys():
            # Count non-empty entries in the column
            non_empty_count = df[col].ne(empty_keywords[col]).sum()
            # Calculate the percentage
            percentage = (non_empty_count / len(df)) * 100
            # Store the result in the dictionary
            column_counts[col] = percentage
            # Append the result to the coverage_data
            coverage_data[col].append(percentage)

        coverage_data['Year'].append(year)

        full_date_info_count, month_year_count, year_count, no_date_count = get_publication_date_granualarity(df)

        # Append the publication date granularity data
        date_granularity_data['Full Date'].append((full_date_info_count / len(df)) * 100)
        date_granularity_data['Month and Year'].append((month_year_count / len(df)) * 100)
        date_granularity_data['Year Only'].append((year_count / len(df)) * 100)
        date_granularity_data['No Date'].append((no_date_count / len(df)) * 100)
        date_granularity_data['Year'].append(year)

    # Print the results
    print(f"Metadata coverage stats for {csv_file}")
    for col, percentage in column_counts.items():
        print(f"{col}: {percentage:.2f}% articles have this field")

    # # Print date results
    # print(f"{((full_date_info_count / len(df)) * 100):.2f}% articles have a full date (##/##/####)")
    # print(f"{((month_year_count / len(df)) * 100):.2f}% articles only have month and year (##/NA/####)")
    # print(f"{((year_count / len(df)) * 100):.2f}% articles only have year (NA/NA/####)")
    # print(f"{((no_date_count / len(df)) * 100):.2f}% articles have no publication date (NA/NA/None)")

    create_plot(coverage_data, date_granularity_data)

    print("total articles " + str(total_articles))
if __name__ == "__main__":
    main()