# PubMed_Scraper

## Overview
This Python script `pubmed_scraper` retrieves PubMed article metadata and abstracts based on the provided query and year range, then saves them to CSV files. A list of "search_terms" are used to query the PubMed database.

This repo also contains several other scripts to analyze/handle the CSV's that result from the PubMed scraping.

## Features
- Searches PubMed database for articles matching a given query within a specified year range.
- Supports testing mode, limiting the number of articles retrieved per year to 20.
- Fetches article metadata and abstracts from PubMed.
- Generates a CSV file containing metadata for each year's articles.

### Example Usage
Specify search terms in "searchterms.txt
```bash
python3 pubmed_scraper.py --start-year 2000 --end-year 2020 --testing