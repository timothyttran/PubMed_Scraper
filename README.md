# PubMed_Scraper

## Overview
This Python script retrieves PubMed article metadata and abstracts based on the provided query and year range, then saves them to CSV files.

## Features
- Searches PubMed database for articles matching a given query within a specified year range.
- Supports testing mode, limiting the number of articles retrieved per year to 20.
- Fetches article metadata and abstracts from PubMed.
- Generates a CSV file containing metadata for each year's articles.

### Example Usage
```bash
python3 pubmed_scraper.py --query "mental health" --start-year 2000 --end-year 2020 --testing