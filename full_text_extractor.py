# full_text_extractor.py

import os
import pandas as pd
from scidownl import scihub_download 
import fitz
import signal

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException("Timed out")

def fetch_full_text(pmid, timeout=10):
    try:
        pdf_path = f"{pmid}.pdf"
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(timeout)  # Set the alarm
        scihub_download(pmid, out=pdf_path)
        signal.alarm(0)
        return pdf_path
    except TimeoutException:
        print(f"Timeout occurred while downloading full text for PMID {pmid}")
        return None
    except Exception as e:
        print(f"Error downloading full text for PMID {pmid}: {e}")
        return None
def extract_text_from_pdf(pdf_path):
    try:
        doc = fitz.open(pdf_path)
        text = ""
        for page in doc:
            text += page.get_text()
        return text
    except Exception as e:
        print(f"Error extracting text from PDF: {e}")
        return ""

def delete_downloaded_pdfs(pdf_path):
    if os.path.exists(pdf_path):
        os.remove(pdf_path)
        print(f"Deleted PDF: {pdf_path}")
    else:
        print(f"PDF does not exist: {pdf_path}")

def create_csv_with_full_text(pmid_list, abstract_list, year, get_full_text=True):
    if not get_full_text:
        return
    # Check if the CSV file already exists for the specific year
    csv_file_path = f'full_text_{year}.csv'
    if os.path.exists(csv_file_path):
        # Read the existing CSV file to check for existing PMIDs
        existing_df = pd.read_csv(csv_file_path)
        existing_pmids = set(existing_df['PMID'])
    else:
        existing_pmids = set()

    full_text_data = []
    for pmid, abstract in zip(pmid_list, abstract_list):
        # Check if the PMID already exists in the CSV file
        if pmid in existing_pmids:
            print(f"PMID {pmid} already exists in the CSV file for {year}. Skipping...")
            continue

        pdf_path = fetch_full_text(pmid)
        full_text_available = pdf_path is not None 
        has_abstract = False if abstract == 'No Abstract' else True
        if pdf_path:
            full_text = extract_text_from_pdf(pdf_path)
            delete_downloaded_pdfs(pdf_path)
            full_text_data.append((pmid, full_text_available, has_abstract, abstract, full_text))
        else:
            full_text_data.append((pmid, full_text_available, has_abstract, abstract, None))

    # Append new data to the existing CSV file or create a new one if it doesn't exist
    df = pd.DataFrame(full_text_data, columns=['PMID','Full_Text_Available', 'Has_Abstract','Abstract', 'Full_Text'])
    df.to_csv(csv_file_path, mode='a', header=not os.path.exists(csv_file_path), index=False)
