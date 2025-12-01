#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pubmed_query.py — Lab03
Interoghează PubMed cu "TP53 AND cancer" și salvează max 5 articole.
"""

import requests
from pathlib import Path

def query_pubmed(term="TP53 AND cancer", max_results=5):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": term,
        "retmax": max_results,
        "retmode": "json"
    }
    r = requests.get(url, params=params)
    data = r.json()
    return data["esearchresult"]["idlist"]

def fetch_details(pmid_list):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pubmed",
        "id": ",".join(pmid_list),
        "retmode": "xml"
    }
    r = requests.get(url, params=params)
    return r.text

if __name__ == "__main__":
    handle = "mariapaula"
    outdir = Path(".")

    print("[*] Interoghez PubMed…")
    pmids = query_pubmed()

    print(f"[*] PMIDs găsite: {pmids}")

    xml_data = fetch_details(pmids)

    outfile = outdir / f"pubmed_{handle}.txt"
    outfile.write_text(xml_data)

    print(f"[OK] Rezultatul a fost salvat în {outfile}")
