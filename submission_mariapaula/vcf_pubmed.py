import requests
from pathlib import Path

def parse_vcf(vcf_path):
    variants = []

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")

            # Protecție împotriva liniilor ciudate
            if len(fields) < 5:
                continue

            chrom, pos, vid, ref, alt = fields[:5]

            if vid == ".":
                query = f"chr{chrom}:{pos} AND TP53"
            else:
                query = vid

            variants.append((chrom, pos, vid, query))

    return variants


def search_pubmed(term, retmax=3):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": term,
        "retmax": retmax,
        "retmode": "json",
    }
    r = requests.get(url, params=params)
    data = r.json()
    return data.get("esearchresult", {}).get("idlist", [])


if __name__ == "__main__":
    handle = "mariapaula"
    vcf = Path("tp53_demo.vcf")
    out = Path(f"variants_{handle}.txt")

    vars = parse_vcf(vcf)

    with open(out, "w") as f:
        for chrom, pos, vid, query in vars:
            f.write(f"=== Variantă: {vid} ({chrom}:{pos}) ===\n")
            f.write(f"Căutare PubMed: {query}\n")
            pmids = search_pubmed(query)
            f.write("PMIDs: " + ", ".join(pmids) + "\n\n")

    print(f"[OK] Am salvat rezultatele în {out}")
