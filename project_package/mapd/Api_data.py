from Bio import Entrez

from time import sleep
from datetime import datetime
from mapd import DATA_DIR, engine, DB_PATH, PUBMED_DIR
from Bio import Medline
from datetime import datetime
from mapd.models import Base,Abstract,Entity
from pathlib import Path
from typing import Optional, Union
from sqlalchemy import select, inspect
from sqlalchemy.orm import Session
from sqlalchemy_utils import database_exists

Entrez.email = "parinishtha.bhalla@gmail.com"
from tqdm import tqdm

class Utilapi:

    def __init__(self, search_query: str):
        self.sleep_time = 0.2
        self.path_outfile = None
        self.search_query = search_query

    def search(self):
        search_info = Entrez.esearch(db="pubmed", term=self.search_query, usehistory='y', retmax=10000)
        record = Entrez.read(search_info)
        identifiers = record['IdList']

        total_abstract_count = int(record["Count"])
        fetch_webenv = record['WebEnv']
        fetch_querykey = record['QueryKey']
        self.get_abstracts(fetch_webenv, fetch_querykey, total_abstract_count)

    def get_abstracts(self, fetch_webenv, fetch_querykey, total_abstract_count):
        start = 0
        batch_size = 1
        with tqdm(total=total_abstract_count, desc="Downloading abstracts") as pbar:
            for start in range(0, total_abstract_count, batch_size):
                end = min(total_abstract_count, start + batch_size)
                fetch_handle = Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=start,
                                             retmax=batch_size, webenv=fetch_webenv, query_key=fetch_querykey, )
                data = fetch_handle.read()
                path_outfile = PUBMED_DIR.joinpath(f'{start + 1}.xml')
                with open(path_outfile, "w", encoding='utf-8') as f:
                    f.write(data)
                sleep(self.sleep_time)
                pbar.update(batch_size)

util = Utilapi('biological AND (clinicaltrial[Filter]) AND (2020:2022[PDAT]) AND (english[Filter])')
util.search()

with engine.begin() as conn:
    for file in tqdm(PUBMED_DIR.glob("*.xml"), desc="Processing Files"):
        with open(file) as handle:
            records = Medline.parse(handle)
            for record in records:
                pmid = record['PMID']
                abstract = record['AB']
                try:
                    date = datetime.strptime(record['DP'], '%Y %b %d').strftime('%Y/%m/%d')
                except (KeyError, ValueError):
                    date = record['DP']

                conn.execute(Abstract.insert().values(pmid=pmid, abstract=abstract, date=date))

