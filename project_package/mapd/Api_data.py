from Bio import Entrez

from time import sleep

from mapd import DATA_DIR, engine, DB_PATH, PUBMED_DIR

from urllib.error import HTTPError
from tqdm import tqdm


Entrez.email = "parinishtha.bhalla@gmail.com"
class Utilapi:
    def __init__(self, search_query: str):
        self.sleep_time = 0.1  # reduce sleep time
        self.batch_size = 100  # increase batch size
        self.search_query = search_query

    def search(self):
        try:
            search_info = Entrez.esearch(db="pubmed", term=self.search_query, usehistory='y', retmax=100)
            record = Entrez.read(search_info)
            identifiers = record['IdList']

            total_abstract_count = int(record["Count"])
            fetch_webenv = record['WebEnv']
            fetch_querykey = record['QueryKey']
            self.get_abstracts(fetch_webenv, fetch_querykey, total_abstract_count)
        except HTTPError as e:
            print(f"Error while fetching data from pubmed: {e}")

    def get_abstracts(self, fetch_webenv, fetch_querykey, total_abstract_count):
        start = 0
        batch_size = 100
        with tqdm(total=total_abstract_count, desc="Downloading abstracts") as pbar:
            for start in range(0, total_abstract_count, batch_size):
                end = min(total_abstract_count, start + batch_size)
                try:
                    fetch_handle = Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=start,
                                                 retmax=batch_size, webenv=fetch_webenv, query_key=fetch_querykey, )
                    data = fetch_handle.read()
                    path_outfile = PUBMED_DIR.joinpath(f'{start + 1}-{end}.xml')
                    with open(path_outfile, "w", encoding='utf-8') as f:
                        f.write(data)
                except HTTPError as e:
                    print(f"Error while fetching data from pubmed: {e}")
                sleep(self.sleep_time)
                pbar.update(batch_size)

util = Utilapi('biological AND (clinicaltrial[Filter]) AND (2020:2022[PDAT]) AND (english[Filter])')
util.search()

