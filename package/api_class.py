from Bio import Entrez

from time import sleep

from startup import DATA_DIR, engine, DB_PATH, PUBMED_DIR
from Bio import Medline
from datetime import datetime


from pathlib import Path
from typing import Optional, Union
from sqlalchemy import select, inspect
from sqlalchemy.orm import Session
from sqlalchemy_utils import database_exists
# from models import Base, Pubmed

Entrez.email = "parinishtha.bhalla@gmail.com"



class Utilapi:

    def __init__(self, search_query: str):
        self.sleep_time = 0.2
        self.path_outfile = None
        self.search_query = search_query
        # self.records = None

    # @staticmethod
    # def cache_file_exists(path: Union[str, Path]) -> bool:
    #     """Checks if cache file exists."""
    #     return path.is_file() if isinstance(path, Path) else Path.is_file(path)

    def search(self) -> None:
        search_info = Entrez.esearch(db="pubmed", term=self.search_query, usehistory='y', retmax=50)
        record = Entrez.read(search_info)
        # print(record)
        identifiers = record['IdList']
        # print(identifiers)

        total_abstract_count = int(record["Count"])
        # print(total_abstract_count)
        self.fetch_webenv = record['WebEnv']
        self.fetch_querykey = record['QueryKey']
        # print(self.fetch_webenv)
        # for uid in record['IdList']:
        # self.get_abstracts
        # return fetch_webenv, fetch_querykey

    def get_abstracts(self):

        start = 0
        batch_size = 1
        for start in range(0, 6, batch_size):
            path_outfile = PUBMED_DIR.joinpath(f'{start + 1}.xml')
            end = min(3, start + batch_size)
            print("Going to download record %i to %i" % (start + 1, end))
            fetch_handle = Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=start,
                                         retmax=batch_size, webenv=self.fetch_webenv, query_key=self.fetch_querykey, )
            data = fetch_handle.read()
            # print(data)
            # if self.cache_file_exists(path_outfile):
            with open(path_outfile, "w") as f:
                f.write(data)

            with open(path_outfile) as handle:
                records = Medline.parse(handle)
                # print(records)
                # print(type(records))
                dic = {}
                for record in records:
                    # print(record)
                    dic['pmid'] = record['PMID']
                    dic['title'] = record['TI']
                    dic['abstract'] = record['AB']
                    try:
                        # x = datetime.strptime(record['DP'], '%Y %b %d').strftime('%Y/%m/%d')
                        x = datetime.strptime(record['DP'], '%Y %b %d').date()
                        print(type(x))
                        dic['pub_date'] = x
                    except (KeyError, ValueError):
                        dic['pub_date'] = record['DP']
                print(dic)
            sleep(self.sleep_time)

        return dic



util = Utilapi('biological AND (clinicaltrial[Filter]) AND (2020:2023[PDAT]) AND (english[Filter])')
util.search()
