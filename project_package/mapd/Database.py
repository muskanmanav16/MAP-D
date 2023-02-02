import json
import os
import time
import logging
import requests
import xmltodict
from tqdm import tqdm
from pathlib import Path
from typing import Optional, Union
from sqlalchemy import select, inspect,create_engine, Column, Integer, String, Table, MetaData,ForeignKey, DATE,text,select, inspect,func
from sqlalchemy.orm import Session
from sqlalchemy_utils import database_exists
from datetime import datetime
from Bio import Medline, Entrez
from os import listdir

from time import sleep
from urllib.error import HTTPError
from sqlalchemy.orm import sessionmaker
from mapd import DATA_DIR, engine, DB_PATH, PUBMED_DIR

from mapd.models import Base, Abstract,Entity
from mapd.constant import ABSTRACT, ENTITY,QUERY_STRING
# from mapd.Api_data import Utilapi
from mapd.NER import EntityPrediction
Entrez.email = "mapd@gmx.net"
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Database:
    """For interfacing with the relational database."""

    def __init__(self, db_engine=engine):
        self.engine = db_engine
        self.session = Session(bind=self.engine)

        # When DB doesn't exist
        if not database_exists(self.engine.url):
            self.build_database()

        # When DB exists, but tables not made yet
        tables = inspect(self.engine).get_table_names()
        if not tables:
            self.build_database()

    def rebuild_database(self) -> None:
        """Burn everything and builds the database."""
        self.drop_database()
        self.build_database()

    def build_database(self) -> None:
        """Build the tables of the database."""
        logger.info("Building database...")
        Base.metadata.create_all(bind=self.engine)

    def drop_database(self) -> None:
        """Drop all of the associated tables in the database."""
        logger.warning("Dropping database...")
        Base.metadata.drop_all(bind=self.engine)

    def add_abstract_to_databse(self):

        if len(listdir(PUBMED_DIR)) < 50:
            Utilapi(QUERY_STRING).search()
        # Utilapi(QUERY_STRING).search()
        if self.session.query(Abstract).count() < 1500:
            with self.engine.begin() as conn:
                for file in tqdm(PUBMED_DIR.glob("*.xml"), desc="Processing Files"):
                    with open(file, encoding='utf-8') as handle:
                        records = Medline.parse(handle)
                        for record in records:
                            pmid = record['PMID']
                            title = record['TI']
                            try:
                                abstract = record['AB']
                            except (KeyError, ValueError):
                                abstract = None
                            date_str = record['DP']
                            try:
                                date = datetime.strptime(date_str, '%Y-%m-%d')
                            except ValueError:
                                try:
                                    date = datetime.strptime(date_str, '%Y %b %d')
                                except ValueError:
                                    try:
                                        date = datetime.strptime(date_str, '%Y %B %d')
                                    except ValueError:
                                        date = None
                            if date in [None, ""] or abstract in [None, ""]:
                                continue
                            else:
                                abstract_entry = Abstract(pubmed_id=pmid, Title=title, abstract_text=abstract, date=date)
                                self.session.add(abstract_entry)
                        self.session.commit()
    def add_entity_data(self):
        self.add_abstract_to_databse()
        # Get all abstracts in the database from Abstract table
        entity_predictor = EntityPrediction(self.session)
        abstracts = self.session.query(Abstract).all()

        # Loop through each abstract and predict entities
        for abstract in tqdm(abstracts, desc="Predicting entities for abstracts"):
            entities = entity_predictor.predict_entities(abstract.abstract_text)
            entity_predictor.insert_entities(abstract.id, entities)

    def get_abstract_info(self,pubmed_id: int) -> Optional[dict]:
        """Get abstract info  for a given pubmedid from the relational database."""
        # self.add_entity_data()
        stmt = select(
            Abstract.pubmed_id,
            Abstract.Title,
            Abstract.date,
            Abstract.abstract_text,
            Entity.entity,
            Entity.labels
        ).filter_by(pubmed_id=pubmed_id).join(Entity, isouter=True)

        entries = self.session.execute(stmt).fetchall()
        ent = {}
        for pubmed_id, title, date, abstract_text, entity, labels in entries:
            ent[labels] = entity
        for pubmed_id, title, date, abstract_text, entity, labels in entries:
            entries_dict = {
                "pubmed_id": pubmed_id,
                "Title": title,
                "date": date.strftime("%Y-%m-%d"),
                "abstract_text": abstract_text,
                "entities": ent
            }
        return entries_dict

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
        except HTTPError:
            pass

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
                except HTTPError:
                    pass
                sleep(self.sleep_time)
                pbar.update(batch_size)
if __name__ == '__main__':
    Database().add_entity_data()
    # Database().get_abstract_info(36680181)