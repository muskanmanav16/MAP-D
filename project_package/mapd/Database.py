
import logging
from tqdm import tqdm
from typing import Optional
from sqlalchemy import select, inspect
from sqlalchemy.orm import Session

from sqlalchemy_utils import database_exists
from datetime import datetime
from Bio import Medline, Entrez

from time import sleep
from urllib.error import HTTPError
from http.client import IncompleteRead
from mapd import engine, PUBMED_DIR

from mapd.models import Base, Abstract, Entity
from mapd.constant import QUERY_STRING
from mapd.NER import EntityPrediction

Entrez.email = "mapd@gmx.net"
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Database:
    """For interfacing with the relational database."""

    def __init__(self, db_engine=engine, cache_dir=PUBMED_DIR, query=QUERY_STRING):

        self.engine = db_engine
        self.session = Session(bind=self.engine)
        self.PUBMED_DIR = cache_dir
        self.QUERY_STRING = query
        self.raw_entity_data = dict()

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

    def add_entity_data(self):
        """Populating the Entity table it also checks if record exists or not to aviod any duplicates"""

        self.add_abstract_to_database()
        # Get all abstracts in the database from Abstract table
        abstracts = self.session.query(Abstract).all()
        existing_entities = self.session.query(Entity).all()
        # Check if the number of entries in the Entity table is less than 900
        if len(existing_entities) >= 800:
            print("Entities have already been predicted.")
            return
        if len(abstracts) > 0:
            # Loop through each abstract and predict entities
            for abstract in tqdm(abstracts, desc="Predicting entities for abstracts"):
                existing_entities = self.session.query(Entity).filter(Entity.abstract_id == abstract.id).all()
                if existing_entities:
                    continue
                else:
                    entity_predictor = EntityPrediction(self.session)
                    entities = entity_predictor.predict_entities(abstract.abstract_text)
                    entity_predictor.insert_entities(abstract.id, entities)

    def add_abstract_to_database(self):
        """First check if files exists in pubmed_dir if not it will download the Abstract in the cache folder
        Also checks the Abstract Table if it filled or not then add the data to the database """

        if len(list(self.PUBMED_DIR.glob("*.xml"))) > 0:
            print('Abstract Downloaded')
        else:
            Utilapi(QUERY_STRING).search()
        for file in tqdm(self.PUBMED_DIR.glob("*.xml"), desc="Processing Files"):
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
                        abstract_entry = self.session.query(Abstract).filter_by(pubmed_id=pmid).first()
                        if abstract_entry is not None:
                            continue
                        else:
                            abstract_entry = Abstract(pubmed_id=pmid, Title=title, abstract_text=abstract,
                                                  date=date)
                            self.session.add(abstract_entry)
                self.session.commit()

    def get_entity_dict(self):
        """Populate the raw_entity_data with the Entity and labels stored in the Entity Table

        returns: dict, with key Entity and Labels as value"""

        self.add_entity_data()
        entities = self.session.query(Entity).all()
        for entity in entities:
            self.raw_entity_data[entity.entity] = entity.labels
        self.raw_entity_data = {k: v for k, v in self.raw_entity_data.items() if v not in self.raw_entity_data.keys()}
        return self.raw_entity_data

    def get_abstract_info(self, pubmed_id: int) -> Optional[dict]:
        """Get abstract info  for a given pubmedid from the relational database.
        param:
        pubmed_id: int

        returns: dict, each dict represents a records in the abstract table along
        with dictionary of entities containing entity and labels"""

        self.add_entity_data()
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
    """For interfacing with the NCBI Entrez API"""

    def __init__(self, search_query: str, cache_dir=PUBMED_DIR):
        self.sleep_time = 0.1  # reduce sleep time
        self.batch_size = 100  # increase batch size
        self.search_query = search_query
        self.PUBMED_DIR = cache_dir

    def search(self):
        """Search pubmed abstract using NCBI Entrez API - esearch"""

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
        return record

    def get_abstracts(self, fetch_webenv, fetch_querykey, total_abstract_count):
        """Batch download of PUBMED ABSTRACTS using NCBI Entrez API - eftech each batch file contains 100 abstracts"""

        start = 0
        batch_size = 100
        with tqdm(total=total_abstract_count, desc="Downloading abstracts") as pbar:
            for start in range(0, total_abstract_count, batch_size):
                end = min(total_abstract_count, start + batch_size)
                try:
                    fetch_handle = Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=start,
                                                 retmax=batch_size, webenv=fetch_webenv, query_key=fetch_querykey, )
                    data = fetch_handle.read()
                    path_outfile = self.PUBMED_DIR.joinpath(f'{start + 1}-{end}.xml')
                    with open(path_outfile, "w", encoding='utf-8') as f:
                        f.write(data)
                except (HTTPError, IncompleteRead):
                    pass
                sleep(self.sleep_time)
                pbar.update(batch_size)

#run below command to download the abstract and populate the db
# if __name__== '__main__':
#     db=Database()
#     db.add_entity_data()