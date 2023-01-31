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
from models import Base, Pubmed
from api_class import Utilapi

class Database:

    def __init__(self, db_engine=engine):
        self.engine = db_engine
        self.session = Session(bind=self.engine)

        if not database_exists(self.engine.url):
            self.build_database()

        tables = inspect(self.engine).get_table_names()
        if not tables:
            self.build_database()

    def rebuild_database(self) -> None:
        self.drop_database()
        self.build_database()

    def build_database(self) -> None:
        Base.metadata.create_all(bind=self.engine)

    def drop_database(self) -> None:
        Base.metadata.drop_all(bind=self.engine)

    def add_data(self, search_query: str) -> None:
        """To add data to the database."""
        profiler = Utilapi(search_query=search_query)
        profiler.search()
        identifiers = profiler.get_abstracts()
        pubmed_entry = Pubmed(pubmed_id=identifiers["pmid"], title=identifiers["title"], abstract_text=identifiers["abstract"], date=identifiers["pub_date"])
        self.session.add(pubmed_entry)
        self.session.commit()


x = Database()
x.rebuild_database()
x.add_data('biological AND (clinicaltrial[Filter]) AND (2020:2023[PDAT]) AND (english[Filter])')
