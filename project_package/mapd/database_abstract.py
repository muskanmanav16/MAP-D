
import os
from mapd import DATA_DIR, engine, DB_PATH, PUBMED_DIR
from Bio import Medline
from tqdm import tqdm
from sqlalchemy import create_engine, Column, Integer, String, Table, MetaData,ForeignKey, DATE,text,select, inspect,func
from mapd.models import Base,Abstract,Entity
from datetime import datetime

from sqlalchemy.orm import sessionmaker

Base.metadata.create_all(engine)
Session = sessionmaker(bind=engine)
session = Session()

with engine.begin() as conn:
    for file in tqdm(PUBMED_DIR.glob("*.xml"), desc="Processing Files"):
        with open(file,encoding='utf-8') as handle:
            records = Medline.parse(handle)
            for record in records:
                pmid = record['PMID']
                title = record['TI']
                try:
                    abstract = record['AB']
                except (KeyError, ValueError):
                    abstract = None
                date_str=record['DP']
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
                if date in [None, ""] or abstract in [None, ""] :
                    print("Date  or abstract is empty or None, not adding to database")
                else:
                    abstract_entry = Abstract(pubmed_id=pmid, Title=title, abstract_text=abstract, date=date)
                    session.add(abstract_entry)
            session.commit()