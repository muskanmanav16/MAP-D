"""Collection of utility methods."""

import pandas as pd
import sqlite3
import pandas as pd
from pathlib import Path
from typing import List
from mapd.models import Base, Abstract,Entity
from sqlalchemy import select, inspect, create_engine,distinct
from sqlalchemy.orm import Session
import re
import os

TOP_FOLDER = Path(__file__).parent.parent.parent
DB_PATH= TOP_FOLDER.joinpath("data/gp2_plab2.db")
# CONN_STRING = f"sqlite:///{DB_PATH}"

# CONN_STRING = 'mysql+pymysql://root:mapdrocks@db:3306/mapddb'
# #CONN_STRING = "sqlite:////app/data/gp2_plab2.db"
# engine = create_engine(CONN_STRING)

# CONN_STRING = 'mysql+pymysql://root:mapdrocks@db:3306/mapddb'
# engine = create_engine(CONN_STRING)

if 'DATABASE_URL' in os.environ:
    CONN_STRING = os.environ['DATABASE_URL']
else:
    CONN_STRING = f"sqlite:///{DB_PATH}"


def get_abstract_info(pubmed_id: int):
    engine = create_engine(CONN_STRING)
    test_session = Session(bind=engine)
    stmt = select(
        Abstract.pubmed_id,
        Abstract.Title,
        Abstract.date,
        Abstract.abstract_text,
        Entity.entity,
        Entity.labels
    ).filter_by(pubmed_id=pubmed_id).join(Entity, isouter=True)

    entries = test_session.execute(stmt).fetchall()
    ent = {}
    for pubmed_id, title, date, abstract_text, entity, labels in entries:
        ent[labels] = entity
    for pubmed_id, title, date, abstract_text, entity, labels in entries:
        entries_dict = {
            "pubmed_id": pubmed_id,
            "Title": title,
            "date": date.strftime("%Y-%m-%d"),
            "abstract_text": abstract_text,
            "entities": ', '.join(ent.values()),
            "labels": ', '.join(ent.keys())
        }
    return entries_dict

def query_database(keyword: str, start_date=None, end_date=None):
    engine = create_engine(CONN_STRING)
    test_session = Session(bind=engine)
    query_ = test_session.query(Abstract.pubmed_id).join(Entity, Abstract.id == Entity.abstract_id)
    query_ = query_.filter(Entity.entity.like('%' + keyword + '%'))
    if start_date and end_date:
        query_ = query_.filter(Abstract.date >= start_date, Abstract.date <= end_date)
    pubmed_id_set = set([row[0] for row in query_.all()])
    pubmed_id_set = list(pubmed_id_set)
    test_list=[]
    for pub_id in pubmed_id_set:
        test_list.append(get_abstract_info(pub_id))
    df = pd.DataFrame(test_list)
    if 'abstract_text' in df.columns:
        df["abstract_text"] = df["abstract_text"].apply(
            lambda x: re.sub(f'({keyword})', r'<mark>\1</mark>', x, flags=re.IGNORECASE))
    return df.to_dict(orient='records')

# if __name__ == '__main__':
#     print(query_database(keyword="amino acid",start_date='2022-01-01',end_date='2022-12-31'))