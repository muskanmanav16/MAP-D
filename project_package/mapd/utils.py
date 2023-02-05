"""Collection of utility methods."""

import pandas as pd
import sqlite3
import pandas as pd
from typing import List
from mapd.models import Base, Abstract,Entity
from sqlalchemy import select, inspect, create_engine
from sqlalchemy.orm import Session
import re
# DB_PATH = 'E:\Desktop\LSISem3\Plab2\Projects\group project plab2\data\gp2_plab2.db'
DB_PATH = 'D:\ClonedRepoSheet1\group2nlp\data\gp2_plab2.db'
CONN_STRING = f"sqlite:///{DB_PATH}"
engine = create_engine(CONN_STRING)
test_session = Session(bind=engine)
try:
    engine = create_engine(CONN_STRING)
    test_session = Session(bind=engine)
except Exception as e:
    print(f"Error connecting to database: {e}")


# def query_database(keyword: str, start_date=None, end_date=None):
#     """Returns all abstracts in the database that have been tagged with the queried keyword, and were published
#     between the start and end dates (if date parameters are input).
#
#     Returns:
#     A list of dict, each dict represents a row in the table, where keys are the column names"""
#
#     query_ = test_session.query(Abstract).outerjoin(Entity, Abstract.id == Entity.abstract_id)
#     query_ = query_.filter(Entity.entity.like('%' + keyword + '%'))
#     if start_date and end_date:
#         query_ = query_.filter(Abstract.date >= start_date, Abstract.date <= end_date)
#
#     df = pd.read_sql(query_.statement, query_.session.bind)
#     # df['date']=
#     # highlight keyword in the abstract text
#     df["abstract_text"] = df["abstract_text"].apply(
#         lambda x: re.sub(f'({keyword})', r'<mark>\1</mark>', x, flags=re.IGNORECASE))
#     return df.to_dict(orient='records')
#     #return df.to_dict(orient='records')

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


if __name__ == '__main__':
    # print(query_database("amino acid","2022-01-01", "2022-12-31"))
    print(query_database("penile"))


