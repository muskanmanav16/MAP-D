"""Collection of utility methods."""

import pandas as pd
import sqlite3
import pandas as pd
from typing import List
from mapd.models import Base, Abstract,Entity
from sqlalchemy import select, inspect, create_engine
from sqlalchemy.orm import Session
import re
DB_PATH = 'E:\Desktop\LSISem3\Plab2\Projects\group project plab2\data\gp2_plab2.db'
CONN_STRING = f"sqlite:///{DB_PATH}"
engine = create_engine(CONN_STRING)
test_session = Session(bind=engine)
try:
    engine = create_engine(CONN_STRING)
    test_session = Session(bind=engine)
except Exception as e:
    print(f"Error connecting to database: {e}")

def query_database(self, keyword: str, start_date=None, end_date=None):
    """Returns all abstracts in the database that have been tagged with the queried keyword, and were published
    between the start and end dates (if date parameters are input).

    Returns:
    A list of dict, each dict represents a row in the table, where keys are the column names"""

    query_ = test_session.query(Abstract).outerjoin(Entity, Abstract.id == Entity.abstract_id)
    query_ = query_.filter(Entity.labels.like('%' + keyword + '%'))
    if start_date and end_date:
        query_ = query_.filter(Abstract.date >= start_date, Abstract.date <= end_date)

    df = pd.read_sql(query_.statement, query_.session.bind)

    # highlight keyword in the abstract text
    df["abstract_text"] = df["abstract_text"].apply(
        lambda x: re.sub(f'({keyword})', r'<mark>\1</mark>', x, flags=re.IGNORECASE))
    return df.to_dict(orient='records')
if __name__ == '__main__':
    query_database("19","2021-01-01", "2021-12-31")
