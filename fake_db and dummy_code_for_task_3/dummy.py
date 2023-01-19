import sqlite3
import pandas as pd
from typing import List

# DB_PATH = r"D:\ClonedRepoSheet1\Group2\front_end\MyDB.db"


def get_info(keyword: str, start_date: str, end_date: str) -> List[dict]:
    """ Queries the DB and retrieves row values where the abstract column contains the keyword and the date of
    publication is between the start_date and end_date entered by the user.
    :param keyword: keyword or phrase entered
    :param start_date: from yyyy-mm-dd
    :param end_date: to yyyy-mm-dd
    :returns a list of dict, each dict represents a row in the table, where keys are the column names."""
    conn = sqlite3.connect('MyDB.db')
    query = "SELECT * FROM PUBMED WHERE abstract LIKE '%" + keyword + "%' AND date >= '" + start_date + \
            "' AND date <= '" + end_date + "'"
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df.to_dict(orient='records')


print(get_info("19", "2021-01-01", "2021-12-31"))


