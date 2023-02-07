"""Utils module unit tests."""

from sqlalchemy import select, inspect, create_engine,distinct
from mapd.utils import get_abstract_info, query_database

DB_PATH = 'project_package/tests/data/Test_DB.db'
CONN_STRING = f"sqlite:///{DB_PATH}"
engine = create_engine(CONN_STRING)


def test_query_database():
    """Checks the method to query the database in utils.py"""

    keyword = "amino acid"
    # start_date = '2022-01-01'
    # end_date = '2022-06-31'
    result = query_database(keyword=keyword, start_date=None, end_date=None)
    print(result)

    assert isinstance(result, list)
    assert all(isinstance(i, dict) for i in result)
    assert all('pubmed_id' in i for i in result)
    assert all('Title' in i for i in result)
    assert all('date' in i for i in result)
    assert all('abstract_text' in i for i in result)
    assert all('entities' in i for i in result)
    assert all('labels' in i for i in result)
    assert all(keyword.lower() in i['abstract_text'].lower() for i in result)

    # assert all('<mark>' + keyword + '</mark>' in i['abstract_text'] for i in result)

    start_date = '2022-01-01'
    end_date = '2022-06-31'
    result = query_database(keyword, start_date, end_date)
    assert all(start_date <= i['date'] <= end_date for i in result)
