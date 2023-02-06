"""Database module unit tests."""

import pytest
from pathlib import Path
import requests
from mapd.Database import Database, Utilapi
from mapd.models import Abstract, Entity
from ncbiutils.ncbiutils import PubMedFetch
from Bio import Medline, Entrez
import unittest
from sqlalchemy import select, inspect, create_engine
from sqlalchemy.orm import Session
from mapd import DATA_DIR, engine, DB_PATH, PUBMED_DIR

query = '("hyaluronan receptors"[MeSH Terms] OR ("hyaluronan"[All Fields] AND "receptors"[All Fields]) OR "hyaluronan receptors"[All Fields] OR "cd44"[All Fields]) AND ((ffrft[Filter]) AND (medline[Filter]) AND (review[Filter]) AND (english[Filter]) AND (2020:2020[pdat]))'

# TEST_DB_PATH = 'project_package/tests/data/Test_DB.db'
TEST_FOLDER = Path(__file__).parent
TEST_DB_PATH = TEST_FOLDER.joinpath("data/Test_DB.db")
TEST_CONN_STRING = f"sqlite:///{TEST_DB_PATH}"
test_engine = create_engine(TEST_CONN_STRING)
test_session = Session(bind=test_engine)

ENTITY_DICT_PATH = '/project_package/tests/data/test_data_entities.txt'


class TestDatabase:
    """Unit tests for Database class in Database.py"""


    # @pytest.fixture(scope='module')
    # @pytest.fixture
    def test_db(self):
        """Create test DB and drop after."""
        db = Database(db_engine=test_engine)
        inspector = inspect(test_engine)

        print('Build Test Database')
        db.build_database()
        assert TEST_DB_PATH.is_file()  # DB is created

        # Check tables and columns
        tables = inspector.get_table_names()
        abstract_cols = ("id", "Title", "pubmed_id", "date", "abstract_text")
        entity_cols = ("id", "entity", "labels", "abstract_id")
        assert all([x in tables for x in ("abstract", "entity")])  # Correct tables
        assert all([x in Abstract.__table__.columns for x in abstract_cols])  # Correct Abstract table columns
        assert all([x in Entity.__table__.columns for x in entity_cols])  # Correct Entity table columns columns


        # yield db

        # print('Delete Test Database')
        # db.session.close()
        # Path.unlink(TEST_DB_PATH)
        # assert not TEST_DB_PATH.is_file()

    def test_row_number(self):
        # x = Database()
        # y = x.session
        abstract_row = test_session.query(Abstract).count()
        entity_row = test_session.query(Entity).count()
        assert abstract_row == 59
        assert entity_row == 930

    def test_get_abstract_info(self):
        db = Database(db_engine=test_engine)
        result = db.get_abstract_info(36298759)

        # Assert the returned result
        assert result == {
            "pubmed_id": '36298759',
            "Title": "Current In Vitro and In Vivo Models to Study MCPyV-Associated MCC.",
            "date": "2022-10-07",
            "abstract_text": "Merkel cell polyomavirus (MCPyV) is the only human polyomavirus currently known to cause human cancer. MCPyV is believed to be an etiological factor in at least 80% of cases of the rare but aggressive skin malignancy Merkel cell carcinoma (MCC). In these MCPyV+ MCC tumors, clonal integration of the viral genome results in the continued expression of two viral proteins: the viral small T antigen (ST) and a truncated form of the viral large T antigen. The oncogenic potential of MCPyV and the functional properties of the viral T antigens that contribute to neoplasia are becoming increasingly well-characterized with the recent development of model systems that recapitulate the biology of MCPyV+ MCC. In this review, we summarize our understanding of MCPyV and its role in MCC, followed by the current state of both in vitro and in vivo model systems used to study MCPyV and its contribution to carcinogenesis. We also highlight the remaining challenges within the field and the major considerations related to the ongoing development of in vitro and in vivo models of MCPyV+ MCC.",
            "entities": {'CANCER': 'skin malignancy Merkel cell carcinoma',
              'CELL': 'Merkel cell polyomavirus',
              'ORGANISM': 'human polyomavirus'},
        }


    def test_get_entity_dict(self):

        with open(ENTITY_DICT_PATH, 'r') as file:
            expected_entity_dict = file.read()

        # expected_entity_dict = ENTITY_DICT_PATH

        # Test the get_entity_dict method
        database = Database(db_engine=test_engine)
        entity_dict = database.get_entity_dict()

        assert entity_dict == expected_entity_dict


class TestApi(unittest.TestCase):
    def test_entrez_search_query(self):
        """Test to check if the NCBI EUtils esearch query works as expected"""

        # Test valid input
        result = Utilapi(search_query=query).search()
        self.assertIsNotNone(result)

        # Test invalid input
        with self.assertRaises(Exception):
            Utilapi(search_query="").search()

        assert result['Count'] == '40', 'Incorrect number of search results'
        # assert result['Retmax'] == '40', 'Incorrect number of search results returned'
        assert result['QueryKey'] is not None, 'Querykey not found'
        assert result['WebEnv'] is not None, 'Webenv not found'

        assert Path(PUBMED_DIR).is_dir()
        assert (any(Path(PUBMED_DIR).iterdir())) is True

    def test_entrez_fetch_query(self):
        """Test to check if the NCBI EUtils efetch query works as expected"""

        #List of PubMed identifiers for those records we wish to retrieve
        uids = ['16186693']

        pubmed_fetch = PubMedFetch(retmax=10)

        record = pubmed_fetch.get_citations(uids)
        self.assertIsNotNone(record)

        # Test invalid input
        with self.assertRaises(Exception):
            pubmed_fetch.get_citations("", id="")



