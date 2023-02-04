"""Database module tests."""

import pytest
from pathlib import Path
import requests
from mapd.Database import Utilapi
from ncbiutils.ncbiutils import PubMedFetch
from Bio import Medline, Entrez
import unittest
from mapd import DATA_DIR, engine, DB_PATH, PUBMED_DIR

query = '("hyaluronan receptors"[MeSH Terms] OR ("hyaluronan"[All Fields] AND "receptors"[All Fields]) OR "hyaluronan receptors"[All Fields] OR "cd44"[All Fields]) AND ((ffrft[Filter]) AND (medline[Filter]) AND (review[Filter]) AND (english[Filter]) AND (2020:2020[pdat]))'


class TestApi(unittest.TestCase):
    def test_entrez_search_query(self):

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
        #List of PubMed identifiers for those records we wish to retrieve
        uids = ['16186693']

        pubmed_fetch = PubMedFetch(retmax=10)

        record = pubmed_fetch.get_citations(uids)
        self.assertIsNotNone(record)

        # Test invalid input
        with self.assertRaises(Exception):
            pubmed_fetch.get_citations("", id="")



