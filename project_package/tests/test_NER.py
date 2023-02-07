"""NER module unit tests."""

# from pathlib import Path
# import spacy
from mapd.NER import EntityPrediction
from sqlalchemy import select, inspect, create_engine
from sqlalchemy.orm import Session
import unittest


DB_PATH = 'project_package/tests/data/Test_DB.db'
CONN_STRING = f"sqlite:///{DB_PATH}"
engine = create_engine(CONN_STRING)
test_session = Session(bind=engine)


class TestModel(unittest.TestCase):
    """Unit tests to check model prediction"""

    def test_entity_prediction(self):
        """Checks if the model prediction is as expected"""

        text = "Merkel cell polyomavirus (MCPyV) is the only human polyomavirus currently known to cause human cancer."
        # expected_output = [("Merkel cell polyomavirus", "CELL"), ("human polyomavirus", "ORGANISM"), ("human cancer", "ORGANISM")]
        expected_output = [{'entity': 'Merkel cell polyomavirus', 'labels': 'CELL'},
                           {'entity': 'human polyomavirus', 'labels': 'ORGANISM'},
                           {'entity': 'human cancer', 'labels': 'ORGANISM'}]
        message = "Predicted entities and expected output does not match"

        test_model = EntityPrediction(test_session)
        entities = test_model.predict_entities(text)

        self.assertEqual(entities, expected_output, message)

    def test_entity_prediction_empty_text(self):
        """Check model prediction for an empty string"""

        text = ""
        expected_output = []
        test_model = EntityPrediction(test_session)
        entities = test_model.predict_entities(text)
        message = "Predicted entities and expected output does not match"

        self.assertEqual(entities, expected_output, message)
