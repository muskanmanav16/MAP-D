import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' #to ignore GPU error

import spacy
import logging
from sqlalchemy import create_engine, Column, Integer, String, Table, MetaData,ForeignKey, DATE,text,select, inspect,func
from sqlalchemy_utils import database_exists
from sqlalchemy.orm import declarative_base, sessionmaker
from typing import Optional, Union
import en_ner_bionlp13cg_md    #The model we are going to use
from mapd import DATA_DIR, engine, DB_PATH, PUBMED_DIR
from mapd.models import Base, Abstract, Entity

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

model_name = "en_ner_bionlp13cg_md"
nlp = spacy.load(model_name)


class EntityPrediction:
    def __init__(self, session):
        self.session = session

    def predict_entities(self, abstract_text):
        doc = nlp(abstract_text)
        entities = []
        for ent in doc.ents:
            entities.append({"entity": ent.text, "labels": ent.label_})
        return entities

    def insert_entities(self, abstract_id, entities):
        for entity in entities:
            e = Entity(entity=entity["entity"], labels=entity["labels"], abstract_id=abstract_id)
            self.session.add(e)
        self.session.commit()
