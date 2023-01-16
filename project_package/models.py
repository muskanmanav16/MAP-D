from sqlalchemy import Column, Integer, String, ForeignKey,DATE
from sqlalchemy.orm import declarative_base

Base = declarative_base()

# Define tables using class definitions
class Abstract(Base):
    __tablename__ = 'Abstract'

    id = Column(Integer, primary_key=True, autoincrement=True)
    Title= Column(String)
    pubmed_id = Column(String)
    date = Column(DATE)
    abstract_text = Column(String)

class Entity(Base):
    __tablename__ = 'entity'

    id = Column(Integer, primary_key=True, autoincrement=True)
    accession = Column(String)
    entity = Column(String)
    labels=Column(String)
    pubmed_id = Column(Integer, ForeignKey(Abstract.id))