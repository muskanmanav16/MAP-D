from sqlalchemy import Column, Integer, String, ForeignKey, Date
from sqlalchemy.orm import declarative_base

Base = declarative_base()

# Define tables using class definitions
class Pubmed(Base):
    __tablename__ = 'pubmed'

    id = Column(Integer, primary_key=True, autoincrement=True)
    pubmed_id = Column(String)
    title = Column(String)
    abstract_text = Column(String)
    date = Column(Date)
    keywords = Column(String) #keywords can be fetched from MESH headings--??

# class Entity(Base):
#     __tablename__ = 'entity'
#
#     id = Column(Integer, primary_key=True, autoincrement=True)
#     accession = Column(String)
#     entity = Column(String)
#     labels = Column(String)
#     pubmed_id = Column(Integer, ForeignKey(Abstract.id))


