from sqlalchemy import Column, Integer, String, ForeignKey, DATE, TEXT
from sqlalchemy.orm import declarative_base

Base = declarative_base()

# Define tables using class definitions
# class Abstract(Base):
#     __tablename__ = 'Abstract'
#
#     id = Column(Integer, primary_key=True, autoincrement=True)
#     Title= Column(String)
#     pubmed_id = Column(String)
#     date = Column(DATE)
#     abstract_text = Column(String)
#     # keywords=Column(String) #keywords can be fetched from MESH headings--??
#
# class Entity(Base):
#     __tablename__ = 'entity'
#
#     id = Column(Integer, primary_key=True, autoincrement=True)
#     entity = Column(String)
#     labels=Column(String)
#     abstract_id = Column(Integer, ForeignKey(Abstract.id))

###################################

# The Text data type in SQLAlchemy is mapped to the TEXT data type in MySQL, which has a maximum length of 65,535 characters.
      
class Abstract(Base):
    __tablename__ = 'abstract'

    id = Column(Integer, primary_key=True, autoincrement=True)
    Title = Column(String(255))
    pubmed_id = Column(String(255))
    date = Column(DATE)
    abstract_text = Column(TEXT)
    # abstract_text = Column(LONGTEXT)
    # keywords=Column(String) #keywords can be fetched from MESH headings--??

class Entity(Base):
    __tablename__ = 'entity'

    id = Column(Integer, primary_key=True, autoincrement=True)
    entity = Column(String(255))
    labels = Column(String(255))
    abstract_id = Column(Integer, ForeignKey(Abstract.id))
    

