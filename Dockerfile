FROM python:3.9


# All members of group should be in the image metadata using LABEL
# metadata
LABEL maintainer="MAPD <mapd@gmx.net>"
LABEL contributors="Astha Anand <s0asanan@uni-bonn.de>, Deepika Pradeep <s0deprad@uni-bonn.de>, Muskan Manav <s0mumana@uni-bonn.de>, Parinishtha Bhalla <s0pabhal@uni-bonn.de>"
LABEL project-name="Natural Language Processing"
LABEL Description="Building a search engine that utilizes NER"
LABEL Usage="To find relevant biomedical papers from a large corpus"
LABEL project-version="1.0"
LABEL created-date="07-02-2023"


ENV FLASK_PORT=5000
# Expose the port defined by the FLASK_PORT environment variable
EXPOSE ${FLASK_PORT}
EXPOSE $FLASK_PORT


COPY frontend/Frontend_progress/run.py /app/
COPY frontend/Frontend_progress/templates /app/templates/
COPY frontend/Frontend_progress/static /app/static/

# setup???
COPY project_package/setup.py /app/
COPY project_package/mapd /app/mapd

# copy db from cache folder
COPY data/gp2_plab2.db /app/data/gp2_plab2.db
COPY data /app/

# copy task 3 func
COPY project_package/mapd/utils.py /app/mapd/utils.py


WORKDIR /app

# Install required packages
RUN pip install -e .
# scispacy model link


ENTRYPOINT ["python", "run.py"]

