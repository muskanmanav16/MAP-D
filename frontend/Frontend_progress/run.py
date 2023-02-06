import csv
import os
from flask import Flask, render_template, request, make_response, send_file
import pandas as pd
import re
from flask import render_template_string, session
from flask_paginate import Pagination, get_page_parameter
from sqlalchemy import select, inspect,create_engine
from mapd.utils import query_database


app = Flask(__name__, static_folder='static')


@app.route('/')
def home():
    """ renders the template index.html """
    return render_template('index.html')


@app.route('/test')
def ss_form():
    """ renders the template combi.html """
    return render_template('combi.html')


@app.route('/results', methods=['POST'])
def search_results():
    """ retrieves the keyword, start and end dates from the form. Calls the query_database()
    with these parameters to retrieve relevant abstracts from the DB.
    Renders the template combi.html with the results passed as a parameter"""

    keyword = request.form['keyword']
    start_date = request.form['start_date']
    end_date = request.form['end_date']
    results = query_database(keyword, start_date, end_date)
    return render_template('combi.html', results=results)


@app.route('/download', methods=['POST'])
def download_results():
    """ converts the results to a dataframe to generate a CSV file, which
    the user can download."""
    keyword = request.form['keyword']
    start_date = request.form['start_date']
    end_date = request.form['end_date']
    results = query_database(keyword, start_date, end_date)  
    # add pubmed page hyperlink in the downloaded csv file and remove html <mark> tags
    for result in results:
        result['pubmed_id'] = f'=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/{result["pubmed_id"]}","{result["pubmed_id"]}")'
        result['abstract_text'] = result['abstract_text'].replace('<mark>', '').replace('</mark>', '')
    df = pd.DataFrame(results)
    df.to_csv('results.csv', index=False)
    return send_file('results.csv', as_attachment=True)


if __name__ == '__main__':
    # Check for the FLASK_PORT environment variable
    flask_port = os.getenv('FLASK_PORT')
    print(flask_port)
    # If the FLASK_PORT environment variable is set, use it as the port for the app
    if flask_port:
        app.run(debug=True, port=flask_port, host="0.0.0.0")
    else:
        # If the FLASK_PORT environment variable is not set, use the default port of 5000
        app.run(debug=True, port=5000, host="0.0.0.0")


