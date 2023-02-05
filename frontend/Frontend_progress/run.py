import csv
import os

from flask import Flask, render_template, request, make_response, send_file
import pandas as pd
import re
from flask import render_template_string, session
from flask_paginate import Pagination, get_page_parameter
from sqlalchemy import select, inspect,create_engine
# from dummy import get_info
#from mapd.Database import Database
from mapd.utils import query_database

app = Flask(__name__, static_folder='static')
app.secret_key = 'ma'

@app.route('/')
def home():
    """ renders the template search_form.html """
    #return render_template('search_form.html')
    return render_template('index.html')


@app.route('/test')
def ss_form():
    """ renders the template search_form.html """
    #return render_template('search_form.html')
    return render_template('combi.html')


# search_form.html contains:
# a form with a text input field for the keyword
# and two date inputs for the start and end date
# and a submit button that triggers the search_results()


# @app.route('/search_results', methods=['GET', 'POST'])
# def search_results():
#     if request.method == 'POST':
#         if 'keyword' in request.form:
#             keyword = request.form['keyword']
#             start_date = request.form['start_date']
#             end_date = request.form['end_date']
#             results = query_database(keyword, start_date, end_date)  # db_path might be a parameter in get_info() - modify acc.
#             # Store the keyword and start_date/end_date in the session
#             session['keyword'] = keyword
#             session['start_date'] = start_date
#             session['end_date'] = end_date
#         else:
#             # handle the case where the keyword field is missing
#             return 'The keyword field is missing from the form. Please try again.'
#         page = request.args.get(get_page_parameter(), type=int, default=1)
#         per_page = 5
#         pagination = Pagination(page=page, total=len(results), per_page=per_page)
#         session['searched_results'] = results
#         return render_template('combi.html',
#                                results=results[(page - 1) * per_page:page * per_page],
#                                pagination=pagination)
#     else:
#         # Retrieve the keyword and start_date/end_date from the session
#         keyword = session.get('keyword', None)
#         start_date = session.get('start_date', None)
#         end_date = session.get('end_date', None)
#         if keyword:
#             results = query_database(keyword, start_date, end_date)
#             page = request.args.get(get_page_parameter(), type=int, default=1)
#             per_page = 5
#             pagination = Pagination(page=page, total=len(results), per_page=per_page)
#             session['searched_results'] = results
#             return render_template('combi.html',
#                                    results=results[(page - 1) * per_page:page * per_page],
#                                    pagination=pagination)
#         else:
#             return 'No keyword found in the session. Please try again.'

@app.route('/results', methods=['POST'])
def search_results():
    """ retrieves the keyword, start and end dates from the form. Calls the get_info()
    with these parameters to retrieve relevant abstracts from the DB.
    Renders the template search_results.html with the results passed as a parameter"""

    keyword = request.form['keyword']
    start_date = request.form['start_date']
    end_date = request.form['end_date']
    results = query_database(keyword, start_date, end_date)
    # if len(results) == 0:
    #     return render_template('combi.html', message="No results found, try again")
    # else:
    return render_template('combi.html', results=results)

# search_results.html
# displays the search results in a table format
# has the hyperlink to pubmed pg by mapping to the PMID
# includes a form with a submit button to download the results
# and a script that uses javascript to toggle the visibility of the abstracts


# search_form.html has a submit button to download file
# @app.route('/download', methods=['POST'])
# def download_results():
#     """ converts the results to a dataframe to generate a CSV file, which
#     the user can download."""
#     keyword = request.form['keyword']
#     results = session.get('searched_results', [])
#     for result in results:
#         result['pubmed_id'] = f'=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/{result["pubmed_id"]}","{result["pubmed_id"]}")'
#         result['abstract_text'] = result['abstract_text'].replace('<mark>', '').replace('</mark>', '')
#         # result['abstract_text'] = result['abstract'].replace('<mark>'+keyword+'</mark>', keyword)
#     df = pd.DataFrame(results)
#     df.to_csv('results.csv', index=False)
#     return send_file('results.csv', as_attachment=True)

@app.route('/download', methods=['POST'])
def download_results():
    """ converts the results to a dataframe to generate a CSV file, which
    the user can download."""
    keyword = request.form['keyword']
    start_date = request.form['start_date']
    end_date = request.form['end_date']
    results = query_database(keyword, start_date, end_date)  # modify when task 3 is complete
    for result in results:
        result['pubmed_id'] = f'=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/{result["pubmed_id"]}","{result["pubmed_id"]}")'
        # result['abstract'] = result['abstract'].replace('<mark>'+keyword+'</mark>', keyword)
        #result['abstract_text'] = re.sub("<mark>.*?</mark>", keyword.lower(), result['abstract_text'])
        result['abstract_text'] = result['abstract_text'].replace('<mark>', '').replace('</mark>', '')
    df = pd.DataFrame(results)
    df.to_csv('results.csv', index=False)
    return send_file('results.csv', as_attachment=True)


# if __name__ == '__main__':
#     app.run(debug=True)

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






