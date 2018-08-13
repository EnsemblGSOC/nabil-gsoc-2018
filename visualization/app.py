"""

    Routes for the Flask server

"""


import os
from flask import Flask, render_template, send_from_directory, request, jsonify
from backend import *
import importlib.util


app = Flask(__name__)


@app.route('/')
def home():
    return render_template("index.html")

@app.route('/favicon.ico') 
def favicon(): 
    return send_from_directory(os.path.join(app.root_path, 'static'), 'favicon.ico', mimetype='image/vnd.microsoft.icon')

@app.route('/orthologs')
def orthologs():
    return render_template("orthologs.html",transcript=request.args.get("transcript"),protein_coding=request.args.get("protein_coding"),gencode_basic=request.args.get("gencode_basic"), match_score= request.args.get("match_score"), mismatch_penalty= request.args.get("mismatch_penalty"), gap_open= request.args.get("gap_open"), gap_extend= request.args.get("gap_extend"), skip_penalty= request.args.get("skip_penalty"), weight_mode=request.args.get("weight_mode"))

@app.route('/extract')
def extract_from_db():
    
    return jsonify(load_from_db(data_type=request.args.get("data_type"), species=request.args.get("species"), gene_id=request.args.get("gene_id"), protein_coding=request.args.get("protein_coding"), gencode_basic=request.args.get("gencode_basic")))

@app.route('/pair_exons')
def pair_exons():
    
    return jsonify(make_pairs(request.args.get("transcript1_id"), request.args.get("transcript2_id"), request.args.get("weight_mode") , float(request.args.get("match_score")) , float(request.args.get("mismatch_penalty")), float(request.args.get("gap_open")), float(request.args.get("gap_extend")), float(request.args.get("skip_penalty")) ))

@app.route('/get_exons')
def get_exons():

    return jsonify(load_exon_sequence(request.args.get("transcript_id"), request.args.get("exon_number")))

@app.route('/get_orthologs')
def get_orthologs():

    return jsonify(load_orthologs(transcript_id=request.args.get("transcript_id"), protein_coding=request.args.get("protein_coding"), gencode_basic=request.args.get("gencode_basic")))
    
@app.route('/visualize')
def process_pair_of_transcript():
    
    return render_template("visualize.html",transcript1=request.args.get("transcript1"), transcript2=request.args.get("transcript2"), match_score=request.args.get("match_score"), mismatch_penalty=request.args.get("mismatch_penalty"), gap_open=request.args.get("gap_open"), gap_extend=request.args.get("gap_extend"), skip_penalty=request.args.get("skip_penalty"), blosum=request.args.get('blosum'), weight_mode=request.args.get("weight_mode") )

@app.route('/scoring')
def scoring():

    return jsonify(get_scoring_metrics(request.args.get('match_score'), request.args.get('mismatch_penalty'), request.args.get('gap_open'), request.args.get('blosum') ))

@app.route('/get_splice_site')
def get_splice_site():

    return jsonify(get_splice_site_info(request.args.get('transcript_id'), request.args.get('exon_id')))

@app.route('/get_splice_site_similarity')
def splice_site_similarity():

    return jsonify(get_transcript_similarity(request.args.get('transcript1_id'), request.args.get('transcript2_id'), request.args.get('exon1_id'), request.args.get('exon2_id'), float(request.args.get('match_score')), float(request.args.get('mismatch_penalty')), float(request.args.get('gap_open')), float(request.args.get('gap_extend')) , 'gaussian'))

@app.route('/get_transcript_similarity')
def transcript_similarity():

    return jsonify(get_transcript_similarity(request.args.get('transcript1_id'), request.args.get('transcript2_id'), request.args.get('exon1_id'), request.args.get('exon2_id'), float(request.args.get('match_score')), float(request.args.get('mismatch_penalty')), float(request.args.get('gap_open')), float(request.args.get('gap_extend')) , 'uniform'))

@app.route('/get_protein_similarity')
def protein_similarity():

    return jsonify(get_protein_similarity(request.args.get('transcript1_id'), request.args.get('transcript2_id'), request.args.get('exon1_id'), request.args.get('exon2_id')))

@app.route('/valid_gene')
def valid_gene():

    return jsonify(extract_from_gene(request.args.get('gene_name'),request.args.get('gene_id'), request.args.get('species'),request.args.get('protein_coding'), request.args.get('gencode_basic')))

@app.route('/valid_transcript')
def valid_transcript():

    return jsonify(extract_from_transcript(request.args.get('transcript_name'),request.args.get('transcript_id')))

@app.route('/get_species')
def get_species():

    return jsonify(extract_species(request.args.get('transcript_id')))


if __name__ == '__main__':
    
    # running the app

    app.debug = False
    app.run(host="127.0.0.1",port="5000")
    