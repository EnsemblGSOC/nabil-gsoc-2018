import os
from flask import Flask, render_template, request, jsonify
from backend import load_from_db, load_orthologs, load_exon_sequence, make_pairs, get_splice_site_info, get_protein_similarity
import importlib.util



app = Flask(__name__)

@app.route('/')
def home():
    return render_template("index.html")

@app.route('/orthologs')
def orthologs():
    return render_template("orthologs.html",transcript=request.args.get("transcript"),protein_coding=request.args.get("protein_coding"),gencode_basic=request.args.get("gencode_basic"))

@app.route('/extract')
def extract_from_db():
    
    return jsonify(load_from_db(data_type=request.args.get("data_type"), species=request.args.get("species"), gene_id=request.args.get("gene_id"), protein_coding=request.args.get("protein_coding"), gencode_basic=request.args.get("gencode_basic")))

@app.route('/pair_exons')
def pair_exons():
    
    return jsonify(make_pairs(request.args.get("transcript1_id"),request.args.get("transcript2_id")))

    '''
    This is the format of the output

    return jsonify({'transcript1_total_exons' : 8,
                   'transcript2_total_exons' : 10,
                  'best_score': 0.81625, 
                  'best_pairing': [{'transcript1_exon_index': 0, 'transcript2_exon_index': 0}, 
                                    {'transcript1_exon_index': 1, 'transcript2_exon_index': 1}, 
                                    {'transcript1_exon_index': 2, 'transcript2_exon_index': 2}, 
                                    {'transcript1_exon_index': 3, 'transcript2_exon_index': 3}, 
                                    {'transcript1_exon_index': 4, 'transcript2_exon_index': 4}, 
                                    {'transcript1_exon_index': 5, 'transcript2_exon_index': 5}, 
                                    {'transcript1_exon_index': 6, 'transcript2_exon_index': 7}, 
                                    {'transcript1_exon_index': 7, 'transcript2_exon_index': 8}]})
    '''

@app.route('/get_exons')
def get_exons():

    return jsonify(load_exon_sequence(request.args.get("transcript_id"), request.args.get("exon_number")))

@app.route('/get_orthologs')
def get_orthologs():
    print('********')
    print(request.args.get("transcript_id"), request.args.get("protein_coding"), request.args.get("gencode_basic"))
    print('********')
    return jsonify(load_orthologs(transcript_id=request.args.get("transcript_id"), protein_coding=request.args.get("protein_coding"), gencode_basic=request.args.get("gencode_basic")))
    

@app.route('/visualize')
def process_pair_of_transcript():
    
    return render_template("visualize.html",transcript1=request.args.get("transcript1"),transcript2=request.args.get("transcript2"))


@app.route('/get_splice_site')
def get_splice_site():

    return jsonify(get_splice_site_info(request.args.get('transcript_id'), request.args.get('exon_id')))


@app.route('/get_protein_similarity')
def protein_similarity():

    return jsonify(get_protein_similarity(request.args.get('transcript1_id'), request.args.get('transcript2_id'), request.args.get('exon1_id'), request.args.get('exon2_id')))


if __name__ == '__main__':
    app.debug = True
    app.run(host="127.0.0.1",port="5000")
    