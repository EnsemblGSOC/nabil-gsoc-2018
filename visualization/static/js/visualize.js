var transcript1 = [];
var transcript2 = [];
var lines = [];
var line_colors=["#99ff99","#3333ff","#cc00ff","#00ffff","#ff6600","#00ff00","#669999"];
var transcript1_id;
var transcript2_id;
var exon1_id;
var exon2_id;

function init(){

    transcript1_id = document.getElementById("transcript1_id").innerHTML;
    transcript2_id = document.getElementById("transcript2_id").innerHTML;

    console.log(transcript1_id);
    console.log(transcript2_id);
    

    $.get(`http://127.0.0.1:5000/pair_exons?transcript1_id=${transcript1_id}&transcript2_id=${transcript2_id}`, function(data, status){
        
        console.log(data);
        console.log(status);
        

        var draw = document.getElementById('draw');
        var height = draw.clientHeight;
        var width = draw.clientWidth;

        var max_len = Math.max(data["transcript1_total_exons"],data["transcript2_total_exons"]);
        
        
        
        var unit_height = Math.floor(height*0.1);
        var unit_width = Math.floor((width/max_len)*0.1);
        var exon_width = Math.floor((width - unit_width * (max_len+1)) / max_len);

        console.log(unit_width);
        console.log(exon_width);
        
        
        var svg_elements =  "";

        for(var i=0; i<data["transcript1_total_exons"]; i++){
            var j = i.toString();
            transcript1['exon'+j] = [];
            transcript1['exon'+j]["x"] = i*(exon_width + unit_width)+unit_width;
            transcript1['exon'+j]["y"] = unit_height;
            transcript1['exon'+j]["mid"] = transcript1['exon'+j]["x"] + Math.floor(exon_width/2) ;
        }

        for(var i=0; i<data["transcript1_total_exons"]; i++){
            var j = i.toString();
            svg_elements += `<rect x="${transcript1['exon'+j]["x"]}" y="${transcript1['exon'+j]["y"]}" width="${exon_width}" height="${unit_height}" style="fill:blue;stroke:pink;stroke-width:5;fill-opacity:0.1;stroke-opacity:0.9" onclick="get_exon_sequence('${transcript1_id}',${j})" />`;
        }

        for(var i=0; i<data["transcript2_total_exons"]; i++){
            var j = i.toString();
            transcript2['exon'+j] = [];
            transcript2['exon'+j]["x"] = i*(exon_width + unit_width)+unit_width;
            transcript2['exon'+j]["y"] = unit_height*8;
            transcript2['exon'+j]["mid"] = transcript2['exon'+j]["x"] + Math.floor(exon_width/2) ;
        }

        for(var i=0; i<data["transcript2_total_exons"]; i++){
            var j = i.toString();
            svg_elements += `<rect x="${transcript2['exon'+j]["x"]}" y="${transcript2['exon'+j]["y"]}" width="${exon_width}" height="${unit_height}" style="fill:blue;stroke:pink;stroke-width:5;fill-opacity:0.1;stroke-opacity:0.9" onclick="get_exon_sequence('${transcript2_id}',${j})" />`;
        }

        for(var i=0;i<data["best_pairing"].length;i++){
            var j = i.toString();
            lines[j] = [];
            lines[j]["exon1"] = data["best_pairing"][i]["transcript1_exon_index"].toString();
            lines[j]["exon2"] = data["best_pairing"][i]["transcript2_exon_index"].toString();
            lines[j]["point1"] = { "x" : transcript1['exon'+ lines[j]["exon1"]]["mid"],
                                   "y" : 2*unit_height };

            lines[j]["point2"] = { "x" : transcript1['exon'+ lines[j]["exon1"]]["mid"],
                                   "y" : 4*unit_height };

            lines[j]["point3"] = { "x" : transcript2['exon'+ lines[j]["exon2"]]["mid"],
                                   "y" : 6*unit_height };

            lines[j]["point4"] = { "x" : transcript2['exon'+ lines[j]["exon2"]]["mid"],
                                   "y" : 8*unit_height };
            
            lines[j]["color"] = line_colors[i%7];
        }

        for(var i=0;i<data["best_pairing"].length;i++){
            var j = i.toString();

            svg_elements += `<polyline points="${lines[j]['point1']['x']},${lines[j]['point1']['y']} ${lines[j]['point2']['x']},${lines[j]['point2']['y']} ${lines[j]['point3']['x']},${lines[j]['point3']['y']} ${lines[j]['point4']['x']},${lines[j]['point4']['y']}" style="fill:none;stroke:${lines[j]['color']};stroke-width:3" onmouseover="on_focus('line${j}')" onmouseout="out_focus('line${j}')" onclick="pop_up_similarity( '${transcript1_id}', '${transcript2_id}', '${lines[j]["exon1"]}', '${lines[j]["exon2"]}')" id="line${j}" />`;
                                
        }

        console.log(transcript1);
        console.log(transcript2);
        console.log(lines);
        

        draw.innerHTML = svg_elements;
        
    })
    

}


function on_focus(id){

    var ele = document.getElementById(id);
    ele.style["stroke-width"] = "10px";
}

function out_focus(id){

    var ele = document.getElementById(id);
    ele.style["stroke-width"] = "3px";
}

function pop_up_similarity(transcript1, transcript2, exon1, exon2){


    clear_modal_line();

    ele = document.getElementById("tab_splice_site");

    $.get(`http://127.0.0.1:5000/get_exons?transcript_id=${transcript1}&exon_number=${exon1}`, function(data, status){

        exon1_id = data["exon_id"];


        $.get(`http://127.0.0.1:5000/get_exons?transcript_id=${transcript2}&exon_number=${exon2}`, function(data, status){


            exon2_id = data["exon_id"];

            $('#modal_line').modal('open');

            ele.innerHTML = `<div style="margin-top:5%">
                                <div class="row">
                                    <div class="col s6">
                                        <div class="row">
                                            <div class="col s12">
                                                <h5>Transcript ID : ${transcript1}</h5>
                                            </div>
                                            <div class="col s12">
                                                <h5>Exon ID : ${exon1_id}</h5>
                                            </div>
                                        </div>
                                    </div>
                                    <div class="col s6">
                                        <div class="row">
                                            <div class="col s12">
                                                <h5>Transcript : ${transcript2}</h5>
                                            </div>
                                            <div class="col s12">
                                                <h5>Exon ID : ${exon2_id}</h5>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>`;

        });

    });

}


function get_exon_sequence(transcript_id,exon_number){
        
    $.get(`http://127.0.0.1:5000/get_exons?transcript_id=${transcript_id}&exon_number=${exon_number}`, function(data, status){
        
        document.getElementById("modal_block_transcript_id").innerHTML = "Transcript ID : " + data["transcript_id"];
        document.getElementById("modal_block_exon_id").innerHTML = "Exon ID : " + data["exon_id"];
        document.getElementById("modal_block_exon_sequence").innerHTML = data["exon_sequence"];

        $('#modal_block').modal('open');
        
    })

}

function get_splice_site_similarity(){

    $.get(`http://127.0.0.1:5000/get_splice_site?transcript_id=${transcript1_id}&exon_id=${exon1_id}`, function(data, status){
        
        $.get(`http://127.0.0.1:5000/get_splice_site?transcript_id=${transcript2_id}&exon_id=${exon2_id}`, function(data2, status){

            var ele = document.getElementById("tab_splice_site");

            ele.innerHTML =`<div style="margin-top:5%">
                                <div class="row">
                                    <div class="col s6">
                                        <div class="row">
                                            <div class="col s12">
                                                <h5>Absolute Start : ${data["absolute_start"]}</h5>
                                            </div>
                                            <div class="col s12">
                                                <h5>Absolute End : ${data["absolute_end"]}</h5>
                                            </div>
                                            <div class="col s12">
                                                <h5>Relative Start : ${data["relative_start"]}</h5>
                                            </div>
                                            <div class="col s12">
                                                <h5>Relative End : ${data["relative_end"]}</h5>
                                            </div>
                                            <div class="col s12">
                                                <h5>Strand : ${data["strand"]}</h5>
                                            </div>
                                        </div>
                                    </div>
                                    <div class="col s6">
                                        <div class="row">
                                            <div class="col s12">
                                                <h5>Absolute Start : ${data2["absolute_start"]}</h5>
                                            </div>
                                            <div class="col s12">
                                                <h5>Absolute End : ${data2["absolute_end"]}</h5>
                                            </div>
                                            <div class="col s12">
                                                <h5>Relative Start : ${data2["relative_start"]}</h5>
                                            </div>
                                            <div class="col s12">
                                                <h5>Relative End : ${data2["relative_end"]}</h5>
                                            </div>
                                            <div class="col s12">
                                                <h5>Strand : ${data2["strand"]}</h5>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>`;
                            
        });
    });
}

function get_protein_similarity(){

    $.get(`http://127.0.0.1:5000/get_protein_similarity?transcript1_id=${transcript1_id}&transcript2_id=${transcript2_id}&exon1_id=${exon1_id}&exon2_id=${exon2_id}`, function(data, status){

        
        
        alignment_1 = data[0];
        alignment_2 = data[1];
        
        var ele = document.getElementById("tab_protein");

        ele.innerHTML = `<p>${alignment_1}</p><p>${alignment_2}</p>`;
        
    }); 

}


function clear_modal_line(){

    var ele = document.getElementById("tabs");
    
    ele.innerHTML = "";
    
    ele.innerHTML = `   <ul class="tabs">
                            <li class="tab col s4" onclick="get_splice_site_similarity()"><a href="#tab_splice_site">Splice Site Similarity</a></li>
                            <li class="tab col s4" onshow=""><a href="#tab_transcript">Transcript Similarity</a></li>
                            <li class="tab col s4" onclick="get_protein_similarity()"><a href="#tab_protein">Protein Similarity</a></li>                    
                        </ul>`;
}