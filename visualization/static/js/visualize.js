var transcript1 = [];
var transcript2 = [];
var lines = [];
var line_colors=["#99ff99","#3333ff","#cc00ff","#00ffff","#ff6600","#00ff00","#669999"];
var transcript1_id;
var transcript2_id;
var exon1_id;
var exon2_id;

var transcript_scoring;
var protein_scoring;

var mode;

var match_score;
var mismatch_penalty;
var gap_open;
var gap_extend;
var skip_penalty;
var blosum_name;

function init(inp_match_score, inp_mismatch_penalty, inp_gap_open, inp_gap_extend, inp_skip_penalty, blosum, weight_mode){

    match_score = inp_match_score;
    mismatch_penalty = inp_mismatch_penalty;
    gap_open = inp_gap_open;
    gap_extend = inp_gap_extend;
    skip_penalty = inp_skip_penalty;
    mode = weight_mode;
    blosum_name = blosum;


    console.log(mode);
    

    $.get(`http://127.0.0.1:5000/scoring?match_score=${match_score}&mismatch_penalty=${mismatch_penalty}&gap_open=${gap_open}&blosum=${blosum_name}`, function(matrices, status){

        transcript_scoring = matrices["transcript_scores"];
        protein_scoring = matrices["protein_scores"];

        transcript_scoring["maxx"] = matrices["transcript_scores_max"];
        transcript_scoring["minn"] = matrices["transcript_scores_min"];

        protein_scoring["maxx"] = matrices["protein_scores_max"];
        protein_scoring["minn"] = matrices["protein_scores_min"];

        
    })

    transcript1_id = document.getElementById("transcript1_id").innerHTML;
    transcript2_id = document.getElementById("transcript2_id").innerHTML;

    $.get(`http://127.0.0.1:5000/get_species?transcript_id=${transcript1_id}`, function(species1, status){
        
        var link1 = `http://www.ensembl.org/${species1["species"]}/Transcript/Exons?t=${transcript1_id}`;
        document.getElementById("link1").href = link1;
    })

    $.get(`http://127.0.0.1:5000/get_species?transcript_id=${transcript2_id}`, function(species2, status){
        
        var link2 = `http://www.ensembl.org/${species2["species"]}/Transcript/Exons?t=${transcript2_id}`;
        document.getElementById("link2").href = link2;
    })

    $.get(`http://127.0.0.1:5000/pair_exons?transcript1_id=${transcript1_id}&transcript2_id=${transcript2_id}&match_score=${match_score}&mismatch_penalty=${mismatch_penalty}&gap_open=${gap_open}&gap_extend=${gap_extend}&skip_penalty=${skip_penalty}&weight_mode=${mode}`, function(data, status){
        
        console.log(data);
        

        var draw = document.getElementById('draw');
        draw.innerHTML = '';
        var height = draw.clientHeight;
        var width = draw.clientWidth;


        var max_len = Math.max(data["transcript1_total_exons"],data["transcript2_total_exons"]);
       
        
        var unit_height = Math.floor(height*0.1);
        var unit_width = Math.floor((width/max_len)*0.1);
        var exon_width = Math.floor((width - unit_width * (max_len+1)) / max_len);     
        

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

        document.getElementById("loading").style.visibility="hidden";
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


    var ele;
    var tab = document.getElementById("tab3");
    if(tab.classList.contains("active")){
        ele = document.getElementById("tab_protein");
        tab.classList.remove("active");
    }

    var tab = document.getElementById("tab2");
    if(tab.classList.contains("active")){
        ele = document.getElementById("tab_transcript");
        tab.classList.remove("active");
    }

    var tab = document.getElementById("tab1");
    if(tab.classList.contains("active")){
        ele = document.getElementById("tab_splice_site");
        tab.classList.remove("active");
    }
    

    
    $('.tabs').tabs();

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
        console.log(data);
            
        document.getElementById("modal_block_transcript_id").innerHTML = "Transcript ID : " + data["transcript_id"];
        document.getElementById("modal_block_exon_id").innerHTML = "Exon ID : " + data["exon_id"];
        if(data["utr_pos"]==-1){
            document.getElementById("modal_block_exon_sequence").innerHTML = `<span class="orange-text">${data["utr_sequence"]}</span>` + data["exon_sequence"];
        }
        else if(data["utr_pos"]==0){
            document.getElementById("modal_block_exon_sequence").innerHTML =    data["exon_sequence"];
        }

        else{
            document.getElementById("modal_block_exon_sequence").innerHTML = data["exon_sequence"] + `<span class="orange-text">${data["utr_sequence"]}</span>`;
        }
        

        $('#modal_block').modal('open');
        
    })

}

function get_splice_coordinates(){

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

    var ele = document.getElementById("tab_protein");
    
    ele.innerHTML = '';

    $.get(`http://127.0.0.1:5000/get_protein_similarity?transcript1_id=${transcript1_id}&transcript2_id=${transcript2_id}&exon1_id=${exon1_id}&exon2_id=${exon2_id}`, function(data, status){

        
                
        if(typeof(data.error)==="undefined"){        
            alignment_1 = data[0];
            alignment_2 = data[1];
            
            
    
            //ele.innerHTML = `<p>${alignment_1}</p><p>${alignment_2}</p>`;
    
            ele.innerHTML =  pretty_format(alignment_1, alignment_2, protein_scoring);
    
            $('.tooltipped').tooltip();
        }

        else{

            alert(data.error);
        }


        
    }); 

}


function get_transcript_similarity(){

    var ele = document.getElementById("tab_transcript");

    ele.innerHTML = '';
    console.log(`http://127.0.0.1:5000/get_transcript_similarity?transcript1_id=${transcript1_id}&transcript2_id=${transcript2_id}&exon1_id=${exon1_id}&exon2_id=${exon2_id}&match_score=${match_score}&mismatch_penalty=${mismatch_penalty}&gap_open=${gap_open}&gap_extend=${gap_extend}&weight_mode=${mode}`);
    
    $.get(`http://127.0.0.1:5000/get_transcript_similarity?transcript1_id=${transcript1_id}&transcript2_id=${transcript2_id}&exon1_id=${exon1_id}&exon2_id=${exon2_id}&match_score=${match_score}&mismatch_penalty=${mismatch_penalty}&gap_open=${gap_open}&gap_extend=${gap_extend}&weight_mode=${mode}`, function(data, status){

        console.log(data);

        var score = data.score;
        alignment_1 = data.alignment[0];
        alignment_2 = data.alignment[1];
        
        

        ele.innerHTML = ele.innerHTML = `<h4 style="text-align:center; margin-top:5%;">Trascript Similarity Score : ${data.score}</h4>` + pretty_format(alignment_1, alignment_2, transcript_scoring);
        
        $('.tooltipped').tooltip();

    }); 

}

function get_splice_site_similarity(){

    var ele = document.getElementById("tab_splice_site");

    ele.innerHTML = '';    
    
    $.get(`http://127.0.0.1:5000/get_splice_site_similarity?transcript1_id=${transcript1_id}&transcript2_id=${transcript2_id}&exon1_id=${exon1_id}&exon2_id=${exon2_id}&match_score=${match_score}&mismatch_penalty=${mismatch_penalty}&gap_open=${gap_open}&gap_extend=${gap_extend}`, function(data, status){

        console.log(data);

        alignment_1 = data.alignment[0];
        alignment_2 = data.alignment[1];
        
        

        //ele.innerHTML = `<p>${alignment_1}</p><p>${alignment_2}</p>`;

        ele.innerHTML = `<h4 style="text-align:center; margin-top:5%;">Splice Site Similarity Score : ${data.score}</h4>` + pretty_format(alignment_1, alignment_2, transcript_scoring);
        
        $('.tooltipped').tooltip();

    }); 

}



function clear_modal_line(){

    return;
    var ele = document.getElementById("tabs");
    
    
    ele.innerHTML = `<ul class="tabs">
                            <li class="tab col s4" onclick="get_splice_site_similarity()"><a href="#tab_splice_site">Splice Site Similarity</a></li>
                            <li class="tab col s4" onshow=""><a href="#tab_transcript">Transcript Similarity</a></li>
                            <li class="tab col s4" onclick="get_protein_similarity()"><a href="#tab_protein">Protein Similarity</a></li>                    
                     </ul>`;
}

function color_map(val,maxx,minn){

    var delta = maxx-minn;

    var red = Math.round(240 + ( ( (102-240) * (val-minn) ) / delta )) ;
    var green = Math.round(83 + ( ( (231-83) * (val-minn) ) / delta )) ;
    var blue = Math.round(79 + ( ( (120-79) * (val-minn) ) / delta )) ;


    return `rgb(${red},${green},${blue})`;
    
}

function pretty_format(line1, line2, scoring){

    var maxx = scoring["maxx"];
    var minn = scoring["minn"];

    var len = line1.length;

    var chars_in_line = 12;

    var table = '<table style="margin-top:5%"><tbody>'

    var tr = "";

    for(var i=0;i<Math.ceil(len/chars_in_line); i++){
        
        
        tr = "<tr>";
        
        for(var j=0; j < chars_in_line; j++){
            
            if(i*chars_in_line+j>=len){
                break;
            }

            
            var score;

            try{
                score = scoring[line1[i*chars_in_line+j]][line2[i*chars_in_line+j]];
            } catch(err) {                
                score = scoring[line2[i*chars_in_line+j]][line1[i*chars_in_line+j]];
                
            }            

            tr +=  `<td class="col s1 tooltipped" data-position="top" data-tooltip="${score}" style="background-color:${color_map(score,maxx,minn)}">${line1[i*chars_in_line+j]}</td>`
        }

        tr += "</tr>";

        tr += "<tr>";
        
        for(var j=0; j <chars_in_line; j++){

            if(i*chars_in_line+j>=len){
                break;
            }

            var score;

            try{
                score = scoring[line1[i*chars_in_line+j]][line2[i*chars_in_line+j]];
            } catch(err) {
                score = scoring[line2[i*chars_in_line+j]][line1[i*chars_in_line+j]];
            }


            tr +=  `<td class="col s1 tooltipped" data-position="top" data-tooltip="${score}" style="background-color:${color_map(score,maxx,minn)}">${line2[i*chars_in_line+j]}</td>`
        }

        tr += "</tr>";

        tr += "<tr>";
        
        for(var j=0; j <chars_in_line; j++){

            if(i*chars_in_line+j>=len){
                break;
            }

            tr +=  `<td class="col s1"> <span style="color:white">_</span> </td>`
        }

        tr += "</tr>";

        table += tr;
    }

    table += "</tbody></table>";

    
    


    return table;

}

