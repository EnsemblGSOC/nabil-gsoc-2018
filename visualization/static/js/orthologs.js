var transcript_id = '';
var protein_coding;
var gencode_basic;
var scores = {};
var mode;
var match_score;
var mismatch_penalty;
var gap_open;
var gap_extend;
var skip_penalty;
var blosum_name;


function init(inp_transcript_id, inp_protein_coding, inp_gencode_basic, inp_match_score, inp_mismatch_penalty, inp_gap_open, inp_gap_extend, inp_skip_penalty, blosum ,weight_mode){

    transcript_id = inp_transcript_id;
    protein_coding = inp_protein_coding;
    gencode_basic = inp_gencode_basic;
    mode = weight_mode;
    match_score = inp_match_score;
    mismatch_penalty = inp_mismatch_penalty;
    gap_open = inp_gap_open;
    gap_extend = inp_gap_extend;
    skip_penalty = inp_skip_penalty;
    blosum_name =blosum;

    document.getElementById('in_match_score').value = match_score;
    document.getElementById('in_mismatch_penalty').value = mismatch_penalty;
    document.getElementById('in_gap_open').value = gap_open;
    document.getElementById('in_gap_extend').value = gap_extend;
    document.getElementById('in_skip_penalty').value = skip_penalty;

    $(document).ready(function() {
        M.updateTextFields();
    });
            
    refresh();
}


function refresh(){

	$.get(`http://127.0.0.1:5000/get_orthologs?transcript_id=${transcript_id}&protein_coding=${protein_coding}&gencode_basic=${gencode_basic}`, function(data, status){
            
        document.getElementById(mode).checked = true;

        var ele = document.getElementById("table_body");
        ele.innerHTML= '';

        var table_data = '';

        scores = {};

        for(var i=0;i<data["data"].length;i++){

            extractScore(transcript_id ,data["data"][i]["transcript_id"])
        }

        var pooling = setInterval( ()=>{            

            var progress = document.getElementById("progress");
            progress.style["width"] = (Math.round((Object.keys(scores).length/data["data"].length)*100)).toString()+'%';
            
            if(Object.keys(scores).length==data["data"].length){
                
                clearInterval(pooling);                

                var scoreSet = new Set();

                for(var i=0;i<data["data"].length;i++){
                    scoreSet.add(scores[data["data"][i]["transcript_id"]]);
                }
                
                sortedScores = Array.from(scoreSet);
                
                sortedScores.sort();
                sortedScores.reverse();                
                
                for(var j=0;j<sortedScores.length;j++){
                    
                    for(var i=0;i<data["data"].length;i++){

                        if(scores[data["data"][i]["transcript_id"]] == sortedScores[j]){                        

                            table_data += `<tr><td>${data["data"][i]["transcript_name"]}</td>
                                    <td>${data["data"][i]["transcript_id"]}</td>
                                    <td>${data["data"][i]["species"]}</td>`;                     
                        
                            if(data["data"][i]["biotype"]==='protein_coding'){
                                table_data += `<td><span class="new badge green accent-4" data-badge-caption="${data["data"][i]["biotype"]}"></span></td>`;
                            }
                            else{
                                table_data += `<td></td>`;
                            }

                            if(data["data"][i]["gencode_basic"]==='GENCODE basic'){
                                table_data += `<td><span class="new badge black" data-badge-caption="${data["data"][i]["gencode_basic"]}"></span></td>`;
                            }
                            else{
                                table_data += `<td></td>`;
                            }
                            
                            // This suffices as we have only 2 species, otherwise for cosidering more species we should use dictionaries and mappings
                            if(data["data"][i]["species"]==='human'){
                                table_data +=`<td><a class="waves-effect waves-light btn btn-flat" href="http://www.ensembl.org/Homo_sapiens/Transcript/Exons?t=${data["data"][i]["transcript_id"]}" target="_blank"><i class="material-icons">open_in_new</i></a></td>`;
                            }

                            else if(data["data"][i]["species"]==='mouse'){
                                table_data +=`<td><a class="waves-effect waves-light btn btn-flat" href="http://www.ensembl.org/Mus_musculus/Transcript/Exons?t=${data["data"][i]["transcript_id"]}" target="_blank"><i class="material-icons">open_in_new</i></a></td>`;
                            }

                            table_data +=`<td><a class="waves-effect waves-light pink btn" href="visualize?transcript1=${transcript_id}&transcript2=${data["data"][i]["transcript_id"]}&match_score=${match_score}&mismatch_penalty=${mismatch_penalty}&gap_open=${gap_open}&gap_extend=${gap_extend}&skip_penalty=${skip_penalty}&blosum=blosum62&weight_mode=${mode}" target="_blank">Similarity</a></td>`;

                            table_data +=`<td>${scores[data["data"][i]["transcript_id"]]}</td></tr>`;
                        }
                    }
                }                        

                ele.innerHTML = table_data;
            }

        }, 500);
                
    });

}


function extractScore(transcript1, transcript2){
        
    $.get(`http://127.0.0.1:5000/pair_exons?transcript1_id=${transcript1}&transcript2_id=${transcript2}&match_score=${match_score}&mismatch_penalty=${mismatch_penalty}&gap_open=${gap_open}&gap_extend=${gap_extend}&skip_penalty=${skip_penalty}&weight_mode=${mode}`, function(data, status){
                       
        scores[transcript2] = data["best_score"];
        
    });

}


function change_mode(new_mode){

    mode = new_mode;
            
}


function update(){

    match_score = document.getElementById('in_match_score').value;
    mismatch_penalty = document.getElementById('in_mismatch_penalty').value;
    gap_open = document.getElementById('in_gap_open').value;
    gap_extend = document.getElementById('in_gap_extend').value;
    skip_penalty = document.getElementById('in_skip_penalty').value;

    document.location.href = `/orthologs?transcript=${transcript_id}&protein_coding=${protein_coding}&gencode_basic=${gencode_basic}&match_score=${match_score}&mismatch_penalty=${mismatch_penalty}&gap_open=${gap_open}&gap_extend=${gap_extend}&skip_penalty=${skip_penalty}&weight_mode=${mode}`;
}