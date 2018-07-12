var _species ;

function init(species){


	if(species=='old')
	{
		species = _species;
	}

	else
	{
		_species = species;
	}

	var done = 0;

	var ele = document.getElementById("loading");
	ele.style.visibility = "visible";

	var protein_coding = document.getElementById('protein_coding').checked;
	var gencode_basic = document.getElementById('gencode_basic').checked;
	
	ele = document.getElementById("gene_id");
	ele.value = "";
	ele = document.getElementById("gene_name");
	ele.value = "";
	ele = document.getElementById("transcript_id");
	ele.value = "";
	ele = document.getElementById("transcript_name");
	ele.value = "";

	var loop = setInterval( ()=>{

		if(done==4){
			var ele = document.getElementById("loading");
			ele.style.visibility = "hidden";

			clearInterval(loop);
		}
	}, 500 )
	
	$.get("http://127.0.0.1:5000/extract?data_type=gene_id&species="+species, function(data, status){
		
		var all_options = ""

		for(var i=0; i<data["data"].length; i++){
			all_options += `<option>${data["data"][i]}</option>`
		}
		
		var ele = document.getElementById("gene_id_data");
		ele.innerHTML = all_options;

		done++;

    });

	$.get("http://127.0.0.1:5000/extract?data_type=gene_name&species="+species, function(data, status){
		
		var all_options = ""

		for(var i=0; i<data["data"].length; i++){
			all_options += `<option>${data["data"][i]}</option>`
		}
		
		var ele = document.getElementById("gene_name_data");
		ele.innerHTML = all_options;

		done++;

	});
	
	$.get("http://127.0.0.1:5000/extract?data_type=transcript_id&species="+species+"&gene_id=all&protein_coding="+protein_coding+"&gencode_basic="+gencode_basic, function(data, status){
		
		var all_options = ""

		for(var i=0; i<data["data"].length; i++){
			all_options += `<option>${data["data"][i]}</option>`
		}
		
		var ele = document.getElementById("transcript_id_data");
		ele.innerHTML = all_options;

		done++;

	});
	
	$.get("http://127.0.0.1:5000/extract?data_type=transcript_name&gene_id=all&species="+species+"&protein_coding="+protein_coding+"&gencode_basic="+gencode_basic, function(data, status){
		
		var all_options = ""

		for(var i=0; i<data["data"].length; i++){
			all_options += `<option>${data["data"][i]}</option>`
		}
		
		var ele = document.getElementById("transcript_name_data");
		ele.innerHTML = all_options;

		done++;
		
    });

	
		
}


function submit(){
	
	var transcript_id = document.getElementById('transcript_id').value;
	var protein_coding = document.getElementById('protein_coding').checked;
	var gencode_basic = document.getElementById('gencode_basic').checked;


	document.location.href = `/orthologs?transcript=${transcript_id}&protein_coding=${protein_coding}&gencode_basic=${gencode_basic}`;

}

