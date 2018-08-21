var _species ;
var gene_modified = 0;


/**
 * Initialize the page
 * @param {str} species : name of the species ('human' or 'mouse' or 'old' which refers to unchanged)
 */
function init(species){

	gene_modified = 0;

	if(species=='old')
	{
		species = _species;
	}

	else
	{
		_species = species;
	}

	var done = 0;  // flag for progress bar, increases as we load each of the data types like gene id, 
				   // gene name, transcript id, transcript name.
				   // when done equals 4, it means all the data have been loaded

	var ele = document.getElementById("loading");  // control the progressbar
	ele.style.visibility = "visible";

	var protein_coding = document.getElementById('protein_coding').checked;		// check protein coding selected or not
	var gencode_basic = document.getElementById('gencode_basic').checked;		// check gencode basic selected or not
	
	ele = document.getElementById("gene_id");			// clearing the fields 
	ele.value = "";
	ele = document.getElementById("gene_name");
	ele.value = "";
	ele = document.getElementById("transcript_id");
	ele.value = "";
	ele = document.getElementById("transcript_name");
	ele.value = "";

	var loop = setInterval( ()=>{			
							// check if all the data has been loaded
		if(done==4){
			var ele = document.getElementById("loading");
			ele.style.visibility = "hidden";

			clearInterval(loop);
		}
	}, 500 )
	
	$.get("http://127.0.0.1:5000/extract?data_type=gene_id&species="+species, function(data, status){
									// request to the server for gene_ids
		
		var all_options = ""

		for(var i=0; i<data["data"].length; i++){					// adding to options
			all_options += `<option>${data["data"][i]}</option>`
		}
		
		var ele = document.getElementById("gene_id_data");
		ele.innerHTML = all_options;

		done++;

    });

	$.get("http://127.0.0.1:5000/extract?data_type=gene_name&species="+species, function(data, status){
									// request to the server for gene_names
		var all_options = ""

		for(var i=0; i<data["data"].length; i++){					// adding to options
			all_options += `<option>${data["data"][i]}</option>`
		}
		
		var ele = document.getElementById("gene_name_data");
		ele.innerHTML = all_options;

		done++;

	});
	
	$.get("http://127.0.0.1:5000/extract?data_type=transcript_id&species="+species+"&gene_id=all&protein_coding="+protein_coding+"&gencode_basic="+gencode_basic, function(data, status){
									// request to the server for transcript_ids
		var all_options = ""

		for(var i=0; i<data["data"].length; i++){					// adding to options
			all_options += `<option>${data["data"][i]}</option>`
		}
		
		var ele = document.getElementById("transcript_id_data");
		ele.innerHTML = all_options;

		done++;

	});
	
	$.get("http://127.0.0.1:5000/extract?data_type=transcript_name&gene_id=all&species="+species+"&protein_coding="+protein_coding+"&gencode_basic="+gencode_basic, function(data, status){
								// request to the server for transcript_names
		var all_options = ""

		for(var i=0; i<data["data"].length; i++){					// adding to options
			all_options += `<option>${data["data"][i]}</option>`
		}
		
		var ele = document.getElementById("transcript_name_data");
		ele.innerHTML = all_options;

		done++;
		
	});
			
}


/**
 * Action of the submit button
 */
function submit(){
	
	var transcript_id = document.getElementById('transcript_id').value;
	var protein_coding = document.getElementById('protein_coding').checked;
	var gencode_basic = document.getElementById('gencode_basic').checked;

			// load the orthologs page with default parameters
	document.location.href = `/orthologs?transcript=${transcript_id}&protein_coding=${protein_coding}&gencode_basic=${gencode_basic}&match_score=1.0&mismatch_penalty=-1.0&gap_open=-0.5&gap_extend=-0.3&skip_penalty=-1&blosum=blosum62&weight_mode=gaussian`;

}


/***
 * Action of changing gene name
 */
$('#gene_name').on('input', function() { 
	var val =$(this).val(); // get the current value of the input field.

	if((val==='') && (gene_modified == 1)){
		init('old');
	}

	else{

		var protein_coding = document.getElementById('protein_coding').checked;
		var gencode_basic = document.getElementById('gencode_basic').checked;

		$.get(`http://127.0.0.1:5000/valid_gene?gene_name=${val}&species=${_species}&protein_coding=${protein_coding}&gencode_basic=${gencode_basic}`, function(data, status){
															// requesting server for validation
			if(data["valid"]){				

				var loading = document.getElementById("loading");
				loading.style.visibility = "visible";

				var ele = document.getElementById("gene_id");
				ele.value = data["gene_id"];
				ele = document.getElementById("transcript_id");
				ele.value = "";
				ele = document.getElementById("transcript_name");
				ele.value = "";				

				var all_options = ""

				for(var i=0; i<data["transcript_ids"].length; i++){
					all_options += `<option>${data["transcript_ids"][i]}</option>`
				}

				var ele = document.getElementById("transcript_id_data");
				ele.innerHTML = all_options;

				var all_options = ""

				for(var i=0; i<data["transcript_names"].length; i++){
					all_options += `<option>${data["transcript_names"][i]}</option>`
				}
			
				var ele = document.getElementById("transcript_name_data");
				ele.innerHTML = all_options;

				loading.style.visibility = "hidden";

				gene_modified = 1;

				$(document).ready(function() {
					M.updateTextFields();
				  });
			}									
		});
	}	  
});


/***
 * Action of changing gene id
 */

$('#gene_id').on('input', function() { 
	var val =$(this).val(); // get the current value of the input field.

	if((val==='') && (gene_modified == 1)){
		init('old');
	}

	else{

		var protein_coding = document.getElementById('protein_coding').checked;
		var gencode_basic = document.getElementById('gencode_basic').checked;

		$.get(`http://127.0.0.1:5000/valid_gene?gene_id=${val}&species=${_species}&protein_coding=${protein_coding}&gencode_basic=${gencode_basic}`, function(data, status){
																// requesting server for validation
			if(data["valid"]){				

				var loading = document.getElementById("loading");
				loading.style.visibility = "visible";

				var ele = document.getElementById("gene_name");
				ele.value = data["gene_name"];
				ele = document.getElementById("transcript_id");
				ele.value = "";
				ele = document.getElementById("transcript_name");
				ele.value = "";				

				var all_options = ""

				for(var i=0; i<data["transcript_ids"].length; i++){
					all_options += `<option>${data["transcript_ids"][i]}</option>`
				}

				var ele = document.getElementById("transcript_id_data");
				ele.innerHTML = all_options;

				var all_options = ""

				for(var i=0; i<data["transcript_names"].length; i++){
					all_options += `<option>${data["transcript_names"][i]}</option>`
				}
			
				var ele = document.getElementById("transcript_name_data");
				ele.innerHTML = all_options;

				loading.style.visibility = "hidden";

				gene_modified = 1;

				$(document).ready(function() {
					M.updateTextFields();
				  });
			}									
		});
	}
});


/***
 * Action of changing transcript name
 */

$('#transcript_name').on('input', function() { 
	var val =$(this).val(); // get the current value of the input field.
	
	if(val===''){
		return;
	}

	else{

		$.get(`http://127.0.0.1:5000/valid_transcript?transcript_name=${val}`, function(data, status){
												// requesting server for validation
			if(data["valid"]){			

				var loading = document.getElementById("loading");
				loading.style.visibility = "visible";

				var ele = document.getElementById("gene_name");
				ele.value = data["gene_name"];
				ele = document.getElementById("gene_id");
				ele.value = data["gene_id"];;
				ele = document.getElementById("transcript_id");
				ele.value = data["transcript_id"];

				loading.style.visibility = "hidden";

				$(document).ready(function() {
					M.updateTextFields();
				  });
			}									
		});
	}
});


/***
 * Action of changing transcript id
 */
$('#transcript_id').on('input', function() { 
	var val =$(this).val(); // get the current value of the input field.
	
	if(val===''){
		return;
	}

	else{

		$.get(`http://127.0.0.1:5000/valid_transcript?transcript_id=${val}`, function(data, status){
													// requesting server for validation
			if(data["valid"]){					

				var loading = document.getElementById("loading");
				loading.style.visibility = "visible";

				var ele = document.getElementById("gene_name");
				ele.value = data["gene_name"];
				ele = document.getElementById("gene_id");
				ele.value = data["gene_id"];;
				ele = document.getElementById("transcript_name");
				ele.value = data["transcript_name"];

				loading.style.visibility = "hidden";

				$(document).ready(function() {
					M.updateTextFields();
				});
			}									
		});
	}	 
});