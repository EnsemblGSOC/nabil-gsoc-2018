var transcript_id = '';

function init(transcript_id, protein_coding, gencode_basic){

	
	$.get(`http://127.0.0.1:5000/get_orthologs?transcript_id=${transcript_id}&protein_coding=${protein_coding}&gencode_basic=${gencode_basic}`, function(data, status){
		
        console.log(data);
        
        var table_data = '';

        for(var i=0;i<data["data"].length;i++){

            table_data += `<tr><td>${data["data"][i]["transcript_name"]}</td>
                            <td>${data["data"][i]["transcript_id"]}</td>
                            <td>${data["data"][i]["species"]}</td>
                            <td>${data["data"][i]["biotype"]}</td>
                            <td>${data["data"][i]["gencode_basic"]}</td>
                            <td><a class="waves-effect waves-light pink btn" href="visualize?transcript1=${transcript_id}&transcript2=${data["data"][i]["transcript_id"]}" target="_blank">Similarity</a></td></tr>`

        }

        var ele = document.getElementById("table_body");

        console.log(ele);
        

        ele.innerHTML = table_data;
        

    });

}