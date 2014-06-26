
//   snp_viewer_get_fasta_newick_data_external_file.js

//  This gets the fasta and possibly newick data from external files
//  and then calls createSNPViewerDataAndCallBacksProvided in the file snp_viewer_loader_data_object_provided.js
//  to build the SNP viewer


	// JavaScript directive:   all variables have to be declared with "var", maybe other things

	"use strict";


///////////////////////////////

//   Params:




var getDataCreateSNPViewer = function( config, inputFilenames ) {

	var rootDivId = config.rootDivId;
	var requestType = config.requestType;
	var successfullyLoadedCallback = config.successfullyLoadedCallback; //   successfullyLoadedCallback - the method called when the SNP viewer is successfully loaded and created
	var failToLoadCallBack = config.failToLoadCallBack;                 //   failToLoadCallBack - the method called when the SNP viewer fails to load or fails to be created
	var noDataCallBack = config.noDataCallBack;                 		//   noDataCallBack - the method called when the server returns no data

	var delayBeforeLoad = config.delayBeforeLoad;

	try {

//  This is the data structure that is required.
//  It can be returned in one piece from the server or assembled in Javascript.


		var protocol =  window.location.protocol;

//		protocol: "file:" // When loaded directly from a file

		if ( protocol === "file:" ) {

			alert("This script must be run from a web server.  It is unable to read the newick and fasta files from the local filesystem.  This script will now exit.");

			throw "Cannot run this script directly from a file";
		}

		var strainSequences = { };
		var sequenceLabelsArray = [ ];

		var errorLoadingFile = false;

		var newickFileData = "";
		var fastaFileData = "";


		//  ajax calls set up as not async to simply knowing when everything is available to process

		if ( inputFilenames.newickFilename !== undefined && inputFilenames.newickFilename !== null ) {

			var url = inputFilenames.newickFilename;

			$.ajax({
				type: "GET",
				url: url,
				async: false,
				dataType: "text",
				statusCode: {
				404: function() {
				alert("data file not found: " + url );
				failToLoadCallBack();
				errorLoadingFile = true;
			}
			}
			}).done(function( newickFileDataParam ) {
				newickFileData = newickFileDataParam;
			});

			if ( errorLoadingFile ) {

				return;
			}

		}

		var fastaURL = inputFilenames.fastaFilename;

		$.ajax({
			type: "GET",
			url: fastaURL,
			async: false,
			dataType: "text",
			statusCode: {
				404: function() {
						alert("data file not found: " + url );
						failToLoadCallBack();
						errorLoadingFile = true;
				}
			}
		}).done(function( fastaFileDataParam ) {
			fastaFileData = fastaFileDataParam;
		});

		if ( errorLoadingFile ) {

			return;
		}


		//  Process newick data


		//		Javascript to convert all line endings to \n

		newickFileData = newickFileData.replace(/(\r\n|\r|\n)/g, '\n');

		var newickFileDataLines = newickFileData.split(/\n/g);


		var newickData = "";

		for ( var newickFileDataLinesIndex = 0; newickFileDataLinesIndex < newickFileDataLines.length; newickFileDataLinesIndex++ ) {

			var fastaFileDataLine = newickFileDataLines[ newickFileDataLinesIndex ];

			newickData += fastaFileDataLine;
		}


		//  Process fasta data


		//		Javascript to convert all line endings to \n

		fastaFileData = fastaFileData.replace(/(\r\n|\r|\n)/g, '\n');

		var fastaFileDataLines = fastaFileData.split(/\n/g);

		var prevLineHeaderLine = false;

		var strainName = null;

		var sequenceData = "";

		for ( var fastaFileDataLinesIndex = 0; fastaFileDataLinesIndex < fastaFileDataLines.length; fastaFileDataLinesIndex++ ) {

			var fastaFileDataLine = fastaFileDataLines[ fastaFileDataLinesIndex ];

			//   bypass empty lines
			if ( fastaFileDataLine === null || fastaFileDataLine.length === 0 ) {

				continue;
			}

			if ( fastaFileDataLine.charAt( 0 ) === ">" ) {

				if ( prevLineHeaderLine ) {

					throw "ERROR: Two header lines in a row found, unable to process fasta file, line count = " + ( fastaFileDataLinesIndex + 1);
				}

				prevLineHeaderLine = true;


				//  Save off previous sequence data, if have a strain name

				if ( strainName !== null ) {

					sequenceLabelsArray.push( strainName );

					strainSequences[ strainName ] = sequenceData;

					sequenceData = "";
				}

				//  collect new strain name

				var fastaFileDataLineSplitWhitespace = fastaFileDataLine.split(/\s/g);

				//  remove leading ">"

				strainName = fastaFileDataLineSplitWhitespace[ 0 ].substring( 1 );


			} else {

				prevLineHeaderLine = false;

				sequenceData += fastaFileDataLine;
			}

		}

		//  save off last sequence

		if ( strainName !== null ) {

			sequenceLabelsArray.push( strainName );

			strainSequences[ strainName ] = sequenceData;
		}



		//  build object for SNP Viewer

		var data = { strainSequences: strainSequences };


		if ( inputFilenames.newickFilename !== undefined && inputFilenames.newickFilename !== null ) {

			data.newick = newickData;

		} else {
			data.sequenceLabelsArray = sequenceLabelsArray;
		}

		var snpVizCreateData = { result:"success", data: data };


		createSNPViewerDataAndCallBacksProvided( config, snpVizCreateData );




	} catch ( t ) {


		var stackTrace = t.stack;

		failToLoadCallBack();

		throw t;


	}







};
