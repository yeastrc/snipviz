
//   snp_viewer_get_data_external_file.js

//  This gets the data from external files
//  and then calls createSNPViewerLocalData
//  to build the SNP viewer


	// JavaScript directive:   all variables have to be declared with "var", maybe other things

	"use strict";


///////////////////////////////

//   Params:




var getDataCreateSNPViewer = function( config, geneName ) {

	var rootDivId = config.rootDivId;
	var requestType = config.requestType;
	var successfullyLoadedCallback = config.successfullyLoadedCallback; //   successfullyLoadedCallback - the method called when the SNP viewer is successfully loaded and created
	var failToLoadCallBack = config.failToLoadCallBack;                 //   failToLoadCallBack - the method called when the SNP viewer fails to load or fails to be created
	var noDataCallBack = config.noDataCallBack;                 		//   noDataCallBack - the method called when the server returns no data

	var delayBeforeLoad = config.delayBeforeLoad;

	try {

//  This is the data structure that is required.
//  It can be returned in one piece from the server or assembled in Javascript.


		var strainSequences = { };
		var errorLoadingFile = false;

		var newickFileData = "";
		var fastaFileData = "";


		//  ajax calls set up as not async to simply knowing when everything is available to process

		var url = "newick_files/" + geneName + ".newick";

		$.ajax({
			type: "GET",
			url: url,
			async: false,
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

		url = "fasta_files/" + geneName + ".fa";

		$.ajax({
			type: "GET",
			url: url,
			async: false,
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

			strainSequences[ strainName ] = sequenceData;
		}



		//  build object for SNP Viewer

		var data = {result:"success",
			data:
				{
				newick: newickData,
				strainSequences: strainSequences
				}
			};



		createSNPViewerLocalData( config, data );




	} catch ( t ) {


		var stackTrace = t.stack;

		failToLoadCallBack();

		throw t;


	}







};