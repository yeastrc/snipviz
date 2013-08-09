

//   snp_viewer_loader_data_object_provided.js

//  This does final setup before calling
//  SNPViewer.createSNPViewer in the file snp_viewer.js
//  to actually create the Snip Viz

//  This function createSNPViewerDataAndCallBacksProvided can be used directly when
//  writing custom code to create the required Javascript data object and config object.

//  This function is also called by the code in the files
//  snp_viewer_get_fasta_newick_data_external_file.js
//  and
//  snip_viz_loader_data_object_provided.js



	// JavaScript directive:   all variables have to be declared with "var", maybe other things

	"use strict";


///////////////////////////////

//   Params:

//  config = {  rootDivId:  //  A string, the HTML id of the div to create the Snip Viz in
//				requestType:  //  either SNPViewer.getRequestTypeDNA() or SNPViewer.getRequestTypeProtein()
//				successfullyLoadedCallback:  // A function to call after the Snip Viz is successfully created.
//				failToLoadCallBack:  // A function to call after the Snip Viz fails to be created.
//				noDataCallBack:  // A function to call when there is no data.
//				delayBeforeLoad:  // Used in the setTimeout() used to decouple the main Snip Viz creation from the page load execution.
//									//  Set to Zero to run immediately.

//	data = { 	nodata:  //  true if there is no data, if true the noDataCallBack() is called.
//				data: {	newick:  //  newick string
//						sequenceLabelsArray: //  Array of sequence labels
//						strainSequences:  //  object containing the sequences.  The keys are the labels/strain names.  The values are the sequences.

//						***  Either the newick or the sequenceLabelsArray must be populated but not both

var createSNPViewerDataAndCallBacksProvided = function( config, data ) {

	var rootDivId = config.rootDivId;
	var requestType = config.requestType;
	var successfullyLoadedCallback = config.successfullyLoadedCallback; //   successfullyLoadedCallback - the method called when the SNP viewer is successfully loaded and created
	var failToLoadCallBack = config.failToLoadCallBack;                 //   failToLoadCallBack - the method called when the SNP viewer fails to load or fails to be created
	var noDataCallBack = config.noDataCallBack;                 		//   noDataCallBack - the method called when the server returns no data

	var delayBeforeLoad = config.delayBeforeLoad;



	// The test for the existence of the root div on the page also is run in the main Javascript  snp_viewer.js
	//  but is also run here so an error can be thrown without waiting for the response from the server

	var $rootDivId = $("#" + rootDivId );

	if ( $rootDivId.size() === 0 ) {

		throw "Unable to find root div id passed in, rootDivID = " + rootDivId;
	}


	var requestTypeLabel = "";

	if ( requestType === SNPViewer.getRequestTypeDNA() ) {

		requestTypeLabel = "bases";

	} else 	if ( requestType === SNPViewer.getRequestTypeProtein() ) {

		requestTypeLabel = "residues";

	} else {

		throw "requestType must be '" + SNPViewer.getRequestTypeDNA() + "' or '" + SNPViewer.getRequestTypeProtein() + "'.";
	}



	//  display loading message

	$("#snp-viewer-loading-message").show();


	var mainLoadFunction = function() {


			try {

				if ( data.nodata === true ) {

					noDataCallBack();

				} else if ( data.result !== "success" ) {

					failToLoadCallBack();

					throw 'data.result !== "success"';


				} else if( data.data !== undefined &&
						( ( data.data.newick !== undefined  && data.data.sequenceLabelsArray === undefined ) ||
							( data.data.newick === undefined  && data.data.sequenceLabelsArray !== undefined ) 	) &&
						data.data.strainSequences !== undefined ) {

					//  data.data.newick or data.data.sequenceLabelsArray must be provided but not both

					SNPViewer.createSNPViewer( rootDivId,
						{ newick: data.data.newick, sequenceLabelsArray: data.data.sequenceLabelsArray, sequences: data.data.strainSequences,
								requestType: requestType, requestTypeLabel: requestTypeLabel },
						{  },
						{ successfullyLoadedCallback: successfullyLoadedCallback, failToLoadCallBack: failToLoadCallBack }  );



				} else {

					failToLoadCallBack();

					throw 'data provided not in valid or recognized format';
				}

			} catch ( t ) {

				var stackTrace = t.stack;

				failToLoadCallBack();

				throw t; // goes to browser since this is a call back function
			}


	};



	if ( delayBeforeLoad > 0 ) {

		// !!! Run the rest of this on a very short delay so that it can run as long as necessary without causing a mesage to the user that it is running too long

		setTimeout( mainLoadFunction, delayBeforeLoad );  //  run in delayBeforeLoad milliseconds

	} else {

		//  No delay so run immediately

		mainLoadFunction();
	}




}


