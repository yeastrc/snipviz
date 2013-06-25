

//   snp_viewer_loader_data_object_provided.js

//  This runs the SNP viewer

//  This is for when the data is passed in



	// JavaScript directive:   all variables have to be declared with "var", maybe other things

	"use strict";


///////////////////////////////

//   Params:




var createSNPViewerLocalData = function( config, data ) {

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



	// !!! Run the rest of this on a very short delay so that it can run as long as necessary without causing a mesage to the user that it is running too long

	setTimeout( function() {


			try {

				if ( data.nodata === true ) {

					noDataCallBack();

				} else if ( data.result !== "success" ) {

					failToLoadCallBack();


				} else if( data.data !== undefined && data.data.newick !== undefined  && data.data.strainSequences !== undefined ) {

					SNPViewer.createSNPViewer( rootDivId, { newick: data.data.newick, sequences: data.data.strainSequences, requestType: requestType, requestTypeLabel: requestTypeLabel },
						{  },
						{ successfullyLoadedCallback: successfullyLoadedCallback, failToLoadCallBack: failToLoadCallBack }  );



				} else {

					failToLoadCallBack();
				}

			} catch ( t ) {

				var stackTrace = t.stack;

				failToLoadCallBack();

				throw t; // goes to browser since this is a call back function
			}


	}, delayBeforeLoad );  //  run in delayBeforeLoad milliseconds






}


