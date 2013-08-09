
//   snip_viz_HTML_config_loader.js

//  This gets the configuration data from the HTML attributes.
//  does initial HTML creation and then calls
//  getDataCreateSNPViewer in the file snip_viz_get_fasta_newick_data_external_file.js
//  to build the SNP viewer


	// JavaScript directive:   all variables have to be declared with "var", maybe other things

	"use strict";




$(document).ready(function(){

	var snpViewerRootDivCounter = 0;

	var snpViewerCounter = 0;

	/////////////////////////////////////////////////////

	//

	var createSnpViewersForARootDiv = function( rootElement ) {

		snpViewerRootDivCounter++;

		var $rootElement = $(rootElement);

		var haveDNAFastaFile = false;
		var haveProteinFastaFile = false;

//		<div class="snp-viewer-create-here" snp-viewer-dna-fasta-file=""  snp-viewer-dna-newick-file="" ></div>

		var snp_viewer_dna_fasta_file = $rootElement.attr("snp-viewer-dna-fasta-file");

		var snp_viewer_dna_newick_file = $rootElement.attr("snp-viewer-dna-newick-file");


		var snp_viewer_protein_fasta_file = $rootElement.attr("snp-viewer-protein-fasta-file");

		var snp_viewer_protein_newick_file = $rootElement.attr("snp-viewer-protein-newick-file");

		//  Create the needed internal HTML

		var htmlMainBlock = '<div id="snp-viewer-root-block-' + snpViewerRootDivCounter + '" ></div>';

		var $htmlMainBlock = $( htmlMainBlock ).appendTo( $rootElement );

		if ( snp_viewer_dna_fasta_file !== undefined && snp_viewer_dna_fasta_file !== '' ) {

			haveDNAFastaFile = true;
		}


		if ( snp_viewer_protein_fasta_file !== undefined && snp_viewer_protein_fasta_file !== '' ) {

			haveProteinFastaFile = true;
		}


		if ( ! haveDNAFastaFile && ! haveProteinFastaFile ) {


			$rootElement.html("Neither SNP viewer attribute ('snp-viewer-dna-fasta-file' and 'snp-viewer-protein-fasta-file') found on this div.  Unable to create SNP viewer.");

			return;  //  EXITTING FUNCTION
		}

		//  for divs

		var snpViewerRootBlockDNA_Id = "snp-viewer-root-block-dna-"  + snpViewerRootDivCounter;

		var snpViewerRootBlockProtein_Id = "snp-viewer-root-block-protein-"  + snpViewerRootDivCounter;

		//  for links

		var snpViewerSelectRootBlockDNA_Id = "snp-viewer-select-root-block-dna-" + snpViewerRootDivCounter;

		var snpViewerSelectRootBlockProtein_Id = "snp-viewer-select-root-block-protein-" + snpViewerRootDivCounter;


		if ( haveDNAFastaFile && haveProteinFastaFile ) {



		}



//				'<h1>SNP viewer</h1>' +
//
//		'<h2>SNP viewer ' + snpViewerTypeLabel + '</h2>' +

		if ( haveDNAFastaFile ) {

			var snpViewerRootBlockDNA = '<div id=' + snpViewerRootBlockDNA_Id + ' >';

			if ( haveProteinFastaFile ) {

				snpViewerRootBlockDNA += '[DNA] <a href="javascript:" id="' + snpViewerSelectRootBlockProtein_Id + '" >Protein</a>';
			} else {

				snpViewerRootBlockDNA += "DNA";
			}

			snpViewerRootBlockDNA += "</div>";

			var $dnaSnpViewerRootBlock = $( snpViewerRootBlockDNA ).appendTo( $htmlMainBlock );

			if ( haveProteinFastaFile ) {

				$("#" + snpViewerSelectRootBlockProtein_Id ).click( function(  ) {

					$("#" + snpViewerRootBlockDNA_Id ).hide();
					$("#" + snpViewerRootBlockProtein_Id ).show();
				});
			}




			var snpViewerInsertionPointDNA_Id = "snp-viewer-insertion-point-dna-" + snpViewerRootDivCounter;

			var snpViewerInsertionPointDNA = '<div id=' + snpViewerInsertionPointDNA_Id + ' ></div>';

			var $dnaSnpViewerInsertionPoint = $( snpViewerInsertionPointDNA ).appendTo( $dnaSnpViewerRootBlock );

			var dnaSnpViewerConfig = { requestType: SNPViewer.getRequestTypeDNA() /* requestType */,
				fastaFile: snp_viewer_dna_fasta_file,
				newickFile: snp_viewer_dna_newick_file };

			createASnpViewer( $dnaSnpViewerInsertionPoint, dnaSnpViewerConfig );
		}

		if ( haveProteinFastaFile ) {

			var snpViewerRootBlockProtein = '<div id=' + snpViewerRootBlockProtein_Id + ' ';

			if ( haveDNAFastaFile ) {

				snpViewerRootBlockProtein += ' style="visibility: hidden;" ><a href="javascript:" id="' + snpViewerSelectRootBlockDNA_Id + '" >DNA</a> [Protein] ';
			} else {

				snpViewerRootBlockProtein += " >Protein";
			}

			snpViewerRootBlockProtein += "</div>";

			var $proteinSnpViewerRootBlock = $( snpViewerRootBlockProtein ).appendTo( $htmlMainBlock );

			var $snpViewerRootBlockProtein_Id = $("#" + snpViewerRootBlockProtein_Id );


			$snpViewerRootBlockProtein_Id.css("visibility", "hidden");


			if ( haveDNAFastaFile ) {

				$("#" + snpViewerSelectRootBlockDNA_Id ).click( function(  ) {

					$("#" + snpViewerRootBlockProtein_Id ).hide();
					$("#" + snpViewerRootBlockDNA_Id ).show();
				});
			}





			var snpViewerInsertionPointProtein_Id = "snp-viewer-insertion-point-protein-" + snpViewerRootDivCounter;

			var snpViewerInsertionPointProtein = '<div id=' + snpViewerInsertionPointProtein_Id + ' ></div>';

			var $proteinSnpViewerRoot = $( snpViewerInsertionPointProtein ).appendTo( $proteinSnpViewerRootBlock )

			var proteinSnpViewerConfig = { requestType: SNPViewer.getRequestTypeProtein() /* requestType */,
				fastaFile: snp_viewer_protein_fasta_file,
				newickFile: snp_viewer_protein_newick_file,
				successCallback: function() { $("#" + snpViewerRootBlockProtein_Id ).hide(); } };

			createASnpViewer( $proteinSnpViewerRoot, proteinSnpViewerConfig );



			$("#" + snpViewerRootBlockProtein_Id ).css("visibility","visible");


		}



	};

	/////////////////////////////////////////////////////

	//

	var createASnpViewer = function( $snpViewerRootElement, snpViewerConfig ) {



		//  Increment the counter
		snpViewerCounter++;


		//  Create the needed internal HTML

		var htmlMainBlock = '<div id="snp-viewer-main-block-' + snpViewerCounter + '" >' +


		'<div id="snp-viewer-failed-to-load-' + snpViewerCounter + '" style="display: none; " >' +

		'  <h2>The SNP viewer has failed to load</h2>' +

		'</div>' +

		'<div id="snp-viewer-no-data-' + snpViewerCounter + '" style="display: none; " >' +

		'  <h2>No Data for the SNP viewer</h2>' +

		'</div>' +


		'<div id="snp-viewer-loading-message-' + snpViewerCounter + '" style="display: none; " >' +

		'  Loading SNP viewer' +

		'</div>' +

		'<div>' +

		'<!--' +


		'<!--   The div where the SNP viewer will be created. -->' +
		'<div id="snp-viewer-' + snpViewerCounter + '" >' +


		'</div>' +
		'</div>' +

		'</div>';


		$snpViewerRootElement.append( htmlMainBlock );


		//////////////////////////////

		//  A function that will be called when the SNP viewer has been successfully created

		var successfullyLoadedCallback = function( ) {

			$("#snp-viewer-loading-message-" + snpViewerCounter ).hide();


			$("#snp-viewer-" + snpViewerCounter ).show();

			if ( snpViewerConfig.successCallback !== undefined ) {

				snpViewerConfig.successCallback();
			}
		};

		//  A function that will be called when the SNP viewer fails to be created

		var failToLoadCallBack = function( ) {

			$("#snp-viewer-loading-message-" + snpViewerCounter ).hide();

			$("#snp-viewer-failed-to-load-" + snpViewerCounter ).show();
			$("#snp-viewer-" + snpViewerCounter ).hide();
		};

		//  A function that will be called when the SNP viewer does not find any data

		var noDataCallBack = function( ) {

			$("#snp-viewer-loading-message-" + snpViewerCounter ).hide();

			$("#snp-viewer-no-data-" + snpViewerCounter ).show();
			$("#snp-viewer-" + snpViewerCounter ).hide();
		};


		var config = { rootDivId: "snp-viewer-" + snpViewerCounter,
			genomeClustering: false,
			requestType: snpViewerConfig.requestType /* requestType */,
			successfullyLoadedCallback: function() {
					successfullyLoadedCallback( );

			},
			failToLoadCallBack: function() { failToLoadCallBack();  },
			noDataCallBack: function() { noDataCallBack( );  },
			delayBeforeLoad: 10 };

		var inputFilenames = { fastaFilename: snpViewerConfig.fastaFile };

		if ( snpViewerConfig.newickFile !== undefined && snpViewerConfig.newickFile !== null && snpViewerConfig.newickFile !== '' ) {

			 inputFilenames.newickFilename = snpViewerConfig.newickFile;
		}


		try {

			getDataCreateSNPViewer( config, inputFilenames );


		} catch ( t ) {

			//  Most of the code inside the function called runs inside a "setTimeout( )" so most exceptions will be caught inside that and never make it here


			var stackTrace = t.stack;

			failToLoadCallBack( );

			throw t;
		}


	};





	var protocol =  window.location.protocol;

	//		protocol: "file:" // When loaded directly from a file

	if ( protocol === "file:" ) {


		//			$("#main-block").hide();
		//			$("#page-loaded-from-file-system-error").show();
		//
		//
		//			$("#snp-viewer-loading-message" ).hide();

		var protocolFileMessage = "The script to create the SNP viewer must be run from a web server. " +
		"It is unable to read the fasta and possibly newick files from the local filesystem.  The script will now exit.";

		$(".snp-viewer-create-here").html( protocolFileMessage );

		alert( protocolFileMessage );

		throw "Cannot run this script directly from a file";
	}



	$(".snp-viewer-create-here").each( function( index, Element ) {

		var rootElement = this;

		createSnpViewersForARootDiv( rootElement, Element );

	});





});