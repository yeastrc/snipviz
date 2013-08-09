

//   snp_viewer.js              Primary Snip Viz code

//////////////////////////////////

// JavaScript directive:   all variables have to be declared with "var", maybe other things

"use strict";


////////////////////////////////////////////

//  To create a SNP Viewer instance:
//      call SNPViewer.createSNPViewer( rootDivIdCreateRoot, requestParamsCreateRoot, configParamsCreateRoot, callbackFunctionsObjCreateRoot )

//      rootDivIdCreateRoot = String with the id of the div to build the SNP viewer in

//				In requestParamsCreateRoot, newick or sequenceLabelsArray must be present but not both

//      requestParamsCreateRoot = { newick:  , sequenceLabelsArray:  , sequences: , requestType: , requestTypeLabel: }

//			newick = String containing newick representation of clustering data

//          sequenceLabelsArray = array with the list of labels for the provided sequences

//          sequences = object containing labels and sequences = { "label1": "seq1", "label2": "seq2", ... }

//          requestType = String containing "dna" or "protein".
//                           These strings are available by calling SNPViewer.getRequestTypeDNA() and SNPViewer.getRequestTypeProtein().

//          requestTypeLabel = String containing label for the given request type = "bases" or "residues"

//      configParamsCreateRoot is not currently used

//	    callbackFunctionsObjCreateRoot = { successfullyLoadedCallback:  , failToLoadCallBack:  }


//			successfullyLoadedCallback = Method to call when the SNP viewer is completed loading and building the DOM and listeners

//          failToLoadCallBack = Method to call when the SNP viewer is failed loading and building the DOM and listeners
//               The code may throw an exception rather than call this method.


////////////////////////////////////////////

//  To delete a SNP Viewer Instance:
//      call SNPViewer.deleteSNPViewerFromRootDivId( rootDivIdCreateRoot )

//      rootDivIdCreateRoot = String with the id of the div to delete the SNP viewer



////////////////////////////////////////////


//    Constant values that drive the code

//    Search for "var constants" for the object with properties that drives much of the code




////////////////////////////////////////////


//  Constructor

function SNPViewer() {


	//  Get the "label" used to store the SNPViewer object in the DOM object

	this.getDataElement = function() {
		return "SNPViewer";
	};

	//  Get the "class" used to store the SNPViewer object in the DOM object of the root div

	var getSNPObjectStorageClassName = function() {
		return "SNPObjectStorageClass";
	};



	var snpObjectCounter = 0;



	//  This Javascript object is for computing the score for proteins

	//  BLOSUM80 matrix
	var proteinScoreStructure =
	{
		A:{A:  5, R: -2, N: -2, D: -2, C: -1, Q: -1, E: -1, G:  0, H: -2, I: -2, L: -2, K: -1, M: -1, F: -3, P: -1, S:  1, T:  0, W: -3, Y: -2, V:  0, B: -2, J: -2, Z: -1, X: -1, '*': -6, '-': -10},
		R:{A: -2, R:  6, N: -1, D: -2, C: -4, Q:  1, E: -1, G: -3, H:  0, I: -3, L: -3, K:  2, M: -2, F: -4, P: -2, S: -1, T: -1, W: -4, Y: -3, V: -3, B: -1, J: -3, Z:  0, X: -1, '*': -6, '-': -10},
		N:{A: -2, R: -1, N:  6, D:  1, C: -3, Q:  0, E: -1, G: -1, H:  0, I: -4, L: -4, K:  0, M: -3, F: -4, P: -3, S:  0, T:  0, W: -4, Y: -3, V: -4, B:  5, J: -4, Z:  0, X: -1, '*': -6, '-': -10},
		D:{A: -2, R: -2, N:  1, D:  6, C: -4, Q: -1, E:  1, G: -2, H: -2, I: -4, L: -5, K: -1, M: -4, F: -4, P: -2, S: -1, T: -1, W: -6, Y: -4, V: -4, B:  5, J: -5, Z:  1, X: -1, '*': -6, '-': -10},
		C:{A: -1, R: -4, N: -3, D: -4, C:  9, Q: -4, E: -5, G: -4, H: -4, I: -2, L: -2, K: -4, M: -2, F: -3, P: -4, S: -2, T: -1, W: -3, Y: -3, V: -1, B: -4, J: -2, Z: -4, X: -1, '*': -6, '-': -10},
		Q:{A: -1, R:  1, N:  0, D: -1, C: -4, Q:  6, E:  2, G: -2, H:  1, I: -3, L: -3, K:  1, M:  0, F: -4, P: -2, S:  0, T: -1, W: -3, Y: -2, V: -3, B:  0, J: -3, Z:  4, X: -1, '*': -6, '-': -10},
		E:{A: -1, R: -1, N: -1, D:  1, C: -5, Q:  2, E:  6, G: -3, H:  0, I: -4, L: -4, K:  1, M: -2, F: -4, P: -2, S:  0, T: -1, W: -4, Y: -3, V: -3, B:  1, J: -4, Z:  5, X: -1, '*': -6, '-': -10},
		G:{A:  0, R: -3, N: -1, D: -2, C: -4, Q: -2, E: -3, G:  6, H: -3, I: -5, L: -4, K: -2, M: -4, F: -4, P: -3, S: -1, T: -2, W: -4, Y: -4, V: -4, B: -1, J: -5, Z: -3, X: -1, '*': -6, '-': -10},
		H:{A: -2, R:  0, N:  0, D: -2, C: -4, Q:  1, E:  0, G: -3, H:  8, I: -4, L: -3, K: -1, M: -2, F: -2, P: -3, S: -1, T: -2, W: -3, Y:  2, V: -4, B: -1, J: -4, Z:  0, X: -1, '*': -6, '-': -10},
		I:{A: -2, R: -3, N: -4, D: -4, C: -2, Q: -3, E: -4, G: -5, H: -4, I:  5, L:  1, K: -3, M:  1, F: -1, P: -4, S: -3, T: -1, W: -3, Y: -2, V:  3, B: -4, J:  3, Z: -4, X: -1, '*': -6, '-': -10},
		L:{A: -2, R: -3, N: -4, D: -5, C: -2, Q: -3, E: -4, G: -4, H: -3, I:  1, L:  4, K: -3, M:  2, F:  0, P: -3, S: -3, T: -2, W: -2, Y: -2, V:  1, B: -4, J:  3, Z: -3, X: -1, '*': -6, '-': -10},
		K:{A: -1, R:  2, N:  0, D: -1, C: -4, Q:  1, E:  1, G: -2, H: -1, I: -3, L: -3, K:  5, M: -2, F: -4, P: -1, S: -1, T: -1, W: -4, Y: -3, V: -3, B: -1, J: -3, Z:  1, X: -1, '*': -6, '-': -10},
		M:{A: -1, R: -2, N: -3, D: -4, C: -2, Q:  0, E: -2, G: -4, H: -2, I:  1, L:  2, K: -2, M:  6, F:  0, P: -3, S: -2, T: -1, W: -2, Y: -2, V:  1, B: -3, J:  2, Z: -1, X: -1, '*': -6, '-': -10},
		F:{A: -3, R: -4, N: -4, D: -4, C: -3, Q: -4, E: -4, G: -4, H: -2, I: -1, L:  0, K: -4, M:  0, F:  6, P: -4, S: -3, T: -2, W:  0, Y:  3, V: -1, B: -4, J:  0, Z: -4, X: -1, '*': -6, '-': -10},
		P:{A: -1, R: -2, N: -3, D: -2, C: -4, Q: -2, E: -2, G: -3, H: -3, I: -4, L: -3, K: -1, M: -3, F: -4, P:  8, S: -1, T: -2, W: -5, Y: -4, V: -3, B: -2, J: -4, Z: -2, X: -1, '*': -6, '-': -10},
		S:{A:  1, R: -1, N:  0, D: -1, C: -2, Q:  0, E:  0, G: -1, H: -1, I: -3, L: -3, K: -1, M: -2, F: -3, P: -1, S:  5, T:  1, W: -4, Y: -2, V: -2, B:  0, J: -3, Z:  0, X: -1, '*': -6, '-': -10},
		T:{A:  0, R: -1, N:  0, D: -1, C: -1, Q: -1, E: -1, G: -2, H: -2, I: -1, L: -2, K: -1, M: -1, F: -2, P: -2, S:  1, T:  5, W: -4, Y: -2, V:  0, B: -1, J: -1, Z: -1, X: -1, '*': -6, '-': -10},
		W:{A: -3, R: -4, N: -4, D: -6, C: -3, Q: -3, E: -4, G: -4, H: -3, I: -3, L: -2, K: -4, M: -2, F:  0, P: -5, S: -4, T: -4, W: 11, Y:  2, V: -3, B: -5, J: -3, Z: -3, X: -1, '*': -6, '-': -10},
		Y:{A: -2, R: -3, N: -3, D: -4, C: -3, Q: -2, E: -3, G: -4, H:  2, I: -2, L: -2, K: -3, M: -2, F:  3, P: -4, S: -2, T: -2, W:  2, Y:  7, V: -2, B: -3, J: -2, Z: -3, X: -1, '*': -6, '-': -10},
		V:{A:  0, R: -3, N: -4, D: -4, C: -1, Q: -3, E: -3, G: -4, H: -4, I:  3, L:  1, K: -3, M:  1, F: -1, P: -3, S: -2, T:  0, W: -3, Y: -2, V:  4, B: -4, J:  2, Z: -3, X: -1, '*': -6, '-': -10},
		B:{A: -2, R: -1, N:  5, D:  5, C: -4, Q:  0, E:  1, G: -1, H: -1, I: -4, L: -4, K: -1, M: -3, F: -4, P: -2, S:  0, T: -1, W: -5, Y: -3, V: -4, B:  5, J: -4, Z:  0, X: -1, '*': -6, '-': -10},
		J:{A: -2, R: -3, N: -4, D: -5, C: -2, Q: -3, E: -4, G: -5, H: -4, I:  3, L:  3, K: -3, M:  2, F:  0, P: -4, S: -3, T: -1, W: -3, Y: -2, V:  2, B: -4, J:  3, Z: -3, X: -1, '*': -6, '-': -10},
		Z:{A: -1, R:  0, N:  0, D:  1, C: -4, Q:  4, E:  5, G: -3, H:  0, I: -4, L: -3, K:  1, M: -1, F: -4, P: -2, S:  0, T: -1, W: -3, Y: -3, V: -3, B:  0, J: -3, Z:  5, X: -1, '*': -6, '-': -10},
		X:{A: -1, R: -1, N: -1, D: -1, C: -1, Q: -1, E: -1, G: -1, H: -1, I: -1, L: -1, K: -1, M: -1, F: -1, P: -1, S: -1, T: -1, W: -1, Y: -1, V: -1, B: -1, J: -1, Z: -1, X: -1, '*': -6, '-': -10},
		'*':{A: -6, R: -6, N: -6, D: -6, C: -6, Q: -6, E: -6, G: -6, H: -6, I: -6, L: -6, K: -6, M: -6, F: -6, P: -6, S: -6, T: -6, W: -6, Y: -6, V: -6, B: -6, J: -6, Z: -6, X: -6, '*':  1, '-': -10},
		'-':{A: -10, R: -10, N: -10, D: -10, C: -10, Q: -10, E: -10, G: -10, H: -10, I: -10, L: -10, K: -10, M: -10, F: -10, P: -10, S: -10, T: -10, W: -10, Y: -10, V: -10, B: -10, J: -10, Z: -10, X: -10, '*':  -10, '-': 1}

	};



	//  This Javascript object is for DNA for determining which differences are "silent"

	// DNA codonds by the Amino Acid they generate

	var dnaCodonsByGenAminoAcidInitialArray =
		[
			[ 'ATT', 'ATC', 'ATA' ],
			[ 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG' ],
			[ 'GTT', 'GTC', 'GTA', 'GTG' ],
			[ 'TTT', 'TTC' ],
			[ 'ATG' ],
			[ 'TGT', 'TGC' ],
			[ 'GCT', 'GCC', 'GCA', 'GCG' ],
			[ 'GGT', 'GGC', 'GGA', 'GGG' ],
			[ 'CCT', 'CCC', 'CCA', 'CCG' ],
			[ 'ACT', 'ACC', 'ACA', 'ACG' ],
			[ 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC' ],
			[ 'TAT', 'TAC' ],
			[ 'TGG' ],
			[ 'CAA', 'CAG' ],
			[ 'AAT', 'AAC' ],
			[ 'CAT', 'CAC' ],
			[ 'GAA', 'GAG' ],
			[ 'GAT', 'GAC' ],
			[ 'AAA', 'AAG' ],
			[ 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG' ],
			[ 'TAA', 'TAG', 'TGA' ]
		];



	//  function to load the array into an object of objects for easier lookup

	var dnaCodonsByGenAminoAcidSetup = function(  ) {

		var dnaCodonsByGenAminoAcid = {};

		for ( var outerCount = 0; outerCount < dnaCodonsByGenAminoAcidInitialArray.length; outerCount++ ) {

			var dnaCodonOriginalAssocArray = dnaCodonsByGenAminoAcidInitialArray[ outerCount ];

			for ( var dnaCodonAssocArrayCount = 0; dnaCodonAssocArrayCount < dnaCodonOriginalAssocArray.length; dnaCodonAssocArrayCount++ ) {

				var dnaCodonInArray = dnaCodonOriginalAssocArray[ dnaCodonAssocArrayCount ];

				var dnaCodonNewSubAssocObject = {};

				for ( var copyCount = 0; copyCount < dnaCodonOriginalAssocArray.length; copyCount++ ) {

					if ( copyCount !== dnaCodonAssocArrayCount ) {

						dnaCodonNewSubAssocObject[ dnaCodonOriginalAssocArray[ copyCount ] ] = true;
					}
				}

				dnaCodonsByGenAminoAcid[ dnaCodonInArray ] = dnaCodonNewSubAssocObject;
			}
		}

		return dnaCodonsByGenAminoAcid;
	}


	//  dnaCodonsByGenAminoAcid = { 'codon1': { 'assocCodon1a': true, 'assocCodon1b': true }, 'codon2': { 'assocCodon2a': true, 'assocCodon2b': true }, ... }

	var dnaCodonsByGenAminoAcid = dnaCodonsByGenAminoAcidSetup( );





	////////////////////////

	//   getProteinScoreForTwoLetters

	//   compute the score between 2 protein letters, based on the above matrix

	var getProteinScoreForTwoLetters = function ( letter1, letter2 ) {

		var proteinScoreGroup = proteinScoreStructure[ letter1 ];

		if ( proteinScoreGroup === undefined ) {

			throw "getProteinScoreForTwoLetters( letter1, letter2 ): letter1 is not found in proteinScoreStructure, letter1 = |" + letter1 + "|";
		}

		var proteinScore = proteinScoreGroup[ letter2 ];

		if ( proteinScore === undefined ) {

			throw "getProteinScoreForTwoLetters( letter1, letter2 ): letter2 is not found in proteinScoreGroup, letter1 = |" + letter1 + "|, letter2 = |" + letter2 + "|";
		}


		return proteinScore;
	};

	////////////////////////

	//   getDNACodonAssocCodons

	//   get the object of associated DNA Codons that generate the same Amino Acid

	var getDNACodonAssocCodons = function ( dnaCodon ) {

		var DNACodonAssocCodons = dnaCodonsByGenAminoAcid[ dnaCodon ];

		if ( ! DNACodonAssocCodons ) {

			DNACodonAssocCodons = null;
		}

		return DNACodonAssocCodons;
	};



	////////////////////////////////////////////////

	//   static function to delete a SNP viewer

	//    This uses jQuery to removing everything inside the enclosing div, removing all listeners to those HTML elements in the process

	//       The overall SNPViewerInternal object is attached to an object that is a child div so when that div is removed,
	//          the overall SNPViewerInternal object is free to be garbage collected unless other code on the page is holding a reference.

	//   this.deleteSNPViewerFromRootDivId

	this.deleteSNPViewerFromRootDivId = function ( rootDivIdCreateRoot ) {

		var $rootDivIdCreateRoot = $("#" + rootDivIdCreateRoot );

		//  getting the actual SNPViewer object is optional since it will go away when it is no longer referenced.

//		var SNPObjectStorageNodeSelector = "#" + rootDivIdCreateRoot + " ." + getSNPObjectStorageClassName();
//
//		var $SNPObjectStorageNode = $( SNPObjectStorageNodeSelector );
//
//		var internalVariableSNPViewer = $SNPObjectStorageNode.data( SNPViewer.getDataElement() );

		$rootDivIdCreateRoot.empty();

	};


	//   static function to create a new object

	this.createSNPViewer = function ( rootDivIdCreateRoot, requestParamsCreateRoot, configParamsCreateRoot, callbackFunctionsObjCreateRoot ) {


		//  internal "class" constructor to create actual stored object

		var SNPViewerInternal = function ( rootDivId, requestParams, configParams, callbackFunctionsObj ) {


			snpObjectCounter++;

			var SNPIDCounterAppend = "-" + snpObjectCounter;

			//	callbackFunctionsObj
			//  { successfullyLoadedCallback: successfullyLoadedCallback, failToLoadCallBack: failToLoadCallBack }





			var isLoggingToConsole = function( ) {

				return false;

				if ( console ) {

					return true;
				}
			};


			var logToConsole = function( itemToLog ) {

				if ( console ) {

//					console.log( itemToLog );
				}
			};




			var $rootDivId = $("#" + rootDivId );

			if ( $rootDivId.size() === 0 ) {

				throw "Unable to find root div id passed in, rootDivID = " + rootDivId;
			}





			//  Hide the root div so the user doesn't see it being developed

			$( $rootDivId ).css("visibility", "hidden");








			var newickString = requestParams.newick;
			var sequenceLabelsArray = requestParams.sequenceLabelsArray;
			var sequencesObject = requestParams.sequences;

			var requestType = requestParams.requestType;

			var requestTypeLabel = requestParams.requestTypeLabel;  //  put on the page after "Showing "


			if ( newickString !== undefined && sequenceLabelsArray !== undefined ) {

				throw "requestParams.newick and requestParams.sequenceLabelsArray cannot both be set to a value";
			}

			if ( ( newickString === undefined || newickString === "" ) && ( sequenceLabelsArray === undefined || sequenceLabelsArray.length === 0 ) ) {

				throw "requestParams.newick or requestParams.sequenceLabelsArray must be set to a value, but not both";
			}

			if ( sequencesObject === undefined || sequencesObject === null ) {

				throw "requestParams.sequences must be set to a value";
			}
			if ( requestType === undefined || requestType === "" ) {

				throw "requestParams.requestType must be set to a value";
			}
			if ( requestTypeLabel === undefined || requestTypeLabel === "" ) {

				throw "requestParams.requestTypeLabel must be set to a value";
			}



			if ( SNPViewer.getRequestTypeDNA() !== requestType &&
					SNPViewer.getRequestTypeProtein() !== requestType ) {

				throw "requestParams.requestType must be '" + SNPViewer.getRequestTypeDNA() + "' or '" + SNPViewer.getRequestTypeProtein() + "'.";
			}


			if ( callbackFunctionsObj.successfullyLoadedCallback === undefined || callbackFunctionsObj.successfullyLoadedCallback === null ) {

				throw "callbackFunctionsObj.successfullyLoadedCallback must be set to a value";
			}

			if ( callbackFunctionsObj.failToLoadCallBack === undefined || callbackFunctionsObj.failToLoadCallBack === null ) {

				throw "callbackFunctionsObj.failToLoadCallBack must be set to a value";
			}

						//	callbackFunctionsObj
			//  { successfullyLoadedCallback: successfullyLoadedCallback, failToLoadCallBack: failToLoadCallBack }











			//////////////////

			///    Constants

			var constants = {


				standardDNACodonLength: 3,  //  number of characters in one DNA codon


				minimumDraggableBoxInnerWidth: 8,


				//  optional widths of a single bar where the selector box is, ( For 1 bar = one letter )
				optionalWidthsOfAboveSingleLineOneToOne: [ 4, 3, 2 ], // [ 5, 4, 3 ],

				//  optional widths of a single bar where the selector box is
				optionalWidthsOfAboveSingleLine: [ 2 ], // [ 5, 4, 3 ],

				//  minimum height of a single bar where the selector box is
				minHeightOfAboveSingleLine: 4,


				//  opacity value for DNA silent mutations on a scale of zero to one.
				dnaSilentMutationOpacity: 0.5,

				//  minimum opacity value on a scale of zero to one.
				minOpacity: 0.4,

				//  For Sequence block, the number of columns displayed will be a multiple of this

				seqBlockColumnMultiple: 5,

				//  For Sequence block, remove this from available width to help ensure it doesn't wrap onto the next line

				seqBlockPaddingRight: 5,  //  in pixels

				//  Dendrogram constants

				//  width in pixels
				widthDendrogramOneLevel: 10,

				dendrogramLineColor: "#ff0000", // red


				//  ID constants

				sequenceIdPrefix: "sequence",

				barIdPrefix: "barPos",

				strainIdPrefix: "strain",


				//  attribute in the strain that identifies it's position

				arrayIndexAttrLabel: "arrayIndex",

				//  attribute in the sequenceLetter that identifies it's position

				sequencePositionAttrLabel: "seqPos",


				strainPermHighlightingAttrYes: "YES",

				strainPermHighlightingAttrNo: "NO",

				strainHoverHighlightingAttrYes: "YES",

				strainHoverHighlightingAttrNo: "NO",


				upperBarChartEntryIdPrefix: "upperBarChartEntry-",



				//  a single pixel PNG with a transparent clear background

				transparentPixelPNG: "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAAXNSR0IArs4c6QAAAAtJREFUCNdjYGAAAAADAAEg1ZTHAAAAAElFTkSuQmCC"

			};


			///////////////////////

			//	Classes

			/////////////

			//  this is for processing the block that displays the sequences

			var SequenceDisplayBlock = function (  ) {


				this.parentToUpdateOnChange = undefined;  //   set later

				this.currentLeftSequencePosition = 0;


				///////////////////////////////////////////////////

				//    setParentToUpdateOnChange

				//    Set the parent to update


				this.setParentToUpdateOnChange = function ( parentToUpdateOnChangeNewValue ) {

					this.parentToUpdateOnChange = parentToUpdateOnChangeNewValue;
				}




				///////////////////////////////////////////////////

				//    populateSequenceBlock

				//     put the chunk of the sequences that will fit into the block

				this.populateSequenceBlock = function (  ) {

					if ( this.parentToUpdateOnChange === undefined || this.parentToUpdateOnChange === null ) {

						throw "SequenceDisplayBlock object state error:  this.parentToUpdateOnChange must be set before calling this.populateSequenceBlock()";
					}

					this.populateBarChartDivsOfSeqBlock();

					this.populateMainSequenceBlock();
				};






				///////////////////////////////////////////////////

				//    updateSequenceBlockForStrainChange

				//     update the sequence block for a strain being clicked on ( selected or deselected )

				this.updateSequenceBlockForStrainChange = function (  ) {


					if ( this.parentToUpdateOnChange === undefined || this.parentToUpdateOnChange === null ) {

						throw "SequenceDisplayBlock object state error:  this.parentToUpdateOnChange must be set before calling this.updateSequenceBlockForStrainChange()";
					}

					this.populateBarChartDivsOfSeqBlock();

					this.populateMainSequenceBlock();
				};


				///////////////////////////////////////////////////

				//   private

				//    populateBarChartDivsOfSeqBlock

				//     populate the Bar Chart Divs Of the Sequence Block

				this.populateBarChartDivsOfSeqBlock = function (  ) {

					var count;

					var barChartRoot = $("#barChart" + SNPIDCounterAppend);

					// remove existing data in the block

					barChartRoot.empty();


					//  Place space characters ("&nbsp") in the barChart area


					var barChartDivsToAdd = "";

					for ( count = 0; count < globalVariables.numberOfSeqCharsDisplayed; count++ ) {

						var seqPosition = count + this.currentLeftSequencePosition;

						var score = globalVariables.sequenceVarianceScoresArray[ seqPosition ];

						if ( score > 0 ) {

							var scoreIE8 = score * 100; // 0 - 100;


							barChartDivsToAdd += "<div id='" + constants.barIdPrefix + count + SNPIDCounterAppend + "' class='bar bar-" + requestType + " position-" + count + "  bar-color-where-differences '" +
							" style='opacity: " + score +
							"; filter:alpha(opacity=" + scoreIE8 + "); ' " + // For IE8 and earlier

							">&nbsp</div>";

						} else {

							barChartDivsToAdd += "<div id='" + constants.barIdPrefix + count + SNPIDCounterAppend + "' class='bar bar-" + requestType + "  position-" + count + "'>&nbsp</div>";
						}

						if ( SNPViewer.getRequestTypeDNA() === requestType  &&
								count < ( globalVariables.numberOfSeqCharsDisplayed - 1 ) &&
								seqPosition % constants.standardDNACodonLength === ( constants.standardDNACodonLength - 1 ) ) {

							barChartDivsToAdd += "<div id='" + constants.barIdPrefix + count + "-Codon-seperator"  + SNPIDCounterAppend + "' class='bar bar-dna-codon-break '>&nbsp</div>";
						}
					}

					barChartRoot.append( barChartDivsToAdd );

				};




				///////////////////////////////////////////////////

				//    populateMainSequenceBlock

				//     populate the Main Sequence Block

				this.populateMainSequenceBlock = function (  ) {

					var seqNum;

					var seqArr = globalVariables.sequencesInStrainOrder.sequences;



					var sequencesRoot =  $("#SequencesRoot" + SNPIDCounterAppend);

					// remove existing data in the block

					sequencesRoot.empty();

					var lastSequenceNum = seqArr.length - 1;

					var lastPosition = globalVariables.numberOfSeqCharsDisplayed - 1;




					//  loop through sequences
					for( seqNum = 0; seqNum < seqArr.length; seqNum++ ){

						var strainHighlighted = false;

						var classAddForStrainHighlighted = "";

						var strainItem = globalVariables.strainList[ seqNum ];

						if ( ( strainItem.strainPermHighlighting === constants.strainPermHighlightingAttrYes ) ||
								( strainItem.strainHoverHighlighting === constants.strainHoverHighlightingAttrYes ) ) {

							strainHighlighted = true;
						}



						var sequence = seqArr[ seqNum ].sequence;

						sequencesRoot.append( "<div id='" + constants.sequenceIdPrefix + seqNum + SNPIDCounterAppend + "' class='sequence-group' ></div> " );

						var sequenceDiv = $("#" + constants.sequenceIdPrefix + seqNum + SNPIDCounterAppend );

						var letterDivs = "";


						//  loop through positions
						for( var pos = 0; pos < globalVariables.numberOfSeqCharsDisplayed; pos++ ){


							var seqSingle = sequence.charAt( pos + this.currentLeftSequencePosition );

							var seqPos = pos + this.currentLeftSequencePosition;

							var seqId = createSequenceId( seqNum, seqPos );

							var score = globalVariables.sequenceVarianceScoresArray[ seqPos ];


							letterDivs += "<div id='" + seqId + SNPIDCounterAppend + "' class=' seq-letter seq-letter-" + requestType + " seq-letter-" + requestType + "-" + seqSingle + " seq-letter-" + seqSingle + " ";

							if ( score > 0 ) {

								var scoreIE8 = score * 100; // 0 - 100;

								letterDivs += "  sequence-color-where-differences ";

								if ( strainHighlighted ) {

									letterDivs +=  " seq-letter-strain-selected-" + requestType + "-" + seqSingle + "  sequence-color-where-differences-strain-selected ' ";

								} else {


									letterDivs += " ' style='opacity: " + score + "; filter:alpha(opacity=" + scoreIE8 + "); ' "; // For IE8 and earlier
								}

							} else {

								if ( strainHighlighted ) {

									letterDivs +=  " seq-letter-strain-selected-" + requestType + "-" + seqSingle + " ' ";

								} else {

									letterDivs += " ' ";
								}

							}

							letterDivs +=
								" " + constants.sequencePositionAttrLabel + "='" + seqPos + "' " +
								">" + seqSingle + "</div> ";

							if ( SNPViewer.getRequestTypeDNA() === requestType  &&
									pos < ( globalVariables.numberOfSeqCharsDisplayed - 1 ) &&
									seqPos % constants.standardDNACodonLength === ( constants.standardDNACodonLength - 1 ) ) {

								letterDivs += "<div id='" + seqId + "-Codon-seperator"  + SNPIDCounterAppend + "' class='bar seq-letter-dna-codon-break ";

								if ( strainHighlighted ) {

									letterDivs += "seq-letter-dna-codon-break-strain-selected";
								}


								letterDivs += " '>&nbsp</div>";
							}

						}

						sequenceDiv.append( letterDivs );


						// var strainItem = globalVariables.strainList[ seqNum ];
//
						// if ( ( strainItem.strainPermHighlighting === constants.strainPermHighlightingAttrYes ) ||
								// ( strainItem.strainHoverHighlighting === constants.strainHoverHighlightingAttrYes ) ) {
//
							// highlightSequenceLetters( seqNum );
						// }
					}


					//  enable/disable the left and right arrows by showing the appropriate arrows

					if ( this.currentLeftSequencePosition <= 0 ) {

						//  page left disabled

						$( "#left-scroll-arrow" + SNPIDCounterAppend ).hide();
						$( "#left-scroll-arrow-disabled" + SNPIDCounterAppend ).show();

					} else {

						//  page left enabled

						$( "#left-scroll-arrow-disabled" + SNPIDCounterAppend ).hide();
						$( "#left-scroll-arrow" + SNPIDCounterAppend ).show();

					}



					if ( this.currentLeftSequencePosition >= ( globalVariables.primaryBarChartData.maxSequenceLength - globalVariables.numberOfSeqCharsDisplayed ) ) {

						//  page right disabled

						$( "#right-scroll-arrow" + SNPIDCounterAppend ).hide();
						$( "#right-scroll-arrow-disabled" + SNPIDCounterAppend ).show();

					} else {

						//  page right enabled

						$( "#right-scroll-arrow-disabled" + SNPIDCounterAppend ).hide();
						$( "#right-scroll-arrow" + SNPIDCounterAppend ).show();

					}





					//			var dnaStartIdent = "<div class='sequence-group sequence-style '>";
					//
					//			//  loop through positions
					//			for( var pos = 0; pos < globalVariables.numberOfSeqCharsDisplayed; pos++ ){
					//
					//				var seqPos = pos + this.currentLeftSequencePosition;
					//
					//				if ( ( seqPos % 3 ) === 0 ) {
					//
					//					dnaStartIdent += "<div class='seq-letter seq-letter-" + requestType + "'>@</div>";
					//				} else {
					//
					//					dnaStartIdent += "<div class='seq-letter seq-letter-" + requestType + "'>&nbsp;</div>";
					//				}
					//
					//
					//					if ( SNPViewer.getRequestTypeDNA() === requestType  && seqPos % constants.standardDNACodonLength === ( constants.standardDNACodonLength - 1 ) ) {
					//
					//						dnaStartIdent += "<div id='" + seqId + "-Codon-seperator"  + SNPIDCounterAppend + "' class='bar bar-dna-codon-break '>&nbsp</div>";
					//					}
					//			}
					//
					//			dnaStartIdent += "</div>";
					//
					//			sequencesRoot.append( dnaStartIdent );

				};


				///////////////////////////////////////////////////

				//    addNavArrowsBelowSequenceBlock - called once on SNP viewer creation to set up the navigation arrows


				this.addNavArrowsBelowSequenceBlock = function() {


					//  Add nav arrows below the sequence block.

					//  Assume that the sequence block is initially created at the start of the sequences ( start position zero ).


					var htmlAdditionForArrow = "<div id='nav-arrows-container" + SNPIDCounterAppend + "' class='nav-arrows-container' style='width: 100%' >" +


							' <div  class="left-scroll-arrow-position " >' +


								//  page left disabled

							'  <div style=""  id="left-scroll-arrow-disabled' + SNPIDCounterAppend + '" ' +
								'             class="scroll-arrow-box left-scroll-arrow-contents-disabled " >' +
								'  </div>' +

								//  page left enabled

							'  <div style="display:none "  id="left-scroll-arrow' + SNPIDCounterAppend + '" ' +
								'  onmouseout="UnTip()"   onmouseover="Tip(' + "'Page Left'" + ') "' +
								'         class="scroll-arrow-box  left-scroll-arrow-contents " >' +
								'  </div>' +

							'</div>' +


							'<div  class="right-scroll-arrow-position "  >' +

								//  page right disabled

							'  <div style="display:none"  id="right-scroll-arrow-disabled' + SNPIDCounterAppend + '" ' +
							'           class="scroll-arrow-box right-scroll-arrow-contents-disabled" >' +
							'  </div>' +

								//  page right enabled

							'  <div style="" id="right-scroll-arrow' + SNPIDCounterAppend + '" ' +
							'  onmouseout="UnTip()"   onmouseover="Tip(' + "'Page Right'" + ') "' +
							'                class="scroll-arrow-box right-scroll-arrow-contents " >' +
							'  </div>' +

							'</div>';



					$("#sequenceScrollWindowOuter" + SNPIDCounterAppend ).append( htmlAdditionForArrow );

					//  Attach click handlers to the arrows.

					var seqBlockThis = this;


					$("#left-scroll-arrow" + SNPIDCounterAppend ).click( function( event, ui ) {


							UnTip();

							seqBlockThis.pageLeft();

						} );


					var $rightArrow = $("#right-scroll-arrow" + SNPIDCounterAppend );

					var $rightArrowSize = $rightArrow.size();

					$rightArrow.click( function( event, ui ) {


							UnTip();

							seqBlockThis.pageRight();

						} );



					//   Change width of enclosing div so that right edge of right arrow aligns with right edge
					//   of rightmost character in the sequence block.

					var $navArrowsContainer = $("#nav-arrows-container" + SNPIDCounterAppend );

					var navArrowsContainerInitWidth = $navArrowsContainer.outerWidth();

					var eightyPercentOfNavArrowsContainerInitWidth = navArrowsContainerInitWidth * 0.8 ;


					var $lastSequenceLine = $("#SequencesRoot" + SNPIDCounterAppend + " .sequence-group").last();


					var $lastSequenceLetter = $lastSequenceLine.children(".seq-letter").last();


					var lastSequenceLetterLeft = $lastSequenceLetter.position().left;

					var lastSequenceLetterWidth = $lastSequenceLetter.outerWidth();

					var newNavArrowsDivWidth = lastSequenceLetterLeft + lastSequenceLetterWidth;


					//  ensure new width not larger than initial width ( css of width: 100% ) nor smaller than 80% of the initial width

					if ( newNavArrowsDivWidth > navArrowsContainerInitWidth ) {

						newNavArrowsDivWidth = navArrowsContainerInitWidth;

					} else if ( newNavArrowsDivWidth < eightyPercentOfNavArrowsContainerInitWidth ) {

						newNavArrowsDivWidth = eightyPercentOfNavArrowsContainerInitWidth;
					}


					$navArrowsContainer.css( "width", newNavArrowsDivWidth );

				};




				///////////////////////////////////////////////////

				//    updateForSelectorChangeFromParent

				this.updateForSelectorChangeFromParent = function( leftmostLetterPosition ) {

					this.currentLeftSequencePosition = leftmostLetterPosition;

					this.populateSequenceBlock( );

					var letterCountDisplay = leftmostLetterPosition + 1;

					$("#sequenceStartPos" + SNPIDCounterAppend).text( letterCountDisplay );

					var lastLetterCountDisplay = leftmostLetterPosition + globalVariables.numberOfSeqCharsDisplayed;

					$("#sequenceEndPos" + SNPIDCounterAppend).text( lastLetterCountDisplay );
				};



				///////////////////////////////////////////////////

				//    pageLeft - Page the shown sequences to the left

				this.pageLeft = function( )  {

					this.currentLeftSequencePosition = this.currentLeftSequencePosition - globalVariables.numberOfSeqCharsDisplayed;

					if ( this.currentLeftSequencePosition < 0 ) {

						this.currentLeftSequencePosition = 0;
					}


					this.updateForSelectorChangeFromParent( this.currentLeftSequencePosition );

					this.parentToUpdateOnChange.updateForSelectorChangeFromChild( this.currentLeftSequencePosition );

				};

				///////////////////////////////////////////////////

				//    pageRight - Page the shown sequences to the left

				this.pageRight = function( )  {

					this.currentLeftSequencePosition = this.currentLeftSequencePosition + globalVariables.numberOfSeqCharsDisplayed;

					if ( this.currentLeftSequencePosition > ( globalVariables.primaryBarChartData.maxSequenceLength - globalVariables.numberOfSeqCharsDisplayed ) ) {

						this.currentLeftSequencePosition = ( globalVariables.primaryBarChartData.maxSequenceLength - globalVariables.numberOfSeqCharsDisplayed );
					}

					this.updateForSelectorChangeFromParent( this.currentLeftSequencePosition );

					this.parentToUpdateOnChange.updateForSelectorChangeFromChild( this.currentLeftSequencePosition );

				};

			}




			/////////////

			//   class  BarChartData

			// declare a class to hold the upper bar chart info.

			//  param:   config - configuration for this instance.
			//  param:   isPrimary - is this the "primary" bar chart.


			var BarChartData = function ( config )
			{


				$("#BarChart" + config.idAppended + SNPIDCounterAppend ).toggle();


				//    start processing the config

				this.isPrimary = config.isPrimary;

				this.label = config.label;

				this.idAppended = config.idAppended;         // idAppended - appended to all IDs, must be unique between BarChartData instances

				this.childToUpdateOnChange = config.childToUpdateOnChange;

				this.parentToUpdateOnChange = config.parentToUpdateOnChange;

				this.startingSequencePosition = config.startingSequencePosition; //  zero based

				this.endingSequencePosition = config.endingSequencePosition; //  zero based

				this.parentSequenceDraggableBoundingBoxWidth = config.parentSequenceDraggableBoundingBoxWidth;

				this.sequenceLength = this.endingSequencePosition - this.startingSequencePosition + 1;

				this.maxSequenceLength = this.sequenceLength;

				//    end processing the config


				//  Set variables for the "id" values in the DOM objects

				var sequenceDraggableID = "sequenceDraggable" + this.idAppended;

				var sequenceDraggableImageID = "sequenceDraggableImage" + this.idAppended;

				var sequenceDraggableBoundingBoxID = "sequenceDraggableBoundingBox" + this.idAppended;


				var upperBarChartID = "upperBarChart" + this.idAppended;

				var upperBarChartClickableID = "upperBarChartClickable" + this.idAppended;

				var upperBarChartClickableImageID = "upperBarChartClickableImage" + this.idAppended;


				var startingSequencePositionID = "starting-sequence-position" + this.idAppended;

				var endingSequencePositionID = "ending-sequence-position" + this.idAppended;

				var curPosUpprBarID = "curPosUpprBar" + this.idAppended;


				//   set other instance variables



				//  width of a single bar where the selector box is
				this.widthOfAboveSingleLine = undefined;


				this.sequenceDraggableBoundingBoxWidth = undefined;

				this.sequenceDraggableBoundingBoxWidthMinusSelectorWidths = undefined;

				this.sequenceDraggableBoundingBoxWidthMinusSelectorWidths = undefined;

				this.selectorBoxLeftBorderWidth = undefined;
				this.selectorBoxRightBorderWidth = undefined;

				this.selectorBoxOuterWidth = undefined;

				this.selectorBoxWidthInLetters = undefined;

				this.numberOfAboveLines = undefined;

				this.sequenceLettersPerLineAboveDraggable = undefined;

				this.sequenceLettersPerPixelAboveDraggable = undefined;

				this.rightMostLetterCount = undefined;

				this.maxDifferenceRatio = undefined;


				this.repaintBarWithNewInfoAfterSlide = undefined;


				//////////////////

				//  methods




				/////////////////////////////////////////////////////////////////////

				//  Called from parent bar ( for updates for second bar when first bar is moved )

				this.updateForSelectorChangeFromParent = function ( startingSequencePosition ) {


					this.startingSequencePosition = startingSequencePosition; //  zero based

					this.endingSequencePosition = startingSequencePosition + this.maxSequenceLength - 1; //  zero based

					if ( this.endingSequencePosition > globalVariables.maxSeqLength - 1 ) {

						this.endingSequencePosition = globalVariables.maxSeqLength - 1;

						this.startingSequencePosition = this.endingSequencePosition - this.sequenceLength + 1;
					}



					this.setStartingAndEndingSequencePositionsOnPage();



					var newLeftEdgeOfBoxInLetters = Math.floor( this.sequenceLength / 2 ) - Math.floor( this.selectorBoxWidthInLetters / 2 ) + this.startingSequencePosition;

					this.updateCurrPosOnSliderBox( newLeftEdgeOfBoxInLetters );


					this.updateForFrameMove( this.startingSequencePosition, newLeftEdgeOfBoxInLetters, true /* calledFrom_updateForSelectorChangeFromParent */  );

					this.childToUpdateOnChange.updateForSelectorChangeFromParent( newLeftEdgeOfBoxInLetters );

				};


				/////////////////////////////////////////////////////////////////////

				// called from child bar or sequence block

				this.updateForSelectorChangeFromChild = function ( startingSequencePositionFromChild ) {


					if ( this.isPrimary ) {

						this.updateForSelectorChangeFromChildUpdatePrimary( startingSequencePositionFromChild );

					} else {

						this.updateForSelectorChangeFromChildUpdateSecondary( startingSequencePositionFromChild );

					}
				};

				/////////////////////////////////////////////////////////////////////

				// called from child bar or sequence block - update Primary bar

				this.updateForSelectorChangeFromChildUpdatePrimary = function ( startingSequencePositionFromChild ) {

					//  position box correctly

					var newLeftPixel = 0;

					if ( startingSequencePositionFromChild >= this.rightMostLetterCount ) {

						//  largestPossibleLeft

						newLeftPixel = this.sequenceDraggableBoundingBoxWidth - this.selectorBoxOuterWidth;

					} else {

						newLeftPixel = ( startingSequencePositionFromChild - this.startingSequencePosition ) / this.sequenceLettersPerPixelAboveDraggable;
					}
					$( "#" + sequenceDraggableID + SNPIDCounterAppend ).css( "left", newLeftPixel + "px" );

					//  update position number on slider box

					this.updateCurrPosOnSliderBox( startingSequencePositionFromChild );
				};


				/////////////////////////////////////////////////////////////////////

				//  called from sequence block - update Secondary bar

				this.updateForSelectorChangeFromChildUpdateSecondary = function ( startingSequencePosition ) {


					this.updateWithNewLeftLetterCountWhenSecondary( startingSequencePosition, undefined /* newLeftPixel */,
						undefined /* calledFrom_updateForSelectorChangeFromParent */,
						true /* calledFrom_updateForSelectorChangeFromChild */ );

				};



				/////////////////////////////////////////////////////////////////////

				this.setStartingAndEndingSequencePositionsOnPage = function ( ) {

					var startingSequencePositionDisplay = this.startingSequencePosition + 1;
					var endingSequencePositionDisplay = this.endingSequencePosition + 1;

					$( "#" + startingSequencePositionID + SNPIDCounterAppend ).text( startingSequencePositionDisplay );

					$( "#" + endingSequencePositionID + SNPIDCounterAppend ).text( endingSequencePositionDisplay );



				};

				/////////////////////////////////////////////////////////////////////

				this.setGlobalVariablesFromExistingSizes = function ( ) {

					this.setStartingAndEndingSequencePositionsOnPage();

					var startingSequencePositionDisplay = this.startingSequencePosition + 1;

					this.updateCurrPosOnSliderBox( this.startingSequencePosition );


					// get width of scroll bar

					var sequenceDraggableBoundingBox = $( "#" + sequenceDraggableBoundingBoxID + SNPIDCounterAppend );

					//   If this is the second bar, "this.parentSequenceDraggableBoundingBoxWidth" will not be "undefined"
					if ( this.parentSequenceDraggableBoundingBoxWidth ) {

						//   If this is the second bar, set the width to the width of the first bar so that the second bar will
						//        always be smaller than the first bar
						sequenceDraggableBoundingBox.css( "width", this.parentSequenceDraggableBoundingBoxWidth + "px" );
					}

					this.sequenceDraggableBoundingBoxWidth = sequenceDraggableBoundingBox.width();

					if ( this.sequenceDraggableBoundingBoxWidth === undefined || this.sequenceDraggableBoundingBoxWidth === null ) {

						throw "Unable to find id sequenceDraggableBoundingBox";
					}

					var sequenceDraggableDOMObject = $( "#" + sequenceDraggableID + SNPIDCounterAppend );


					var selectorBoxBorderWidthLeftString = sequenceDraggableDOMObject.css("border-left-width");

					var selectorBoxBorderWidthLeft = parseInt( selectorBoxBorderWidthLeftString, 10 );

					var selectorBoxBorderWidthRightString = sequenceDraggableDOMObject.css("border-right-width");

					var selectorBoxBorderWidthRight = parseInt( selectorBoxBorderWidthRightString, 10 );

					this.selectorBoxLeftBorderWidth = selectorBoxBorderWidthLeft;
					this.selectorBoxRightBorderWidth = selectorBoxBorderWidthRight;


					//  subtract  the selector borders since the bars need to be inside them
					this.sequenceDraggableBoundingBoxWidthMinusSelectorWidths  = this.sequenceDraggableBoundingBoxWidth - this.selectorBoxLeftBorderWidth - this.selectorBoxRightBorderWidth;


				};



				////////////////////////////////////////////////////

				//          isDifferenceAtThisPosition

				//          Determine if there are differences in the sequences at this position

				/////////////////////////////////////////////


				this.isDifferenceAtThisPosition = function ( position ) {

					var differenceAtThisPosition = false;

					if ( globalVariables.strainsSelectedCount <= 1 ) {


						if ( ! globalVariables.sequencesAllMatchAtPositionArray[ position ] ) {

							differenceAtThisPosition = true;
						}


					} else {

						if ( ! globalVariables.sequencesAllMatchAtPositionArray[ position ] ) {

							var initValue = "";

							var prevValue = initValue;

							var seqArr = globalVariables.sequencesInStrainOrder.sequences;


							for( var strainArrayIndex = 0; strainArrayIndex < globalVariables.strainList.length; strainArrayIndex++){

								var strainItem = globalVariables.strainList[ strainArrayIndex ];

								if ( strainItem.strainPermHighlighting !== undefined &&
										strainItem.strainPermHighlighting === constants.strainPermHighlightingAttrYes ) {

									var sequence = seqArr[ strainArrayIndex ].sequence;

									var seqSingle = sequence.charAt( position );

									if ( prevValue !== initValue && prevValue !== seqSingle  ) {

										differenceAtThisPosition = true;

										break;
									}

									prevValue = seqSingle;
								}

							}

						}
					}



					return differenceAtThisPosition;
				};



				///////////////////////////////////////////////////

				//     getMaxScoreForUpperBarIndividBars

				this.getMaxScoreForUpperBarIndividBars = function (  ) {


					var barsAboveMaxScoresArray = [];

					var maxScorePerBar = 0;

					var currSeqNode;
					var seqNum;

					var position = 0;

					var lineAboveCounter = 0;  // the line currently working on

					//////////////////////

					//   determine the previous break in the sequence letters per line based on starting at zero

					var prevSeqBreak = Math.floor( this.startingSequencePosition / this.sequenceLettersPerLineAboveDraggable ) * this.sequenceLettersPerLineAboveDraggable;


					//  the break point to process the next line above,
					//  first "- 1" because "position" variable is zero based,
					//  second "- 1" so that "position"  exceeds "lineAboveNextDelimination" at a point where it can be processed.
					var lineAboveNextDelimination = this.sequenceLettersPerLineAboveDraggable + prevSeqBreak - 1 - 1;

					//				WAS		var lineAboveNextDelimination = this.sequenceLettersPerLineAboveDraggable + this.startingSequencePosition - 1 - 1;


					var positionsToprocess = false;

					var loopStart = this.startingSequencePosition;

					var loopEnd = this.endingSequencePosition;


					//  loop through positions

					for( position = loopStart; position <= loopEnd; position++){

						positionsToprocess = true;

						var positionScore = globalVariables.sequenceVarianceScoresArray[ position ];

						if ( positionScore > maxScorePerBar ) {

							maxScorePerBar = positionScore;
						}


						//  determine if have processed all sequences for current bar above

						if ( position > lineAboveNextDelimination ) {

							//  advance the break point to process the next line above
							lineAboveNextDelimination += this.sequenceLettersPerLineAboveDraggable;

							barsAboveMaxScoresArray[ lineAboveCounter ] = maxScorePerBar;

							maxScorePerBar = 0;

							positionsToprocess = false;

							lineAboveCounter++;
						}
					}

					if ( positionsToprocess ) { // if any positions after the last time the max is stored inside the loop

						barsAboveMaxScoresArray[ lineAboveCounter ] = maxScorePerBar;
					}

					return barsAboveMaxScoresArray;
				};





				///////////////////////////////////////////////////

				//     getRatiosForUpperBarIndividBars

				this.getRatiosForUpperBarIndividBars = function ( params ) {

					var computeMaxOnly = false;

					if ( params && params.computeMaxOnly ) {

						computeMaxOnly = true;
					}

					// declare function for processing the data for one line above

					var processForLineAbove = function ( currentBarChart ) {

						if ( numberOfDifferences > 0 ) {


							var differenceRatio = numberOfDifferences / currentBarChart.sequenceLettersPerLineAboveDraggable; // possibleNumberOfDifferences;

							if ( ! computeMaxOnly ) {

								barsAboveDiffRatiosArray[ lineAboveCounter ] = differenceRatio;
							}

							if ( differenceRatio > maxDifferenceRatio ) {

								maxDifferenceRatio = differenceRatio;
							}

						} else {

							if ( ! computeMaxOnly ) {
								barsAboveDiffRatiosArray[ lineAboveCounter ] = 0;
							}
						}

					};  //  end of internally defined function


					// start code in function getRatiosForUpperBarIndividBars()

					var sequences = globalVariables.sequencesInStrainOrder.sequences;

					var barsAboveDiffRatiosArray = [];

					var currSeqNode;
					var seqNum;

					var barPositionsText = "";

					var sequencesCount = sequences.length;

					var position = 0;

					var lineAboveCounter = 0;  // the line currently working on

					//////////////////////

					//   determine the previous break in the sequence letters per line based on starting at zero

					var prevSeqBreak = Math.floor( this.startingSequencePosition / this.sequenceLettersPerLineAboveDraggable ) * this.sequenceLettersPerLineAboveDraggable;


					//  the break point to process the next line above,
					//  first "- 1" because "position" variable is zero based,
					//  second "- 1" so that "position"  exceeds "lineAboveNextDelimination" at a point where it can be processed.
					var lineAboveNextDelimination = this.sequenceLettersPerLineAboveDraggable + prevSeqBreak - 1 - 1;

					var maxDifferenceRatio = 0;

					var numberOfDifferences = 0;

					var possibleNumberOfDifferences = 0;

					var loopStart = this.startingSequencePosition;

					var loopEnd = this.endingSequencePosition;

					if ( computeMaxOnly ) {

						lineAboveNextDelimination = this.sequenceLettersPerLineAboveDraggable - 1 - 1;

						loopStart = 0;
						loopEnd = globalVariables.maxSeqLength - 1;
					}

					//  loop through positions

					for( position = loopStart; position <= loopEnd; position++){

						possibleNumberOfDifferences++;

						if ( this.isDifferenceAtThisPosition( position ) ) {

							numberOfDifferences++;
						}

						//  determine if have processed all sequences for current bar above

						if ( position > lineAboveNextDelimination ) {

							//  advance the break point to process the next line above
							lineAboveNextDelimination += this.sequenceLettersPerLineAboveDraggable;

							processForLineAbove( this );

							numberOfDifferences = 0;
							possibleNumberOfDifferences = 0;

							lineAboveCounter++;
						}
					}

					if ( possibleNumberOfDifferences > 0 ) {

						processForLineAbove( this );
					}

					if ( computeMaxOnly ) {

						this.maxDifferenceRatio = maxDifferenceRatio;

						return maxDifferenceRatio;
					}

					if ( this.maxDifferenceRatio ) {

						maxDifferenceRatio = this.maxDifferenceRatio;
					}

					//  refactor the diffRatio based on the max difference ratio

					var barsAboveDiffRatiosRefactoredArray = [];

					for ( var counter = 0; counter < barsAboveDiffRatiosArray.length; counter++ ) {

						var differenceRatio = barsAboveDiffRatiosArray[ counter ];

						if ( differenceRatio === 0 ) {

							barsAboveDiffRatiosRefactoredArray.push( differenceRatio );
						} else {
							var refactoredDiffRatio = differenceRatio / maxDifferenceRatio;
							barsAboveDiffRatiosRefactoredArray.push( refactoredDiffRatio );
						}

					}

					return barsAboveDiffRatiosRefactoredArray;
				};





				/////////////////////////////////////////////////////////////////////

				this.updateDraggableBoundingBoxForMaxLength = function ( ) {

					var count;


					//   Prep section of page at the top where the scroll bar will be

					var sequenceDraggableBoundingBox = $( "#" + sequenceDraggableBoundingBoxID + SNPIDCounterAppend );


					var numberOfAboveLinesMatchesMaxSeqLength = false;

					//  look for a line width where there is a one to one ratio of lines to characters

					for ( count = 0; count < constants.optionalWidthsOfAboveSingleLineOneToOne.length; count++ ) {

						this.widthOfAboveSingleLine = constants.optionalWidthsOfAboveSingleLineOneToOne[ count ];


						this.numberOfAboveLines = Math.floor( this.sequenceDraggableBoundingBoxWidthMinusSelectorWidths / this.widthOfAboveSingleLine );

						// reduce by one since the last one would wrap onto the next line
						if ( this.sequenceDraggableBoundingBoxWidthMinusSelectorWidths % this.widthOfAboveSingleLine === 0 ) {

							this.numberOfAboveLines -= 1;
						}

						if ( this.isPrimary && this.sequenceLength < this.numberOfAboveLines ) {

							//  handle where the number of characters is less than the number of bars.

							//   If the number of characters is less than the number of bars, shrink the above box to the number of characters

							this.numberOfAboveLines = this.sequenceLength;

							//  update global variable for sequenceLettersPerLineAboveDraggable

							this.sequenceLettersPerLineAboveDraggable = 1;

							this.sequenceLettersPerPixelAboveDraggable = 1 / this.widthOfAboveSingleLine;


							numberOfAboveLinesMatchesMaxSeqLength = true;

							break;
						}
					}


					if ( ! numberOfAboveLinesMatchesMaxSeqLength ) {

						//  Unable to find a 1 to 1 ratio so:

						//  Find the best choice of line width and shrink the bar above so that the number of characters to a bar is an integer.

						var largestNumberOfBars = 0;

						var largestNumberOfBarsMatchingBarWidth = 0;


						for ( count = 0; count < constants.optionalWidthsOfAboveSingleLine.length; count++ ) {

							var optionalWidthOfAboveSingleLine = constants.optionalWidthsOfAboveSingleLine[ count ];


							var numberOfBarsMaxForThisWidth = Math.floor( this.sequenceDraggableBoundingBoxWidthMinusSelectorWidths / optionalWidthOfAboveSingleLine );

							var lettersPerBarBeforeCeil = this.sequenceLength / numberOfBarsMaxForThisWidth;
							var lettersPerBar = Math.ceil( lettersPerBarBeforeCeil );

							var numberOfBarsBeforeCeil = this.sequenceLength / lettersPerBar;

							var numberOfBars = Math.ceil( numberOfBarsBeforeCeil );

							if ( count === 0 ) {

								largestNumberOfBars = numberOfBars;

								largestNumberOfBarsMatchingBarWidth = optionalWidthOfAboveSingleLine;

							} else if ( numberOfBars < largestNumberOfBars ) {

								largestNumberOfBars = numberOfBars;

								largestNumberOfBarsMatchingBarWidth = optionalWidthOfAboveSingleLine;
							}
						}

						this.widthOfAboveSingleLine = largestNumberOfBarsMatchingBarWidth;

						// compute the number of bars based on an integer number of letters per bar

						var numberOfBarsMaxBeforeFloor = this.sequenceDraggableBoundingBoxWidthMinusSelectorWidths / this.widthOfAboveSingleLine;
						var numberOfBarsMax = Math.floor( numberOfBarsMaxBeforeFloor );

						var lettersPerBarFinalBeforeCeil = this.sequenceLength / numberOfBarsMax;
						var lettersPerBarFinal = Math.ceil( lettersPerBarFinalBeforeCeil );

						var numberOfBarsFinalBeforeCeil = this.sequenceLength / lettersPerBarFinal;
						var numberOfBarsFinal = Math.ceil( numberOfBarsFinalBeforeCeil );


						this.numberOfAboveLines = numberOfBarsFinal;

						//				this.numberOfAboveLines = Math.floor( this.sequenceDraggableBoundingBoxWidthMinusSelectorWidths / this.widthOfAboveSingleLine );



						//  update global variable for sequenceLettersPerLineAboveDraggable

						this.sequenceLettersPerLineAboveDraggable = lettersPerBarFinal;


						this.sequenceLettersPerPixelAboveDraggable = lettersPerBarFinal / this.widthOfAboveSingleLine;


					}


					//////////

					//  compute maxDifferenceRatio across the whole sequence for this bar, done here since this is where all the globals needed are finally set


					this.maxDifferenceRatio = this.getRatiosForUpperBarIndividBars( { computeMaxOnly: true } );


					////////////


					//  Fix the width of boxes so they match the new width

					this.sequenceDraggableBoundingBoxWidth = ( this.widthOfAboveSingleLine * this.numberOfAboveLines ) + this.selectorBoxLeftBorderWidth + this.selectorBoxRightBorderWidth;

					//  subtract  the selector borders since the bars need to be inside them
					this.sequenceDraggableBoundingBoxWidthMinusSelectorWidths  = this.sequenceDraggableBoundingBoxWidth - this.selectorBoxLeftBorderWidth - this.selectorBoxRightBorderWidth;


					//  Fix positioning of "upperBarChart" to match left width of selector box

					$( "#" + upperBarChartID + SNPIDCounterAppend ).css("left", this.selectorBoxLeftBorderWidth + "px" );

					sequenceDraggableBoundingBox.css( "width", this.sequenceDraggableBoundingBoxWidth + "px" );

					if ( ! this.isPrimary ) {

						//  Also shrink outer box to width ( outer box only on second box )

						var $outerSequenceDraggableBoundingBox = $("#outerSequenceDraggableBoundingBox-2" + SNPIDCounterAppend);

						$outerSequenceDraggableBoundingBox.css( "width", this.sequenceDraggableBoundingBoxWidth + "px" );

						//  Center second box under the first box

						var boxMarginLeft = Math.floor( ( this.parentSequenceDraggableBoundingBoxWidth - this.sequenceDraggableBoundingBoxWidth ) / 2 );

						$("#BarChart-2" + SNPIDCounterAppend).css( "margin-left", boxMarginLeft + "px" );
					}


					var widthOfClickable = this.numberOfAboveLines * this.widthOfAboveSingleLine;


					var upperBarChartClickableDOMObject = $( "#" + upperBarChartClickableID + SNPIDCounterAppend );

					//  set this since this is the visible bounded box
					upperBarChartClickableDOMObject.css( "width", widthOfClickable + "px" );

					$( "#" + upperBarChartClickableImageID + SNPIDCounterAppend ).css( "width", widthOfClickable + "px" );


					//  slide to right to account for width of selector box

					var upperBarChartClickableBorderWidthLeftString = upperBarChartClickableDOMObject.css("border-left-width");

					var upperBarChartClickableBorderWidthLeft = parseInt( upperBarChartClickableBorderWidthLeftString, 10 );

					var upperBarChartClickableLeft = this.selectorBoxLeftBorderWidth - upperBarChartClickableBorderWidthLeft;

					upperBarChartClickableDOMObject.css( "left", upperBarChartClickableLeft + "px" );

					// set box above bar to same width so it lines up.  This box contains the length of the sequence.

					$( "#" + endingSequencePositionID + SNPIDCounterAppend ).parent().css( "width", widthOfClickable + "px" );

					//  put bars in upper bar

					this.addBarsToUpperBar (  );

					///////////////////////////////

					this.rightMostLetterCount = this.sequenceLength - globalVariables.numberOfSeqCharsDisplayed;


				};


				///////////////////////////////////////////////////

				// addBarsToUpperBar


				//   add the bars to the upper chart
				//  Set all the lines where the selector box is to a specific height based on the number of differences for the letters that go into a single line

				this.addBarsToUpperBar = function (  ) {

					var barsAboveDiffRatiosArray = this.getRatiosForUpperBarIndividBars( );

					var barsAboveMaxScoresArray = this.getMaxScoreForUpperBarIndividBars();

					//  Place characters in the Above barChart area

					var upperBarChartDiv = $( "#" + upperBarChartID + SNPIDCounterAppend );

					if ( ! this.maxPossibleBarHeight ) {

						var upperBarChartHtChkID = "upperBarChartHtChk"  + this.idAppended;

						upperBarChartDiv.append( "<div id='" + upperBarChartHtChkID + SNPIDCounterAppend + "'  class='upperBarChartEntry'  >&nbsp;<div>" );

						var firstBar = $("#" + upperBarChartHtChkID + SNPIDCounterAppend );

						this.maxPossibleBarHeight = firstBar.height();

						var initialTopMarginFromCSSText = $("#" + upperBarChartHtChkID + SNPIDCounterAppend ).css("margin-top");

						this.initialTopMarginFromCSS = parseInt( initialTopMarginFromCSSText, 10 );

						upperBarChartDiv.empty();
					}

					upperBarChartDiv.empty();

					var upperBarChartDivsToAdd = "";

					for ( var counter = 0; counter < barsAboveDiffRatiosArray.length && counter < this.numberOfAboveLines; counter++ ) {

						try {

							var differenceRatio = barsAboveDiffRatiosArray[ counter ];

							var maxScoreForBar = barsAboveMaxScoresArray[ counter ];

							var maxScoreForBarForIE = maxScoreForBar * 100;

							upperBarChartDivsToAdd += "<div id='" + constants.upperBarChartEntryIdPrefix + counter + this.idAppended + SNPIDCounterAppend +
							"' style=' width:" +  this.widthOfAboveSingleLine + "px; " +
							" opacity: " + maxScoreForBar +
							"; filter:alpha(opacity=" + maxScoreForBarForIE + "); "; // For IE8 and earlier

							if ( differenceRatio > 0 ) {

								var barHeightAddOn = Math.floor( ( this.maxPossibleBarHeight - constants.minHeightOfAboveSingleLine ) * differenceRatio );

								var barHeight = constants.minHeightOfAboveSingleLine + barHeightAddOn;

								var marginTop = this.maxPossibleBarHeight - barHeight + this.initialTopMarginFromCSS;

								upperBarChartDivsToAdd += " height:" +  barHeight + "px; " +
								" margin-top:" +  marginTop + "px; " +
								"' class='upperBarChartEntry  upperBarChartEntry-Color'";

							} else {

								upperBarChartDivsToAdd += "' class='upperBarChartEntry' ";

							}

							upperBarChartDivsToAdd += ">&nbsp</div>";

						} catch ( t ) {

							11;

							throw t;
						}

					}

					upperBarChartDiv.append( upperBarChartDivsToAdd );

				};

				///////////////////////////////////////////////////

				//     updateDraggableBox


				this.updateDraggableBox = function (  ) {


					//  This will be replaced with a different value for the primary if a secondary bar is created.
					this.selectorBoxWidthInLetters = globalVariables.numberOfSeqCharsDisplayed;



					var sequenceDraggableDOMObject = $( "#" + sequenceDraggableID + SNPIDCounterAppend );


					// change to force left and right margins to zero

					sequenceDraggableDOMObject.css("margin-left", "0px");
					sequenceDraggableDOMObject.css("margin-right", "0px");

					// border pixels to remove from available width
					var borderAndMarginPixels = this.selectorBoxLeftBorderWidth + this.selectorBoxRightBorderWidth;

					//  Adjust the width of the selector box in the aboveLineDraggable

					var newSelectorBoxWidthInPixelsBeforeFloor = globalVariables.numberOfSeqCharsDisplayed / this.sequenceLettersPerPixelAboveDraggable;

					var newSelectorBoxWidthInPixels = Math.floor( newSelectorBoxWidthInPixelsBeforeFloor ); // - borderAndMarginPixels;  (borderAndMarginPixels already factored in?)

					if ( ( ! this.isPrimary ) || newSelectorBoxWidthInPixels >= constants.minimumDraggableBoxInnerWidth ) {

						this.setSelectorBoxWidth( newSelectorBoxWidthInPixels );

					} else {

						this.fixSelectorBoxTooNarrow();
					}



				};


				//////////////////////////////////


				this.setSelectorBoxWidth = function( newSelectorBoxWidth ) {

					var sequenceDraggableDOMObject = $( "#" + sequenceDraggableID + SNPIDCounterAppend );

					sequenceDraggableDOMObject.css( "width", newSelectorBoxWidth + "px" );

					$( "#" + sequenceDraggableImageID + SNPIDCounterAppend ).css( "width", newSelectorBoxWidth + "px" );

					this.selectorBoxOuterWidth =  sequenceDraggableDOMObject.outerWidth();
				};




				//////////////////////////////////


				this.fixSelectorBoxTooNarrow = function() {


					//			var secondDraggableWidthInLetters = this.sequenceLettersPerPixelAboveDraggable * constants.minimumDraggableBoxInnerWidth;
					//
					//			var newDraggableBoxInnerWidth = constants.minimumDraggableBoxInnerWidth;


					//			if ( secondDraggableWidthInLetters < this.numberOfAboveLines ) {
					//
					//				secondDraggableWidthInLetters = this.numberOfAboveLines;
					//
					//				var newDraggableBoxInnerWidthBeforeFloor = this.numberOfAboveLines / this.sequenceLettersPerPixelAboveDraggable;
					//
					//				newDraggableBoxInnerWidth = Math.floor( newDraggableBoxInnerWidthBeforeFloor ); // - borderAndMarginPixels;  (borderAndMarginPixels already factored in?)
					//			}
					//
					//			if ( globalVariables.maxSeqLength > 9000 ) {
					//
					//				secondDraggableWidthInLetters = this.numberOfAboveLines * 2;
					//
					//				var newDraggableBoxInnerWidthBeforeFloor = ( this.numberOfAboveLines * 2 ) / this.sequenceLettersPerPixelAboveDraggable;
					//
					//				newDraggableBoxInnerWidth = Math.floor( newDraggableBoxInnerWidthBeforeFloor ); // - borderAndMarginPixels;  (borderAndMarginPixels already factored in?)
					//			}

					//  special test:    this.sequenceLength / 14  -  Seems to work for bars of width 900 - doesn't seem to work for narrower bars

					var secondDraggableWidthInLetters = Math.floor( this.sequenceLength / 14 );



					var newDraggableBoxInnerWidthBeforeFloor = secondDraggableWidthInLetters / this.sequenceLettersPerPixelAboveDraggable;



					var newDraggableBoxInnerWidth = Math.floor( newDraggableBoxInnerWidthBeforeFloor ); // - borderAndMarginPixels;  (borderAndMarginPixels already factored in?)



					//			if ( secondDraggableWidthInLetters < 800 ) {
					//
					//				secondDraggableWidthInLetters = 800;
					//
					//				var newDraggableBoxInnerWidthBeforeFloor = 800 / this.sequenceLettersPerPixelAboveDraggable;
					//
					//				newDraggableBoxInnerWidth = Math.floor( newDraggableBoxInnerWidthBeforeFloor ); // - borderAndMarginPixels;  (borderAndMarginPixels already factored in?)
					//			}


					this.selectorBoxWidthInLetters = secondDraggableWidthInLetters;

					this.createSecondDraggableBar( secondDraggableWidthInLetters );

					//  TODO  Need to update this if need to resize the second bar if it is too small

					this.setSelectorBoxWidth( newDraggableBoxInnerWidth );


					this.rightMostLetterCount = this.sequenceLength - secondDraggableWidthInLetters - 1;

				};



				//////////////////////////////////



				this.createSecondDraggableBar = function( secondDraggableWidthInLetters ) {

					globalVariables.secondaryBarChartData =
						new BarChartData( { label: 2,
							idAppended: "-2",
							startingSequencePosition: 0 ,  // zero based
							endingSequencePosition: secondDraggableWidthInLetters,  //  zero based
							parentSequenceDraggableBoundingBoxWidth : this.sequenceDraggableBoundingBoxWidth,
							parentToUpdateOnChange: this,
							childToUpdateOnChange: globalVariables.sequenceDisplayBlock } );



					globalVariables.sequenceDisplayBlock.setParentToUpdateOnChange( globalVariables.secondaryBarChartData );


					this.childToUpdateOnChange = globalVariables.secondaryBarChartData;

				};







				//   This hack is not needed at this time due to other solutions

				//     The following hack is to deal with the width of the draggable box being less than 1.

				//  TODO  This hack will require the overall bounding box to change width and the visible box div "#upperBarChartClickable" to change position
				//                 so for now don't do this

				//			//  TODO  Hack to deal with negative width
				//
				//			if (  selectorBoxMarginLeft === 0 && selectorBoxMarginRight === 0 ) {
				//
				//				if ( newSelectorBoxWidth < 0 ) {
				//
				//					if ( newSelectorBoxWidth < -2 ) {
				//
				//						//  shrink box left border to 1 pixel width and right border to zero pixel width
				//
				//						sequenceDraggableDOMObject.css( {"border-left-width": "1px" , "border-right-width": "0px"  });
				//
				//						newSelectorBoxWidth += 3;
				//
				//					} else {
				//
				//						//  shrink box border to 1 pixel width
				//
				//						sequenceDraggableDOMObject.css( {"border-left-width": "1px" , "border-right-width": "1px"  });
				//
				//						newSelectorBoxWidth += 2;
				//					}
				//				}
				//
				//			}


				//   Update the current position on the slider box

				this.updateCurrPosOnSliderBox = function( letterCount ) {

					var letterCountDisplay = letterCount + 1;

					$("#" + curPosUpprBarID + SNPIDCounterAppend ).text( letterCountDisplay );
				};



				//////////////////////////////////////

				//  Update the "current position" number just under the slider box on the page

				this.updateCurPosUpprBar = function ( leftEdgePos ) {

					var letterCount  = Math.floor( this.sequenceLettersPerPixelAboveDraggable * leftEdgePos );

					if ( letterCount > this.rightMostLetterCount ) {

						letterCount = this.rightMostLetterCount;
					}


					var letterCountDisplay = letterCount + this.startingSequencePosition;


					this.updateCurrPosOnSliderBox( letterCountDisplay );
				};



				//////////////////////////////////////

				//   run when the user stops dragging the box in the upper bar


				this.processEndDragUpperBar = function ( ui ) {

					var newLeft = ui.position.left;

					this.updateWithNewLeft( newLeft );

				};



				//////////////////////////////////////

				//  respond to when the user clicks on the bar at the top to move the box and update the page

				this.upperBarChartClickableClickEvt = function ( eventObject, ui ) {

					if ( ! this.isPrimary ) {

						//  stop the current animation if one is running and run the callback on it immediately

						//  passing jumpToEnd = true causes animation to complete and call backs to run

						$( "#" + sequenceDraggableID + SNPIDCounterAppend ).stop( true, true /* [clearQueue] [, jumpToEnd] */ );
						$( "#" + "upperBarChart-slidable-2" + SNPIDCounterAppend ).stop( true, true /* [clearQueue] [, jumpToEnd] */ );

						//  run the callback if it is not cleared
						if ( this.repaintBarWithNewInfoAfterSlide !== undefined ) {

							this.repaintBarWithNewInfoAfterSlide();
						}

					}


					var originalLeft = $( "#" + sequenceDraggableID + SNPIDCounterAppend ).position().left;

					var offsetX = eventObject.offsetX;


					//  compute from other values if 'eventObject.offsetX' not defined
					if (!(offsetX)) {


						//  The 'X' position of the click on the page
						var eventObjectPageX = eventObject.pageX;

						var eventTarget = $(eventObject.target);

						var eventTargetParent = eventTarget.parent();

						//  The distance of the object from the left edge of the page
						var targetOffsetLeft = $(eventObject.target).offset().left;

						offsetX = eventObjectPageX - targetOffsetLeft;
					}

					//  center the box where the user clicked
					var newLeft = offsetX - ( this.selectorBoxOuterWidth / 2 );


					//  mod newLeft to "this.widthOfAboveSingleLine"


					var newLeftDivSingleLine = newLeft / this.widthOfAboveSingleLine;

					var newLeftDivSingleLineFloor = Math.floor( newLeftDivSingleLine );

					var newLeftFloorSingleLine = newLeftDivSingleLineFloor * this.widthOfAboveSingleLine;


					this.updateWithNewLeft( newLeftFloorSingleLine );


				};


				//////////////////////////////////////

				//  Update sequence display area with new left edge ( from drag end or click )

				this.updateWithNewLeft = function( newLeftPixel ) {


					if ( newLeftPixel < 0 ) {

						newLeftPixel = 0;
					}


					var largestPossibleLeft = this.sequenceDraggableBoundingBoxWidth - this.selectorBoxOuterWidth;

					if ( newLeftPixel > largestPossibleLeft ) {

						newLeftPixel = largestPossibleLeft;
					}

					$( "#" + sequenceDraggableID + SNPIDCounterAppend ).css( "left", newLeftPixel + "px" );

					this.updateFromTopSelector( newLeftPixel );
				};

				//////////////////////////////////////

				//  Update the page when the user relocates the box in the upper bar


				this.updateFromTopSelector = function ( newLeftPixel ) {

					if ( this.isPrimary ) {

						this.updateFromTopSelectorWhenPrimary( newLeftPixel );

					} else {

						this.updateFromTopSelectorWhenSecondary( newLeftPixel );
					}
				};


				//////////////////////////////////////

				//  Update the page when the user relocates the box in the upper bar ( top bar / primary bar )


				this.updateFromTopSelectorWhenPrimary = function  ( newLeftPixel ) {

					var newLeftletterCount  = Math.floor( this.sequenceLettersPerPixelAboveDraggable * newLeftPixel );

					if ( newLeftletterCount > this.rightMostLetterCount ) {

						newLeftletterCount = this.rightMostLetterCount;
					}

					var newLeftletterCountDisplay = newLeftletterCount + this.startingSequencePosition;

					this.updateWithNewLeftLetterCountWhenPrimary( newLeftletterCountDisplay );

				};


				//////////////////////////////////////

				//  Update the page when the user relocates the box in the upper bar ( second bar / secondary bar )


				this.updateFromTopSelectorWhenSecondary = function ( newLeftPixel, calledFrom ) {

					var calledFrom_updateForSelectorChangeFromParent = false;

					if ( calledFrom && calledFrom.updateForSelectorChangeFromParent ) {

						calledFrom_updateForSelectorChangeFromParent = true;
					}


					if ( isLoggingToConsole() ) {

						logToConsole( "In this.updateFromTopSelectorWhenSecondary: newLeftPixel = " + newLeftPixel + ", calledFrom_updateForSelectorChangeFromParent = "
							+ calledFrom_updateForSelectorChangeFromParent );

					}



					var newLeftEdgeLetterCount  = Math.floor( this.sequenceLettersPerPixelAboveDraggable * newLeftPixel ) + this.startingSequencePosition;

					if ( newLeftEdgeLetterCount + this.selectorBoxWidthInLetters > globalVariables.maxSeqLength ) {

						newLeftEdgeLetterCount = globalVariables.maxSeqLength - this.selectorBoxWidthInLetters;
					}

					this.updateWithNewLeftLetterCountWhenSecondary( newLeftEdgeLetterCount, newLeftPixel, calledFrom_updateForSelectorChangeFromParent  );

				};


				//////////////////////////////////////

				//  Update with new left edge Letter Count

				this.updateWithNewLeftLetterCount = function( newLeftEdgeLetterCount  ) {

					if ( this.isPrimary ) {

						this.updateWithNewLeftLetterCountWhenPrimary( newLeftEdgeLetterCount );

					} else {

						this.updateWithNewLeftLetterCountWhenSecondary( newLeftEdgeLetterCount  );
					}


				};

				//////////////////////////////////////

				//  Update with new left edge Letter Count When Primary

				this.updateWithNewLeftLetterCountWhenPrimary = function( newLeftletterCount ) {

					var newLeftletterCountDisplay = newLeftletterCount + this.startingSequencePosition;

					this.updateCurrPosOnSliderBox( newLeftletterCountDisplay );


					var newLeftletterCountForChild = newLeftletterCount + this.startingSequencePosition;

					this.childToUpdateOnChange.updateForSelectorChangeFromParent( newLeftletterCountForChild );


				};

				//////////////////////////////////////

				//  Update with new left edge Letter Count When Secondary

				this.updateWithNewLeftLetterCountWhenSecondary = function( newLeftEdgeLetterCount, newLeftPixel, calledFrom_updateForSelectorChangeFromParent, calledFrom_updateForSelectorChangeFromChild  ) {



					this.updateCurrPosOnSliderBox( newLeftEdgeLetterCount );

					this.childToUpdateOnChange.updateForSelectorChangeFromParent( newLeftEdgeLetterCount );


					//  Compute new start and end letter positions of frame

					var newCenterOfBoxInLetters = newLeftEdgeLetterCount + ( Math.floor( this.selectorBoxWidthInLetters / 2 ) );

					var newFrameStartPosition = newCenterOfBoxInLetters - ( Math.floor( this.sequenceLength / 2 ) );

					var newFrameEndPosition = newFrameStartPosition + this.sequenceLength - 1;

					if ( newFrameEndPosition > globalVariables.maxSeqLength ) {

						newFrameStartPosition = globalVariables.maxSeqLength - this.sequenceLength;

						newFrameEndPosition = this.sequenceLength - 1;

					} else if ( newFrameStartPosition < 0 ) {

						newFrameStartPosition = 0;

						newFrameEndPosition = newFrameStartPosition + this.sequenceLength - 1;
					}

					//  compute box location in new frame

					var leftEdgeLetterCountOffsetFromLeftFrameEdge = newLeftEdgeLetterCount - newFrameStartPosition;

					//			WAS  var leftEdgePixelOffset = leftEdgeLetterCountOffsetFromLeftFrameEdge * this.sequenceLettersPerPixelAboveDraggable;

					var leftEdgePixelOffset = leftEdgeLetterCountOffsetFromLeftFrameEdge / this.sequenceLettersPerPixelAboveDraggable;

					//  compute pixel slide from current box position to new position

					var changeInLeft = 0;

					if ( newLeftPixel ) {

						changeInLeft = leftEdgePixelOffset - newLeftPixel;
					}


					var changeInLeftWithOffset = "+=" + Math.abs( changeInLeft );

					if ( changeInLeft < 0 ) {

						changeInLeftWithOffset = "-=" + Math.abs( changeInLeft );
					}

					if ( calledFrom_updateForSelectorChangeFromParent || calledFrom_updateForSelectorChangeFromChild ) {

						this.updateForFrameMove( newFrameStartPosition, newLeftEdgeLetterCount, calledFrom_updateForSelectorChangeFromParent );

					} else {

						if ( newFrameStartPosition ===  this.startingSequencePosition || changeInLeftWithOffset === 0 ) {

									this.updateForFrameMove( newFrameStartPosition, newLeftEdgeLetterCount, calledFrom_updateForSelectorChangeFromParent );

						} else {

							//  slide frame


							/////////////////////


							var barChartThis = this;


							this.repaintBarWithNewInfoAfterSlide = function() {

								//  This function runs as a callback to the "this" is not the BarChart object

								//  reposition bar to original position
								$( "#" + "upperBarChart-slidable-2" + SNPIDCounterAppend ).css("left", "0px");

								//  draw new frame, position box correctly

								barChartThis.updateForFrameMove( newFrameStartPosition, newLeftEdgeLetterCount, calledFrom_updateForSelectorChangeFromParent );

								barChartThis.repaintBarWithNewInfoAfterSlide = undefined;
							};


							$( "#" + sequenceDraggableID + SNPIDCounterAppend ).animate(
								{
									left: changeInLeftWithOffset
								}, 500 );

							$( "#" + "upperBarChart-slidable-2" + SNPIDCounterAppend ).animate(
								{
									left: changeInLeftWithOffset
								}, 500 , this.repaintBarWithNewInfoAfterSlide
							);


						}

					}
				};



				/////////////////////////////////////////////////////////////////////


				this.updateForFrameMove = function ( startingSequencePosition, newLeftEdgeLetterCount, calledFrom_updateForSelectorChangeFromParent ) {


					this.startingSequencePosition = startingSequencePosition; //  zero based

					this.endingSequencePosition = startingSequencePosition + this.maxSequenceLength - 1; //  zero based

					this.setStartingAndEndingSequencePositionsOnPage();



					//  put bars in upper bar

					//			var barsAboveDiffRatiosRefactoredArray = this.getRatiosForUpperBarIndividBars( );

					//			this.addBarsToUpperBar ( barsAboveDiffRatiosRefactoredArray );

					this.addBarsToUpperBar (  );


					//  position box correctly

					var newLeftPixel = ( newLeftEdgeLetterCount - this.startingSequencePosition ) / this.sequenceLettersPerPixelAboveDraggable;


					$( "#" + sequenceDraggableID + SNPIDCounterAppend ).css( "left", newLeftPixel + "px" );


					if ( this.parentToUpdateOnChange && ! calledFrom_updateForSelectorChangeFromParent ) {

						this.parentToUpdateOnChange.updateForSelectorChangeFromChild( startingSequencePosition );

					}

				};


				//////////////////////////////////////


				this.attachEventHandlers = function() {

					var barChartThis = this;

					//  attach the draggable

					$( "#" + sequenceDraggableID + SNPIDCounterAppend ).draggable(  {
						containment: "parent",

						//  BROKEN:  grid of "2" doesn't work if position of selector all the way to the right is an odd number, which happens sometimes.
						//                  It happens in particular when the selector width is an odd number

						//  If want "grid", also need to fix the "click" to place the selector on the grid ( set a value that is valid for the grid )
						//				grid: [ this.widthOfAboveSingleLine, this.widthOfAboveSingleLine ],

						axis: "x",
						drag: function(event, ui) {

						try {
							barChartThis.updateCurPosUpprBar( ui.position.left );
						} catch ( e ) {

							var stracktrace = e.stack;
						}
					} ,
					stop: function(event, ui) {
						try {
							barChartThis.processEndDragUpperBar( ui );
						} catch ( e ) {

							var stracktrace = e.stack;
						}
					}

					});

					//////////////////////////////////////

					//  attach the click event handler


					$( "#" + upperBarChartClickableImageID + SNPIDCounterAppend ).click( function( event, ui ) {
						barChartThis.upperBarChartClickableClickEvt( event, ui );
					}
					);


				};


			}



			/////////////

			//   class  Strain

			// declare a class to hold the info related to a single strain.

			//

			var Strain = function ( name )
			{
				this.name = name;

				this.strainPermHighlighting = undefined;
			}


			//////////////////////////////////

			//   global variables - within this function

			//   globalVariables.

			var globalVariables = {

				sequenceDisplayBlock: undefined,  //   new SequenceDisplayBlock(  ),

				primaryBarChartData: undefined,  //   new BarChartData( 1, "-1", sequenceDisplayBlock  ),

				secondaryBarChartData: undefined,  //   new BarChartData( 2, "-2", sequenceDisplayBlock  ),

				strainList: undefined,

				strainsSelectedCount: 0,


				widthOfOneChar: undefined,

				numberOfSeqCharsDisplayed: undefined, //  In the scroll window below

				maxSeqLength: undefined,



				//  TODO  not used
				leftMarginUpperBarChart: undefined,




				sequencesInStrainOrder: undefined,   //  array of strings

				//  Whether or not all of the letters of the sequences at a given position match
				sequencesAllMatchAtPositionArray: [], // array of boolean



				//  The variance of the letters of the sequences at a given position
				sequenceVarianceScoresArray: undefined, // array of numbers

				//  The previous variance of the letters of the sequences at a given position
				prevSequenceVarianceScoresArray: undefined // array of numbers

			};

			//////////////////

			//  Start of code:






			/////////////////////



			////////////////////////////////////


			//  start functions


			//  create the id that is assigned for a given sequence number ( index of that sequence ) and position in that sequence
			var createSequenceId = function ( seqNum, position ) {

				return "seq" + seqNum + "Pos" + position;
			};



			////////////////////////////////////

			//  parse the newick string

			var parseNewickString = function ( newickStringToParse ) {

				var ancestors = [];
				var tree = {};

				try {


					var prevToken = "";

					var position = 0;


					var subtree;

					var textString;



					while (  position < newickStringToParse.length ) {

						var token = newickStringToParse.charAt( position );

						switch (token) {
						case '(': // new branchset
							subtree = {};
							tree.branchset = [subtree];
							ancestors.push(tree);
							tree = subtree;
							prevToken = token;
							position++;
							break;
						case ',': // another branch
							subtree = {};
							ancestors[ ancestors.length - 1 ].branchset.push(subtree);
							tree = subtree;
							prevToken = token;
							position++;
							break;
						case ')': // optional name next
							tree = ancestors.pop();
							prevToken = token;
							position++;
							break;
						case ':': // optional length next
							prevToken = token;
							position++;
							break;
						default:

							//  Not found a token so found text

							textString = token;
						position++;

						//  get the text
						while (  position < newickStringToParse.length &&
								newickStringToParse.charAt( position ) !== "(" &&
								newickStringToParse.charAt( position ) !== "," &&
								newickStringToParse.charAt( position ) !== ")" &&
								newickStringToParse.charAt( position ) !== ":" ) {

							textString += newickStringToParse.charAt( position );

							position++;
						}

						//  determine if the text is a name or a length
						if ( prevToken === "(" ||
								prevToken === "," ||
								prevToken === ")" ) {

							if ( textString !== ";" ) {

								//  set name
								tree.name = textString;
							}

						} else if ( prevToken === ":" ) {

							//  set length
							tree.length = parseFloat(textString);
						}
						}
					}

				} catch ( t ) {

					var stackTrace = t.stack;

					throw t;
				}
				return tree;
			};

			///////////////////////////////////////////////////

			//    highlightSequenceLetters

			//   highlight a row of sequence letters - swap the font color with the background color


			var highlightSequenceLetters = function ( strainArrayIndex ) {

				var $sequenceDiv = $("#" + constants.sequenceIdPrefix + strainArrayIndex + SNPIDCounterAppend );

				$sequenceDiv.children( ".seq-letter" ).each(function(){

					try {
						var $seqDiv = $(this);

						highlightSequenceLetter( $seqDiv );

					} catch ( t ) {

						var stackTrace = t.stack;

						throw t;
					}

				});


				$sequenceDiv.children( ".seq-letter-dna-codon-break" ).each(function(){

					try {
						var $seqDiv = $(this);

						$seqDiv.addClass("seq-letter-dna-codon-break-strain-selected");

					} catch ( t ) {

						var stackTrace = t.stack;

						throw t;
					}

				});


			};

			///////////////////////////////////////////////////

			//    highlightSequenceLetter

			//   highlight a sequence letter - swap the font color with the background color


			var highlightSequenceLetter = function( $seqDiv ) {



				var seqDivLetter = $seqDiv.html();


				$seqDiv.addClass( "seq-letter-strain-selected-" + requestType + "-" + seqDivLetter );


				var seqPos = $seqDiv.attr( constants.sequencePositionAttrLabel );

				var score = globalVariables.sequenceVarianceScoresArray[ seqPos ];

				if ( score !== 0 ) {

					$seqDiv.css({ opacity: "" }); //  remove opacity setting


					$seqDiv.addClass( "sequence-color-where-differences-strain-selected" );
				}


			};


			///////////////////////////////////////////////////

			var putStrainsInDOM = function ( strains ) {

				var count;
				var strain;

				var strainsRoot =  $("#StrainsRoot" + SNPIDCounterAppend);

				//  loop through sequences
				for(count = 0; count < strains.length; count++){

					strain = strains[ count ].name;

					strainsRoot.append( "<div id='" + constants.strainIdPrefix + count + SNPIDCounterAppend + "' class='strain strain-instance-" + SNPIDCounterAppend + "  ' " + constants.arrayIndexAttrLabel + "='" + count +  "'>" + strain + "</div>");
				}
			};



			///////////////////////////////////////////////////


			var drawDendrogramLine = function (  dendrogramGraphics, lineStartX, lineStartY, lineEndX, lineEndY ) {

				if ( lineStartX < 0 || lineStartY < 0 ) {

					var z = 0;
				}

				dendrogramGraphics.setColor( constants.dendrogramLineColor );
				dendrogramGraphics.drawLine(lineStartX, lineStartY, lineEndX, lineEndY); // co-ordinates related to surrounding DIVS
			};


			///////////////////////////////////////////////////


			var renderDendrogramPartOfImageRecursive = function ( node, maxLevelsTotal, maxLevelsBelowOfAboveNode, strainHeight, strainHeightHalf, dendrogramGraphics ) {

				var childrenNodes = node.branchset;

				var thisPosition = node.nodePosition;

				var level = node.depthAtNode;

				var maxLevelsBelow = node.maxDepthBelow;

				if ( childrenNodes !== undefined && (  childrenNodes.length !== 0 ) ) {

					for ( var childNodeCounter = 0; childNodeCounter < childrenNodes.length; childNodeCounter++ ) {

						var childNode = childrenNodes[ childNodeCounter ];

						renderDendrogramPartOfImageRecursive( childNode, maxLevelsTotal, maxLevelsBelow, strainHeight, strainHeightHalf, dendrogramGraphics );
					}



					if ( childrenNodes.length !== 1 ) {

						var firstChildPosition = childrenNodes[ 0 ].nodePosition;

						var lastChildIndex = childrenNodes.length - 1;

						var lastChild = childrenNodes[ lastChildIndex ];

						var lastChildPosition =  lastChild.nodePosition;

						if ( lastChildPosition === undefined ) {

							var z = 0;
						}

						// draw horizontal line to connect to the next level up

						var lineStartXV = ( ( maxLevelsTotal - maxLevelsBelowOfAboveNode - 1 ) * constants.widthDendrogramOneLevel );
						var lineStartYV = ( strainHeightHalf + ( strainHeight * thisPosition ) );
						var lineEndXV   = ( ( maxLevelsTotal - maxLevelsBelow - 1 ) * constants.widthDendrogramOneLevel );
						var lineEndYV   = lineStartYV;

						drawDendrogramLine( dendrogramGraphics, lineStartXV, lineStartYV, lineEndXV, lineEndYV );


						// draw vertical  line to connect the levels below

						var lineStartXH = ( ( maxLevelsTotal - maxLevelsBelow - 1 ) * constants.widthDendrogramOneLevel );
						var lineStartYH =  (  strainHeightHalf + ( strainHeight * firstChildPosition ) );

						var lineEndXH   =  lineStartXH;
						var lineEndYH   =  (  strainHeightHalf + ( strainHeight * lastChildPosition ) );

						//  draw the actual line
						drawDendrogramLine( dendrogramGraphics, lineStartXH, lineStartYH, lineEndXH, lineEndYH );

					}



				} else {

					var lineStartXS = ( ( maxLevelsTotal - 1 ) * constants.widthDendrogramOneLevel );
					var lineStartYS = ( strainHeightHalf + ( strainHeight * thisPosition ) );
					var lineEndXS   = ( ( maxLevelsTotal - maxLevelsBelowOfAboveNode - maxLevelsBelow ) * constants.widthDendrogramOneLevel );

					var lineEndYS   = lineStartYS;

					//  draw the actual line
					drawDendrogramLine( dendrogramGraphics, lineStartXS, lineStartYS, lineEndXS, lineEndYS );

				}

			};


			///////////////////////////////////////////////////

			var drawDendrogram = function ( updatedTree ) {

				if ( updatedTree.leafs.length === 0 ) {

					return;
				}

				var strainZeroObj = $("#" + constants.strainIdPrefix + "0" + SNPIDCounterAppend );
				var strainZeroVal = strainZeroObj.first().val();

				if ( strainZeroVal === undefined ) {

					return; //  exit if cannot find first strain text
				}

				//  get height of a single strain block

				var strainHeight = strainZeroObj.height();

				var strainHeightHalf = strainHeight / 2;

				var strainsDendrogram = $("#StrainsDendrogram" + SNPIDCounterAppend);

				//  set height and width;

				var height = updatedTree.leafs.length * ( strainHeight + 1 );

				var width = constants.widthDendrogramOneLevel * ( updatedTree.maxDepth + 1 );

				strainsDendrogram.css( "height", height + "px" );
				strainsDendrogram.css( "width", width + "px" );

				var dendrogramGraphics = new jsGraphics("StrainsDendrogram" + SNPIDCounterAppend);

				var tree = updatedTree.tree;

				var maxLevelsBelowOfAboveNode = tree.maxDepthBelow;

				var maxLevelsTotal = tree.maxDepthBelow + 1;

				renderDendrogramPartOfImageRecursive( tree, maxLevelsTotal, maxLevelsBelowOfAboveNode, strainHeight, strainHeightHalf, dendrogramGraphics );


				dendrogramGraphics.paint();


			};


			//////////////////////////////////////////

			//


			var seqAreaComputeWidthOneChar = function( ) {

				var sequencesRoot =  $("#SequencesRoot" + SNPIDCounterAppend);

				var measurementDivID = "measurementDiv";

				var measurementDivHTML = "<div id='" + measurementDivID + SNPIDCounterAppend + "' class='sequence-group'  style='' >";

				if ( SNPViewer.getRequestTypeDNA() === requestType ) {

					for ( var counter = 0; counter < constants.standardDNACodonLength; counter++ ) {

						measurementDivHTML += "<div class='seq-letter seq-letter-" + requestType + " ' >S</div> ";
					}

					measurementDivHTML += "<div class='bar bar-dna-codon-break '>&nbsp</div>";


				} else if ( SNPViewer.getRequestTypeProtein() === requestType ) {

					measurementDivHTML += "<div class='seq-letter seq-letter-" + requestType + " ' >&nbsp;</div> ";
				}

				measurementDivHTML += "	</div> ";

				sequencesRoot.append( measurementDivHTML );

				var widthOfMeasurementDiv = $("#" + measurementDivID + SNPIDCounterAppend ).outerWidth();


				//  compute the width of one characters

				var widthOfOneChar = widthOfMeasurementDiv;


				if ( SNPViewer.getRequestTypeDNA() === requestType ) {

					widthOfOneChar = widthOfMeasurementDiv / constants.standardDNACodonLength;
				}

				globalVariables.widthOfOneChar = widthOfOneChar;


			}

			/////////////////////////////////////////////////////////////////////

			var setSequenceAreaWidth = function ( ) {

				//  compute the width of one characters

				seqAreaComputeWidthOneChar();

				//  compute the number of characters displayed in the scroll window

				var overallWidth = $rootDivId.width();
				var strainDendrogramWidth = $("#StrainsDendrogramBlock" + SNPIDCounterAppend).outerWidth( true /* [includeMargin] */ );



				var largestPossiblewidthOfSequenceScrollWindow = overallWidth - strainDendrogramWidth - globalVariables.widthOfOneChar - constants.seqBlockPaddingRight;

				globalVariables.numberOfSeqCharsDisplayed = Math.floor( largestPossiblewidthOfSequenceScrollWindow / globalVariables.widthOfOneChar );



				if ( globalVariables.numberOfSeqCharsDisplayed > globalVariables.maxSeqLength ) {

					globalVariables.numberOfSeqCharsDisplayed = globalVariables.maxSeqLength;

				} else {

					// round down to the nearest constants.seqBlockColumnMultiple

					var numberOfSeqCharsDisplayedMultiplier = Math.floor( globalVariables.numberOfSeqCharsDisplayed / constants.seqBlockColumnMultiple );

					globalVariables.numberOfSeqCharsDisplayed = numberOfSeqCharsDisplayedMultiplier * constants.seqBlockColumnMultiple;
				}


				var widthOfSequenceScrollWindow = ( globalVariables.numberOfSeqCharsDisplayed * globalVariables.widthOfOneChar ) + globalVariables.widthOfOneChar;

				var sequenceScrollWindowOuterDOM = $("#sequenceScrollWindowOuter" + SNPIDCounterAppend);

				sequenceScrollWindowOuterDOM.css( "width", widthOfSequenceScrollWindow + "px" );

				//  set last character displayed onto the page.  This may need to be updated if the initial position is not the left edge.
				$("#sequenceEndPos" + SNPIDCounterAppend).text( globalVariables.numberOfSeqCharsDisplayed );


				$("#SequencesShowingBlock" + SNPIDCounterAppend).show();

			};



			/////////////////////////////////////////////////////////////////////

			var determineAndSetSequencesMaxLength = function ( ) {

				var sequences = globalVariables.sequencesInStrainOrder.sequences;

				var currSeqNode;
				var seqNum;

				globalVariables.maxSeqLength = -1;

				//  loop through sequences to get max length
				for(seqNum = 0; seqNum < sequences.length; seqNum++){

					currSeqNode = sequences[ seqNum ];

					if ( currSeqNode !== null ) {

						var sequence = currSeqNode.sequence;

						if ( sequence === undefined ) {

							var z = 0;

							throw "sequence === undefined for currSeqNode = " + currSeqNode ;

						} else if ( sequence === null ) {

							var z = 0;

							throw "sequence === null for currSeqNode = " + currSeqNode ;


						} else {
							if ( sequence.length > globalVariables.maxSeqLength ) {

								globalVariables.maxSeqLength = sequence.length;  // maxSeqLength is a global variable
							}
						}
					}
				}
			};


			///////////////////////////////////////////////////

			//     populate_sequencesAllMatchAtPositionArray

			//  The array sequencesAllMatchAtPositionArray is used for where to put the bars in the upper bar

			var populate_sequencesAllMatchAtPositionArray = function (  ) {

				var sequences = globalVariables.sequencesInStrainOrder.sequences;

				var sequencesArrayLength = sequences.length;

				var sequence;

				var allLettersMatch;

				var positionCounter;
				var sequenceCounter;

				var maxLength = globalVariables.maxSeqLength;

				var firstLetter;
				var currentLetter;

				for ( positionCounter = 0; positionCounter < maxLength; positionCounter++ ) {

					//  process all the letters for a given position

					sequence =  sequences[ 0 ].sequence;

					firstLetter = sequence[ positionCounter ];
					currentLetter = undefined;

					allLettersMatch = true;

					for ( sequenceCounter = 1; sequenceCounter < sequencesArrayLength; sequenceCounter++ ) {

						sequence = sequences[ sequenceCounter ].sequence;

						currentLetter = sequence[ positionCounter ];

						if ( currentLetter !== firstLetter ) {

							allLettersMatch = false;
							break;
						}
					}

					globalVariables.sequencesAllMatchAtPositionArray[ positionCounter ] = allLettersMatch;
				}

			};

			////////////////////////////////////

			//

			var computeScores = function (  ) {

				if ( globalVariables.sequenceVarianceScoresArray ) {

					globalVariables.prevSequenceVarianceScoresArray = globalVariables.sequenceVarianceScoresArray.slice(0); // clone the array
				}

				globalVariables.sequenceVarianceScoresArray = [];


				if ( SNPViewer.getRequestTypeDNA() === requestType ) {


					return computeScoresDNA(  );

				} else if ( SNPViewer.getRequestTypeProtein() === requestType ) {

					return computeScoresProtein(  );

				} else {

					return;
				}

			};




			////////////////////////////////////

			//   computeScoreDNAAllFullCodons

			//   Compute Score for DNA when all passed sequence parts are full codons in length

			var computeScoreDNAAllFullCodons = function ( codonsForPosition, position ) {

				var differenceScore = 0;

				var scores = [];

				var count;

				var firstCodon = codonsForPosition[ 0 ];

				// get associated codons for silent mutation check
				var firstCodonAssocCodonsObj = getDNACodonAssocCodons( firstCodon );

				for ( count = 1; count < codonsForPosition.length; count++ ) {

					var codon = codonsForPosition[ count ];

					if ( firstCodon !== codon ) {

						var silentMutation = false;

						if ( firstCodonAssocCodonsObj !== undefined && firstCodonAssocCodonsObj !== null ) {

							var assocCodonLookup = firstCodonAssocCodonsObj[ codon ];
							if ( assocCodonLookup ) {

								silentMutation = true;
							}
						}

						if ( silentMutation ) {

							if ( differenceScore < constants.dnaSilentMutationOpacity ) {

								differenceScore = constants.dnaSilentMutationOpacity;
							}
						} else {

							if ( differenceScore < 1 ) {

								differenceScore = 1;
							}
						}

						for ( var codonPos = 0; codonPos < firstCodon.length; codonPos++ ) {

							var firstCodonLetter = firstCodon.charAt( codonPos );
							var codonLetter = codon.charAt( codonPos );

							if ( firstCodonLetter !== codonLetter ) {

								if ( ! scores[ codonPos ] ) {

									scores[ codonPos ] = differenceScore;
								} else {
									if ( scores[ codonPos ] && scores[ codonPos ] < differenceScore ) {

										scores[ codonPos ] = differenceScore;
									}
								}
							}
						}
					}
				}

				//  transfer scores into the global array

				for ( count = 0; count < firstCodon.length; count++ ) {

					if ( scores[ count ] ) {
						globalVariables.sequenceVarianceScoresArray[ position + count ] = scores[ count ];
					} else {

						globalVariables.sequenceVarianceScoresArray[ position + count ] = 0;
					}
				}

			};




			////////////////////////////////////

			//   computeScoreDNASomePartialCodons

			//   Compute Score for DNA when some passed sequence parts are partial codons in length

			//   Since some are partials, no checking for silent mutations

			var computeScoreDNASomePartialCodons = function ( codonsForPosition, position ) {

				var scores = [];

				var maxCodonLength = 0;

				var count = 0;
				var letterIndex = 0;

				var codon = null;
				var comparisonLetter = null;

				//  determine max length
				for ( count = 0; count < codonsForPosition.length; count++ ) {

					codon = codonsForPosition[ count ];

					if ( codon.length > maxCodonLength ) {

						maxCodonLength = codon.length;
					}
				}

				for ( letterIndex = 0; letterIndex < maxCodonLength ; letterIndex++  ) {

					comparisonLetter = null;

					for ( count = 0; count < codonsForPosition.length; count++ ) {

						codon = codonsForPosition[ count ];

						var codonLetter = codon.charAt( letterIndex );

						if ( codonLetter ) {
							if ( comparisonLetter === null ) {
								comparisonLetter = codonLetter;

							} else {

								if ( comparisonLetter !== codonLetter ) {

									scores[ letterIndex ] = 1;
								}
							}
						}
					}
				}

				//  transfer scores into the global array

				for ( count = 0; count < maxCodonLength; count++ ) {

					if ( scores[ count ] ) {
						globalVariables.sequenceVarianceScoresArray[ position + count ] = scores[ count ];
					} else {

						globalVariables.sequenceVarianceScoresArray[ position + count ] = 0;
					}
				}

			};


			///////////////////////////////////////////////////

			//     computeScoresDNA

			//  Compute scores for the sequences.

			var computeScoresDNA = function (  ) {


				var sequences = globalVariables.sequencesInStrainOrder.sequences;

				var currSeqNode;
				var seqNum;

				var sequencesCount = sequences.length;

				var shortestCodonLength = null;

				var foundCharacters = true;

				//  loop through positions, constants.standardDNACodonLength at a time
				for ( var position = 0; position < globalVariables.maxSeqLength; position+= constants.standardDNACodonLength ){

					//  optimization that sets the score to zero if all sequences match

					var sequencesAllMatchForAllCodonPositions = true;

					for ( var codonLengthCounter = 0; codonLengthCounter < constants.standardDNACodonLength; codonLengthCounter++ ) {

						if ( ! globalVariables.sequencesAllMatchAtPositionArray[ position + codonLengthCounter ] ) {

							sequencesAllMatchForAllCodonPositions = false;
						}
					}

					if ( sequencesAllMatchForAllCodonPositions ) {

						for ( codonLengthCounter = 0; codonLengthCounter < constants.standardDNACodonLength; codonLengthCounter++ ) {

							globalVariables.sequenceVarianceScoresArray[ position + codonLengthCounter ] = 0;
						}

					} else {

						var codonsForPosition = [];

						foundCharacters = false;  // was text found for these positions for the selected strains?

						//  handle sequences of different lengths

						//  loop through sequences

						for(seqNum = 0; seqNum < sequencesCount; seqNum++){

							var processSequence = true;

							if ( globalVariables.strainsSelectedCount > 1 ) {

								var strainItem = globalVariables.strainList[ seqNum ];

								if ( strainItem.strainPermHighlighting === undefined ||
										strainItem.strainPermHighlighting !== constants.strainPermHighlightingAttrYes ) {

									processSequence = false;
								}
							}

							if ( processSequence ) {

								currSeqNode = sequences[ seqNum ];

								var currSeq = currSeqNode.sequence;

								var singleCodon =  currSeq.substring( position, position + constants.standardDNACodonLength );


								if ( singleCodon === undefined || singleCodon === "" ) {


								} else {

									codonsForPosition.push( singleCodon );

									if ( shortestCodonLength === null ) {

										shortestCodonLength = singleCodon.length;
									} else {

										if ( singleCodon.length < shortestCodonLength  ) {

											shortestCodonLength = singleCodon.length;
										}
									}

									foundCharacters = true;
								}

							}
						}

						if ( foundCharacters ) {

							if ( shortestCodonLength === constants.standardDNACodonLength ) {

								computeScoreDNAAllFullCodons( codonsForPosition, position );

							} else {

								computeScoreDNASomePartialCodons( codonsForPosition, position );
							}

						} else {

						}

					}
				}

			};



			////////////////////////////////////

			//   computeScoreProtein

			//   compute Protein Score for one sequence position

			var computeScoreProtein = function ( charsForPosition ) {

				var score = 0;

				var letter1;
				var letter2;
				var count1, count2;

				var charsLenMinusOne = charsForPosition.length - 1;

				for ( count1 = 0; count1 < charsLenMinusOne; count1++ ) {

					letter1 = charsForPosition[ count1 ];

					for ( count2 = count1 + 1; count2 < charsForPosition.length; count2++ ) {

						letter2 = charsForPosition[ count2 ];

						var singleScore = getProteinScoreForTwoLetters( letter1, letter2 );

						score += singleScore;
					}
				}

				return score;
			};




			///////////////////////////////////////////////////

			//     computeScoresProtein

			//  Compute Protein Scores for the sequences.

			var computeScoresProtein = function (  ) {

				var sequencesAllMatchAtPositionForProcessedStrainsArray = [];

				var scoresBeforeRescalingArray = [];

				var maxScore = null;
				var minScore = null;

				var sequences = globalVariables.sequencesInStrainOrder.sequences;

				var currSeqNode;
				var seqNum;

				var sequencesCount = sequences.length;

				var position = 0;

				var foundSeqLetter = true;

				//  loop through positions
				for ( position = 0; position < globalVariables.maxSeqLength; position++ ){

					//  optimization that sets the score to null if all sequences match

					if ( globalVariables.sequencesAllMatchAtPositionArray[ position ] ) {

						sequencesAllMatchAtPositionForProcessedStrainsArray[ position ] = true;

						scoresBeforeRescalingArray[ position ] = null;

					} else {

						var charsForPosition = [];

						foundSeqLetter = false;  // handle sequences of different lengths

						var comparisonLetter = null;

						var allLettersAtPositionMatch = true;

						//  loop through sequences

						for(seqNum = 0; seqNum < sequencesCount; seqNum++){

							var processSequence = true;

							if ( globalVariables.strainsSelectedCount > 1 ) {

								var strainItem = globalVariables.strainList[ seqNum ];

								if ( strainItem.strainPermHighlighting === undefined ||
										strainItem.strainPermHighlighting !== constants.strainPermHighlightingAttrYes ) {

									processSequence = false;
								}
							}

							if ( processSequence ) {

								currSeqNode = sequences[ seqNum ];

								var currSeq = currSeqNode.sequence;

								var seqSingle =  currSeq.charAt( position );

								if ( seqSingle === undefined || seqSingle === "" ) {

								} else {

									charsForPosition.push( seqSingle );

									foundSeqLetter = true;

									if ( comparisonLetter === null ) {

										comparisonLetter = seqSingle;  //  store off first letter

									} else {

										if ( comparisonLetter !== seqSingle ) {

											allLettersAtPositionMatch = false;
										}

									}
								}
							}
						}

						if ( foundSeqLetter ) {

							if ( allLettersAtPositionMatch ) {

								sequencesAllMatchAtPositionForProcessedStrainsArray[ position ] = true;

								scoresBeforeRescalingArray[ position ] = null;

							} else {

								sequencesAllMatchAtPositionForProcessedStrainsArray[ position ] = false;

								var score = computeScoreProtein( charsForPosition );

								scoresBeforeRescalingArray[ position ] = score;

								if ( maxScore === null ) {

									maxScore = score;
									minScore = score;

								} else if ( score > maxScore ) {

									maxScore = score;

								} else if ( score < minScore ) {

									minScore = score;
								}
							}

						} else {

						}
					}
				}

				//  Scale the scores based on the max and min

				var maxScoreMinScoreDifference = maxScore - minScore;

				var minusMinOpacityLocal = constants.minOpacity;

				var oneMinusMinOpacity = 1 - constants.minOpacity;

				//  loop through positions
				for ( position = 0; position < globalVariables.maxSeqLength; position++ ){

					//  optimization that sets the score to zero if all sequences match

					if ( sequencesAllMatchAtPositionForProcessedStrainsArray[ position ] ) {

						globalVariables.sequenceVarianceScoresArray[ position ] = 0;

					} else {

						if ( maxScoreMinScoreDifference === 0 ) {

							// Max and Min are equal so every score is the same so give all a "1"

							globalVariables.sequenceVarianceScoresArray[ position ] = 1;

						} else {

							//  scale score to zero to one
							var scaledScore = ( scoresBeforeRescalingArray[ position ] - minScore ) / maxScoreMinScoreDifference;

							//  scale scaledScore to constants.minOpacity to one
							var scaledScoreToMinimum = ( ( oneMinusMinOpacity ) * scaledScore ) + minusMinOpacityLocal;

							//  round to nearest hundredth
							var roundedScore = ( Math.round( scaledScoreToMinimum * 100 ) ) / 100;


							globalVariables.sequenceVarianceScoresArray[ position ] = roundedScore;
						}
					}

				}

			};


			///////////////////////////////////////////////////

			// updateTreeWithPositionEtc

			///   Create updated tree representing parsed Newick data
			//     along with additional information needed for the rendering
			//     of the dendrogram.

			var updateTreeWithPositionEtc = function ( tree ) {

				var updatedNewickData = {};

				var leafs = [];

				globalVariables.strainList = [];

				updatedNewickData.tree = tree;

				updatedNewickData.leafs = leafs;

				var maxDepth = 0;

				var leafPosition = -1;

				//  internally defined function/closure

				var updateTreeWithPositionEtcInternal = function ( node, depthAtNodeParam ) {

					node.depthAtNode = depthAtNodeParam;

					maxDepth = Math.max( maxDepth, depthAtNodeParam );

					var maxDepthBelow = 0;

					var depthAtNodeChildren = depthAtNodeParam + 1;

					var childrenNodes = node.branchset;

					if ( childrenNodes !== undefined && (  childrenNodes.length !== 0 ) ) {

						for ( var childNodeCounter = 0; childNodeCounter < childrenNodes.length; childNodeCounter++ ) {

							var childNode = childrenNodes[ childNodeCounter ];

							updateTreeWithPositionEtcInternal( childNode, depthAtNodeChildren );

							maxDepthBelow = Math.max( maxDepthBelow, childNode.maxDepthBelow );
						}

						node.maxDepthBelow = maxDepthBelow + 1;

						if ( childrenNodes.length === 1 ) {

							node.nodePosition =  childrenNodes[0].nodePosition;

						} else {

							var firstChildPosition = childrenNodes[ 0 ].nodePosition;

							var lastChildIndex = childrenNodes.length - 1;

							var lastChild = childrenNodes[ lastChildIndex ];

							var lastChildPosition =  lastChild.nodePosition;

							var thisPosition = ( firstChildPosition + lastChildPosition ) / 2;

							node.nodePosition = thisPosition;
						}

					} else {

						leafPosition++;

						node.nodePosition = leafPosition;

						node.maxDepthBelow = 1;

						leafs.push( node );

						var strain = new Strain( node.name );

						globalVariables.strainList.push( strain );

					}
				};

				updateTreeWithPositionEtcInternal( tree, 1 /* depthAtNodeParam  */ );

				updatedNewickData.leafCount = leafPosition + 1;

				updatedNewickData.maxDepth = maxDepth;

				return updatedNewickData;
			};


			////////////////////////////////////////

			//   getSequencesInStrainOrder

			//   put sequences from data structure into an array in the same order as the strains


			var getSequencesInStrainOrder = function ( leafs /* array of node */, sequencesInObject ) {

				var sequencesInStrainOrderLocal = {};

				var sequences = [];

				var count;

				sequencesInStrainOrderLocal.sequences = sequences;

				//  loop through sequences
//				for(count = 0; count < leafs.length; count++){
//
//					var strainName = leafs[ count ].name;

				//  loop through sequences
				for(count = 0; count < globalVariables.strainList.length; count++){

					var strainName = globalVariables.strainList[ count ].name;


					var sequenceInObject = sequencesInObject[ strainName ];

					if ( sequenceInObject === undefined || sequenceInObject === null ) {

						throw "ERROR: Strain name in newick file or provided strain List not provided in object containing sequences, strain name = '" + strainName + "'";
					}

					var sequenceObj = {};

					sequenceObj.sequence = sequenceInObject;
					sequenceObj.strainName = strainName;

					sequences.push( sequenceObj );
				}

				return sequencesInStrainOrderLocal;
			};




			////////////////////////////////////////


			//  The following code is not in a function.  It runs immediately.


			var mainHTML =

				'	<div class="' + getSNPObjectStorageClassName() + '" style="clear: both;">' +
				'		<div id="starting-sequence-position-1' + SNPIDCounterAppend + '" class="starting-sequence-position">' +
				'			&nbsp;' +
				'		</div>' +
				'		<div id="ending-sequence-position-1' + SNPIDCounterAppend + '" class="ending-sequence-position" >' +
				'			&nbsp;' +
				'		</div>' +
				'	</div>' +
				'	<div style="clear: both;"></div>' +
				'	<div id="sequenceDraggableBoundingBox-1' + SNPIDCounterAppend + '" class="sequenceDraggableBoundingBox">' +

				'		<div id="upperBarChartClickable-1' + SNPIDCounterAppend + '" class="upperBarChartClickable upperBarChartClickableBorder" >' +

				'			<img id="upperBarChartClickableImage-1' + SNPIDCounterAppend + '" class="upperBarChartClickable"' +
				'				src="' + constants.transparentPixelPNG + '" />' +
				'		</div>' +

				'		<div id="upperBarChart-1' + SNPIDCounterAppend + '" class="upperBarChart" ></div>' +



				'		<div id="sequenceDraggable-1' + SNPIDCounterAppend + '" class="upperBarChartSelector" >' +

				'			<img id="sequenceDraggableImage-1' + SNPIDCounterAppend + '" class="upperBarChartDraggableImage"' +
				'				src="' + constants.transparentPixelPNG + '" />' +

				'			<div id="curPosUpprBar-1' + SNPIDCounterAppend + '" class="currPosNumOnUpprBar" >' +
				'				&nbsp;' +
				'			</div>' +

				'		</div>' +
				'	</div>' +



				//	<%--  Start second bar --%>' +

				'	<div id="BarChart-2' + SNPIDCounterAppend + '" style="display: none;">' +

				'		<br />' +
				'		<br />' +

				'		<div style="clear: both;"></div>' +




				'		<div >' +
				'			<div id="starting-sequence-position-2' + SNPIDCounterAppend + '" class="starting-sequence-position">' +
				'				&nbsp;' +
				'			</div>' +
				'			<div id="ending-sequence-position-2' + SNPIDCounterAppend + '" class="ending-sequence-position" >' +
				'				&nbsp;' +
				'			</div>' +
				'		</div>' +
				'		<div style="clear: both;"></div>' +

				//		<%--  width of 'outerSequenceDraggableBoundingBox-2' is updated in the code to match width of 'sequenceDraggableBoundingBox-2' --%>' +
				//		<%--  padding added to allow the display of the top and bottom of the selector box and the position number below the selector box --%>' +


				'		<div id="sequenceDraggableBoundingBox-2' + SNPIDCounterAppend + '" ' +
				'			class="sequenceDraggableBoundingBox' + '">' +

				'			<div id="outerSequenceDraggableBoundingBox-2' + SNPIDCounterAppend + '"' +
				'				class="sequenceDraggableBoundingBox"' +
				'				style="overflow: hidden;">' +

				'				<div id="upperBarChart-slidable-2' + SNPIDCounterAppend + '"' +
				'					style="position: relative; width: 100%">' +

				'					<div id="upperBarChart-2' + SNPIDCounterAppend + '" class="upperBarChart" ></div>' +




				'				</div>' +

				'			</div>' +
				'			<div id="sequenceDraggable-2' + SNPIDCounterAppend + '" class="upperBarChartSelector" >' +

				'				<img id="sequenceDraggableImage-2' + SNPIDCounterAppend + '" class="upperBarChartDraggableImage"' +
				'					src="' + constants.transparentPixelPNG + '" />' +

				'				<div id="curPosUpprBar-2' + SNPIDCounterAppend + '"  class="currPosNumOnUpprBar" >' +
				'					&nbsp;' +
				'				</div>' +

				'			</div>' +



				'			<div id="upperBarChartClickable-2' + SNPIDCounterAppend + '"' +
				'				class="upperBarChartClickable upperBarChartClickableBorder" >' +

				'				<img id="upperBarChartClickableImage-2' + SNPIDCounterAppend + '" class="upperBarChartClickable"' +
				'					src="' + constants.transparentPixelPNG + '" />' +
				'			</div>' +

				'		</div>' +

				'	</div>' +
				'	<br />' +
//				'	<br />' +

				'	<div style="clear: both;"></div>' +


				//	 "Showing" start and end positions data block


				'		<div id="barChartStrainsDendrogramBlock' + SNPIDCounterAppend + '" class="showing-sequence-pos-block" >' +
				'			<div id="SequencesShowingBlock' + SNPIDCounterAppend + '" >' +
				'				Showing' +
				'				<span id="sequenceRequestTypeLabel' + SNPIDCounterAppend + '"></span>' +
				'				<span id="sequenceStartPos' + SNPIDCounterAppend + '">1</span> -' +
				'				<span id="sequenceEndPos' + SNPIDCounterAppend + '">0</span>' +
				'			</div>' +
				'		</div>' +


				//	 Start Dendrogram and sequence data block

				'	<div id="StrainsDendrogramBlock' + SNPIDCounterAppend + '" class="StrainsDendrogramBlock">' +

				'		<div  class="barchart" >&nbsp;</div>' +





				'		<div id="DendrogramBlock' + SNPIDCounterAppend + '" class="DendrogramBlock">' +
				'			<div id="StrainsDendrogram' + SNPIDCounterAppend + '" class="StrainsDendrogram">' +


				'			</div>' +

				'		</div>' +


				'		<div id="StrainsBlock' + SNPIDCounterAppend + '" class="strain-block ">' +

				'			<div id="StrainsRoot' + SNPIDCounterAppend + '" class=" strain-style">' +
				'			</div>' +


				'		</div>' +
				'	</div>' +




				'	<div id="sequenceScrollWindowOuter' + SNPIDCounterAppend + '"' +
				'		class="sequence-block sequence-style"' +
				'		style="overflow: hidden; position: relative;">' +

				'			<div id="barChart' + SNPIDCounterAppend + '" class="barchart" style="">' +

				'			</div>' +

				'			<div id="SequencesRoot' + SNPIDCounterAppend + '">' +


				'			</div>' +


				'	</div>' +
				'	<div style="clear: both;">' ;



			$("#" + rootDivId  ).append( mainHTML );



			$("#sequenceRequestTypeLabel" + SNPIDCounterAppend).text( requestTypeLabel );

			var updatedTree = undefined;

			if ( newickString !== undefined ) {

				var treeParsedNewick = parseNewickString( newickString );

				var updatedTree = updateTreeWithPositionEtc( treeParsedNewick ); //  updates globalVariables.strainList

			} else {

				// populate the "globalVariables.strainList" from sequenceLabelsArray

				globalVariables.strainList = [];

				for ( var sequenceLabelsArrayCounter = 0; sequenceLabelsArrayCounter < sequenceLabelsArray.length; sequenceLabelsArrayCounter++ ) {

					var strainEntry = new Strain( sequenceLabelsArray[ sequenceLabelsArrayCounter ] );
					globalVariables.strainList.push( strainEntry );
				}
			}

			putStrainsInDOM( globalVariables.strainList );


			if ( newickString !== undefined ) {
				drawDendrogram( updatedTree );
			}

			try {

				// setting global variables

				globalVariables.sequencesInStrainOrder = getSequencesInStrainOrder( null /* updatedTree.leafs */, sequencesObject );


				determineAndSetSequencesMaxLength(  );

				populate_sequencesAllMatchAtPositionArray();


				computeScores( );


				setSequenceAreaWidth();


				globalVariables.sequenceDisplayBlock = new SequenceDisplayBlock();



				globalVariables.primaryBarChartData  =

					new BarChartData( { label: 1 ,
						idAppended: "-1" ,
						isPrimary: true ,
						startingSequencePosition: 0 ,  // zero based
						endingSequencePosition: globalVariables.maxSeqLength - 1 ,  // zero based
						childToUpdateOnChange: globalVariables.sequenceDisplayBlock } );


				globalVariables.sequenceDisplayBlock.setParentToUpdateOnChange( globalVariables.primaryBarChartData );


				globalVariables.primaryBarChartData.setGlobalVariablesFromExistingSizes();

				globalVariables.primaryBarChartData.updateDraggableBoundingBoxForMaxLength( );


				globalVariables.primaryBarChartData.updateDraggableBox( );




			} catch ( tt ) {

				var stackTrace = tt.stack;

				throw tt;
			}

			//  seperate the following parts into a seperate "run" to avoid the timeout problem

			//	setTimeout( function() {


			try {

				// perform initial population of sequences
				globalVariables.sequenceDisplayBlock.populateSequenceBlock(  );

				globalVariables.sequenceDisplayBlock.addNavArrowsBelowSequenceBlock();

				//   More functions related to user manipulaton of the page


				globalVariables.primaryBarChartData.attachEventHandlers();



				//  seperate the following parts into a seperate "run" to avoid the timeout problem

				//			setTimeout( function() {


				try {


					try {

						//   set up second bar chart if created
						if ( globalVariables.secondaryBarChartData ) {

							globalVariables.secondaryBarChartData.setGlobalVariablesFromExistingSizes();

							globalVariables.secondaryBarChartData.updateDraggableBoundingBoxForMaxLength( );


							globalVariables.secondaryBarChartData.updateDraggableBox( );

							globalVariables.secondaryBarChartData.attachEventHandlers();

						}

					} catch ( t ) {

						var stackTrace = t.stack;

						throw t;
					}


					////////////////////////////////////

					//   deselectStrain


					var deselectStrain = function( strainDiv ) {


						var strainArrayIndex = strainDiv.attr( constants.arrayIndexAttrLabel  );

						var strainItem = globalVariables.strainList[ strainArrayIndex ];

						strainItem.strainHoverHighlighting = constants.strainHoverHighlightingAttrNo;

						if ( strainItem.strainPermHighlighting === undefined ||
								strainItem.strainPermHighlighting !== constants.strainPermHighlightingAttrYes ) {

							var $sequenceDiv = $("#" + constants.sequenceIdPrefix + strainArrayIndex + SNPIDCounterAppend );

							$sequenceDiv.children( ".seq-letter" ).each(function(){

								try {
									var $seqDiv = $(this);


									var seqDivLetter = $seqDiv.html();


									$seqDiv.removeClass( "seq-letter-strain-selected-" + requestType + "-" + seqDivLetter );


									var seqPos = $seqDiv.attr( constants.sequencePositionAttrLabel );

									var score = globalVariables.sequenceVarianceScoresArray[ seqPos ];

									if ( score !== 0 ) {

										$seqDiv.css( "opacity", score );


										$seqDiv.removeClass( "sequence-color-where-differences-strain-selected" );

										//									} else {
										//
										//										seqDiv.css( "opacity", "" );

									}
								} catch ( t ) {

										var stackTrace = t.stack;

										var z = 0;

										throw t;
								}
							});


							$sequenceDiv.children( ".seq-letter-dna-codon-break" ).each(function(){

								try {
									var $seqDiv = $(this);

									$seqDiv.removeClass("seq-letter-dna-codon-break-strain-selected");

								} catch ( t ) {

									var stackTrace = t.stack;

									throw t;
								}

							});



						} else {


						}

					};



					//   attach handlers to the strain labels

					$( "div.strain-instance-" + SNPIDCounterAppend ).click( function( event, ui ) {

						var strainDiv = $(this);

						var strainArrayIndex = strainDiv.attr( constants.arrayIndexAttrLabel  );

						var strainItem = globalVariables.strainList[ strainArrayIndex ];

						if ( strainItem.strainPermHighlighting && strainItem.strainPermHighlighting === constants.strainPermHighlightingAttrYes ) {

							strainItem.strainPermHighlighting = constants.strainPermHighlightingAttrNo;

							globalVariables.strainsSelectedCount--;

							if ( globalVariables.strainsSelectedCount < 0 ) {

								globalVariables.strainsSelectedCount = 0;
							}

							strainDiv.removeClass("strain-selected");

						} else {

							strainItem.strainPermHighlighting = constants.strainPermHighlightingAttrYes;

							globalVariables.strainsSelectedCount++;

							strainDiv.addClass("strain-selected");

						}

						computeScores( );

						globalVariables.primaryBarChartData.getRatiosForUpperBarIndividBars( { computeMaxOnly: true } );

						globalVariables.primaryBarChartData.addBarsToUpperBar();

						if ( globalVariables.secondaryBarChartData ) {

							globalVariables.secondaryBarChartData.getRatiosForUpperBarIndividBars( { computeMaxOnly: true } );

							globalVariables.secondaryBarChartData.addBarsToUpperBar();
						}

						// perform update of sequences
						globalVariables.sequenceDisplayBlock.updateSequenceBlockForStrainChange(  );


					}
					);




					$("div.strain-instance-" + SNPIDCounterAppend).hover(    //       hover( enterFunction, leaveFunction )

						function(){  //  enterFunction

							var strainDiv = $(this);

							var strainArrayIndex = strainDiv.attr( constants.arrayIndexAttrLabel  );

							var strainItem = globalVariables.strainList[ strainArrayIndex ];

							strainItem.strainHoverHighlighting = constants.strainHoverHighlightingAttrYes;

							if ( strainItem.strainPermHighlighting === undefined ||
									strainItem.strainPermHighlighting !== constants.strainPermHighlightingAttrYes ) {

								highlightSequenceLetters( strainArrayIndex );

							} else {


							}


						}, function(){  //  leaveFunction

							var strainDiv = $(this);

							deselectStrain( strainDiv );
						});



					//  show the root div when it is done being developed

					$rootDivId.css("visibility", "visible");

					callbackFunctionsObj.successfullyLoadedCallback();

				} catch ( t ) {


					var stackTrace = t.stack;

					//  Only "throw t" if NOT inside a timer
					throw t;

					//	callbackFunctionsObj
					//  { successfullyLoadedCallback: successfullyLoadedCallback, failToLoadCallBack: failToLoadCallBack }

					//  Only "callbackFunctionsObj.failToLoadCallBack();" if  inside a timer
					//  callbackFunctionsObj.failToLoadCallBack();

				}


				//			} ,50);



			} catch ( t ) {


				var stackTrace = t.stack;

				//	callbackFunctionsObj
				//  { successfullyLoadedCallback: successfullyLoadedCallback, failToLoadCallBack: failToLoadCallBack }

				//  Only "throw t" if NOT inside a timer
				throw t;

				//  Only "callbackFunctionsObj.failToLoadCallBack();" if  inside a timer
				//			callbackFunctionsObj.failToLoadCallBack();

			}


			//	} ,50);


		};

		////////////////////////

		//  Code after the constructor for SNPViewerInternal

		var internalVariableSNPViewer = new SNPViewerInternal ( rootDivIdCreateRoot, requestParamsCreateRoot, configParamsCreateRoot, callbackFunctionsObjCreateRoot );

		var SNPObjectStorageNodeSelector = "#" + rootDivIdCreateRoot + " ." + getSNPObjectStorageClassName();

		var $SNPObjectStorageNode = $( SNPObjectStorageNodeSelector );

		//  attach the newly created object to the DOM object that is a child of the id provided

		if ( $SNPObjectStorageNode.size() === 0 ) {

			throw "Unable to find SNPObjectStorageNode under rootDivID, rootDivID = " + rootDivId + ", SNPObjectStorageNodeSelector = " + SNPObjectStorageNodeSelector;
		}


		$SNPObjectStorageNode.data( SNPViewer.getDataElement(), internalVariableSNPViewer );


		return internalVariableSNPViewer;

	};

};



var SNPViewer = new SNPViewer();  // Declare SNPViewer namespace

SNPViewer.getRequestTypeDNA = function() {
	return "dna";
};


SNPViewer.getRequestTypeProtein = function() {
	return "protein";
};




