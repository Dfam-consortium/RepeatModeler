/////////////////////////////////////////////////////////////////////////
/*
 *
 *  Alignment Viewer Class
 *    - Robert Hubley 11/2013
 *
 */
////////////////////////////////////////////////////////////////////////
function AlignmentViewer(canvas, json, options)
{
  this.json = json;
  this.canvas = canvas;
  this.context = this.canvas.getContext("2d");

  // Defaults
  this.context.font = "16px Courier New, monospace";
  this.viewType = "norm";
  this.highlightCpGs = 1;  // TODO: Hook this up
  this.noteTrans = 1;      // TODO: Hook this up
  this.dotEquivBases = 1;  // TODO: Hook this up
  this.showInsert = 1;     // TODO: Hook this up
  this.lineSpacing = 10;    // px
  this.rulerHeight = 25;    // px
  this.scoreGraphHeight = 20; //px
  this.topBottomMargin = 5; // px
  this.rulerTickDistance = 25; // bases

  // Derived values
  this.fontWidth = this.context.measureText("A").width;
  // NOTE: There is a way to get the font height using
  //       a technique of rendering and measuring.  Unless
  //       we find the need we are going to assume it's
  //       the same as the width.
  this.fontHeight = this.fontWidth;


  // Find max score
  this.maxScore = 0;
  for (i = 0; i < this.json.alignmentScore.length; i++) {
    if ( this.json.alignmentScore[i] > this.maxScore )
      this.maxScore = this.json.alignmentScore[i];
  }

  // Find reference sequence in the json object
  // and the deepest section of the alignment.
  this.referenceSeqIdx = -1;   // Index of sequence with id = "reference" with json.alignment[]
  this.maxIDLen = 0;  // The string length of the longest ID in the json.alignment[] collection
  var coverageArray = []; // A temporary array used in calcuation of max alignment depth
  for (i = 0; i < this.json.alignment.length; i++) {
    if (this.json.alignment[i].id == "reference") {
      this.referenceSeqIdx = i;
    }else
    {
      if ( this.maxIDLen < this.json.alignment[i].id.length )
        this.maxIDLen = this.json.alignment[i].id.length;
    }
    for ( j = 0; j < this.json.alignment[i].sequence.length; j++ ) {
      if ( coverageArray[j+this.json.alignment[i].start] == undefined )
      {
        coverageArray[j + this.json.alignment[i].start ] = 1;
      }else
      {
        coverageArray[j + this.json.alignment[i].start ]++;
      }
    }
  }
  this.maxDepth = 0;
  for ( i = 0; i < coverageArray.length; i++ ) {
    if ( coverageArray[i] > this.maxDepth ) 
      this.maxDepth = coverageArray[i];
  }
  this.referenceSeq = this.json.alignment[this.referenceSeqIdx].sequence;
  this.rulerPosIdx = []; // A sparse array with string position to base position translation for ruler drawing
  var pos = 0;
  for ( i = 0; i < this.referenceSeq.length; i++ ) {
    if ( this.referenceSeq.charAt(i) != "-" ) 
    {
      pos++;
      if ( pos == 1 || (pos % this.rulerTickDistance) == 0 )
        this.rulerPosIdx[i] = pos;
    }
  }

  // Alignment viewport size ( bp and lines )
  this.viewCols = (this.canvas.width - ( this.maxIDLen * this.fontWidth )) / this.fontWidth;
  this.viewLines = (this.canvas.height - this.topBottomMargin - this.rulerHeight -
                     this.lineSpacing - this.fontHeight) / ( this.fontWidth + this.lineSpacing );

  // create Scroller instance
  var that = this;
  this.scroller = new Scroller(function(left, top, zoom) {
          that.render(left, top, zoom);
     }, options);

  // bind events
  this.bindEvents();

  // reflow for the first time
  this.reflow();
}


AlignmentViewer.prototype.setViewType = function(type) {
  this.viewType = type;
  var values = this.scroller.getValues(); // { left, top, zoom }
  this.render( values.left, values.top, values.zoom );
}; 


AlignmentViewer.prototype.reflow = function() {
  // NOTE: Not sure why I need to pad out the Y scale below
  this.scroller.setDimensions(this.viewCols, this.viewLines, this.referenceSeq.length, this.maxDepth + 10 );
};


AlignmentViewer.prototype.render = function(left, top, zoom) {
    // DEBUGING
    //var foo = document.getElementById("status");
    //foo.innerHTML="viewLines = " + this.viewLines;

    // Clear the canvas
    this.context.clearRect(0, 0, this.canvas.width, this.canvas.height);

    // Draw the bounding box
    this.context.beginPath();
    this.context.moveTo(0, 0);
    this.context.lineTo(this.canvas.width, 0);
    this.context.lineTo(this.canvas.width, this.canvas.height);
    this.context.lineTo(0, this.canvas.height);
    this.context.lineTo(0, 0);
    this.context.stroke();
 
    // Save space for IDs
    var curX = this.maxIDLen * this.fontWidth;
    var curY = this.topBottomMargin;

    // Pull out the reference sub-sequence
    var conSeq = "";
    if ( left < 0 )
    {
       conSeq = this.referenceSeq.substr(0, this.viewCols);
    }else 
    {
       conSeq = this.referenceSeq.substr(left, this.viewCols);
    }

    // Allow stretch animation past beginning
    if ( left < 0 )
    {
      curX = curX + (-left * this.fontWidth );
    }

    // Draw the ruler
    this.context.fillStyle = 'black';
    var tmpX = curX;
    for ( i = 0; i < conSeq.length; i++ ) 
    {
      var idx = (left>0)?Math.round(left)+i:i;
      if ( this.rulerPosIdx[idx] != undefined )
      {
        this.context.beginPath();
        this.context.moveTo(tmpX + (this.fontWidth / 2), curY + this.fontHeight + 2);
        this.context.lineTo(tmpX + (this.fontWidth / 2), curY + this.rulerHeight);
        this.context.stroke();
        this.context.fillText(this.rulerPosIdx[idx],
              tmpX + (this.fontWidth / 2) - (( this.rulerPosIdx[idx].toString().length / 2 ) * this.fontWidth), curY + 10);
      }
      tmpX = tmpX + this.fontWidth;
    }

    curY = curY + this.rulerHeight + this.lineSpacing

    // Draw score
    // TODO: Some cleanup work
    //       - refine spacing
    //       - Draw bounding box around sub-graph
    //       - Clean up redundant computation
    tmpX = curX;
    var startBlock = -1;
    for ( i = 0; i < conSeq.length; i++ ) 
    {
      var idx = (left>0)?Math.round(left)+i:i;
      this.context.fillStyle = 'red';
      this.context.fillRect(tmpX, curY + this.scoreGraphHeight - ((this.json.alignmentScore[idx] / this.maxScore) * this.scoreGraphHeight), 
                            this.fontWidth, ((this.json.alignmentScore[idx] / this.maxScore) * this.scoreGraphHeight));
      if ( idx > 0 && this.json.alignmentScore[idx-1] != this.json.alignmentScore[idx] )
      {
        if ( startBlock > -1 )
        {
          // write it out
          var halfIdx = Math.round((idx-startBlock)/2);
          this.context.font = "10px Courier New, monospace";
          var thisX = tmpX - ( halfIdx * this.fontWidth ) - ( (this.json.alignmentScore[idx-1].toString().length / 2 ) * this.fontWidth);
          this.context.fillText(this.json.alignmentScore[idx-1], thisX, curY );
          this.context.font = "16px Courier New, monospace";
        }
        if ( this.json.alignmentScore[idx] >= 1 )
        {
          startBlock = idx;
        }else
        {
          startBlock = -1;
        }
      }
      tmpX = tmpX + this.fontWidth;
    }

    curY = curY + this.scoreGraphHeight + this.lineSpacing + this.fontHeight;
     
    // Draw the reference 
    //    - TODO: Draw reference and/or consensus
    this.context.fillStyle = 'blue';
    this.context.fillText("Reference",0, curY);
    this.context.fillText(conSeq, curX, curY);
    curY = curY + ( this.fontHeight + this.lineSpacing );

    // Allow stretch animation when alignment is pulled past top
    if ( top < 0 )
    {
      curY = curY + (-top * ( this.fontWidth + 10 ) );
    }

    // Lookup table for transitions
    var transitions = { "CT": 1, "TC": 1, "AG": 1, "GA": 1 };

    // 
    // Display the alignments
    //
    var visibleLineCounter = 0;
    var inCpG = 0;
    for (i = 0; i < this.json.alignment.length; i++) {
       if ( i == this.referenceSeqIdx)
         continue;

        var align = this.json.alignment[i];
        var aSeq = align.sequence;
        var aStart = align.start;
        var aEnd = align.start + aSeq.length;

        right = left + this.viewCols;
        if ( this.referenceSeq.length < right )
           right =  this.referenceSeq.length;
 
        bottom = top + this.viewLines;
        
        // Is this alignment within our current sequence range?
        if ( aStart > right || aEnd < left ) continue;
        visibleLineCounter++;

        // Is this sequence beyond our vertical view boundary?
        if ( visibleLineCounter > bottom ) break;

        // Is this sequence inside our vertical biew boundary?
        if ( visibleLineCounter >= top )
        {
            inCpG = 0;
            this.context.fillStyle = 'black';
            var startIndx = left;
            if ( aStart > startIndx ) startIndx = aStart;
            var tmpX = curX + ((startIndx - ((left<0)?0:left)) * this.fontWidth);
            this.context.fillText(align.id, 0, curY);

            if ( this.viewType == "norm" ) {
              // Optimisation: draw full alignment data in one fillText 
              //               if in normal display mode.
              var start = 0;
              if ( aStart < left ) start = left - aStart;
              var subLen = ((left + this.viewCols ) - start + 1);
              if ( subLen > 0 )
              {
                var text = aSeq.substr( start, subLen );
                this.context.fillText(text, tmpX, curY);
              }
            }else{
              for (j = startIndx; j < right; j++) {
                var cBase = this.referenceSeq.charAt(j);
                // NOTE: This can go past the end of the referenceSeq 
                //       but unlike other languages js returns null if
                //       it exceeds the length.  No need to waste time
                //       checking for safety.
                var cBasePlus = this.referenceSeq.charAt((j + 1));
                aBase = aSeq.charAt(j - aStart);
                if (cBase == "-") {
                    if (aBase == "-") {
                        // Draw blank space
                        aBase = " ";
                    } else {
                        // Draw character
                    }
                } else if (aBase == cBase) {
                    aBase = ".";

                } else {
                    if ( transitions[ aBase + cBase ] != undefined )
                       aBase = "i";
                }

                // Check for cons CG coloring
                if (cBase == "C" && cBasePlus == "G") {
                    // Change background color
                    this.context.fillStyle = 'red';
                    inCpG = 1;
                }

                this.context.fillText(aBase, tmpX, curY);
                tmpX = tmpX + this.fontWidth;

                if (inCpG && cBase == "G") {
                    // change back the color
                    this.context.fillStyle = 'black';
                    inCpG = 0;
                }
              } // for j over reference length
           } // else -- diff display mode
           curY = curY + this.fontHeight + this.lineSpacing;
        } // if visibleLineCounter >= top
    } // for i over alignments
};


AlignmentViewer.prototype.bindEvents = function() {

	var that = this;

	// reflow handling
	window.addEventListener("resize", function() {
		that.reflow();
	}, false);

	// touch devices bind touch events
	if ('ontouchstart' in window) {

		this.canvas.addEventListener("touchstart", function(e) {

			// Don't react if initial down happens on a form element
			if (e.touches[0] && e.touches[0].target && e.touches[0].target.tagName.match(/input|textarea|select/i)) {
				return;
			}

			that.scroller.doTouchStart(e.touches, e.timeStamp);
			e.preventDefault();

		}, false);

		document.addEventListener("touchmove", function(e) {
			that.scroller.doTouchMove(e.touches, e.timeStamp, e.scale);
		}, false);

		document.addEventListener("touchend", function(e) {
			that.scroller.doTouchEnd(e.timeStamp);
		}, false);

		document.addEventListener("touchcancel", function(e) {
			that.scroller.doTouchEnd(e.timeStamp);
		}, false);

	// non-touch bind mouse events
	} else {
		
		var mousedown = false;

		this.canvas.addEventListener("mousedown", function(e) {

			if (e.target.tagName.match(/input|textarea|select/i)) {
				return;
			}
		
			that.scroller.doTouchStart([{
				pageX: e.pageX,
				pageY: e.pageY
			}], e.timeStamp);

			mousedown = true;
			e.preventDefault();

		}, false);

		document.addEventListener("mousemove", function(e) {

			if (!mousedown) {
				return;
			}
			
			that.scroller.doTouchMove([{
				pageX: e.pageX,
				pageY: e.pageY
			}], e.timeStamp);

			mousedown = true;

		}, false);

		document.addEventListener("mouseup", function(e) {

			if (!mousedown) {
				return;
			}
			
			that.scroller.doTouchEnd(e.timeStamp);

			mousedown = false;

		}, false);

		this.canvas.addEventListener("mousewheel", function(e) {
			if(that.options.zooming) {
				that.scroller.doMouseZoom(e.wheelDelta, e.timeStamp, e.pageX, e.pageY);	
			}
		}, false);

	}
};
