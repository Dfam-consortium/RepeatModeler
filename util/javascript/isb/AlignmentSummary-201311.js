/////////////////////////////////////////////////////////////////////////
/*
 *
 *  Alignment Summary
 *    - Robert Hubley 11/2013
 *
 */
////////////////////////////////////////////////////////////////////////
function AlignmentSummary(align_canvas, guide_canvas, json, options)
{
  this.json = json;
  this.align_canvas = align_canvas;
  this.guide_canvas = guide_canvas;

  // Adjust height/width of canvas
  this.align_canvas.height = ( json.num_alignments * 2 ) + 10 + 8 + 10 + 20;
  this.guide_canvas.height = align_canvas.height;
  this.cdiv = this.align_canvas.parentNode;
  this.cdiv.style.height = this.align_canvas.height;

  // Get drawing contexts
  this.guide_context = this.guide_canvas.getContext("2d");
  this.align_context = this.align_canvas.getContext("2d");

  // Heatmap Colors
  this.qualColor = ["#ff6600", "#ffcc00", "#ccff00", "#66ff00", "#00ff00",
    "#00ff66", "#00ffcc", "#00ccff", "#0066ff", "#0000ff"];

  // Constants to reduce lookup(?) in event listener
  this.WIDTH = this.align_canvas.width;
  this.HEIGHT = this.align_canvas.height;
  this.pixelToBP = this.json.length / this.WIDTH;
  this.currRulerY = 0;

  var that = this;
  this.guide_canvas.addEventListener( "mousemove", function( evt ) {
       that.mouseMoveHndlr( evt );
     }, false );

  this.render("norm");
}


AlignmentSummary.prototype.getMousePos = function(canvas, evt) {
    var rect = canvas.getBoundingClientRect();
    return {
        x: evt.clientX - rect.left,
        y: evt.clientY - rect.top
    };
};

AlignmentSummary.prototype.mouseMoveHndlr = function(evt){
    var mousePos = this.getMousePos(this.guide_canvas, evt);
    if (mousePos.x >= 10) {
        this.guide_context.clearRect(0, 0, this.WIDTH, this.HEIGHT);
        this.guide_context.strokeStyle = "#ff0000";
        this.guide_context.beginPath();
        this.guide_context.moveTo(mousePos.x, 0);
        this.guide_context.lineTo(mousePos.x, this.HEIGHT);
        this.guide_context.stroke();

        this.guide_context.font = "italic 11pt Calibri";
        var txt = "" + Math.round(((mousePos.x - 10) * this.pixelToBP) + 1);
        var text_width = this.guide_context.measureText(txt).width;
        var text_height = 12; //Estimated based on font ( no height call in HTML5 )
        var textXPos = mousePos.x - (text_width / 2);
        // TODO: Use coordinates of ruler....get somehow
        if (textXPos < 10) {
            textXPos = 10;
        }
        if (textXPos + text_width > this.WIDTH) {
            textXPos = this.WIDTH - this.text_width - 10;
        }

        this.guide_context.fillStyle = "#FAF7F8";
        this.guide_context.fillRect(textXPos, this.currRulerY, text_width, text_height);
        this.guide_context.fillStyle = "#000000";

        this.guide_context.fillText(txt, textXPos, this.currRulerY + 11);

    }
};


AlignmentSummary.prototype.ruler = function (x, y, width, height, minVal, maxVal, minorTickInterval, majorTickInterval) {
    this.align_context.beginPath();
    this.align_context.moveTo(x, y);
    this.align_context.lineTo(x + width, y);
    this.align_context.stroke();

    // Translations
    var pixelsPerUnit = width / (maxVal - minVal + 1);
    var pixelsPerMajorTick = majorTickInterval * pixelsPerUnit;
    var pixelsPerMinorTick = minorTickInterval * pixelsPerUnit;

    for (var i = 0; i < width; i += pixelsPerMajorTick) {
        this.align_context.beginPath();
        this.align_context.moveTo(x + i, y);
        this.align_context.lineTo(x + i, y + height);
        this.align_context.stroke();
    }

    for (var i = 0; i < width; i += pixelsPerMinorTick) {
        this.align_context.beginPath();
        this.align_context.moveTo(x + i, y);
        this.align_context.lineTo(x + i, y + (height / 2));
        this.align_context.stroke();
    }
};


//
//
//
AlignmentSummary.prototype.render = function(order) {
    // Visual Constants
    var divMargin = 10; // Left margin in div block in pixels
    var alignmentGlyphHeight = 1;
    var alignmentSpacing = alignmentGlyphHeight + 1;
    var rulerHeight = 8;
    var rulerVerticalMargin = 10;
    
    var viewWidth = this.align_canvas.width - divMargin; // Width of reference sequence in pixels
    var xScale = viewWidth / this.json.length;
    var alignments = this.json.alignments;
    var qualWidthBP = this.json.qualityBlockLen;
    
    // Clear overlayed canvases
    this.align_context.clearRect(0, 0, this.align_canvas.width,
    this.align_canvas.height);
    this.guide_context.clearRect(0, 0, this.guide_canvas.width,
    this.guide_canvas.height);

    // Select ordering
    if (order == "orient") {
        alignments.sort(function (a, b) {
            if (a[4] === b[4]) {
                if (a[4] === "R") {
                    if (a[1] === b[1]) {
                        return (b[2] - a[2]);
                    } else {
                        return (a[1] - b[1]);
                    }
                } else {
                    if (a[1] === b[1]) {
                        return (a[2] - b[2]);
                    } else {
                        return (b[1] - a[1]);
                    }
                }
            } else {
                return a[4] < b[4] ? -1 : a[4] > b[4] ? 1 : 0;
            }
        });
    } else if ( order == "end" )
    {
        alignments.sort(function (a, b) {
            if ((a[1] + a[2]) == (b[1] + b[2])) {
                return (a[1]  - b[1]);
            } else {
                return ((b[1] + b[2]) - (a[1] + a[2]));
            }
        });        
    } else if ( order == "div" )
    {
        alignments.sort(function (a, b) {
            return (a[5] - b[5]);
        });        
     }else {
        alignments.sort(function (a, b) {
            if (a[1] === b[1]) {
                return (b[2] - a[2]);
            } else {
                return (a[1] - b[1]);
            }
        });
    }


    var curY = 0;
    var referenceDrawn = 0;
    if (order != "orient") {
        this.ruler(divMargin, 0, viewWidth, rulerHeight, 1, 946, 10, 100);
        this.currRulerY = 0;
        curY = rulerVerticalMargin + rulerHeight;        
        referenceDrawn = 1;
    }
    for (var i = 0; i < alignments.length; i += 1) {
        if (referenceDrawn == 0 && alignments[i][4] == "R") {
            curY = curY + rulerVerticalMargin;
            this.ruler(divMargin, curY + (i * alignmentSpacing),
            viewWidth, rulerHeight, 1, 946, 10, 100);
            this.currRulerY = curY + (i * alignmentSpacing);
            curY = curY + rulerHeight + rulerVerticalMargin;
            referenceDrawn = 1;
        }

        var xOffset = alignments[i][1];
        var qualities = alignments[i][3];
        var qualIdx = 0;
        for (var j = 0; j < alignments[i][2]; j += qualWidthBP) {

            // TODO fix this indexing error
            if (qualIdx < qualities.length) {

                var grd = this.align_context.createLinearGradient(
                            0, 0, (xScale * qualWidthBP), 0 );

                var grad = this.qualColor[qualities[qualIdx] - 1] + " - ";
                grd.addColorStop(0, this.qualColor[qualities[qualIdx] - 1]);
                if (qualIdx == qualities.length - 1) {
                    grd.addColorStop(1, 
                                   this.qualColor[qualities[qualIdx] - 1]);
                    grad = grad + this.qualColor[qualities[qualIdx] - 1];
                } else {
                    grd.addColorStop(1, 
                                 this.qualColor[qualities[qualIdx + 1] - 1]);
                    grad = grad + this.qualColor[qualities[qualIdx + 1] - 1];
                }
                this.align_context.fillStyle = grd;
                this.align_context.fillRect(divMargin + (xOffset * xScale) +
                                       (j * xScale),
                curY + (i * alignmentSpacing), (xScale * qualWidthBP),
                alignmentGlyphHeight);
            }

            qualIdx++;
        }
    }
};
