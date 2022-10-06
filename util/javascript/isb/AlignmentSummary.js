//
// HTML5 Canvas Sequence Alignment Colleciton Summary Viewer
// Robert Hubley 2013-2015
//
//   It is common practice in bioinformatics to generate
//  a set of alignments given a single query sequence/model.
//  This shouldn't be confused with Multiple Sequence Alignment
//  , which attempts to find a mutually optimal alignment among all 
//  members of a collection.  
// 
//  The summary view is an extension of a whisker plot which depicts
//  the relative position and extents of each alignment as a single
//  horizontal bar anchored on the reference ( query ) sequence.
//  Here we also color the bar in non-overlapping windows using a 
//  quality metric from 1-10. The normal mode of the visualization
//  displays the alignments in: start position and length sorted order.
//  The "orient" mode displays the forward strand alignments on top
//  of the reference ( ruler ) bar and reverse strand hits below.
//
//  2/2/2021
//    Added an extension to visualize an MSA.  The only change was to 
//    allow a score of 0.  If 0 is used if the sequence is not present
//    ( e.g gap in MSA ), and 1 for the sequences for which there are
//    insertion bases present.
//  
//
//
//  Example invocation:
//  -------------------
//  HTML:
//    <div id="canvasesdiv" style="position:relative">
//        <canvas id="alignment_canvas" width="800" height="1600" style="z-index:1;position:absolute;left:0px;top:0px;">Canvas not supported</canvas>
//        <canvas id="guideline_canvas" width="800" height="1600" style="z-index:3;position:absolute;left:0px;top:0px;">Canvas not supported</canvas>
//        <canvas id="detail_canvas" width="1200" height="2000" style="z-index:2;position:absolute;left:0px;top:0px;">Canvas not supported</canvas>
//    </div>
//
// Javascript:
//   var summaryData = {
//       "qualityBlockLen": 10,
//       "length": 3302,
//       "seedStart": 2325,
//       "seedEnd": 2386,
//        "alignments": [
//           ["JH827946_1_1_1639", 2201, 169, [8, 7, 8, 5, 9, 2, 3, 8, 8, 9, 6, 8, 6, 7, 8, 9, 9], "F", "0.23", "1363", "1510"],
//           ["JH829156_0_2916_5235_R", 2136, 128, [8, 7, 7, 6, 8, 9, 9, 7, 6, 10, 6, 5, 9], "F", "0.24", "1", "120"],
//           ["JH829156_0_2916_5235_R", 2256, 109, [9, 7, 9, 9, 3, 5, 8, 8, 8, 7], "R", "0.16", "298", "396"],
//      ],
//      "num_alignments": 1252
//    };
//
//   var mySummary = new AlignmentSummary(
//                 document.getElementById('alignment_canvas'),
//                 document.getElementById('guideline_canvas'),
//                 document.getElementById('detail_canvas'),
//                 summaryData, {});
//
//  Example JSFIDDLE: http://jsfiddle.net/4wGm8/101/
//
function AlignmentSummary(align_canvas, guide_canvas, detail_canvas, json, options) {
    this.json = json;
    this.align_canvas = align_canvas;
    this.guide_canvas = guide_canvas;
    this.detail_canvas = detail_canvas;

    this.alignDetailVisible = false;
    this.alignDetailXPos = 0;
    this.alignDetailYPos = 0;
    this.alignDetailWidth = 400;
    this.alignDetailHeight = 150;

    // TODO: Use above general layout variables to determine canvas height
    this.align_canvas.height = (this.json.num_alignments * 2) + 10 + 8 + 10 + 10;
    this.guide_canvas.height = this.align_canvas.height;
    this.detail_canvas.height = this.align_canvas.height + this.alignDetailHeight;
    this.cdiv = this.align_canvas.parentNode;
    this.cdiv.style.height = this.align_canvas.height;

    // Get drawing contexts
    this.guide_context = this.guide_canvas.getContext("2d");
    this.align_context = this.align_canvas.getContext("2d");
    this.detail_context = this.detail_canvas.getContext("2d");

    // Heatmap Colors
    this.qualColor = ["#ff6600", "#ffcc00", "#ccff00", "#66ff00", "#00ff00",
        "#00ff66", "#00ffcc", "#00ccff", "#0066ff", "#0000ff"];

    // Constants to reduce lookup(?) in event listener
    this.WIDTH = this.align_canvas.width;
    this.HEIGHT = this.align_canvas.height;
    this.pixelToBP = this.json.length / this.WIDTH;
    this.currRulerY = 0;

    var that = this;
    this.guide_canvas.addEventListener("mousemove", function (evt) {
        that.mouseMoveHndlr(evt);
    }, false);
    this.guide_canvas.addEventListener("mousedown", function (evt) {
        that.mouseDownHndlr(evt);
    }, false);

    this.render("norm");
}


AlignmentSummary.prototype.getMousePos = function (canvas, evt) {
    var rect = canvas.getBoundingClientRect();
    return {
        x: evt.clientX - rect.left,
        y: evt.clientY - rect.top
    };
};


// Detect clicks on individual alignment lines and 
// draw a "alignDetail" box alongside the mouse pointer. 
// Subsequent clicks on the box will make it disappear.
AlignmentSummary.prototype.mouseDownHndlr = function (evt) {
    var mousePos = this.getMousePos(this.guide_canvas, evt);

    if (this.alignDetailVisible) {
        this.detail_context.beginPath();
        this.detail_context.rect(this.alignDetailXPos, this.alignDetailYPos,
        this.alignDetailWidth, this.alignDetailHeight);
        if (this.detail_context.isPointInPath(mousePos.x, mousePos.y)) {
            this.detail_context.clearRect(0, 0, this.detail_canvas.width,
            this.detail_canvas.height);
            this.alignDetailVisible = false;
            return;
        }
    }
    this.alignDetailXPos = mousePos.x;
    if ( this.alignDetailXPos + this.alignDetailWidth > this.detail_canvas.width )
    { 
      this.alignDetailXPos = this.alignDetailXPos -
                             ((this.alignDetailXPos + this.alignDetailWidth) - this.detail_canvas.width);
    }
    this.alignDetailYPos = mousePos.y;
    this.alignDetailVisible = true;
    var alignIdx = parseInt((mousePos.y - (8 + 10)) / (1 + 1));
    if ( alignIdx >= 0 )
    {
      this.drawAlignDetail2(this.alignDetailXPos, this.alignDetailYPos,
                           this.alignDetailWidth, this.alignDetailHeight, alignIdx);
    }
};

// Draw a popup box on the "detail_canvas" containing the alignment 
// details a given sequence from the alignment collection. 
// TODO: Testing variant.  This one displays the whole instance sequence rather than
//       discrete portions of the instance sequence.
AlignmentSummary.prototype.drawAlignDetail2 = function (x, y, width, height, alignIdx) {
    var popupWidth = width;
    var popupHeight = height;
    var titleHeight = 25;
    var margin = 10;
    var alignViewWidth = popupWidth - (2 * margin);
    var alignSpacing = 20;

    // Sequence clicked on by the user (alignIdx) may also used in other
    // alignments.  We are only interested in showing the alignments
    // that are nearby (maxGroupingDist) or overlaping the sequence the
    // user clicked on. 
    var name = this.json.alignments[alignIdx][0];
    var refStart = this.json.alignments[alignIdx][6];
    var refEnd = this.json.alignments[alignIdx][7];
    var maxGroupingDist = 10;

    // Identify the length of the genomic sequence that covers
    // all alignments we will be displaying.
    var idxs = [];
    var alignLevels = [];
    var minAlignPos = -1;
    var maxAlignPos = -1;
    for (var j = 0; j < this.json.alignments.length; j += 1) {
        var instStart = this.json.alignments[j][6];
        var instEnd = this.json.alignments[j][7];
        if (this.json.alignments[j][0] === name &&
            ((instStart > (refStart-maxGroupingDist) && instStart < (refEnd+maxGroupingDist))  || 
             (instEnd > (refStart-maxGroupingDist) && instEnd < (refEnd+maxGroupingDist)))) 
        {
            idxs[idxs.length] = this.json.alignments[j];
            if (minAlignPos == -1 || instStart < minAlignPos) minAlignPos = instStart;
            if (maxAlignPos == -1 || instEnd > maxAlignPos) maxAlignPos = instEnd;
        }
    }
    idxs.sort(function (a, b) {
        if (a[6] === b[6]) {
            return ((b[7] - b[6]) - (a[7] - a[6]));
        } else {
            return (a[6] - b[6]);
        }
    });

    var referenceXOffset;
    var alignedLen = maxAlignPos - minAlignPos + 1;
    var xsc;
    if (alignedLen > this.json.length) {
        xsc = alignViewWidth / alignedLen;
        referenceXOffset = parseInt(x + margin + (((alignedLen - this.json.length) / 2) * xsc));
    } else {
        referenceXOffset = parseInt(x + margin);
        xsc = alignViewWidth / this.json.length;
    }
    
    // Draw the popup frame
    var calcHeight = (2*margin) + titleHeight + ((idxs.length + 2) * alignSpacing);
    if ( calcHeight > height ){
      height = calcHeight;
    }
    this.detail_context.clearRect(0, 0, 1200, 2000);
    this.detail_context.fillStyle = "rgba(255, 255, 255, 1.0)";
    this.detail_context.fillRect(x, y, width, height);
    this.detail_context.beginPath();
    this.detail_context.rect(x, y, width, height);
    this.detail_context.lineWidth = 2;
    this.detail_context.strokeStyle = 'black';
    this.detail_context.stroke();

    // Write detail header
    this.detail_context.font = "15px Georgia";
    this.detail_context.fillStyle = 'black';
    this.detail_context.fillText(alignIdx + " : " + this.json.alignments[alignIdx][0], x + margin, y + margin + 5);

    // Draw forward strand reference line
    this.detail_context.beginPath();
    this.detail_context.lineWidth = 2;
    this.detail_context.strokeStyle = 'green';
    this.detail_context.moveTo(referenceXOffset, y + margin + titleHeight);
    this.detail_context.lineTo(referenceXOffset + parseInt(this.json.length * xsc), y + margin + titleHeight);
    this.detail_context.stroke();

    var levels = [];
    for (var j = 0; j < idxs.length; j += 1) {
        for (var k = 0; k <= levels.length; k += 1) {
            if (k == levels.length) {
                levels[k] = [];
                levels[k][0] = idxs[j];
                break;
            }

            var prevEle = levels[k][levels[k].length - 1];
            var prevEnd = prevEle[7];
            var curStart = idxs[j][6];
            if (parseInt(prevEnd) <= parseInt(curStart)) {
                levels[k][levels[k].length] = idxs[j];
                break;
            }
        }
    }


    // Draw main instance sequence line [gray]
    this.detail_context.strokeStyle = 'gray';
    this.detail_context.lineWidth = 2;
    var offset = parseInt(x + margin + ((alignViewWidth/2) - ((alignedLen * xsc)/2)));
    var alignYOffset = y + height - margin - (levels.length * (alignSpacing + 2));

    this.detail_context.beginPath();
    this.detail_context.moveTo(offset, alignYOffset + (levels.length * alignSpacing));
    this.detail_context.lineTo(offset+(alignedLen*xsc), alignYOffset + (levels.length * alignSpacing));
    this.detail_context.stroke();

    for (var j = 0; j < levels.length; j += 1) {
        for (var k = 0; k < levels[j].length; k += 1) {
            this.detail_context.beginPath();
            this.detail_context.strokeStyle = 'black';
            var start = offset + ((parseInt(levels[j][k][6]) - minAlignPos) * xsc);
            //var start = offset + ((levels[j][k][6]) * xsc);
            var end = start + ((parseInt(levels[j][k][7]) - parseInt(levels[j][k][6]) + 1) * xsc);

            // Draw alignment line
            this.detail_context.moveTo(start, alignYOffset + (j * alignSpacing));
            this.detail_context.lineWidth = 2;
            this.detail_context.lineTo(end, alignYOffset + (j * alignSpacing));
            if (levels[j][k][4] === "F") this.detail_context.strokeStyle = 'black';
            else this.detail_context.strokeStyle = 'red';
            this.detail_context.stroke();
            // write divergence
            this.detail_context.font = "8px Georgia";
            this.detail_context.fillStyle = 'black';
            this.detail_context.fillText(levels[j][k][5],
                                         start,alignYOffset + (j * alignSpacing) + 8);

            // Draw start connector
            this.detail_context.strokeStyle = 'black';
            this.detail_context.lineWidth = 1;
            this.detail_context.setLineDash([2, 3]);
            this.detail_context.beginPath();
            this.detail_context.moveTo(start, alignYOffset + (j * alignSpacing));
            this.detail_context.lineTo(parseInt(referenceXOffset + (levels[j][k][1] * xsc)), y + margin + titleHeight);
            this.detail_context.stroke();
            // Draw end connector
            this.detail_context.beginPath();
            this.detail_context.moveTo(end, alignYOffset + (j * alignSpacing));
            this.detail_context.lineTo(parseInt(referenceXOffset + ((levels[j][k][1] + levels[j][k][2]) * xsc)), y + margin + titleHeight);
            this.detail_context.stroke();

            this.detail_context.setLineDash([]);
        }
    }
};

// Draw a popup box on the "detail_canvas" containing the alignment 
// details a given sequence from the alignment collection. 
AlignmentSummary.prototype.drawAlignDetail = function (x, y, width, height, alignIdx) {
    var popupWidth = width;
    var popupHeight = height;
    var titleHeight = 25;
    var margin = 10;
    var alignViewWidth = popupWidth - (2 * margin);
    var alignSpacing = 20;

    this.detail_context.clearRect(0, 0, 1200, 2000);
    this.detail_context.fillStyle = "rgba(255, 255, 255, 1.0)";
    this.detail_context.fillRect(x, y, width, height);
    //
    this.detail_context.beginPath();
    this.detail_context.rect(x, y, width, height);
    this.detail_context.lineWidth = 2;
    this.detail_context.strokeStyle = 'black';
    this.detail_context.stroke();
    // Write detail header
    this.detail_context.font = "15px Georgia";
    this.detail_context.fillStyle = 'black';
    this.detail_context.fillText(this.json.alignments[alignIdx][0] + "    [ " + alignIdx + " ]",
    x + margin, y + margin + 5);

    var name = this.json.alignments[alignIdx][0];
    var idxs = [];
    var alignLevels = [];
    var minAlignPos = 99999999999999999999999999;
    var maxAlignPos = 0;
    for (var j = 0; j < this.json.alignments.length; j += 1) {
        if (this.json.alignments[j][0] === name) {
            idxs[idxs.length] = this.json.alignments[j];
            if (parseInt(this.json.alignments[j][6]) < minAlignPos) minAlignPos = this.json.alignments[j][6];
            if (parseInt(this.json.alignments[j][7]) > maxAlignPos) maxAlignPos = this.json.alignments[j][7];
        }
    }
    idxs.sort(function (a, b) {
        if (a[6] === b[6]) {
            return ((b[7] - b[6]) - (a[7] - a[6]));
        } else {
            return (a[6] - b[6]);
        }
    });
    var referenceXOffset;
    var alignedLen = maxAlignPos - minAlignPos + 1;
    var xsc;
    if (alignedLen > this.json.length) {
        xsc = alignViewWidth / alignedLen;
        referenceXOffset = parseInt(x + margin + (((alignedLen - this.json.length) / 2) * xsc));
    } else {
        referenceXOffset = parseInt(x + margin);
        xsc = alignViewWidth / this.json.length;
    }

    // Draw forward strand reference line
    this.detail_context.beginPath();
    this.detail_context.lineWidth = 2;
    this.detail_context.strokeStyle = 'green';
    this.detail_context.moveTo(referenceXOffset, y + margin + titleHeight);
    this.detail_context.lineTo(referenceXOffset + parseInt(this.json.length * xsc), y + margin + titleHeight);
    this.detail_context.stroke();

    var offset = parseInt(x + margin + ((alignViewWidth - (alignedLen * xsc)) / 2));

    var levels = [];
    for (var j = 0; j < idxs.length; j += 1) {
        for (var k = 0; k <= levels.length; k += 1) {
            if (k == levels.length) {
                levels[k] = [];
                levels[k][0] = idxs[j];
                break;
            }

            var prevEle = levels[k][levels[k].length - 1];
            var prevEnd = prevEle[7];
            var curStart = idxs[j][6];
            if (parseInt(prevEnd) <= parseInt(curStart)) {
                levels[k][levels[k].length] = idxs[j];
                break;
            }
        }
    }

    var alignYOffset = y + popupHeight - margin - (levels.length * (alignSpacing + 2));
    for (var j = 0; j < levels.length; j += 1) {
        for (var k = 0; k < levels[j].length; k += 1) {
            this.detail_context.beginPath();
            this.detail_context.strokeStyle = 'black';
            var start = offset + ((levels[j][k][6] - minAlignPos) * xsc);
            var end = start + ((levels[j][k][7] - levels[j][k][6] + 1) * xsc);

            // Draw alignment line
            this.detail_context.moveTo(start, alignYOffset + (j * alignSpacing));
            this.detail_context.lineWidth = 2;
            this.detail_context.lineTo(end, alignYOffset + (j * alignSpacing));
            if (levels[j][k][4] === "F") this.detail_context.strokeStyle = 'black';
            else this.detail_context.strokeStyle = 'red';
            this.detail_context.stroke();
            // write divergence
            this.detail_context.font = "8px Georgia";
            this.detail_context.fillStyle = 'black';
            this.detail_context.fillText(levels[j][k][5],
                                         start,alignYOffset + (j * alignSpacing) + 8);

            // Draw start connector
            this.detail_context.strokeStyle = 'black';
            this.detail_context.lineWidth = 1;
            this.detail_context.setLineDash([2, 3]);
            this.detail_context.beginPath();
            this.detail_context.moveTo(start, alignYOffset + (j * alignSpacing));
            this.detail_context.lineTo(parseInt(referenceXOffset + (levels[j][k][1] * xsc)), y + margin + titleHeight);
            this.detail_context.stroke();
            // Draw end connector
            this.detail_context.beginPath();
            this.detail_context.moveTo(end, alignYOffset + (j * alignSpacing));
            this.detail_context.lineTo(parseInt(referenceXOffset + ((levels[j][k][1] + levels[j][k][2]) * xsc)), y + margin + titleHeight);
            this.detail_context.stroke();

            this.detail_context.setLineDash([]);
        }
    }
};

AlignmentSummary.prototype.mouseMoveHndlr = function (evt) {
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
            textXPos = this.WIDTH - text_width;
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
AlignmentSummary.prototype.render = function (order) {
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
    this.detail_context.clearRect(0, 0, this.detail_canvas.width,
    this.detail_canvas.height);

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
    } else if (order == "end") {
        alignments.sort(function (a, b) {
            if ((a[1] + a[2]) == (b[1] + b[2])) {
                return (a[1] - b[1]);
            } else {
                return ((b[1] + b[2]) - (a[1] + a[2]));
            }
        });
    } else if (order == "div") {
        alignments.sort(function (a, b) {
            return (a[5] - b[5]);
        });
        //}else if (order == "groupById") {
        //    alignments.sort(function (a, b) {
        //
        //        if (a[0] < b[0]) {
        //            return -1;
        //        } else if (a[0] > b[0]) {
        //            return 1;
        //        } else { // nothing to split them
        //            return 0;
        //        }
        //    });
    } else {
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
                0, 0, (xScale * qualWidthBP), 0);

                //var grad = this.qualColor[qualities[qualIdx] - 1] + " - ";
                if (  qualities[qualIdx] > 0 ) {
                  grd.addColorStop(0, this.qualColor[qualities[qualIdx] - 1]);
                }else { 
                  grd.addColorStop(0,"#ffffff");
                }
                if (qualIdx == qualities.length - 1) {
                    if ( qualities[qualIdx+1] > 0 ) {
                      grd.addColorStop(1,
                      this.qualColor[qualities[qualIdx] - 1]);
                    }else {
                      grd.addColorStop(1,"#ffffff");
                    }
                    //grad = grad + this.qualColor[qualities[qualIdx] - 1];
                } else {
                    if ( qualities[qualIdx+1] > 0 ) {
                      grd.addColorStop(1,
                      this.qualColor[qualities[qualIdx + 1] - 1]);
                    }else {
                      grd.addColorStop(1,"#ffffff");
                    }
                    //grad = grad + this.qualColor[qualities[qualIdx + 1] - 1];
                }
                this.align_context.fillStyle = grd;
                this.align_context.fillRect(divMargin + (xOffset * xScale) + (j * xScale),
                curY + (i * alignmentSpacing), (xScale * qualWidthBP),
                alignmentGlyphHeight);
            }

            qualIdx++;
        }
    }
    if (this.json.seedStart) {
        this.align_context.fillStyle = "rgba(10, 10, 10, 0.25)";
        this.align_context.fillRect((divMargin + (this.json.seedStart * xScale)),
        rulerHeight + rulerVerticalMargin, ((this.json.seedEnd - this.json.seedStart + 1) * xScale), ((alignmentGlyphHeight + alignmentSpacing) * this.json.alignments.length));
    }

};

