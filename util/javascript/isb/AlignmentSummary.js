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

    // Shared horizontal scale.  When several MSAs are displayed together the
    // caller passes options.maxLength ( the longest MSA among them ) so that
    // every panel is rendered at the same bp-per-pixel scale and their rulers
    // line up.  When absent we fall back to this MSA's own length.
    this.maxLength = (options && options.maxLength) ? options.maxLength : this.json.length;

    this.alignDetailVisible = false;
    this.alignDetailXPos = 0;
    this.alignDetailYPos = 0;
    this.alignDetailWidth = 400;
    this.alignDetailHeight = 150;

    // TODO: Use above general layout variables to determine canvas height
    this.align_canvas.height = (this.json.num_alignments * 2) + 10 + 8 + 10 + 10 + 12;
    this.guide_canvas.height = this.align_canvas.height;
    this.detail_canvas.height = this.align_canvas.height + this.alignDetailHeight;
    this.cdiv = this.align_canvas.parentNode;
    this.cdiv.style.height = this.align_canvas.height + "px";
    // The detail canvas extends alignDetailHeight px below the others to draw
    // click-popups.  We deliberately do NOT reserve that space ( it would add a
    // large gap between stacked MSA panels ).  Because it is a transparent
    // overlay that overflows onto the next panel, disable its pointer events so
    // it cannot swallow clicks meant for that panel's buttons -- all mouse
    // handling lives on the guide canvas.
    this.detail_canvas.style.pointerEvents = "none";

    // Get drawing contexts
    this.guide_context = this.guide_canvas.getContext("2d");
    this.align_context = this.align_canvas.getContext("2d");
    this.detail_context = this.detail_canvas.getContext("2d");

    // 
    this.maxGroupingDist = 2000;

    // Heatmap Colors
    this.qualColor = ["#ff6600", "#ffcc00", "#ccff00", "#66ff00", "#00ff00",
        "#00ff66", "#00ffcc", "#00ccff", "#0066ff", "#0000ff"];

    // Constants to reduce lookup(?) in event listener
    this.WIDTH = this.align_canvas.width;
    this.HEIGHT = this.align_canvas.height;
    this.pixelToBP = this.maxLength / (this.WIDTH - 10);
    this.currRulerY = 0;
    // Right edge ( in pixels ) of the drawn, scaled alignment.  Updated in
    // render(); the guide overlay is confined to [10, rulerRightX].
    this.rulerRightX = this.WIDTH;

    var that = this;
    this.guide_canvas.addEventListener("mousemove", function (evt) {
        that.mouseMoveHndlr(evt);
    }, false);
    this.guide_canvas.addEventListener("mousedown", function (evt) {
        that.mouseDownHndlr(evt);
    }, false);

    this.render("norm", this.maxGroupingDist);
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
    // First alignment row sits below the ruler labels + ruler + margin
    // ( rulerLabelHeight 11 + rulerVerticalMargin 10 + rulerHeight 8 ); rows are
    // alignmentSpacing ( alignmentGlyphHeight 1 + 1 ) px tall.
    var alignIdx = parseInt((mousePos.y - (11 + 10 + 8)) / (1 + 1));
    if ( alignIdx >= 0 )
    {
      this.drawAlignDetail2(this.alignDetailXPos, this.alignDetailYPos,
                           this.alignDetailWidth, this.alignDetailHeight, alignIdx, this.maxGroupingDist);
    }
};

// Draw a popup box on the "detail_canvas" containing the alignment 
// details a given sequence from the alignment collection. 
// TODO: Testing variant.  This one displays the whole instance sequence rather than
//       discrete portions of the instance sequence.
AlignmentSummary.prototype.drawAlignDetail2 = function (x, y, width, height, alignIdx, maxGroupingDist) {
    var popupWidth = width;
    var popupHeight = height;
    var titleHeight = 40;
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
    this.detail_context.fillText("line " + alignIdx + " : " + this.json.alignments[alignIdx][0] + " : " + this.json.alignments[alignIdx][6] + "-" + this.json.alignments[alignIdx][7] , x + margin, y + margin + 5);

    // Draw forward strand reference line
    var refY = y + margin + titleHeight;
    var refLeftX = referenceXOffset;
    var refRightX = referenceXOffset + parseInt(this.json.length * xsc);
    this.detail_context.beginPath();
    this.detail_context.lineWidth = 2;
    this.detail_context.strokeStyle = 'green';
    this.detail_context.moveTo(refLeftX, refY);
    this.detail_context.lineTo(refRightX, refY);
    this.detail_context.stroke();

    // "MSA Reference: # bp" caption, one font-line above the coordinate-label
    // row, centered over the green bar.
    this.detail_context.font = "8px Georgia";
    this.detail_context.fillStyle = 'green';
    var refCaption = "MSA Reference: " + this.json.length + " bp";
    var refCaptionW = this.detail_context.measureText(refCaption).width;
    this.detail_context.fillText(refCaption,
        refLeftX + ((refRightX - refLeftX) / 2) - (refCaptionW / 2), refY - 13);

    // The connector-landing coordinates ( drawn in the instance loop below )
    // share a row just above the green bar; track occupied spans so they don't
    // overprint one another.
    var refLabelSpans = [];

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

    // Instance (grey) line near the bottom of the popup.
    var offset = parseInt(x + margin + ((alignViewWidth / 2) - ((alignedLen * xsc) / 2)));
    var greyY = y + height - margin;
    var greyRightX = offset + (alignedLen * xsc);
    this.detail_context.strokeStyle = 'gray';
    this.detail_context.lineWidth = 2;
    this.detail_context.beginPath();
    this.detail_context.moveTo(offset, greyY);
    this.detail_context.lineTo(greyRightX, greyY);
    this.detail_context.stroke();

    // Grey caption above the instance line: sequence id and its length.
    this.detail_context.font = "8px Georgia";
    this.detail_context.fillStyle = 'gray';
    this.detail_context.fillText(name + "  " + alignedLen + " bp", offset, greyY - 3);

    // Black fragment lines are centered vertically between the green reference
    // line ( top ) and the grey instance line ( bottom ).
    var blackBandHeight = (levels.length - 1) * alignSpacing;
    var blackBandTop = ((refY + greyY) / 2) - (blackBandHeight / 2);

    for (var j = 0; j < levels.length; j += 1) {
        var blackY = blackBandTop + (j * alignSpacing);
        for (var k = 0; k < levels[j].length; k += 1) {
            var start = offset + ((parseInt(levels[j][k][6]) - minAlignPos) * xsc);
            var end = start + ((parseInt(levels[j][k][7]) - parseInt(levels[j][k][6]) + 1) * xsc);

            // Draw alignment line
            this.detail_context.beginPath();
            this.detail_context.moveTo(start, blackY);
            this.detail_context.lineWidth = 2;
            this.detail_context.lineTo(end, blackY);
            if (levels[j][k][4] === "F") this.detail_context.strokeStyle = 'black';
            else this.detail_context.strokeStyle = 'red';
            this.detail_context.stroke();

            // write divergence, centered under the black fragment line
            this.detail_context.font = "8px Georgia";
            this.detail_context.fillStyle = 'black';
            var divLbl = "Div: " + levels[j][k][5];
            var divW = this.detail_context.measureText(divLbl).width;
            this.detail_context.fillText(divLbl,
                ((start + end) / 2) - (divW / 2), blackY + 8);

            // Dotted connectors: from where they land on the green reference bar,
            // through the black fragment endpoint, continuing straight down to
            // the grey instance line.
            var startLandX = parseInt(referenceXOffset + (levels[j][k][1] * xsc));
            var endLandX = parseInt(referenceXOffset +
                ((levels[j][k][1] + levels[j][k][2]) * xsc));
            this.detail_context.strokeStyle = 'black';
            this.detail_context.lineWidth = 1;
            this.detail_context.setLineDash([2, 3]);
            this.detail_context.beginPath();
            this.detail_context.moveTo(startLandX, refY);
            this.detail_context.lineTo(start, blackY);
            this.detail_context.lineTo(start, greyY);
            this.detail_context.stroke();
            this.detail_context.beginPath();
            this.detail_context.moveTo(endLandX, refY);
            this.detail_context.lineTo(end, blackY);
            this.detail_context.lineTo(end, greyY);
            this.detail_context.stroke();
            this.detail_context.setLineDash([]);

            // Label the reference coordinates where the connectors land on the
            // green bar, just above it, when there is horizontal room for them.
            this.detail_context.font = "8px Georgia";
            this.detail_context.fillStyle = 'green';
            var sLbl = "" + parseInt(levels[j][k][1]);
            var sLblW = this.detail_context.measureText(sLbl).width;
            if (this._reserveLabelSpan(refLabelSpans,
                    startLandX - (sLblW / 2), startLandX + (sLblW / 2))) {
                this.detail_context.fillText(sLbl, startLandX - (sLblW / 2), refY - 3);
            }
            var eLbl = "" + (parseInt(levels[j][k][1]) + parseInt(levels[j][k][2]));
            var eLblW = this.detail_context.measureText(eLbl).width;
            if (this._reserveLabelSpan(refLabelSpans,
                    endLandX - (eLblW / 2), endLandX + (eLblW / 2))) {
                this.detail_context.fillText(eLbl, endLandX - (eLblW / 2), refY - 3);
            }
        }
    }
};


// Reserve [lo,hi] on a list of occupied horizontal spans.  Returns false ( and
// reserves nothing ) when the range overlaps one already taken, so callers can
// skip a label that would overprint another.
AlignmentSummary.prototype._reserveLabelSpan = function (spans, lo, hi) {
    for (var s = 0; s < spans.length; s += 1) {
        if (lo < spans[s][1] && hi > spans[s][0]) {
            return false;
        }
    }
    spans.push([lo, hi]);
    return true;
};

AlignmentSummary.prototype.mouseMoveHndlr = function (evt) {
    var mousePos = this.getMousePos(this.guide_canvas, evt);
    this.guide_context.clearRect(0, 0, this.WIDTH, this.HEIGHT);
    // Confine the column marker to the drawn, scaled alignment region.
    if (mousePos.x >= 10 && mousePos.x <= this.rulerRightX) {
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
        if (textXPos + text_width > this.rulerRightX) {
            textXPos = this.rulerRightX - text_width;
        }

        this.guide_context.fillStyle = "#FAF7F8";
        this.guide_context.fillRect(textXPos, this.currRulerY, text_width, text_height);
        this.guide_context.fillStyle = "#000000";

        this.guide_context.fillText(txt, textXPos, this.currRulerY + 11);

    }
};


// Choose a "nice" tick interval ( 1/2/5 x 10^n ) so that roughly
// targetTicks major ticks span the given range.
AlignmentSummary.prototype.niceTickInterval = function (range, targetTicks) {
    if (range <= 0) {
        return 1;
    }
    var raw = range / targetTicks;
    var mag = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10));
    var norm = raw / mag;
    // Round the interval *up* to the next nice value so the resulting tick
    // spacing is at least what the caller asked for ( never tighter ).
    var nice;
    if (norm <= 1) {
        nice = 1;
    } else if (norm <= 2) {
        nice = 2;
    } else if (norm <= 5) {
        nice = 5;
    } else {
        nice = 10;
    }
    return nice * mag;
};


// Draw a horizontal ruler from value minVal ( at x ) to maxVal ( at x+width ).
// Major ticks run the full height, minor ticks half height.  The endpoints
// ( minVal and maxVal ) are always labeled, and interior major ticks carry
// their value as well.
AlignmentSummary.prototype.ruler = function (x, y, width, height, minVal, maxVal, minorTickInterval, majorTickInterval) {
    var ctx = this.align_context;
    ctx.strokeStyle = "#000000";
    ctx.fillStyle = "#000000";
    ctx.lineWidth = 1;

    // Baseline
    ctx.beginPath();
    ctx.moveTo(x, y);
    ctx.lineTo(x + width, y);
    ctx.stroke();

    // Map a value in [minVal,maxVal] linearly onto [x, x+width].
    var range = maxVal - minVal;
    if (range < 1) {
        range = 1;
    }
    var pixelsPerUnit = width / range;

    // Minor ticks
    for (var v = minVal; v <= maxVal; v += minorTickInterval) {
        var px = x + (v - minVal) * pixelsPerUnit;
        ctx.beginPath();
        ctx.moveTo(px, y);
        ctx.lineTo(px, y + (height / 2));
        ctx.stroke();
    }

    // Major ticks.  The loop starts on minVal ( the start tick ); maxVal rarely
    // lands on a major step, so draw an explicit full-height end tick to match.
    for (var v = minVal; v <= maxVal + 0.0001; v += majorTickInterval) {
        var px = x + (v - minVal) * pixelsPerUnit;
        ctx.beginPath();
        ctx.moveTo(px, y);
        ctx.lineTo(px, y + height);
        ctx.stroke();
    }
    ctx.beginPath();
    ctx.moveTo(x + width, y);
    ctx.lineTo(x + width, y + height);
    ctx.stroke();

    // Labels: only the start/end points of the MSA are labeled.
    ctx.font = "9px Calibri";
    ctx.textBaseline = "bottom";
    ctx.textAlign = "left";
    ctx.fillText(minVal, x, y - 1);
    ctx.textAlign = "right";
    ctx.fillText(maxVal, x + width, y - 1);

    // Restore text defaults
    ctx.textAlign = "left";
    ctx.textBaseline = "alphabetic";
};


//
//
//
AlignmentSummary.prototype.render = function (order, maxGroupingDist) {
    // Visual Constants
    var divMargin = 10; // Left margin in div block in pixels
    var alignmentGlyphHeight = 1;
    var alignmentSpacing = alignmentGlyphHeight + 1;
    var rulerHeight = 8;
    var rulerVerticalMargin = 10;

    var viewWidth = this.align_canvas.width - divMargin; // Width of reference sequence in pixels
    // Normalize: scale against the longest MSA so that all panels share one
    // bp-per-pixel scale.  This panel's ruler only spans its own length.
    var xScale = viewWidth / this.maxLength;
    var rulerWidth = this.json.length * xScale;
    this.rulerRightX = divMargin + rulerWidth;
    // Derive tick spacing from the ruler's on-screen width so that ticks never
    // crowd closer than these pixel thresholds ( and thin out automatically
    // when this panel is short relative to the shared scale ).
    var minMajorPx = 70; // minimum pixels between full-height ticks
    var minMinorPx = 10; // minimum pixels between minor ticks
    var majorTick = this.niceTickInterval(this.json.length,
        Math.max(1, Math.floor(rulerWidth / minMajorPx)));
    var minorTick = this.niceTickInterval(this.json.length,
        Math.max(1, Math.floor(rulerWidth / minMinorPx)));
    if (minorTick > majorTick) {
        minorTick = majorTick;
    }
    var rulerLabelHeight = 11; // vertical room reserved above the ruler for labels
    var alignments = this.json.alignments;
    var qualWidthBP = this.json.qualityBlockLen;

    // Clear overlayed canvases
    this.align_context.clearRect(0, 0, this.align_canvas.width,
    this.align_canvas.height);
    this.guide_context.clearRect(0, 0, this.guide_canvas.width,
    this.guide_canvas.height);
    this.detail_context.clearRect(0, 0, this.detail_canvas.width,
    this.detail_canvas.height);

    // Reset the max grouping dist, but only when the caller supplied one.  The
    // sort buttons call render(order) with no second argument; clobbering this
    // with undefined would break drawAlignDetail2's proximity filter ( every
    // comparison becomes NaN ) so the detail popup loses all instance lines.
    if (maxGroupingDist !== undefined) {
        this.maxGroupingDist = maxGroupingDist;
    }

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
    }else if (order == "groupById") {
        alignments.sort(function (a, b) {
               if (a[0] < b[0]) {
                  return -1;
                } else if (a[0] > b[0]) {
                  return 1;
                }else if (a[1] === b[1]) {
                  return (b[2] - a[2]);
                } else {
                  return (a[1] - b[1]);
                }
            });
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

    // Orientation sort normally splits the ruler between forward ( above ) and
    // reverse ( below ) strands.  With no reverse strands there is nothing to
    // place below the ruler, so draw it at the top like the other sorts.
    var hasReverse = false;
    for (var r = 0; r < alignments.length; r += 1) {
        if (alignments[r][4] == "R") {
            hasReverse = true;
            break;
        }
    }
    if (order != "orient" || !hasReverse) {
        this.ruler(divMargin, rulerLabelHeight, rulerWidth, rulerHeight,
        1, this.json.length, minorTick, majorTick);
        this.currRulerY = rulerLabelHeight;
        curY = rulerLabelHeight + rulerVerticalMargin + rulerHeight;
        referenceDrawn = 1;
    }
    for (var i = 0; i < alignments.length; i += 1) {
        if (referenceDrawn == 0 && alignments[i][4] == "R") {
            curY = curY + rulerVerticalMargin + rulerLabelHeight;
            var rulerY = curY + (i * alignmentSpacing);
            this.ruler(divMargin, rulerY,
            rulerWidth, rulerHeight, 1, this.json.length, minorTick, majorTick);
            this.currRulerY = rulerY;
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

