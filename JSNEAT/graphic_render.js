/* 
 * The MIT License
 *
 * Copyright 2015 yl9.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


var speed = 4;
var render_count = 0;
var last_count = 0;
var last = 0;

var foodCanvas;
var foodContext;

var cellCanvas;
var cellContext;

function initialise(container) {
    var new_canvas = d3.select("#" + container)
            .append("canvas")
            .attr("width", poolWidth)
            .attr("height", poolHeight)
            .style("position", "absolute")
            .style("z-index", "2");

    cellCanvas = new_canvas;
    cellContext = new_canvas.node().getContext("2d");

    var new_canvas2 = d3.select("#" + container)
            .append("canvas")
            .attr("width", poolWidth)
            .attr("height", poolHeight)
            .style("position", "absolute")
            .style("z-index", "1");

    foodCanvas = new_canvas2;
    foodContext = new_canvas2.node().getContext("2d");


}

function render_eraseLastFrame2() {
    cellContext.clearRect(0, 0, poolWidth, poolHeight);
    foodContext.clearRect(0, 0, poolWidth, poolHeight);
}

function render_eraseLastFrame(cellData) {
    //cellContext.clearRect(0, 0, poolWidth, poolHeight);
    cellData.each(function (d, i) {
        //console.log(d);
        var x = Math.floor(d.x);
        var y = Math.floor(d.y);
        var lastx = Math.floor(d.lastx);
        var lasty = Math.floor(d.lasty);

        if (x === lastx && y === lasty) {
            return;
        }
        var r = Math.floor(d.radius + 1);
        var sx = lastx - r - 1;
        var sy = lasty - r - 1;
        var ex = lastx + r + 1;
        var ey = lasty + r + 1;
        if (sx < 0)
            sx = 0;
        if (sy < 0)
            sy = 0;
        if (ex > poolWidth)
            ex = poolWidth;
        if (ey > poolHeight)
            ey = poolHeight;
        cellContext.clearRect(sx, sy, ex, ey);

    });
}

function render_cellbody(cellData)
{
    cellData.each(function (d, i) {
        //console.log(d);
        cellContext.beginPath();
        var x = Math.floor(d.x);
        var y = Math.floor(d.y);

        cellContext.arc(x, y, d.radius, 0, 2 * Math.PI);
        cellContext.fillStyle = "rgb(" + d.color[0] + "," + d.color[1] + "," + d.color[2] + ")";
        if (cells.length == 1) {
            cellContext.strokeStyle = "red";
            cellContext.stroke();
        }
        cellContext.fill();

        cellContext.closePath();
    });
}

function render_foods(foodData) {
    foodData.each(function (d, i) {
        if (d.is_eaten)
            return;
        foodContext.beginPath();
        var x = Math.floor(d.x);
        var y = Math.floor(d.y);

        foodContext.arc(x, y, d.radius, 0, 2 * Math.PI);
        foodContext.fillStyle = "rgb(" + d.color[0] + "," + d.color[1] + "," + d.color[2] + ")";
        foodContext.fill();

        foodContext.closePath();
    });
}
var debugid = 0;
function render_cellfeature(cellData) {

    debugid = cells[0].id;

    cellData.each(function (d, i) {
        var headx = Math.floor(d.x + d.radius * Math.cos(d.r));
        var heady = Math.floor(d.y + d.radius * Math.sin(d.r));
        cellContext.beginPath();
        cellContext.lineWidth = 2;
        // set line color
        if (d.id === debugid) {
            cellContext.strokeStyle = "white";
            cellContext.arc(d.x, d.y, 107, 0, 2 * Math.PI);
            cellContext.stroke();
            cellContext.closePath();
            var a = d.sense();

            if (a[a.length - 1] === 1) {
                cellContext.beginPath();

                cellContext.strokeStyle = "yellow";
                cellContext.arc(d.x, d.y, 15, 0, 2 * Math.PI);
                cellContext.stroke();
                cellContext.closePath();
            }

            //console.log(a);
          

            if (a[0] < 107) {
                //console.log(a);
                cellContext.beginPath();
                cellContext.strokeStyle = "red";
                cellContext.moveTo(d.x, d.y);
                cellContext.lineTo(d.x + a[0] * Math.cos(a[1] + d.r), d.y + a[0] * Math.sin(a[1] + d.r));
                cellContext.stroke();
                cellContext.closePath();

            }
            if (a[2] < 107) {
                //console.log(a);
                cellContext.beginPath();
                cellContext.strokeStyle = "yellow";
                cellContext.moveTo(d.x, d.y);
                cellContext.lineTo(d.x + a[2] * Math.cos(a[3] + d.r), d.y + a[2] * Math.sin(a[3] + d.r));
                cellContext.stroke();
                cellContext.closePath();

            }
            var pos = foodgrid.xyTocr(d.x, d.y);
            var col = pos[0];
            var row = pos[1];
            var grids = foodgrid.getNeighborGrids(col, row);

            for (var i = 0; i < grids.length; i++) {

                for (var j = 0; j < grids[i].length; j++) {
                    cellContext.beginPath();
                    cellContext.strokeStyle = "yellow";
                    cellContext.rect(grids[i][j].x - 10, grids[i][j].y - 10, 20, 20);
                    cellContext.stroke();
                    cellContext.closePath();
                }
            }
            var pos = cellgrid.xyTocr(d.x, d.y);
            var col = pos[0];
            var row = pos[1];
            var grids = cellgrid.getNeighborGrids(col, row);
            //console.log(neignbors.length);

            for (var i = 0; i < grids.length; i++) {

                for (var j = 0; j < grids[i].length; j++) {
                    cellContext.beginPath();
                    cellContext.strokeStyle = "red";
                    cellContext.rect(grids[i][j].x - 12, grids[i][j].y - 12, 24, 24);
                    cellContext.stroke();
                    cellContext.closePath();
                }
            }
        } else {

        }

        cellContext.beginPath();

        cellContext.strokeStyle = '#0f3af0';
        cellContext.moveTo(d.x, d.y);
        cellContext.lineTo(headx, heady);
        cellContext.stroke();
        cellContext.closePath();
    });
}


function render(cellData, foodData) {
    render_eraseLastFrame2();
    //console.log("renderbody");
    render_cellbody(cellData);

    render_foods(foodData);
    //console.log("renderfeature");
    render_cellfeature(cellData);
}