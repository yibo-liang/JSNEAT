<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8">
        <title>D3: Drawing divs with data</title>

        <script src="http://d3js.org/d3.v3.min.js"></script>

        <script type="text/javascript" src="debug.js"></script>.
        <script type="text/javascript" src="JSNEAT.js"></script>
        <script type="text/javascript" src="core.js"></script>
        <script type="text/javascript" src="graphic_render.js"></script>
        <script type="text/javascript" src="controller.js"></script>
        <style type="text/css">

            div.fullDiv {
                height: 100%;
                width: 100%;
                left: 0;
                top: 0;
                overflow: hidden;
                position: fixed;
                background-color: black;
                z-index: 0;
            }

            div.cell {
                position: absolute;
                background-color: red;
                border-radius: 50%;
                border-top: 5px solid white;
                border-left: 3.5px solid transparent;
                border-right: 3.5px solid transparent;


            }

            div.fps{
                left: 0;
                top: 0;
                position: fixed;
                background-color: white;
                z-index: 4;
            }

            div.canvas_container{
                left:0px;
                top:0;
                position: fixed;
                background-color: black;
                width: 1400px;
                height: 800px
            }

            div.panel{
                left:1400px;
                top:0px;
                position: fixed;
            }

            div.control{
                top: 800px;
                position: fixed;
            }

        </style>
    </head>
    <body>
        <div class="fps" id="fps"></div>
        <div class="canvas_container" id="canvas_container"></div>
        <div class="control"><button type="button" onclick="speedToggle();">speed mode toggle</button></div>
        <script type="text/javascript">


            initialise("canvas_container");
            var info = new Info(0, cellpopulation, foodamount);
            var pa = d3.select("body").append("div");
            pa.classed("panel", true);
            panel = d3.select("body").select(".panel").selectAll("p").data([info]);
            panel.enter().append("p").text(function (d) {
                return  "gen:" + d.gen
                        + ", food n=" + d.foodnum
                        + ", cell n=" + d.cellnum;
            });

            var cells = createCells(cellpopulation);
            var genomes = initialisePool(cellDefault.inputNum, cellDefault.outputNum, cellpopulation);
            for (var i = 0; i < cellpopulation; i++) {
                cells[i].genome = genomes[i];
            }

            initFoods();

            //console.log(cells);
            var detachedContainer = document.createElement("container");
            var dataContainer = d3.select(detachedContainer);





            for (var i = 0; i < cellpopulation; i++) {
                cellgrid.updateObjPos(cells[i]);
            }

            var fstart = new Date().getTime();
            var framecount = 0;
            var current = 0;


            d3.timer(function () {
                stepframe();


                framecount++;
                if (framecount % 18 === 0) {
                    var fend = new Date().getTime();
                    var dt = (fend - fstart) / 1000;
                    var df = framecount - current;
                    var fps = df / dt;
                    document.getElementById("fps").innerHTML = "fps:" + Math.floor(fps);
                    current = framecount;
                    fstart = fend;
                }
            });

        </script>



    </body>
</html>