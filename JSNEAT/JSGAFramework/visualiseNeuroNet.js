/* global d3 */

function getNeuronColor(layer) {
    var c = {
        B: "red",
        I: "green",
        H: "grey",
        O: "yellow"
    }
    return c[layer];
}

function getNeuronLayerName(name) {
    //console.log("name", name)
    var s = name.split(".");
    return s[0];
}

function getNeuronIndex(name) {
    var l = getNeuronLayerName(name);
    var s = name.split(".");
    //console.log(l);
    if (l === "I" || l === "O") {
        return +s[1];
    } else if (l === "B") {
        return 6;
    } else if (l === "H") {
        return +s[2] - 1;
    }
}

function getNeuronLayerIndex(name, maxLayer) {
    var l = getNeuronLayerName(name);
    if (l === "I" || l === "B") {
        return 0;
    } else if (l === "H") {
        var s = name.split(".");
        return +s[1];
    } else if (l === "O") {
        return maxLayer;
    }

}

function getNeuronPosition(name, maxLayer) {
    var x = getNeuronLayerIndex(name, maxLayer);
    var y = getNeuronIndex(name);
    return {x: x, y: y};
}

function renderMLPNEATNetwork(genome, network, maxLayer, Container) {
    var height = 500;
    var width = 700;
    var margin = {top: 50, bottom: 50, left: 70, right: 80};
    var innerWidth = width - margin.left - margin.right;
    var innerHeight = height - margin.top - margin.bottom;
    var maxLayer = genome.properties["hiddenLayers"] + 1;
    var radius = 18;

    var xscale = d3.scale.linear()
            .domain([0, maxLayer])
            .range([margin.left, innerWidth + margin.left]);
    var yscale = d3.scale.linear()
            .domain([0, 6])
            .range([margin.top, innerHeight + margin.top]);


    var container = d3.select("#" + Container);
    var svg;
    if (container.selectAll(".bestNet")[0].length === 0) {
        svg = container.append("svg").classed("bestNet", true).attr("height", height).attr("width", width);
    } else {
        svg = container.select(".bestNet");
    }
    svg.selectAll("*").remove();
    var data = [];
    var rawdata = network.neurons;
    for (var i in rawdata) {
        if (rawdata.hasOwnProperty(i)) {
            rawdata[i].nname = i;
            data.push(rawdata[i]);
        }
    }


    //add neurons
    var selection = svg.selectAll(".neuron").data(data);
    var enter = selection.enter().append("g")
            .classed("neuron", true);

    var group1 = enter.append("g").classed("group1", true).attr("transform", function (d) {
        var pos = getNeuronPosition(d.nname, maxLayer);
        var x = xscale(pos.x);
        var y = yscale(pos.y);
        d.pos = {x: x, y: y};
        return "translate(" + x + "," + y + ")";
    });
    ;
    var group2 = enter.append("g").classed("group2", true);

    var neuronCircle = group1.append("circle");
    neuronCircle.each(function (d) {
        //console.log(d.nname, pos.x, pos.y);
        //console.log(xscale(pos.x), yscale(pos.y))
        var layer = getNeuronLayerName(d.nname);
        d3.select(this)
                .attr("r", radius)
                .attr("stroke", "black")
                .attr("fill", getNeuronColor(layer))
                .attr("stroke-width", 1);

    });
    group1.append("text").each(
            function (d) {
                var layer = getNeuronLayerName(d.nname);
                if (layer !== "I" && layer !== "B") {
                    var name = d.activationFunction.name;
                    d3.select(this).text(name)
                            .attr("x", radius + 5)
                            .attr("fill", "black");
                }
            });


    var rawlinks = genome.chromosomes["NEAT"].genes;
    var links = [];
    for (var i = 0; i < rawlinks.length; i++) {
        var inN = rawdata[rawlinks[i].in].pos;
        var outN = rawdata[rawlinks[i].out].pos;
        links.push([inN, outN, rawlinks[i].weight]);
    }

   // console.log(links)
    var colorScale = d3.scale.linear()
            .domain([
                d3.min(rawlinks, function (d) {
                    return d.weight;
                }),
                d3.max(rawlinks, function (d) {
                    return d.weight;
                })
            ])
            .range([0, 255]);

    var linkGroup = svg.append("g").classed("linkgroup", true);
    var select2 = linkGroup.selectAll(".neuroLink").data(links);
    var enter = select2.enter();
    enter.append("line")
            .classed("neuroLink", true)
            .each(function (d, i) {
                var p1 = d[0];
                var p2 = d[1];
                var c = Math.floor(255 - colorScale(+d[2]));
                //console.log(c)
                d3.select(this)
                        .attr("x1", p1.x)
                        .attr("y1", p1.y)
                        .attr("x2", p2.x)
                        .attr("y2", p2.y)
                        .attr("style", "stroke:rgb(" + 0 + "," + c + "," + c + ")");

            });

}