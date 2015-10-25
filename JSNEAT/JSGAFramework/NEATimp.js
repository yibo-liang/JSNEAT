var MaxNodes = 10;


var useBias = true;

var isRecurrent = false;


var CrossoverChance = 0.75;
var AdoptionChanche = 0.01;

var PointMutationChance = 0.15;
var LinkMutationChance = 0.5;
var NodeMutationChance = 0.3;
var BiasMutationChance = 0.10;

var DisableMutationChance = 0.4;
var EnableMutationChance = 0.2;

var PerturbStep = 0.1;
var PerturbChance = 0.80;

var PruningLinkChance = 0.003;
var PruningNodeChance = 0.005;


var NeuronGenotype = function (type, index) {
    this.type = type;
    this.index = index;
};
NeuronGenotype.prototype.copy = function () {
    var result = new NeuronGenotype(this.type, this.index);
    return result;
};

function newInnovation(pool) {
    var i = pool.innovation;
    pool.innovation++;
    return i;
}

var Gene = function () {
    this.in = null;
    this.out = null;
    this.weight = 0;
    this.enable = true;
    this.innovation = 0;
};

Gene.prototype.copy = function () {
    var copy = new Gene();
    copy.in = this.in;
    copy.out = this.out;
    copy.weight = this.weight;
    copy.enable = this.enable;
    copy.innovation = this.innovation;
    return copy;
};

function sigmoid(x) {
    return 2 / (1 + Math.exp(-4.9 * x)) - 1;
}

var Neuron = function () {
    this.incomingLinks = [];
    this.outcomingLinks = [];
    this.cValue = 0;
    this.nValue = 0;
    this.thresholdFunction = sigmoid;
};

var Network = function (chromesome) {
    this.chromesome = chromesome;
    this.neurons = {};
    this.indices = [];
};
Network.prototype.copy = function () {
    return generateNN(this.chromesome);
};

function generateNN(chromesome) {
    var network = new Network(chromesome);

    //console.log(chromesome);
    var inputNeuronNum = chromesome.genome.properties["inputs"];
    var outputNeuronNum = chromesome.genome.properties["outputs"];
    for (var i = 0; i < inputNeuronNum; i++) {
        var n = new Neuron();
        network.neurons["I" + i] = n;
        network.indices.push("I" + i);
    }
    for (var i = 0; i < outputNeuronNum; i++) {
        var n = new Neuron();

        network.neurons["O" + i] = n;
        network.indices.push("O" + i);
    }
    if (useBias === true) {
        var n = new Neuron();
        n.cValue = 1;
        n.nValue = 1;
        network.neurons["B"] = n;
        network.indices.push("B");
    }


    for (var i = 0; i < chromesome.genes.length; i++) {
        var gene = chromesome.genes[i];
        if (gene.enable === true) {
            if (typeof network.neurons[gene.out] === "undefined") {
                network.neurons[gene.out] = new Neuron();
                network.indices.push(gene.out);
            }
            var neuron_i = gene.out;
            var neuron = network.neurons[neuron_i];

            neuron.incomingLinks.push(gene);
            if (typeof network.neurons[gene.in] === "undefined") {
                network.neurons[gene.in] = new Neuron();
                network.indices.push(gene.in);
            }

            for (var k = 0; k < chromesome.genes.length; k++) {
                var genek = chromesome.genes[k];
                if (genek.in === neuron_i) {
                    neuron.outcomingLinks.push(genek);
                }
            }
        }
    }

    chromesome.genome.properties["NeuroNetwork"] = network;
    //console.log("generate nn=", network);
    return network;


}

function evaluateNeuroNetwork(genome, network, inputArray, isDebuging) {
    var inputVec = [];
    for (var i = 0; i < inputArray.length; i++) {
        inputVec.push(inputArray[i]);
    }
    var inputs = genome.properties["inputs"];
    if (inputVec.length !== inputs) {
        console.log("ERROR: INPUT NUMBER IS DIFFERENT FROM NETWORK", inputVec);
        return;
    }


    for (var i = 0; i < inputs; i++) {
        network.neurons["I" + i].cValue = inputVec[i];
    }
    var debugstr = "";
    for (var i in network.neurons) {
        if (network.neurons.hasOwnProperty(i)) {
            if (i === "B")
                continue
            var neuron = network.neurons[i];

            debugstr += "Debug neuron i='" + i + "';\n";

            //console.log(i, neuron);
            var sum = 0;
            for (var j = 0; j < neuron.incomingLinks.length; j++) {
                var incoming = neuron.incomingLinks[j];
                var other = network.neurons[incoming.in];
                //console.log(incoming, other);
                sum += incoming.weight * other.cValue;
                debugstr += "value += n[" + incoming.in + "](" + other.cValue + ")*" + incoming.weight + "\n";
                //console.log("sum=",sum);
            }
            debugstr += "value=" + sum + ", sigmoid(value)=" + neuron.thresholdFunction(sum) + "\n";
            if (neuron.incomingLinks.length > 0) {
                if (typeof isDebuging !== "undefined" && isDebuging) {
                    console.log(debugstr);
                    console.log(network);

                }
                var temp = neuron.thresholdFunction(sum);
                //if recurrent network is used, then save the new value for next iteration,
                //otherwise, save the new value for current iteration
                if (isRecurrent === true) {
                    neuron.nValue = temp;
                } else {
                    neuron.cValue = temp;
                }
                //console.log("v=",neuron.v);
            }

            //console.log("i=", i, "  v=", network.neurons[i].v);
        }
    }

    if (isRecurrent === true)
        for (var i in network.neurons) {
            if (network.neurons.hasOwnProperty(i)) {
                var neuron = network.neurons[i];
                neuron.cValue = neuron.nValue;
            }
        }

    var outputArray = [];

    var outputs = genome.properties["outputs"];
    for (var i = 0; i < outputs; i++) {
        var v = network.neurons["O" + i].cValue;
        //if (v > 0) {
        //    outputArray.push(true);
        //} else {
        //    outputArray.push(false);
        //}
        outputArray.push(v);
    }
    return outputArray;

}

function disjoint(geneslist1, geneslist2) {
    var l1 = {};
    var l2 = {};

    for (var i = 0; i < geneslist1.length; i++) {
        var gene = geneslist1[i];
        l1[gene.innovation] = true;
    }

    for (var i = 0; i < geneslist2.length; i++) {
        var gene = geneslist2[i];
        l2[gene.innovation] = true;
    }
    //console.log(l1, l2);
    var disjointGenes = 0;
    for (var i = 0; i < geneslist1.length; i++) {

        var gene = geneslist1[i];
        //console.log(i, gene, genes1);
        if (typeof l2[gene.innovation] === "undefined" || l2[gene.innovation] !== true) {
            disjointGenes = disjointGenes + 1;
        }

    }
    for (var i = 0; i < geneslist2.length; i++) {

        var gene = geneslist2[i];
        if (typeof l1[gene.innovation] === "undefined" || l1[gene.innovation] !== true) {
            disjointGenes = disjointGenes + 1;
        }

    }

    var n = Math.max(geneslist1.length, geneslist2.length);
    //console.log(disjointGenes,n);
    if (n === 0)
        return 0;
    return disjointGenes / n;
}


function weights(geneslist1, geneslist2) {
    //console.log("weights:",genes1, genes2);
    var i2 = {};
    for (var i = 1; i < geneslist2.length; i++) {
        var gene = geneslist2[i];
        i2[+gene.innovation] = gene;
    }
    //console.log(i2);
    var sum = 0;
    var coincident = 0;
    for (var i = 1; i < geneslist1.length; i++) {
        var gene = geneslist1[i];
        if (typeof i2[+gene.innovation] !== "undefined") {
            var gene2 = i2[+gene.innovation];
            sum = sum + Math.abs(gene.weight - gene2.weight);
            coincident = coincident + 1;
        }
    }
    //console.log("sum, coinci", sum, coincident, sum/coincident);
    if (coincident === 0) {
        return 1;
    }
    return sum / coincident;
}

var NEATcrossover = function (chromesome1, chromesome2) {


    var child = new Chromesome(
            chromesome1.genome,
            "NEAT",
            chromesome1.evaluator,
            chromesome1.chromesomeOperator,
            chromesome1.deltaFunction,
            chromesome1.deltaWeight
            );
    var innovations2 = {};
    for (var i = 0; i < chromesome2.genes.length; i++) {
        var gene = chromesome2.genes[i];
        innovations2[gene.innovation] = gene;
    }


    for (var i = 0; i < chromesome1.genes.length; i++) {
        var gene1 = chromesome1.genes[i];
        var gene2 = innovations2[gene1.innovation];
        if (typeof gene2 !== "undefined" && gene2.enable === true) {
            child.genes.push(gene2.copy());
        } else {
            child.genes.push(gene1.copy());
        }
    }

    return child;
};

function isInputNeuron(index) {
    if (index.substring(0, 1) === "I") {
        return true;
    }
    return false;
}

function isBiasNeuron(index) {
    return index === "B";
}

function randomNeuron(genome, NoInputs, NoBias) {
    var candidates = [];

    //console.log(genome.properties["NeuroNetwork"]);
    var indices = genome.properties["NeuroNetwork"].indices;
    for (var i = 0; i < indices.length; i++) {
        var index = indices[i];
        if (NoInputs && isInputNeuron(index) || NoBias && isBiasNeuron(index)) {

        } else {
            candidates.push(index);
        }
    }
    var i = RandomIntInclusive(0, candidates.length - 1);
    //console.log(candidates)
    return candidates[i];
}

function containsLink(genes, link) {
    for (var i = 0; i < genes.length; i++) {
        var gene = genes[i];
        if (gene.in === link.in && gene.out === link.out) {
            return true;
        }
    }
    return false;
}


function pointMutate(chromesome) {
    //console.log("pointMutation");
    var step = chromesome.genome.properties["PertubStep"];
    var chance = chromesome.genome.properties["PertubChance"];
    for (var i = 0; i < chromesome.genes.length; i++) {
        var gene = chromesome.genes[i];
        if (Math.random() < chance) {
            gene.weight = gene.weight + Math.random() * step * 2 - step;
        } else {
            gene.weight = Math.random() * step * 4 - 2;
        }
    }

}



function linkMutation(chromesome, forceBias) {
    //console.log("linkMutation");

    var neuron1 = randomNeuron(chromesome.genome, false, false);
    var neuron2 = randomNeuron(chromesome.genome, true, true);
    var inputs = chromesome.genome.properties["inputs"];

    if (isInputNeuron(neuron1) && isInputNeuron(neuron2)) {
        // both inputs neuron
        return;
    }

    if (neuron1 === neuron2 && !isRecurrent) {
        //if input is output but recurrent is not allowed
        return;
    }

    if (isInputNeuron(neuron2)) {
        //swap 1 and 2
        var temp = neuron1;
        neuron1 = neuron2;
        neuron2 = temp;
    }

    var newLink = new Gene();
    newLink.in = neuron1;
    newLink.out = neuron2;
    if (forceBias) {
        newLink.in = "B";
    }
    if (containsLink(chromesome.genes, newLink)) {
        return;
    }
    newLink.innovation = newInnovation(chromesome.genome.pool);
    //console.log("new link=", newLink)
    newLink.weight = Math.random() * 4 - 2;
    chromesome.genes.push(newLink);
}


function nodeMutation(chromesome) {
    var genome = chromesome.genome;
    if (genome.properties["hiddenNeurons"] >= MaxNodes) {
        return;
    }

    if (chromesome.genes.length === 0) {
        return;
    }

    var gene = chromesome.genes[RandomIntInclusive(0, chromesome.genes.length - 1)];
    if (gene.enable === false) {
        return;
    }

    genome.properties["hiddenNeurons"] += 1;

    gene.enable = false;
    var gene1 = gene.copy();
    gene1.out = "H" + genome.properties["hiddenNeurons"];
    gene1.weight = 1.0;
    gene1.innovation = newInnovation(genome.pool);
    gene1.enable = true;
    chromesome.genes.push(gene1);
    var gene2 = gene.copy();
    gene2.in = "H" + genome.properties["hiddenNeurons"];
    gene2.innovation = newInnovation(genome.pool);
    gene2.enable = true;
    chromesome.genes.push(gene2);
    //console.log(gene, ",", gene1, ",", gene2);
    //console.log();
}



function enableDisableMutate(chromesome, enable) {
    var candidates = [];
    for (var i = 0; i < chromesome.genes.length; i++) {
        if (chromesome.genes[i].enable !== enable) {
            candidates.push(chromesome.genes[i]);
        }
    }
    if (candidates.length === 0) {
        return;
    }
    var gene = candidates[RandomIntInclusive(0, candidates.length - 1)];
    gene.enable = enable;
}


function deleteOneLink(chromesome, gene) {
    var len1 = chromesome.genes.length;
    for (var i = 0; i < chromesome.genes.length; i++) {
        var g = chromesome.genes[i];
        if (gene.in === g.in && g.out === gene.out)
        {

            var deleted = chromesome.genes.splice(i, 1);
            var len2 = chromesome.genes.length;
            console.log("delete link", gene.in, gene.out, "lens=", len1, len2);
            return;
        }
    }
}

function deleteDisabledLinkMutate(chromesome) {

    var candidates = [];
    for (var i = 0; i < chromesome.genes.length; i++) {
        var g = chromesome.genes[i];
        if (g.enable === false)
            candidates.push(g);

    }

    if (candidates.length === 0)
        return;

    var i = getRandomIntInclusive(0, candidates.length - 1);
    console.log("link found..");
    deleteOneLink(chromesome, candidates[i]);
    generateNN(chromesome);


}

var NEATmutate = function (chromesome) {
    if (Math.random() < chromesome.genome.properties["PointMutationChance"]) {
        pointMutate(chromesome);
    }
    var genome = chromesome.genome;

    var p = genome.properties["LinkMutationChance"];
    while (p > 0) {
        if (Math.random() < p) {
            linkMutation(chromesome, false);
        }
        p -= 1;
    }

    var p = genome.properties["BiasMutationChance"];
    while (p > 0) {
        if (Math.random() < p) {
            linkMutation(chromesome, true);
        }
        p -= 1;
    }

    var p = genome.properties["NodeMutationChance"];
    while (p > 0) {
        if (Math.random() < p) {
            nodeMutation(chromesome);
        }
        p -= 1;
    }

    var p = genome.properties["EnableMutationChance"];
    while (p > 0) {
        if (Math.random() < p) {
            enableDisableMutate(chromesome, true);
        }
        p -= 1;
    }
    var p = genome.properties["DisableMutationChance"];
    while (p > 0) {
        if (Math.random() < p) {
            enableDisableMutate(chromesome, false);
        }
        p -= 1;
    }
    return chromesome;


};

var NEAToperator = new chromesomeOperator(NEATcrossover, NEATmutate);
var NEATdeltaWeight = 1;


function NEATdeltaFunction(chromesome1, chromesome2) {
    var DeltaDisjoint = 0.4;
    var DeltaWeights = 0.4;
    var dd = DeltaDisjoint * disjoint(chromesome1.genes, chromesome2.genes);
    var dw = DeltaWeights * weights(chromesome1.genes, chromesome2.genes);
    //console.log(dd, dw);
    return dd + dw;
}

function baseNEATChromesome(genome, initInputs, initOutputs) {
    genome.properties["inputs"] = initInputs;
    genome.properties["outputs"] = initOutputs;
    genome.properties["hiddenNeurons"] = 0;

    genome.properties["PertubStep"] = PerturbStep;
    genome.properties["PertubChance"] = PerturbChance;

    genome.properties["PointMutationChance"] = PointMutationChance;
    genome.properties["LinkMutationChance"] = LinkMutationChance;
    genome.properties["NodeMutationChance"] = NodeMutationChance;
    genome.properties["BiasMutationChance"] = BiasMutationChance;
    genome.properties["DisableMutationChance"] = DisableMutationChance;
    genome.properties["EnableMutationChance"] = EnableMutationChance;



    //var inputs = genome.properties["inputs"];
    //var outputs = genome.properties["outputs"];

    var result = new Chromesome(
            genome,
            "NEAT",
            generateNN,
            NEAToperator,
            NEATdeltaFunction,
            NEATdeltaWeight
            );

    return result;
}
;


function totalComplexity(pool) {
    //console.log(pool)
    var sum = 0;
    var population = 0;
    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        for (var j in species.genomes) {
            if (species.genomes.hasOwnProperty(j)) {
                population++;

                if (typeof species.genomes[j].properties["NeuroNetwork"] === "undefined")
                    generateNN(species.genomes[j].chromesomes["NEAT"]);

                var network = species.genomes[j].properties["NeuroNetwork"];
                //console.log(network);
                for (var ni in network.neurons) {
                    if (network.neurons.hasOwnProperty(ni)) {
                        sum++;
                        for (var l = 0; l < network.neurons[ni].incomingLinks.length; l++) {
                            sum++;
                        }
                    }
                }
            }
        }
    }

    var result = sum / population;
    //console.log("complexity", sum, population, result);
    return result;
}

function basicNEATGenome(pool, inputn, outputn) {
    var basicGenome = new Genome(pool);
    var chromesome = baseNEATChromesome(basicGenome, inputn, outputn);
    basicGenome.chromesomes["NEAT"] = chromesome;
    basicGenome.properties["NeuroNetwork"] = generateNN(chromesome);
    //console.log(chromesome);
    chromesome.chromesomeOperator.mutateOperator(chromesome);
    return basicGenome;
}

function initNEATPool(inputn, outputn, population) {
    var init_genomes = [];
    var pool = new Pool(population, true);

    for (var i = 0; i < population / 2; i++) {
        var basicGenome = new Genome(pool);
        var chromesome = baseNEATChromesome(basicGenome, inputn, outputn);
        basicGenome.chromesomes["NEAT"] = chromesome;
        basicGenome.properties["NeuroNetwork"] = generateNN(chromesome);
        chromesome.chromesomeOperator.mutateOperator(chromesome);
        addToSpecies(basicGenome, pool);
        init_genomes.push(basicGenome);
    }

    pool.complexity = totalComplexity(pool);
    return {pool: pool, genomes: init_genomes};
}