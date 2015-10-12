

var inputsize;
var inputs;
var outputs;
var pool;
var innovation = 0;
var Population;
var DeltaDisjoint = 1.5;
var DeltaWeights = 0.4;
var DeltaThreshold = 1.3;
var StaleSpecies = 30;
var dropOffAge = 30;
var ageSignificance = 1.2;

var MutateConnectionsChance = 0.25;
var PerturbChance = 0.90;
var CrossoverChance = 0.75;

var AdoptionChanche = 0.01;
var LinkMutationChance = 2;
var NodeMutationChance = 0.5;
var BiasMutationChance = 0.30;
var StepSize = 0.1;
var DisableMutationChance = 0.4;
var EnableMutationChance = 0.2;
var MaxNodes = 100000;
function sigmoid(x) {
    return 2 / (1 + Math.exp(-4.9 * x)) - 1;
}

function newInnovation() {
    pool.innovation = pool.innovation + 1;
    return pool.innovation;
}

var Pool = function () {
    this.species = [];
    this.generation = 0;
    this.speciesCount = 0;
    this.innovation = outputs;
    this.currentSpecies = 1;
    this.currentGenome = 1;
    this.currentFrame = 0;
    this.maxFitness = 0;
    this.isPruning = false;
};
//pool = new Pool();

var Species = function () {
    this.id = pool.speciesCount++;
    this.topFitness = 0;
    this.topAverageFitness = 0;
    this.staleness = 0;
    this.genomes = [];
    this.averageFitness = 0;
    this.age = 0;
};
var Genome = function () {

    this.genes = [];
    this.fitness = 0;
    this.adJustedFitness = 0;
    this.network = [];
    this.maxNeuron = 0;
    this.globalRank = 0;
    this.mutationRates = {};
    this.mutationRates["connections"] = MutateConnectionsChance;
    this.mutationRates["link"] = LinkMutationChance;
    this.mutationRates["bias"] = BiasMutationChance;
    this.mutationRates["node"] = NodeMutationChance;
    this.mutationRates["enable"] = EnableMutationChance;
    this.mutationRates["disable"] = DisableMutationChance;
    this.mutationRates["adopt"] = AdoptionChanche;
    this.mutationRates["step"] = StepSize;
};
Genome.prototype.copy = function () {
    var copy = new Genome();
    for (var i = 0; i < this.genes.length; i++) {
        copy.genes.push(this.genes[i]);
    }

    copy.maxNeuron = this.maxNeuron;
    copy.mutationRates["connections"] = this.mutationRates["connections"];
    copy.mutationRates["link"] = this.mutationRates["link"];
    copy.mutationRates["bias"] = this.mutationRates["bias"];
    copy.mutationRates["node"] = this.mutationRates["node"];
    copy.mutationRates["enable"] = this.mutationRates["enable"];
    copy.mutationRates["disable"] = this.mutationRates["disable"];
    copy.mutationRates["adopt"] = this.mutationRates["AdoptionChanche"];

    copy.mutationRates["step"] = this.mutationRates["step"];
    return copy;
};
var Gene = function () {
    this.in = 0;
    this.out = 0;
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
function basicGenome() {
    var result = new Genome();
    var genes = result.genes;
    var basic_innovatation = 0;
    result.maxNeuron = inputs;
    for (var i = 0; i < inputs; i++) {
        for (var j = 0; j < outputs; j++) {
            var newGene = new Gene();
            newGene.in = i;
            newGene.out = MaxNodes + j;
            newGene.enable = true;
            newGene.innovation = basic_innovatation++;
            newGene.weight = 0; //Math.random() * 4 - 2;
            //genes.push(newGene);
        }
    }

    return result;
}
;

var neuron_count = 0;
var Neuron = function () {
    this.incomingLinks = [];
    this.v = 0;
    this.id = neuron_count++;
};

function generateNeuroNetwork(genome) {
    //console.log("start generate NN - - - - ");
    var network = {neurons: {}};
    for (var i = 0; i < inputs; i++) {
        network.neurons[i] = new Neuron();
    }

    for (var i = 0; i < outputs; i++) {
        network.neurons[i + MaxNodes] = new Neuron();
    }

    genome.genes.sort(function (a, b) {
        if (a.out < b.out) {
            return -1;
        } else {
            return 1;
        }
    });
    //console.log(genome.genes);
    for (var i = 0; i < genome.genes.length; i++) {
        var gene = genome.genes[i];
        if (gene.enable === true) {
            if (typeof network.neurons[gene.out] === "undefined") {
                network.neurons[gene.out] = new Neuron();
            }
            var neuron = network.neurons[gene.out];
            //console.log(network);
            neuron.incomingLinks.push(gene);
            if (typeof network.neurons[gene.in] === "undefined") {
                network.neurons[gene.in] = new Neuron();
            }
        }

    }
    genome.network = network;
    //console.log("generate network", network);
    //return network;
}

function evaluateNeuroNetwork(network, inputArray) {
    //console.log("start eval nn -   ------")
    var inputVec = [];
    for (var i = 0; i < inputArray.length; i++) {
        inputVec.push(inputArray[i]);
    }
    inputVec.push(1);

    if (inputVec.length !== inputs) {
        console.log("ERROR: INPUT NUMBER IS DIFFERENT FROM NETWORK", inputVec);
        return;
    }

    for (var i = 0; i < inputs; i++) {
        network.neurons[i].v = inputVec[i];

    }

    for (var i in network.neurons) {

        if (network.neurons.hasOwnProperty(i)) {

        }
    }

    //console.log(inputVec);
    //console.log("netwrok=", network);
    for (var i in network.neurons) {
        if (network.neurons.hasOwnProperty(i)) {
            var neuron = network.neurons[i];
            //console.log(i, neuron);
            var sum = 0;
            for (var j = 0; j < neuron.incomingLinks.length; j++) {
                var incoming = neuron.incomingLinks[j];
                var other = network.neurons[incoming.in];
                //console.log(incoming, other);
                sum += incoming.weight * other.v;
                //console.log("sum=",sum);
            }
            if (neuron.incomingLinks.length > 0) {
                neuron.v = sigmoid(sum);
                //console.log("v=",neuron.v);
            }

            //console.log("i=", i, "  v=", network.neurons[i].v);
        }
    }

    var outputArray = [];
    for (var i = 0; i < outputs; i++) {
        //console.log("i=", i + MaxNodes, "v=", network.neurons[i + MaxNodes].v)
        if (network.neurons[i + MaxNodes].v > 0) {
            outputArray.push(true);
        } else {
            outputArray.push(false);
        }
    }
    return outputArray;
}

function crossover(g1, g2) {
    if (g2.fitness > g1.fitness) {
        var temp = g2;
        g2 = g1;
        g1 = temp;
    }



    var child = new Genome();
    var innovations2 = {};
    for (var i = 0; i < g2.genes.length; i++) {
        var gene = g2.genes[i];
        innovations2[gene.innovation] = gene;
    }


    for (var i = 0; i < g1.genes.length; i++) {
        var gene1 = g1.genes[i];
        var gene2 = innovations2[gene1.innovation];
        if (typeof gene2 !== "undefined" && gene2.enable === true) {
            child.genes.push(gene2.copy());
        } else {
            child.genes.push(gene1.copy());
        }
    }

    child.maxNeuron = Math.max(g1.maxNeuron, g2.maxNeuron);
    for (var mutation in g1.mutationRates) {
        child.mutationRates[mutation] = g1.mutationRates[mutation];
    }
    return child;
}

function getRandomIntInclusive(min, max) {
    return Math.floor(Math.random() * (max - min + 1)) + min;
}

function randomNeuron(genes, NoInputs) {
    var neurons = {};
    if (NoInputs !== true) {
        for (var i = 0; i < inputs; i++) {
            neurons[i] = true;
        }
    }

    for (var i = 0; i < outputs; i++) {
        neurons[MaxNodes + i] = true;
    }

    for (var i = 0; i < genes.length; i++) {
        if (NoInputs === false || genes[i].in + 1 > inputs) {
            neurons[genes[i].in] = true;
        }
        if (NoInputs === false || genes[i].out + 1 > inputs) {
            neurons[genes[i].out] = true;
        }
    }

    var count = 0;
    for (var n in neurons) {
        count++;
    }
    var n = getRandomIntInclusive(1, count);
    for (var k in neurons) {
        if (neurons.hasOwnProperty(k)) {
//alert("Key is " + k + ", value is" + target[k]);
            n--;
            if (n === 0) {
                return +k;
            }
        }
    }
    return 0;
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

function pointMutate(genome) {
    //console.log("pointMutation");
    var step = genome.mutationRates["step"];
    for (var i = 0; i < genome.genes.length; i++) {
        var gene = genome.genes[i];
        if (Math.random() < PerturbChance) {
            gene.weight = gene.weight + Math.random() * step * 2 - step;
        } else {
            gene.weight = Math.random() * step * 4 - 2;
        }
    }

}

function linkMutation(genome, forceBias) {
    //console.log("linkMutation");
    var neuron1 = +randomNeuron(genome, false);
    var neuron2 = +randomNeuron(genome, true);

    if (neuron1 <= inputs - 1 && neuron2 <= inputs - 1 || +neuron1 === +neuron2) {
// both inputs neuron
        return;
    }
    if (neuron2 <= inputs - 1) {
//swap 1 and 2
        var temp = neuron1;
        neuron1 = neuron2;
        neuron2 = temp;
    }

    var newLink = new Gene();
    newLink.in = +neuron1;
    newLink.out = +neuron2;
    if (forceBias === true) {
        newLink.in = inputs;
    }
    if (containsLink(genome.genes, newLink) === true) {
        return;
    }
    newLink.innovation = newInnovation();
    newLink.weight = Math.random() * 4 - 2;
    genome.genes.push(newLink);
}

function nodeMutation(genome) {
    if (genome.genes.length === 0) {
        return;
    }

    var gene = genome.genes[getRandomIntInclusive(0, genome.genes.length - 1)];
    if (gene.enable === false) {
        return;
    }

    genome.maxNeuron += 1;

    gene.enable = false;
    var gene1 = gene.copy();
    gene1.out = +genome.maxNeuron;
    gene1.weight = 1.0;
    gene1.innovation = newInnovation();
    gene1.enable = true;
    genome.genes.push(gene1);
    var gene2 = gene.copy();
    gene2.in = +genome.maxNeuron;
    gene2.innovation = newInnovation();
    gene2.enable = true;
    genome.genes.push(gene2);
}


function enableDisableMutate(genome, enable) {
    var candidates = [];
    for (var gene in genome.genes) {
        if (gene.enable === false) {
            candidates.push(gene);
        }
    }

    if (candidates.length === 0) {
        return;
    }

    var gene = candidates[getRandomIntInclusive(0, candidates.length - 1)];
    gene.enable = true;
}


function mutate(genome) {
    for (var mutation in genome.mutationRates) {
        if (genome.mutationRates.hasOwnProperty(mutation))
            if (Math.random() > 0.5) {
                genome.mutationRates[mutation] *= 0.95;
            } else {
                genome.mutationRates[mutation] *= 1.05263157895;
            }
    }

    if (Math.random() < genome.mutationRates["connections"]) {
        pointMutate(genome);
    }

    var p = genome.mutationRates["link"];
    while (p > 0) {
        if (Math.random() < p) {
            linkMutation(genome, false);
        }
        p -= 1;
    }

    var p = genome.mutationRates["bias"];
    while (p > 0) {
        if (Math.random() < p) {
            linkMutation(genome, true);
        }
        p -= 1;
    }


    var p = genome.mutationRates["node"];
    while (p > 0 && pool.isPruning === false) {
        if (Math.random() < p) {
            nodeMutation(genome);
        }
        p -= 1;
    }

    var p = genome.mutationRates["enable"];
    while (p > 0 && pool.isPruning === false) {
        if (Math.random() < p) {
            enableDisableMutate(genome, true);
        }
        p -= 1;
    }
    var p = genome.mutationRates["disable"];
    while (p > 0 && pool.isPruning === false) {
        if (Math.random() < p) {
            enableDisableMutate(genome, false);
        }
        p -= 1;
    }
}

function disjoint(genes1, genes2) {
    var l1 = {};
    var l2 = {};

    for (var i = 0; i < genes1.length; i++) {
        var gene = genes1[i];
        l1[gene.innovation] = true;
    }

    for (var i = 0; i < genes2.length; i++) {
        var gene = genes2[i];
        l2[gene.innovation] = true;
    }
    //console.log(l1, l2);
    var disjointGenes = 0;
    for (var i = 0; i < genes1.length; i++) {

        var gene = genes1[i];
        //console.log(i, gene, genes1);
        if (typeof l2[gene.innovation] === "undefined" || l2[gene.innovation] !== true) {
            disjointGenes = disjointGenes + 1;
        }

    }
    for (var i = 0; i < genes2.length; i++) {

        var gene = genes2[i];
        if (typeof l1[gene.innovation] === "undefined" || l1[gene.innovation] !== true) {
            disjointGenes = disjointGenes + 1;
        }

    }

    var n = Math.max(genes1.length, genes2.length);
    //console.log(disjointGenes,n);
    return disjointGenes / n;
}


function weights(genes1, genes2) {
    //console.log("weights:",genes1, genes2);
    var i2 = {};
    for (var i = 1; i < genes2.length; i++) {
        var gene = genes2[i];
        i2[+gene.innovation] = gene;
    }
    //console.log(i2);
    var sum = 0;
    var coincident = 0;
    for (var i = 1; i < genes1.length; i++) {
        var gene = genes1[i];
        if (typeof i2[+gene.innovation] !== "undefined") {
            var gene2 = i2[+gene.innovation];
            sum = sum + Math.abs(gene.weight - gene2.weight);
            coincident = coincident + 1;
        }
    }
    //console.log("sum, coinci", sum, coincident, sum/coincident);
    return sum / coincident;
}


function sameSpecies(genome1, genome2) {
    var dd = DeltaDisjoint * disjoint(genome1.genes, genome2.genes);
    var dw = DeltaWeights * weights(genome1.genes, genome2.genes);
    return dd + dw < DeltaThreshold;
}


function rankGlobally() {
    var global = [];
    for (var s = 0; s < pool.species.length; s++) {
        var species = pool.species[s];
        for (var g = 1; g < species.genomes.length; g++) {
            global.push(species.genomes[g]);
        }
    }
    global.sort(function (a, b) {
        return (b.fitness - a.fitness);
    });
    for (var g = 0; g < global.length; g++) {
        global[g].globalRank = g;
    }
}

function calculateAverageFitness(species) {
    var total = 0;
    var a = 1;
    if (species.age <= 10) {
        a = ageSignificance;
    } else {
        a = 1;
    }

    for (var g = 0; g < species.genomes.length; g++) {
        var genome = species.genomes[g];
        genome.adJustedFitness = a * genome.fitness / Math.log(10 * species.genomes.length);
        //console.log("adjfit:",genome.adJustedFitness);
        if (genome.fitness > pool.maxFitness) {
            pool.maxFitness = genome.fitness;
        }
        total = total + genome.adJustedFitness;//genome.globalRank;
    }


    species.averageFitness = total / species.genomes.length;
}


function totalAverageFitness() {
    var total = 0;
    for (var s = 0; s < pool.species.length; s++) {
        var species = pool.species[s];
        total = total + species.averageFitness;
    }

    return total;
}

function cullSpecies(cutToOne) {
    for (var s = 0; s < pool.species.length; s++) {
        var species = pool.species[s];
        species.genomes.sort(function (a, b) {
            return (b.fitness - a.fitness);
        });

        var remaining = Math.ceil(species.genomes.length * 0.5);
        if (cutToOne === true) {
            remaining = 1;
        }
        //console.log("before cull", species.genomes.length);
        //for (var i = 0; i < species.genomes.length; i++) {
        //    console.log("before g fit=", species.genomes[i].fitness);
        //}
        species.genomes.splice(remaining);

        //console.log("after cull", species.genomes.length);
        //for (var i = 0; i < species.genomes.length; i++) {
        //    console.log("remain g fit=", species.genomes[i].fitness);
        //}
    }
}


function tourment(species) {
    var n = 2;
    var t = [];
    //console.log("test s=", species);
    for (var i = 0; i < n; i++) {
        var g = species.genomes[getRandomIntInclusive(0, species.genomes.length - 1)];
        //console.log(g);
        t.push(g);
    }

    var maxg;
    var maxfit = 0;
    for (var i = 0; i < n; i++) {
        if (t[i].fitness > maxfit) {
            maxfit = t[i].fitness;
            maxg = t[i];
        }
    }
    return maxg;

}

function breedChild(species, potentialAdoptSpecies) {
    var child;
    var g1;
    var g2;
    var g;
    //console.log("breeding:", species);
    if (Math.random() < CrossoverChance && pool.isPruning === false) {
        if (Math.random() < AdoptionChanche) {
            g1 = tourment(species);//species.genomes[getRandomIntInclusive(0, species.genomes.length - 1)];
            g2 = tourment(potentialAdoptSpecies);//species.genomes[getRandomIntInclusive(0, species.genomes.length - 1)];
            //console.log("breed from g1 fit=", g1.fitness, " and g2 fit=", g2.fitness);
            child = crossover(g1, g2);
        } else {
            g1 = tourment(species);//species.genomes[getRandomIntInclusive(0, species.genomes.length - 1)];
            g2 = tourment(species);//species.genomes[getRandomIntInclusive(0, species.genomes.length - 1)];
            //console.log("breed from g1 fit=", g1.fitness, " and g2 fit=", g2.fitness);
            child = crossover(g1, g2);
        }

    } else {
        if (Math.random() < AdoptionChanche) {
            g = tourment(potentialAdoptSpecies);//species.genomes[getRandomIntInclusive(0, species.genomes.length - 1)];
            //console.log("breed from g fit=", g.fitness);
            child = g.copy();
        } else {
            g = tourment(species);//species.genomes[getRandomIntInclusive(0, species.genomes.length - 1)];
            //console.log("breed from g fit=", g.fitness);
            child = g.copy();
        }
    }

    mutate(child);

    return child;
}

function adjustAllFitness() {
    for (var s = 0; s < pool.species.length; s++) {
        var species = pool.species[s];

        var numOfParents;
        var count;
        var ageDebt = species.staleness - dropOffAge;
        if (ageDebt === 0)
            ageDebt = 1;
        for (var i = 0; i < species.genomes.length; i++) {
            var genome = species.genomes[i];
            if (ageDebt >= 1) {
                genome.adjustedFitness = genome.fitness * 0.01;
            }

            if (genome.age <= 10) {
                genome.adjustedFitness = genome.fitness * ageSignificance;
            }

            genome.adjustedFitness = genome.adjustedFitness / species.genomes.length;
        }
    }
}

function removeStaleSpecies() {
    var survived = [];
    for (var s = 0; s < pool.species.length; s++) {
        var species = pool.species[s];

        species.genomes.sort(function (a, b) {
            return (b.fitness - a.fitness);
        });
        if (species.genomes[0].fitness > species.topFitness) {
            species.topFitness = species.genomes[0].fitness;
            species.staleness = 0;
        } else if (species.averageFitness > species.topAverageFitness) {
            species.topAverageFitness = species.averageFitness;
            species.staleness = 0;
        } else {
            species.staleness++;
        }
        if (species.staleness < StaleSpecies || species.topFitness >= pool.maxFitness) {// 
            survived.push(species);

        }
    }

    pool.species = survived;
}

function removeWeakSpecies() {
    var survided = [];

    var sum = totalAverageFitness();
    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        var br = species.averageFitness / sum * Population;
        var breed = Math.floor(br + 0.01);
        console.log("removeWeakSp, breed=", breed, ", br=", br);
        if (breed >= 1) {
            survided.push(species);
        }
    }
    //console.log("surv num=" + survided.length, ", origin sp=", pool.species);
    pool.species = survided;
}

function addToSpecies(child) {
    var foundSpecies = false;
    //console.log(child);
    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        if (foundSpecies === false && sameSpecies(child, species.genomes[0])) {
            foundSpecies = true;
            species.genomes.push(child);
            break;
        }
    }

    if (foundSpecies === false) {
        var childSpecies = new Species();
        childSpecies.genomes.push(child);
        //console.log("pushed", childSpecies);
        pool.species.push(childSpecies);
    }
}

function newGeneration() {
    console.log();
    cullSpecies(false);
    //rankGlobally();
    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        calculateAverageFitness(species);
    }
    if (pool.species.length >= 4)
        removeStaleSpecies();

    //rankGlobally();
    //rankGlobally();
    //adadjustAllFitness();
    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        calculateAverageFitness(species);
    }

    removeWeakSpecies();

    var survived = pool.species;
    var str = "surv=[";
    for (var i = 0; i < survived.length; i++) {
        str += "(id=" + survived[i].id
                + ", n=" + survived[i].genomes.length
                + ", af=" + survived[i].averageFitness
                + ", taf" + survived[i].topAverageFitness
                + ", tf=" + survived[i].topFitness + ")";
        if (i < survived.length - 1) {
            str += ", \n";
        }
    }
    console.log(str + "], n=" + survived.length);

    var sum = totalAverageFitness();
    var children = [];

    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        //console.log("t2 s=", species);
        species.age++;
        var breed = Math.floor(species.averageFitness / (sum) * (Population)) - species.genomes.length;
        //console.log(i, " breed=" + breed);
        for (var j = 0; j < breed; j++) {
            var randomSpecies = pool.species[getRandomIntInclusive(0, pool.species.length - 1)];

            children.push(breedChild(species, randomSpecies));
        }

    }
    console.log("newgen child len=", children.length, "pool spe len=", pool.species.length);
    cullSpecies(true);
    //each species now only has 1 child. so the number of children + number of species < population
    //number of species= number of remained genome of last generation
    console.log("sp n=", pool.species.length);
    while (children.length + pool.species.length < Population) {
        var ri = getRandomIntInclusive(0, pool.species.length - 1);
        var species = pool.species[ri];

        var randomSpecies = pool.species[getRandomIntInclusive(0, pool.species.length - 1)];
        //console.log("ri=", ri, "species=", species)
        var child = breedChild(species, randomSpecies);
        children.push(child);
        //addToSpecies(child);
    }
    //console.log("child len!!!!!=", children.length);
    for (var i = 0; i < children.length; i++) {
        var child = children[i];
        addToSpecies(child);
    }

    pool.generation++;


}

function initialisePool(inputn, outputn, population) {
    var init_genomes = [];
    inputsize = inputn;
    outputs = outputn;
    inputs = inputsize + 1;
    innovation = 0;
    Population = population;
    pool = new Pool();

    for (var i = 0; i < Population; i++) {
        var basic = basicGenome();
        mutate(basic);
        var newg = basic.copy();
        addToSpecies(newg);
        init_genomes.push(newg);
    }
    return init_genomes;
}