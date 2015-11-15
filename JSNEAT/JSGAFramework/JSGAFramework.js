//this is a framework for genetic algorithm to support not only NEAT, but also any other kind of genome structure

var ageSignificance = 1.2;

var speciationThreshold = 1;

var maxStaleness = 15;

var minSpecies = 18;
//helper function, 
function RandomIntInclusive(min, max) {
    return Math.floor(Math.random() * (max - min + 1)) + min;
}

//tourment selection implementation,
//given list of genomes, returns one genome.
function tourment(genomes) {
    var n = 2;
    var t = [];
    for (var i = 0; i < n; i++) {
        var g = genomes[RandomIntInclusive(0, genomes.length - 1)];
        t.push(g);
    }
    var maxg = genomes[0];
    var maxfit = 0;
    for (var i = 0; i < n; i++) {
        if (t[i].fitness > maxfit) {
            maxfit = t[i].fitness;
            maxg = t[i];
        }
    }
    if (typeof maxg === "undefined") {
        console.log(genomes)
    }
    return maxg;
}




// 
var Pool = function (population, isSpeciating, selectionOperator) {
    this.population = population;
    this.species = [];
    this.generation = 0;
    this.speciesCount = 0;
    this.maxFitness = 0;
    this.innovation = 0;
    this.isPruning = false;

    //wether use speciation concept described in NEAT paper
    this.isSpeciating = typeof isSpeciating === "undefined" ? false : isSpeciating;


    this.crossoverEnabled = true;
    this.mutateEnabled = true;

    this.staleness = 0;
    this.complexity = 0;
    this.minComplexity = 0;

    this.selectionOperator =
            typeof selectionOperator === "undefined"
            ? tourment : selectionOperator;

    this.CrossoverChance = 0.75;
    this.AdoptionChanche = 0.01;
};

var Species = function (pool) {
    this.pool = pool;
    // console.log(pool)
    this.id = this.pool.speciesCount++;
    this.topFitness = 0;
    this.topAverageFitness = 0;
    this.staleness = 0;
    this.genomes = [];
    this.averageFitness = 0;
    this.age = 0;
};

var Genome = function (pool) {
    this.pool = pool;
    this.properties = {};
    this.chromosomes = {};
    this.fitness = 0;
    this.adJustedFitness = 0;
    this.log = "LOG: \n";
};

var copyFunc = function (old) {
    return JSON.parse(JSON.stringify(old));

};

Genome.prototype.copy = function () {
    var copy = new Genome();
    copy.pool = this.pool;
    for (var i in this.chromosomes) {
        if (this.chromosomes.hasOwnProperty(i))
            copy.chromosomes[i] = this.chromosomes[i].copy();
    }
    for (var i in this.properties) {
        if (this.properties.hasOwnProperty(i)) {
            if (typeof this.properties[i].copy === "undefined") {
                copy.properties[i] = copyFunc(this.properties[i]);
            } else {
                copy.properties[i] = this.properties[i].copy();

            }
        }
    }

    return copy;
};
var chromosomeOperator = function (crossoverOperator, mutationOperator) {
    this.crossoverOperator = crossoverOperator;
    this.mutateOperator = mutationOperator;
};

var csid = 0;
//evaluator : evaluate this single chromosome to the phenotype corresponding to it
var Chromesome = function (genome, name, evaluator, chromosomeOperator, deltaFunction, deltaWeight) {
    this.genome = genome;
    this.genes = [];
    this.id = csid++;
    this.name = name;
    this.chromosomeOperator = chromosomeOperator;

    this.evaluator = evaluator;
    this.deltaFunction = deltaFunction;
    this.deltaWeight = deltaWeight;
};
Chromesome.prototype.copy = function () {
    var copy = new Chromesome();
    for (var i = 0; i < this.genes.length; i++) {
        copy.genes.push(this.genes[i].copy());
    }
    copy.chromosomeOperator = this.chromosomeOperator;
    copy.evaluator = this.evaluator;

    copy.deltaFunction = this.deltaFunction;
    copy.deltaWeight = this.deltaWeight;
    copy.genome = this.genome;
    copy.id = csid++;
    copy.genome.log = this.genome.log + "chromosome id=" + copy.id + " copied once from " + this.id + "\n";
    copy.name = this.name;

    //console.log("copy:", copy, "this", this);
    return copy;
};



function isSameSpecies(pool, genome1, genome2, threshold) {

    //if not using speciation, then always return true
    if (pool.isSpeciating === false) {
        return true;
    }

    //if has different number of chromosomes, then false


    //if chromosomes have different deltaFunction(comparing function)

    for (var i in genome1.chromosomes) {
        if (genome1.chromosomes.hasOwnProperty(i)) {
            // console.log(i, genome1, genome2);
            if (typeof genome1.chromosomes[i] === "undefined" ||
                    typeof genome2.chromosomes[i] === "undefined" ||
                    genome1.chromosomes[i].deltaFunction !== genome2.chromosomes[i].deltaFunction) {
                return false;
            }
        }
    }
    var sumdf = 0;
    for (var i in genome1.chromosomes) {
        if (genome1.chromosomes.hasOwnProperty(i)) {
            var c1 = genome1.chromosomes[i];
            var c2 = genome2.chromosomes[i];
            var f = c1.deltaFunction;
            // console.log("df=", f, "df(c1c2)=",f(c1, c2));
            var df = f(c1, c2) * (c1.deltaWeight + c2.deltaWeight) / 2;

            sumdf += df;
        }
    }
    //console.log(sumdf);
    return sumdf < threshold;

}


function calculateAverageFitness(species) {
    var total = 0;
    var a = 1;
    if (species.age <= 20) {
        a = ageSignificance;
    } else {
        a = 1;
    }

    for (var g = 0; g < species.genomes.length; g++) {
        var genome = species.genomes[g];
        genome.adJustedFitness = a * genome.fitness / Math.log(10 * species.genomes.length);
        //console.log("adjfit:",genome.adJustedFitness);
        if (genome.fitness > species.pool.maxFitness) {
            species.pool.maxFitness = genome.fitness;
            species.pool.staleness = 0;
        }
        total = total + genome.adJustedFitness;//genome.globalRank;
    }


    species.averageFitness = total / species.genomes.length;
}

function totalAverageFitness(pool) {
    var total = 0;
    for (var s = 0; s < pool.species.length; s++) {
        var species = pool.species[s];
        total = total + species.averageFitness;
    }

    return total;
}

function cullSpecies(pool, cutToOne) {
    //console.log(pool)
    for (var s = 0; s < pool.species.length; s++) {
        var species = pool.species[s];
        species.genomes.sort(function (a, b) {
            return (b.fitness - a.fitness);
        });

        var remaining = Math.ceil(species.genomes.length * 0.5);
        if (cutToOne === true) {
            remaining = 1;
        }
        species.genomes.splice(remaining);
    }
}

function removeStaleSpecies(pool) {

    if (pool.species.length <= minSpecies) {
        return;
    }
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
        if (species.staleness < maxStaleness || species.topFitness >= pool.maxFitness) {// 
            survived.push(species);

        }
    }

    pool.species = survived;
}

function removeWeakSpecies(pool) {
    if (pool.species.length <= minSpecies) {
        return;
    }
    var survided = [];

    var sum = totalAverageFitness(pool);
    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        var br = species.averageFitness / sum * pool.population;
        var breed = Math.floor(br + 0.01);
        //console.log("removeWeakSp, breed=", breed, ", br=", br);
        if (breed >= 1) {
            survided.push(species);
        }
    }
    //console.log("surv num=" + survided.length, ", origin sp=", pool.species);
    pool.species = survided;
}


function addToSpecies(child, pool) {
    var foundSpecies = false;
    //console.log(child);
    //console.log(pool)
    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        if (foundSpecies === false && isSameSpecies(pool, child, species.genomes[0], speciationThreshold)) {
            foundSpecies = true;
            species.genomes.push(child);
            break;
        }
    }

    if (foundSpecies === false) {
        var childSpecies = new Species(pool);
        childSpecies.genomes.push(child);
        //console.log("pushed", childSpecies);
        pool.species.push(childSpecies);
    }
}



function crossoverChromesomes(g1, g2, genome) {
    //console.log("crossover")
    var newChromesomes = {};

    for (var chromosome in g1.chromosomes) {
        if (g1.chromosomes.hasOwnProperty(chromosome) &&
                g2.chromosomes.hasOwnProperty(chromosome)) {
            var c1 = g1.chromosomes[chromosome];
            var c2 = g2.chromosomes[chromosome];

            var crossover1 = c1.chromosomeOperator.crossoverOperator;
            var crossover2 = c2.chromosomeOperator.crossoverOperator;

            if (crossover1 === crossover2) {
                var result = crossover1(c1, c2);
                genome.log += "created from crossover, p1=" + g1.chromosomes["NEAT"].id + ", p2=" + g2.chromosomes["NEAT"].id + "\n";

                result.genome = genome;
                newChromesomes[chromosome] = result;
            } else {
                return null;
            }
        }
    }
    return newChromesomes;

}

function mutateChromesomes(g) {
    var newChromesomes = {};
    for (var chromosome in g.chromosomes) {
        if (g.chromosomes.hasOwnProperty(chromosome)) {
            var c = g.chromosomes[chromosome];
            //console.log(c);
            var mutate = c.chromosomeOperator.mutateOperator;
            // console.log(g, c, mutate);
            var result = mutate(c);
            //result.genome = g;
            newChromesomes[chromosome] = result;
        }
    }

    return newChromesomes;
}


function breedChild(pool, species, potentialAdoptSpecies) {
    var child;
    var g1;
    var g2;
    var g;

    var select = pool.selectionOperator;
    if (Math.random() < pool.CrossoverChance && pool.crossoverEnabled) {
        if (Math.random() < pool.AdoptionChanche) {
            g1 = select(species.genomes);
            g2 = select(potentialAdoptSpecies.genomes);
            // console.log(g1);

            child = g1.copy();
            child.fitness = 0;
            child.adJustedFitness = 0;
            child.chromosomes = crossoverChromesomes(g1, g2, child);


        } else {
            g1 = select(species.genomes);
            g2 = select(species.genomes);
            child = g1.copy();
            child.fitness = 0;
            child.adJustedFitness = 0;
            child.chromosomes = crossoverChromesomes(g1, g2, child);

        }
    } else {
        if (Math.random() < pool.AdoptionChanche) {
            g = select(potentialAdoptSpecies.genomes);

            child = g.copy();
        } else {
            //console.log(select, pool);
            g = select(species.genomes);
            child = g.copy();
        }
    }
    //console.log(child);
    for (var i in child.chromosomes) {
        if (child.chromosomes.hasOwnProperty(i)) {
            child.chromosomes[i].genome = child;
        }
    }
    if (pool.mutateEnabled)
        child.chromosomes = mutateChromesomes(child);

    //console.log("newchild .c =", child.chromosomes);
    return child;
}

function nextGeneration(pool) {
    //pool.staleness++;
    cullSpecies(pool, false);
    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        calculateAverageFitness(species);
    }
    if (pool.species.length > minSpecies) {
        removeStaleSpecies(pool);
        for (var i = 0; i < pool.species.length; i++) {
            var species = pool.species[i];
            calculateAverageFitness(species);
        }
        removeWeakSpecies(pool);
        for (var i = 0; i < pool.species.length; i++) {
            var species = pool.species[i];
            calculateAverageFitness(species);
        }
    }

    var sum = totalAverageFitness(pool);
    var children = [];

    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        //console.log("t2 s=", species);
        species.age++;
        var breed = Math.floor(species.averageFitness / (sum) * (pool.population));// - species.genomes.length;
        //console.log(i, " breed=" + breed);
        for (var j = 0; j < breed; j++) {
            var randomSpecies = pool.species[RandomIntInclusive(0, pool.species.length - 1)];

            children.push(breedChild(pool, species, randomSpecies));
        }

    }
    cullSpecies(pool, true);

    while (children.length + pool.species.length < pool.population) {
        var ri = RandomIntInclusive(0, pool.species.length - 1);
        var species = pool.species[ri];

        var randomSpecies = pool.species[RandomIntInclusive(0, pool.species.length - 1)];
        //console.log("ri=", ri, "species=", species)
        var child = breedChild(pool, species, randomSpecies);
        children.push(child);
        //addToSpecies(child);
    }
    //console.log("child len!!!!!=", children.length);
    for (var i = 0; i < children.length; i++) {
        var child = children[i];
        addToSpecies(child, pool);
    }


    pool.generation++;
}