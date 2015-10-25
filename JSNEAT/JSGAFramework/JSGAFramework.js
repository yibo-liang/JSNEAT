//this is a framework for genetic algorithm to support not only NEAT, but also any other kind of genome structure

var ageSignificance = 1.5;

var speciationThreshold = 1;


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

    this.staleness = 0;
    this.complexity = 0;
    this.minComplexity = 0;

    this.selectionOperator =
            typeof selectionOperator === "undefined"
            ? tourment : selectionOperator;

    this.CrossoverChance = 0.75;
    this.AdoptionChanche = 0.1;
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
    this.chromesomes = {};
    this.fitness = 0;
    this.adJustedFitness = 0;
};
Genome.prototype.copy = function () {
    var copy = new Genome();
    copy.pool = this.pool;
    for (var i in this.chromesomes) {
        if (this.chromesomes.hasOwnProperty(i))
            copy.chromesomes[i] = this.chromesomes[i].copy();
    }

    for (var i in this.properties) {
        if (this.properties.hasOwnProperty(i)) {
            if (typeof this.properties[i].copy === "undefined") {
                copy.properties[i] = this.properties[i];
            } else {
                copy.properties[i] = this.properties[i].copy();

            }
        }
    }

    return copy;
};
var chromesomeOperator = function (crossoverOperator, mutationOperator) {
    this.crossoverOperator = crossoverOperator;
    this.mutateOperator = mutationOperator;
};

//evaluator : evaluate this single chromesome to the phenotype corresponding to it
var Chromesome = function (genome, name, evaluator, chromesomeOperator, deltaFunction, deltaWeight) {
    this.genome = genome;
    this.genes = [];
    this.name = name;
    this.chromesomeOperator = chromesomeOperator;

    this.evaluator = evaluator;
    this.deltaFunction = deltaFunction;
    this.deltaWeight = deltaWeight;
};
Chromesome.prototype.copy = function () {
    var copy = new Chromesome();
    for (var i = 0; i < this.genes.length; i++) {
        copy.genes.push(this.genes[i].copy());
    }
    copy.chromesomeOperator = this.chromesomeOperator;
    copy.evaluator = this.evaluator;

    copy.deltaFunction = this.deltaFunction;
    copy.deltaWeight = this.deltaWeight;
    copy.genome = this.genome;
    copy.name = this.name;

    //console.log("copy:", copy, "this", this);
    return copy;
};



function isSameSpecies(pool, genome1, genome2, threshold) {

    //if not using speciation, then always return true
    if (pool.isSpeciating === false) {
        return true;
    }

    //if has different number of chromesomes, then false


    //if chromesomes have different deltaFunction(comparing function)

    for (var i in genome1.chromesomes) {
        if (genome1.chromesomes.hasOwnProperty(i)) {
            // console.log(i, genome1, genome2);
            if (typeof genome1.chromesomes[i] === "undefined" ||
                    typeof genome2.chromesomes[i] === "undefined" ||
                    genome1.chromesomes[i].deltaFunction !== genome2.chromesomes[i].deltaFunction) {
                return false;
            }
        }
    }
    var sumdf = 0;
    for (var i in genome1.chromesomes) {
        if (genome1.chromesomes.hasOwnProperty(i)) {
            var c1 = genome1.chromesomes[i];
            var c2 = genome2.chromesomes[i];
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
    if (species.age <= 10) {
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

var StaleSpecies = 15;
function removeStaleSpecies(pool) {
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

function removeWeakSpecies(pool) {
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
    var newChromesomes = {};
    for (var chromesome in g1.chromesomes) {
        if (g1.chromesomes.hasOwnProperty(chromesome) &&
                g2.chromesomes.hasOwnProperty(chromesome)) {
            var c1 = g1.chromesomes[chromesome];
            var c2 = g2.chromesomes[chromesome];

            var crossover1 = c1.chromesomeOperator.crossoverOperator;
            var crossover2 = c2.chromesomeOperator.crossoverOperator;

            if (crossover1 === crossover2) {
                var result = crossover1(c1, c2);
                result.genome = genome;
                newChromesomes[chromesome] = result;
            } else {
                return null;
            }
        }
    }
    return newChromesomes;

}

function mutateChromesomes(g) {
    var newChromesomes = {};
    for (var chromesome in g.chromesomes) {
        if (g.chromesomes.hasOwnProperty(chromesome)) {
            var c = g.chromesomes[chromesome];
            //console.log(c);
            var mutate = c.chromesomeOperator.mutateOperator;
            // console.log(g, c, mutate);
            var result = mutate(c);
            result.genome = g;
            newChromesomes[chromesome] = result;
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

    if (Math.random() < pool.CrossoverChance) {
        if (Math.random() < pool.AdoptionChanche) {
            g1 = select(species.genomes);
            g2 = select(potentialAdoptSpecies.genomes);
            // console.log(g1);

            child = g1.copy();
            child.fitness = 0;
            child.adJustedFitness = 0;
            child.chromesomes = crossoverChromesomes(g1, g2, child);


        } else {
            g1 = select(species.genomes);
            g2 = select(species.genomes);
            child = g1.copy();
            child.fitness = 0;
            child.adJustedFitness = 0;
            child.chromesomes = crossoverChromesomes(g1, g2, child);

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
    child.chromesomes = mutateChromesomes(child);
    //console.log("newchild .c =", child.chromesomes);
    return child;
}

function nextGeneration(pool) {
    //pool.staleness++;
    cullSpecies(pool, false);
    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        calculateAverageFitness(species);
    }
    if (pool.species.length >= 4) {
        removeStaleSpecies(pool);
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