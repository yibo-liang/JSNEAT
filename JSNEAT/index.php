<!DOCTYPE html>

<html>
    <head>
        <meta charset="UTF-8">
        <title></title> 
        <script type="text/javascript" src="JSNEAT/JSNEAT.js"></script>

        <script type="text/javascript" src="JSNEAT/box2dweb.js"></script>

        <script type="text/javascript" src="JSNEAT/core.js"></script>
    </head>
    <body>
        <script>
            initialisePool(2, 1, 20);

            function xor(a, b) {
                var l1;
                var l2;
                if (a === 0) {
                    l1 = false;
                } else {
                    l1 = true;
                }

                if (b === 0) {
                    l2 = false;
                } else {
                    l2 = true;
                }
                return (l1 || l2) && !(l1 && l2);
            }
            var inputvec = [[0, 0], [0, 1], [1, 0], [1, 1]];
            var done = false;
            var resultGenome;
            var genc = 0;
            //die();
            do {
                //for each species
                for (var i = 0; i < pool.species.length; i++) {
                    var species = pool.species[i];

                    //for each genome in this species
                    for (var j = 0; j < species.genomes.length; j++) {
                        var genome = species.genomes[j];
                        generateNeuroNetwork(genome);

                        //for each inputvec
                        var c = 0;
                        for (var k = 0; k < inputvec.length; k++) {
                            //console.log("input=", inputvec[k]);
                            var output = evaluateNeuroNetwork(genome.network, inputvec[k]);

                            var answer = output[0];


                            var correctAnswer = xor(inputvec[k][0], inputvec[k][1]);
                            //console.log("output=", answer, "correct=", correctAnswer);
                            if (answer === correctAnswer) {
                                c++;
                            }
                        }
                        genome.fitness = c / 4;
                        if (genome.fitness > pool.maxFitness) {
                            pool.maxFitness = genome.fitness;
                            console.log("new max fitness = " + pool.maxFitness);
                        }
                        if (c === 4) {
                            done = true;
                            resultGenome = genome;
                            break;
                        }


                    }
                    if (done === true) {
                        break;
                    }
                    species.age++;
                }
                newGeneration();
                genc++;
                console.log("new gen g", pool);
            } while (done !== true && genc < 30);
            console.log(resultGenome);
        </script>
    </body>
</html>
